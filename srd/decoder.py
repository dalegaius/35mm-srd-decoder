"""EC11 mantissa decoding (allocation + quantizer) and audio synthesis.

Allocation rule and field order: read from the ZR38000 code by Mike (see README).
Synthesis: evenly-stacked TDAC (Princen-Bradley), as described in US 5,583,962:
half-block 0 = modified DCT, half-block 1 = modified DST, N = 512, hop 256, KBD window alpha 4.
"""
import numpy as np
from .ec11_parse import parse
from .tables import H

WIDTHS = [1, 1, 1, 1, 1, 2, 2, 3, 3, 4, 6, 7, 8, 10, 13, 16, 20, 26, 36, 56]   # coefficients per band (217)
START = np.cumsum([0] + WIDTHS)[:-1]
BTAB = [20, 20, 20, 20, 20, 20, 17, 16, 16, 16, 16, 16, 13, 12, 12, 12, 12, 9, 2, 0]
N, M = 512, 256

def _s20(v):
    v &= 0xFFFFF; return (v - (1 << 20) if v & 0x80000 else v) / (1 << 19)
def _s12(v):
    v &= 0xFFF; return (v - 4096 if v & 0x800 else v) / 2048

_TRI = _BOOK = _B3 = _STEP = None
def _init():
    global _TRI, _BOOK, _B3, _STEP
    if _TRI is not None: return
    _TRI = []
    for a in range(37):
        p = H(0x2B2F + a) & 0xFFFF
        _TRI.append((H(p) & 0xFFFFF, H(p + 1) & 0xFFFFF, H(p + 2) & 0xFFFF))
    _STEP = [(H(0x2B54 + a) & 0xFFFFF) / 2.0 ** 18 for a in range(37)]
    _BOOK, _B3 = {}, {}

def _book(band, prm):
    p = H((H(0x2AFD + band) & 0xFFFF) + prm) & 0xFFFF
    if p not in _BOOK:
        P = H(p) & 0xFFFFF; ent = [H(p + 2 + k) & 0xFFFFF for k in range(1 << P)]
        vals, lens, k = [], [], 0
        while k < len(ent):
            run = 1
            while run * 2 <= len(ent) - k and k % (run * 2) == 0 and all(e == ent[k] for e in ent[k:k + run * 2]):
                run *= 2
            vals += [_s20(ent[k])] * run; lens += [P - int(np.log2(run))] * run; k += run
        _BOOK[p] = (P, vals, lens)
    return _BOOK[p]

def _book3(band):
    p = H(H(0x2AFD + band) & 0xFFFF) & 0xFFFF
    if p not in _B3:
        _B3[p] = [(_s12(H(p + k) >> 20), _s12(H(p + k) >> 8), H(p + k) & 0xFF) for k in range(16)]
    return _B3[p]

def _peek(bits, p, n):
    v = 0
    for k in range(n): v = v * 2 + (int(bits[p + k]) if p + k < len(bits) else 0)
    return v

def _slots(r):
    s = [(0, b, r['exps'][0][b]) for b in range(14)] + [(0, 14 + k, r['exps'][6][k]) for k in range(6)]
    for c in range(1, 5): s += [(c, b, r['exps'][c][b]) for b in range(14)]
    return s + [(5, b, r['exps'][5][b]) for b in range(6)]

def alloc_index(r, b, x):
    xm = max(sum(r['exps'], []))
    v = xm - x + BTAB[b] + (4 if b >= 14 else 0) + (r['f6'] - 32)
    if v > 0: v = max(v, 8)
    v = min(max(v, 0), 36)
    return int(min(v, max(4 * (20 - x), 0)))

SHRINK = 0.0     # restore: fraction of the half-step by which non-zero values are pulled toward zero

def _shrink(v, T):
    return 0.0 if v == 0 else float(np.sign(v) * max(abs(v) - SHRINK * T, 0.0))

def decode_mantissas(bits, r, rng=None, ser=True):
    """-> list of (group, band, values[]) ; values interleaved (half-block 0, half-block 1) per coefficient."""
    _init(); p = r['mant_start']; out = []
    for g, b, x in _slots(r):
        vals = []
        if b < 2:
            for _ in range(1 if b == 0 else 2):
                v = _peek(bits, p, 9); vals.append((v - 512 if v >= 256 else v) / 256); p += 9
        else:
            a = alloc_index(r, b, x); D7, I3, h = _TRI[a]
            for _ in range(WIDTHS[b]):
                if h == 0x355D:
                    vals += [0.0, 0.0]
                elif h == 0x358F:
                    for _ in range(2):
                        v = _peek(bits, p, D7); vals.append((v - (1 << D7) if v >= 1 << (D7 - 1) else v) / (1 << (D7 - 1))); p += D7
                elif h == 0x3535:
                    P, vv, ll = _book(b, I3)
                    for _ in range(2):
                        c = _peek(bits, p, P); vals.append(vv[c]); p += ll[c]
                else:  # $3569: two 3-level values per code
                    v0, v1, L = _book3(b)[_peek(bits, p, 4)]; vals += [v0, v1]; p += L
            if SHRINK and h in (0x3535, 0x3569, 0x358F):
                T = _STEP[a]; vals = [_shrink(v, T) for v in vals]
            if ser and rng is not None and x < 20:       # statistically enhanced reconstruction
                T = _STEP[a]; vals = [v + rng.uniform(-T, T) for v in vals]
        out.append((g, b, vals))
    return out, p

# Per-band gain of the coupled bands, measured against the DTS soundtrack of the same print
# (same mix, same reel, serial 65414). Below 4.7 kHz the decoder already matches DTS to 0.1 dB;
# these figures correct the only remaining deficit, which grows with frequency.
CPL_EQ_DB = [2.4, 3.6, 6.2, 8.9, 9.9, 11.1]     # bands 14 to 19
COUPLING_EQ = True


def coupling_coord(e, m, j):
    g = 10 ** (CPL_EQ_DB[j] / 20) if COUPLING_EQ else 1.0
    return np.sqrt((m + 16) / 32 * 2.0 ** (-e)) * g

def validate(bits, r=None):
    """Strict check for frames recovered with erasures: in-range exponents and a parse that ends at
    the data end (within 16 bits of the last non-zero bit)."""
    try:
        r = r or parse(bits)
        ex = sum(r['exps'], [])
        if min(ex) < 0 or max(ex) > 30: return False
        _, end = decode_mantissas(bits, r, None, False)
        nz = np.nonzero(bits)[0]; last = int(nz[-1]) + 1 if len(nz) else 0
        return 0 <= end - last <= 16
    except Exception:
        return False

def spectra(bits, rng=None, ser=True):
    """-> spec[6, 2, 256] for channels in stream order (L, R, Ls, Rs, C, LFE)."""
    r = parse(bits); dec, _ = decode_mantissas(bits, r, rng, ser)
    sp = np.zeros((6, 2, M)); comp = {}
    for g, b, vals in dec:
        x = r['exps'][g][b] if (g < 5 and b < 14) else (r['exps'][6][b - 14] if b >= 14 else r['exps'][5][b])
        sc = 2.0 ** (-x); c0 = START[b]
        if b == 0: sp[g, 0, 0] = vals[0] * sc; continue
        for i in range(WIDTHS[b]):
            v0, v1 = vals[2 * i] * sc, vals[2 * i + 1] * sc
            if b >= 14: comp[c0 + i] = (v0, v1)
            else: sp[g, 0, c0 + i] = v0; sp[g, 1, c0 + i] = v1
    for ch in range(5):                                   # coupling (bands 14-19, above ~4.7 kHz)
        for j in range(6):
            e, m = r['r1'][ch * 6 + j], r['r2'][ch * 6 + j]
            cf = coupling_coord(e, m, j); b = 14 + j
            for i in range(WIDTHS[b]):
                k = START[b] + i
                if k in comp: sp[ch, 0, k] = comp[k][0] * cf; sp[ch, 1, k] = comp[k][1] * cf
    return sp

_n, _k = np.arange(N), np.arange(M)
KC = np.cos(2 * np.pi / N * np.outer(_n + 128, _k))
KS = np.sin(2 * np.pi / N * np.outer(_n + 128, _k))
W0 = np.ones(M); W0[0] = 0.5
_kw = np.kaiser(N // 2 + 1, np.pi * 4.0); _c = np.cumsum(_kw[:N // 2]); _h = np.sqrt(_c / np.sum(_kw))
KBD = np.concatenate([_h, _h[::-1]])

def synthesize(spec_list):
    """spec_list: list of spec[6,2,256] (one per film block, 512 samples each) -> (6, n) float audio @ 48 kHz."""
    out = np.zeros((6, len(spec_list) * 512 + N))
    for t, sp in enumerate(spec_list):
        for ch in range(6):
            out[ch, t * 512:t * 512 + N] += (KC @ (sp[ch, 0] * W0)) * KBD
            out[ch, t * 512 + M:t * 512 + M + N] += (KS @ sp[ch, 1]) * KBD
    return out
