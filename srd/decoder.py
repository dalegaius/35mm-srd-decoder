"""EC11 (SR-D revision 11): mantissa decoding (allocation and quantizers), coupling, and audio synthesis.

Synthesis: evenly stacked TDAC (Princen-Bradley), as described in US patent 5,583,962: half-block 0 is a modified
DCT, half-block 1 a modified DST, N = 512, hop 256, KBD window alpha 4.
"""
import numpy as np
from .ec11_parse import parse
from . import codebooks

WIDTHS = [1, 1, 1, 1, 1, 2, 2, 3, 3, 4, 6, 7, 8, 10, 13, 16, 20, 26, 36, 56]   # coefficients per band (217)
START = np.cumsum([0] + WIDTHS)[:-1]
BTAB = [20, 20, 20, 20, 20, 20, 17, 16, 16, 16, 16, 16, 13, 12, 12, 12, 12, 9, 2, 0]   # base allocation per band
N, M = 512, 256


def quantizer(a):
    """allocation index 0..36 -> ('zero',) | ('pair3',) | ('huff', levels) | ('raw', bits)"""
    if a <= 3: return ('zero',)
    if a <= 8: return ('pair3',)
    if a <= 11: return ('huff', 5)
    if a <= 13: return ('huff', 7)
    if a <= 18: return ('huff', 15)
    if a <= 22: return ('huff', 31)
    if a <= 26: return ('huff', 63)
    if a <= 30: return ('raw', 7)
    if a <= 34: return ('raw', 8)
    return ('raw', 9)


def step(a):
    """quantizer step for this allocation index: one over the number of levels"""
    q = quantizer(a)
    return 1.0 if q[0] == 'zero' else 1 / 3 if q[0] == 'pair3' else 1 / q[1] if q[0] == 'huff' else 1 / (1 << q[1])


def _init():          # kept for callers of the former interface
    pass

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
            a = alloc_index(r, b, x); q = quantizer(a)
            for _ in range(WIDTHS[b]):
                if q[0] == 'zero':
                    vals += [0.0, 0.0]
                elif q[0] == 'raw':
                    nb = q[1]
                    for _ in range(2):
                        v = _peek(bits, p, nb); vals.append((v - (1 << nb) if v >= 1 << (nb - 1) else v) / (1 << (nb - 1))); p += nb
                elif q[0] == 'huff':
                    P, vv, ll = codebooks.huffman(b, q[1])
                    for _ in range(2):
                        c = _peek(bits, p, P); vals.append(vv[c]); p += ll[c]
                else:                                    # two three-level values in one code
                    P, vv, ll = codebooks.pairs(b); c = _peek(bits, p, 4); vals += list(vv[c]); p += ll[c]
            if SHRINK and q[0] != 'zero':
                T = step(a); vals = [_shrink(v, T) for v in vals]
            if ser and rng is not None and x < 20:       # statistically enhanced reconstruction
                T = step(a); vals = [v + rng.uniform(-T, T) for v in vals]
        out.append((g, b, vals))
    return out, p

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
    # Coupling (bands 14-19, above about 4.7 kHz). The composite only gives the SHAPE of each band: it is brought to
    # unit RMS over the band (both half-blocks). Each channel then has its own ABSOLUTE level per band,
    # (m + 16) / 32 * 2^-E with E = e + 3 * c2[channel]: e and m are the channel's 4-bit exponent and mantissa for
    # the band, c2 its 2-bit field. At the last exponent (E = 24) the mantissa has no implicit leading one, so
    # that a level of exactly zero can be sent. Checked against the DTS soundtrack of the same reel: band levels
    # agree within a few tenths of a dB on average, 1.3 dB rms per block.
    for j in range(6):
        b = 14 + j; ks = [START[b] + i for i in range(WIDTHS[b]) if START[b] + i in comp]
        if not ks: continue
        c = np.array([comp[k] for k in ks]) / (2.0 ** (-r['exps'][6][j]))           # mantissas, exponent removed
        rms = np.sqrt(np.mean(c ** 2))
        if rms == 0: continue
        shape = c / rms
        for ch in range(5):
            e, m = r['r1'][ch * 6 + j], r['r2'][ch * 6 + j]; E = e + 3 * r['c2'][ch]
            lvl = (m + 16) / 32 * 2.0 ** (-E) if E < 24 else m / 32 * 2.0 ** (-E)
            for (v0, v1), k in zip(shape, ks): sp[ch, 0, k] = v0 * lvl; sp[ch, 1, k] = v1 * lvl
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
