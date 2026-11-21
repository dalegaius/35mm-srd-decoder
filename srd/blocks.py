"""Block level: 19x38 byte grid -> RS correction -> frame words -> de-whitening.

Grid format: 722 bytes per block, row-major 19 rows x 38 columns, as written by
davidferguson/dolby-srd srd-recognize (fixels 2 wide x 4 high per byte, column-major, LSB first).
"""
import numpy as np, reedsolo as rs
from . import layout


class tables:                      # the symbol layout comes from rules (srd/layout.py), not from a table
    @staticmethod
    def load():
        outer, inner = layout.build()
        return {'outer': outer.reshape(-1), 'inner': inner.reshape(-1)}

rs.init_tables(prim=0x187, generator=2, c_exp=8)
_addr = lambda a: ((int(a) - 1) // 38, (int(a) - 1) % 38)

ERASURE_STEP = 2      # erasure counts tried: 2, 4, ... ; set to 4 for speed on large batches

def _decode_cw(cw, nsym, order=None):
    """RS decode; on failure, retry with the k least-confident symbols as erasures.
    A decode with k erasures and e further errors uses k + 2e of the nsym check symbols; what is left over is
    the only protection against a wrong answer (with k = nsym ANY word 'decodes'). So a result is accepted only
    if k + 2e <= nsym - 4, or if it is erasures-only with k <= nsym - 2."""
    try:
        msg, ecc, _ = rs.rs_correct_msg(cw, nsym, fcr=120, generator=2); return bytes(msg) + bytes(ecc), 0
    except rs.ReedSolomonError:
        pass
    if order is not None:
        for k in range(ERASURE_STEP, nsym - 1, ERASURE_STEP):
            er = list(order[:k])
            try:
                msg, ecc, errata = rs.rs_correct_msg(cw, nsym, fcr=120, generator=2, erase_pos=er)
            except (rs.ReedSolomonError, ValueError):
                continue
            full = bytes(msg) + bytes(ecc)
            extra = sum(1 for i in range(len(cw)) if full[i] != cw[i] and i not in er)   # errors outside the erasures
            if k + 2 * extra <= nsym - 4 or (extra == 0 and k <= nsym - 2):
                return full, k
    return None, 0

def correct_single(grid, conf=None):
    """Outer RS(87,67) x 8, then inner RS(177,157) x 3. conf (19x38, optional): per-byte read
    confidence, used to retry failed codewords with erasures (idea and procedure by Mike).
    Returns (grid, outer_ok, inner_ok, erasures_used)."""
    t = tables.load(); g = grid.copy(); ok = 0; eras = 0
    c = None if conf is None else conf
    for k in range(8):
        offs = t['outer'][k * 87:(k + 1) * 87]; pos = [_addr(a) for a in offs]
        cw = bytearray(int(g[p]) for p in pos)
        order = None if c is None else sorted(range(87), key=lambda i: c[pos[i]])
        full, e = _decode_cw(cw, 20, order)
        if full is None: continue
        for p, b in zip(pos, full): g[p] = b
        ok += 1; eras += e
    iok = 0
    for cw_offs in t['inner'].reshape(3, 177):
        idx = [i for i, a in enumerate(cw_offs) if int(a)]
        pos = {i: _addr(cw_offs[i]) for i in idx}
        cw = bytearray(int(g[pos[i]]) if i in pos else 0 for i in range(177))
        if max(rs.rs_calc_syndromes(cw, 20, fcr=120, generator=2)) == 0:
            iok += 1; continue
        order = None if c is None else sorted(idx, key=lambda i: c[pos[i]])
        full, e = _decode_cw(cw, 20, order)
        if full is None: continue
        for i in idx: g[pos[i]] = full[i]
        iok += 1; eras += e
    return g, ok, iok, eras

def _syn0(cw):
    return max(rs.rs_calc_syndromes(cw, 20, fcr=120, generator=2)) == 0

def correct(grid, conf=None, rounds=3):
    """Product-code decoding: outer RS(87,67) x 8 and inner RS(177,157) x 3, repeated in turn. A codeword that
    decodes fixes its bytes (they become fully trusted), which lowers the error count of the codewords crossing
    it, so a block that fails a single outer+inner pass can still be recovered. conf (19x38, optional): per-byte
    read confidence; failed codewords are retried with their least-confident bytes as erasures.
    Returns (grid, outer_ok, inner_ok, erasures_used), where outer_ok / inner_ok count the codewords that are
    consistent ON THE FINAL GRID (see the end of the function): 8 and 3 together mean the whole block is
    self-consistent, which a miscorrected codeword cannot fake."""
    t = tables.load(); g = grid.copy(); eras = 0
    c = None if conf is None else np.array(conf, np.float32).copy()
    outer = [[_addr(a) for a in t['outer'][k * 87:(k + 1) * 87]] for k in range(8)]
    inner = []
    for cw_offs in t['inner'].reshape(3, 177):
        idx = [i for i, a in enumerate(cw_offs) if int(a)]
        inner.append((idx, {i: _addr(cw_offs[i]) for i in idx}))
    def o_cw(k): return bytearray(int(g[p]) for p in outer[k])
    def i_cw(j):
        idx, pos = inner[j]; return bytearray(int(g[pos[i]]) if i in pos else 0 for i in range(177))
    for rnd in range(rounds):
        progress = False
        if rnd == 1 and sum(_syn0(o_cw(k)) for k in range(8)) < 3: break     # hopeless block: do not insist
        for k in range(8):
            cw = o_cw(k)
            if _syn0(cw): continue
            order = None if c is None else sorted(range(87), key=lambda i: c[outer[k][i]])
            full, e = _decode_cw(cw, 20, order)
            if full is None: continue
            for p, b in zip(outer[k], full):
                g[p] = b
                if c is not None: c[p] = np.inf
            eras += e; progress = True
        for j in range(3):
            cw = i_cw(j)
            if _syn0(cw): continue
            idx, pos = inner[j]
            order = None if c is None else sorted(idx, key=lambda i: c[pos[i]])
            full, e = _decode_cw(cw, 20, order)
            if full is None: continue
            for i in idx:
                g[pos[i]] = full[i]
                if c is not None: c[pos[i]] = np.inf
            eras += e; progress = True
        if not progress: break
    # Final consistency. The inner codewords must have zero syndromes. The outer codewords are checked with a
    # tolerance of 2 symbols each: on some prints a few header bytes (e.g. the case of the letters in the
    # 'DOLBY SR*D' text, used as flags) are written over inner-code symbols after inner encoding, so the inner
    # decoder legitimately restores them and the outer code then sees those bytes as differing.
    iok = sum(_syn0(i_cw(j)) for j in range(3)); ok = 0
    for k in range(8):
        cw = o_cw(k)
        if _syn0(cw): ok += 1; continue
        try:
            _, _, errata = rs.rs_correct_msg(cw, 20, fcr=120, generator=2)
            if len(errata) <= 2: ok += 1
        except (rs.ReedSolomonError, ValueError):
            pass
    return g, ok, iok, eras

def frame_words(g, n=217):
    """n 16-bit frame words, taken in turn from the 3 inner codewords (message symbols, pairs)."""
    t = tables.load()
    S = [[None if int(a) == 0 else int(g[_addr(a)]) for a in cw[:157]] for cw in t['inner'].reshape(3, 177)]
    off = (0, 0, 1)
    out = []
    for w in range(n):
        j, p = w % 3, w // 3
        i0, i1 = 2 * p + off[j], 2 * p + 1 + off[j]
        hi = S[j][i0] if i0 < 157 else None; lo = S[j][i1] if i1 < 157 else None
        out.append(((hi or 0) << 8) | (lo or 0))
    return np.array(out, np.uint32)

def dewhiten(words, seed=0xD944):
    key = (seed + np.arange(len(words), dtype=np.uint32) * 0x2777) & 0xFFFF
    return words ^ key

def block_number(g):
    lower, middle, upper = int(g[11, 25]), int(g[11, 27]), int(g[11, 26])
    return ((upper & 3) << 16) | (middle << 8) | lower, upper >> 2

def audio_bits(grid, conf=None):
    """-> dict(block, reel, version, bits) or None for non-audio / uncorrectable blocks.
    Blocks that needed erasures are kept only if all 11 codewords decode and the frame validates."""
    g, ok, iok, eras = correct(grid, conf)
    if ok < 8 or iok < 3: return None           # the whole block must be self-consistent (all 11 codewords)
    w = dewhiten(frame_words(g))
    if (w[1] >> 8) not in (9, 11) or w[0] not in (0x8000, 0x4000, 0x0000):
        return None
    bits = ((w[2:217, None] >> np.arange(15, -1, -1)) & 1).reshape(-1).astype(np.int8)
    if (w[1] >> 8) == 9:
        from . import ec9
        if not ec9.validate(bits): return None
    elif eras:
        from . import decoder
        if not decoder.validate(bits): return None
    bn, reel = block_number(g)
    raw = frame_words(g, 235)
    audio_no = ((int(raw[218]) >> 8 & 3) << 16) | ((int(raw[218]) & 0xFF) << 8) | (int(raw[219]) & 0xFF)
    return dict(block=bn, audio_block=audio_no, reel=reel, recovered=bool(eras), mode=int(w[0]), version=int(w[1] >> 8), bits=bits, inner_ok=iok)
