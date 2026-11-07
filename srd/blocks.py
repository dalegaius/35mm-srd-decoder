"""Block level: 19x38 byte grid -> RS correction -> frame words -> de-whitening.

Grid format: 722 bytes per block, row-major 19 rows x 38 columns, as written by
davidferguson/dolby-srd srd-recognize (fixels 2 wide x 4 high per byte, column-major, LSB first).
"""
import numpy as np, reedsolo as rs
from . import tables

rs.init_tables(prim=0x187, generator=2, c_exp=8)
_addr = lambda a: ((int(a) - 1) // 38, (int(a) - 1) % 38)

def _decode_cw(cw, nsym, order=None):
    """RS decode; on failure, retry with the k least-confident symbols as erasures (k = 2..20)."""
    try:
        msg, ecc, _ = rs.rs_correct_msg(cw, nsym, fcr=120, generator=2); return bytes(msg) + bytes(ecc), 0
    except rs.ReedSolomonError:
        pass
    if order is not None:
        for k in range(2, nsym + 1, 2):
            try:
                msg, ecc, _ = rs.rs_correct_msg(cw, nsym, fcr=120, generator=2, erase_pos=list(order[:k]))
                return bytes(msg) + bytes(ecc), k
            except (rs.ReedSolomonError, ValueError):
                continue
    return None, 0

def correct(grid, conf=None):
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
    if ok < 8: return None
    if eras and iok < 3: return None
    w = dewhiten(frame_words(g))
    if (w[1] >> 8) != 11 or w[0] not in (0x8000, 0x4000):
        return None
    bits = ((w[2:217, None] >> np.arange(15, -1, -1)) & 1).reshape(-1).astype(np.int8)
    if eras:
        from . import decoder
        if not decoder.validate(bits): return None
    bn, reel = block_number(g)
    raw = frame_words(g, 235)
    audio_no = ((int(raw[218]) >> 8 & 3) << 16) | ((int(raw[218]) & 0xFF) << 8) | (int(raw[219]) & 0xFF)
    return dict(block=bn, audio_block=audio_no, reel=reel, recovered=bool(eras), mode=int(w[0]), version=int(w[1] >> 8), bits=bits, inner_ok=iok)
