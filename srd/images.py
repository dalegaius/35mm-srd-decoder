"""Find and read SR-D blocks in arbitrary images (film photos with several blocks, scans, single-block close-ups).

For each image: detect fine-textured, roughly square regions, rectify each, read the 76 x 76 fixels
with a local threshold, pack into 722 bytes, and keep the orientation (8 candidates) that the
Reed-Solomon outer code accepts. RAW files (DNG, CR2, NEF, ARW...) need the optional `rawpy` package.
"""
import os, numpy as np, cv2
from . import blocks

RAW_EXT = {'.dng', '.cr2', '.cr3', '.nef', '.arw', '.raf', '.orf', '.rw2', '.pef'}
IMG_EXT = {'.jpg', '.jpeg', '.png', '.tif', '.tiff', '.bmp'} | RAW_EXT
CELL = 10                        # pixels per fixel after rectification

def read_gray(path):
    ext = os.path.splitext(path)[1].lower()
    if ext in RAW_EXT:
        import rawpy
        with rawpy.imread(path) as r:
            rgb = r.postprocess(output_bps=16, gamma=(1, 1), no_auto_bright=True, use_camera_wb=True)
        g = rgb.astype(np.float32).mean(2)
    else:
        g = cv2.imread(path, cv2.IMREAD_GRAYSCALE)
        if g is None: raise ValueError(f'cannot read {path}')
        g = g.astype(np.float32)
    g -= g.min(); g /= max(g.max(), 1e-6)
    return g

def find_blocks(g):
    """-> list of rotated rectangles ((cx, cy), (w, h), angle) in full-resolution pixel coordinates.

    Blocks are square areas of very high local contrast at the fixel scale. Thin edges (perforations,
    analog track) are removed by a morphological opening; film grain has much lower contrast."""
    H, W = g.shape; s = min(1.0, 3000.0 / max(H, W))
    sm = cv2.resize(g, (int(W * s), int(H * s)), interpolation=cv2.INTER_AREA) if s < 1 else g.copy()
    rects = []
    for kw in (3, 5, 9):
        mu = cv2.blur(sm, (kw, kw)); sd = np.sqrt(np.maximum(cv2.blur(sm * sm, (kw, kw)) - mu * mu, 0))
        sd = sd / (np.percentile(sd, 99.5) + 1e-9)
        m = (sd > 0.35).astype(np.uint8) * 255
        m = cv2.morphologyEx(m, cv2.MORPH_CLOSE, np.ones((kw * 2 + 1, kw * 2 + 1), np.uint8))
        for ko in (int(min(sm.shape) / 40), int(min(sm.shape) / 20), int(min(sm.shape) / 10)):
            ko = max(5, ko)
            mo = cv2.morphologyEx(m, cv2.MORPH_OPEN, np.ones((ko, ko), np.uint8))
            cnts, _ = cv2.findContours(mo, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            for c in cnts:
                (cx, cy), (w, h), a = cv2.minAreaRect(c)
                if min(w, h) < max(2.5 * ko, 76 * 1.2) or max(w, h) / max(min(w, h), 1) > 1.25: continue
                if cv2.contourArea(c) < 0.8 * w * h: continue
                hw = max(w, h) / 2
                if cx - hw < 2 or cy - hw < 2 or cx + hw > sm.shape[1] - 2 or cy + hw > sm.shape[0] - 2: continue
                r = ((cx / s, cy / s), (w / s, h / s), a)
                if all(np.hypot(r[0][0] - q[0][0], r[0][1] - q[0][1]) > 0.5 * min(q[1]) for q in rects):
                    rects.append(r)
    return rects

def _pack(fx):
    g = np.zeros((19, 38), np.uint8)
    for bit in range(8):
        col, row = bit // 4, bit % 4
        g |= (fx[row::4, col::2].astype(np.uint8) << bit)
    return g

_TPL = None
def _template():
    global _TPL
    if _TPL is None:
        z = np.load(os.path.join(os.path.dirname(__file__), 'template.npz'))
        m, v = z['mask'], z['value']
        _TPL = []
        for k in range(4):
            for flip in (False, True):
                mm, vv = np.rot90(m, k), np.rot90(v, k)
                if flip: mm, vv = mm[:, ::-1], vv[:, ::-1]
                _TPL.append(((mm * (2 * vv.astype(np.float32) - 1)).astype(np.float32), k, flip, float(mm.sum())))
    return _TPL

PITCH = 7.0                      # target pixels per fixel in the rectified image

def _rectify(g, rect, margin=1.35):
    (cx, cy), (w, h), a = rect; side = max(w, h) * margin
    S = int(76 * margin * PITCH)
    box = cv2.boxPoints(((cx, cy), (side, side), a)).astype(np.float32)
    o = np.argsort(np.arctan2(box[:, 1] - cy, box[:, 0] - cx)); box = box[o]
    dst = np.float32([[0, 0], [S, 0], [S, S], [0, S]])
    src = cv2.resize(g, None, fx=1, fy=1) if True else g
    return cv2.warpPerspective(src, cv2.getPerspectiveTransform(box, dst), (S, S), flags=cv2.INTER_AREA)

def _lattice(w, pX, pY, px, py):
    n = int((w.shape[0] - py) / pY) - 1; m = int((w.shape[1] - px) / pX) - 1
    ys = (py + (np.arange(n) + 0.5) * pY).astype(int); xs = (px + (np.arange(m) + 0.5) * pX).astype(int)
    I = cv2.blur(w, (3, 3)); L = I[np.ix_(ys, xs)].astype(np.float32)
    return L > cv2.blur(L, (15, 15))

def _coarse(w, hint=None):
    """pitch, phase, offset and orientation from the fixed patterns, on a rough uniform lattice.
    hint = (pitch, k, flip) from a block already read in the same image: narrows the search."""
    best = None
    prange = np.arange(PITCH * 0.82, PITCH * 1.18, PITCH * 0.012) if hint is None else hint[0] * np.arange(0.96, 1.041, 0.012)
    for p in prange:
        for ph in ((0, 0), (0.5, 0), (0, 0.5), (0.5, 0.5)):
            L = _lattice(w, p, p, ph[0] * p, ph[1] * p)
            if L.shape[0] < 76 or L.shape[1] < 76: continue
            Lp = 2 * L.astype(np.float32) - 1
            for tpl, k, flip, nk in _template():
                if hint is not None and (k, flip) != (hint[1], hint[2]): continue
                _, mx, _, loc = cv2.minMaxLoc(cv2.matchTemplate(Lp, tpl, cv2.TM_CCORR))
                if best is None or mx / nk > best[0]: best = (mx / nk, p, ph[0] * p, ph[1] * p, loc, k, flip)
    return best

def _orient(img, k, flip):
    """apply the inverse of the template transform to the image (so fixel (0,0) is top-left)."""
    if flip: img = img[:, ::-1]
    return np.ascontiguousarray(np.rot90(img, -k))

_REGIONS = [(0, 0), (0, 68), (68, 0), (68, 68), (31, 31)]      # 8x8 / 14x14 fixed regions: corners, logo

def _refine(w, p, px, py, ox, oy):
    """locate each fixed region precisely (sub-pixel) and fit a homography fixel -> pixel."""
    z = np.load(os.path.join(os.path.dirname(__file__), 'template.npz')); m, v = z['mask'], z['value']
    wn = (w - cv2.blur(w, (int(6 * p) | 1, int(6 * p) | 1))).astype(np.float32)
    src, dst = [], []
    for (r0, c0) in _REGIONS:
        size = 14 if (r0, c0) == (31, 31) else 8
        mm = m[r0:r0 + size, c0:c0 + size]; vv = v[r0:r0 + size, c0:c0 + size]
        if mm.sum() < 10: continue
        pat = (mm * (2 * vv.astype(np.float32) - 1))
        tpl = cv2.resize(pat, (int(round(size * p)), int(round(size * p))), interpolation=cv2.INTER_NEAREST)
        ex = px + (ox + c0) * p; ey = py + (oy + r0) * p; R = int(2.5 * p)
        x0, y0 = int(max(0, ex - R)), int(max(0, ey - R))
        win = wn[y0:int(ey + R) + tpl.shape[0], x0:int(ex + R) + tpl.shape[1]]
        if win.shape[0] < tpl.shape[0] or win.shape[1] < tpl.shape[1]: continue
        sc = cv2.matchTemplate(win, tpl, cv2.TM_CCORR)
        _, _, _, (bx, by) = cv2.minMaxLoc(sc)
        dx = dy = 0.0                                                      # sub-pixel parabola
        if 0 < bx < sc.shape[1] - 1:
            a, b, c = sc[by, bx - 1], sc[by, bx], sc[by, bx + 1]; dx = 0.5 * (a - c) / (a - 2 * b + c + 1e-9)
        if 0 < by < sc.shape[0] - 1:
            a, b, c = sc[by - 1, bx], sc[by, bx], sc[by + 1, bx]; dy = 0.5 * (a - c) / (a - 2 * b + c + 1e-9)
        src.append([c0, r0]); dst.append([x0 + bx + dx, y0 + by + dy])
    if len(src) < 4: return None
    Hm, _ = cv2.findHomography(np.float32(src), np.float32(dst), 0)
    return Hm

def _sample_h(w, Hm, p, sharpen=(0.6, 3.0), win=15):
    """read the 76x76 fixel centres through the homography; optional unsharp mask against lens blur."""
    jj, ii = np.meshgrid(np.arange(76) + 0.5, np.arange(76) + 0.5)
    pts = cv2.perspectiveTransform(np.stack([jj, ii], -1).reshape(-1, 1, 2).astype(np.float32), Hm).reshape(76, 76, 2)
    I = w.astype(np.float32)
    if sharpen:
        sig, a = sharpen; I = I + a * (I - cv2.GaussianBlur(I, (0, 0), sig * p))
    L = cv2.remap(I, pts[..., 0].astype(np.float32), pts[..., 1].astype(np.float32), cv2.INTER_LINEAR)
    d = L - cv2.blur(L, (win, win))
    return d > 0, np.abs(d) / (d.std() + 1e-9)

def read_block(g, rect, hint=None):
    """-> (grid 19x38, outer_ok, hint)."""
    w = _rectify(g, rect)
    c = _coarse(w, hint)
    if c is None: return None, 0, None, None
    score, p, px, py, (ox, oy), k, flip = c
    newhint = (p, k, flip)
    best = (None, -1)
    # coarse lattice reading (fallback)
    L = _lattice(w, p, p, px, py)[oy:oy + 76, ox:ox + 76]
    f = np.rot90((L[:, ::-1] if flip else L), -k)
    cand = [(f, np.ones(f.shape))]
    # homography refinement in the canonical orientation
    wo = _orient(w, k, flip)
    co = _coarse_fixed(wo, p)
    if co is not None:
        _, p2, px2, py2, (ox2, oy2) = co
        Hm = _refine(wo, p2, px2, py2, ox2, oy2)
        if Hm is not None:
            for sh in ((0.6, 3.0), (0.4, 2.0), (0.8, 2.0), None):
                cand.insert(len(cand) - 1, _sample_h(wo, Hm, p2, sh))
    best = (None, -1, None)
    for fx, mg in cand:
        if fx.shape != (76, 76): continue
        grid = _pack(fx); _, ok, _, _ = blocks.correct(grid)
        if ok > best[1]: best = (grid, ok, _byte_conf(mg))
        if ok == 8: break
    return best[0], best[1], (newhint if best[1] == 8 else None), best[2]

def _byte_conf(margin):
    """per-byte confidence = smallest decision margin among its 8 fixels"""
    c = np.full((19, 38), np.inf)
    for bit in range(8):
        col, row = bit // 4, bit % 4
        c = np.minimum(c, margin[row::4, col::2])
    return c.astype(np.float32)

def _coarse_fixed(w, p0):
    """same as _coarse but orientation already canonical (template index 0)."""
    tpl, _, _, nk = _template()[0]; best = None
    for p in p0 * np.array([0.99, 1.0, 1.01]):
        for ph in ((0, 0), (0.5, 0), (0, 0.5), (0.5, 0.5)):
            L = _lattice(w, p, p, ph[0] * p, ph[1] * p)
            if L.shape[0] < 76 or L.shape[1] < 76: continue
            _, mx, _, loc = cv2.minMaxLoc(cv2.matchTemplate(2 * L.astype(np.float32) - 1, tpl, cv2.TM_CCORR))
            if best is None or mx / nk > best[0]: best = (mx / nk, p, ph[0] * p, ph[1] * p, loc)
    return best

def grids_from_image(path):
    """-> list of (grid 19x38, confidence 19x38) read from one image. Blocks that the outer RS code
    does not fully accept are kept when most codewords pass, for erasure recovery later."""
    g = read_gray(path); out = []; hint = None
    for r in find_blocks(g):
        grid, ok, h, conf = read_block(g, r, hint)
        if ok < 8 and hint is not None:                   # retry without the hint
            grid, ok, h, conf = read_block(g, r, None)
        if ok == 8: hint = hint or h
        if grid is not None and ok >= 4: out.append((grid, conf))
    return out
