#!/usr/bin/env python3
"""Extrait les blocs SR-D d'une vidéo de scan 35 mm en grilles 19x38, avec leur carte de confiance.

Usage : python srd_extract_video.py video.mp4 -o sortie.npz
        python srd_decode.py sortie.npz -o film
Sortie (.npz) : grids (n,19,38) uint8, conf (n,19,38), frame, kind (0 = bloc entier, 1 = bloc recollé), q.

Blocs coupés : quand la pellicule n'avance pas d'exactement 4 blocs entre deux images, un bloc se retrouve
coupé entre le bas d'une image et le haut de la suivante. Le script garde ces morceaux et les recolle d'après
les numéros des blocs entiers voisins (il faut pour cela le module srd du dépôt, à côté du script) ; un
recollage n'est gardé que s'il est corrigé sans effacement douteux et porte le numéro attendu.
Vidéos AV1 : si OpenCV ne sait pas les lire, convertir d'abord, par exemple
        ffmpeg -i video.mp4 -map 0:v:0 -c:v ffv1 video.mkv
Dépendances : numpy, opencv-python.
"""
import argparse, cv2, numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('video'); ap.add_argument('-t', '--template', default=None)
ap.add_argument('-o', '--out', default='srd_grids.npz')
ap.add_argument('--x0', type=float, default=0.03, help='début de la bande SR-D (fraction de la largeur)')
ap.add_argument('--x1', type=float, default=0.17, help='fin de la bande SR-D (fraction de la largeur)')
ap.add_argument('--start', type=int, default=0); ap.add_argument('--count', type=int, default=0)
ap.add_argument('--fast', action='store_true', help='une seule lecture par bloc (plus rapide, moins de blocs récupérés)')
ap.add_argument('--no-stitch', action='store_true', help='ne pas recoller les blocs coupés entre deux images')
a = ap.parse_args()

TEMPLATE = 'eNqrZ4AA/v8foKz//8kVYhgBQEcB6NEPDAIOzAdAfmdg5GCwaAALKTAwnwALMYEVMsPDxR4SMApQIUaEEANYIywc8Qg1GPD4f9kANtbB/d88sJAJ7/t9ICEOZgb/39+Hdqj+p2IqfNK0OYzB9KaEu55HXsZr+e0PPxY9iTwnUdg94fNzxfX2jHvkE6zNeaXUFpk5/a1i+nDP3X+r6QSjR0uLpNbOcWnImqM44QL/0obX3W3a7gxlGvkf1+RvspAyKZ7Rb6cTERTV+X3f93s3WZgSdz1m1JJaWhG2OOGL8JEt7BvXuhevUp5x3qf4zdc9V5vnqk2XYqx8IKCj/DSqXdMnrOGExHPm9sfGR14XN755MTv4ng/f2snO9d0vL3t+r1h2m6/USnhBtaq21NYtBtZiWdJT0qxltzT/OyumuSHxX6ODdt2vv+UT09qXrlrhNf1tFWuj2bO9aqYK+uUqvHeNbrDnvndRimyaxcFiZn1e/7b/0RdzO9o/ljE/+Ox0/LvQsp3saiu3/ltyQPvn1jVtLGacqg1L70fv+/45qHFF6B/ZFVkbd+5nZ2Cd5vC7Rfvv8gAGhvTn3ffl/R8kJF3szamSfln2na8m67/4g+3xMyrflZUaCCupNLvP+7b0Ksuzr0/lKh9Grfl0+NmkPvN54XwTK/StricVvA4KdvgZN9nqGNvcqG7mJEGFAgbWpsD21HVnzzEwKFdINdtFhdTJ3V7SYM6S+Oyt6BVmjZbFcd4Wz2as9Poueijs5J3zy2asUGLyjJ1wMUdsynfmDvPYFZcEZsae36x7dnl03TL5aaE75Zq1boS/vMollPfRoO6jTd51K3XVpe6/fhhrXn0W/efZAgEWs3v2P6rVHxqZfw6/cHah+5Od9QwaRgr2r48avorbl3FW+V7ll3+7XsqaKFfEnr31fk1cyJNIoffJjO+dmD4Hef78+GpKu+2NT4/stv1Vmcz1OnXChKmWURuzqs0/29lfVYuOSNOOWRThEVyxeKHXieurWY9X/35/Sc/QpkmpZh/ThztN3xwYFvI/zdQvVnDg89hk89Tt1YzvRU9SG8UDqm02LdfzYvuozyi0dVeR1HchlvfGaW8ZLi6eF2P3eq5zDO+T4xKbE9pjJNwVl+4QmmyYnqTKCwAPX5nt'
if a.template:
    z = np.load(a.template); M, V = z['mask'].astype(bool), z['value'].astype(np.uint8)
else:
    import base64, zlib
    raw = np.frombuffer(zlib.decompress(base64.b64decode(TEMPLATE)), np.uint8); h = len(raw) // 2
    M = np.unpackbits(raw[:h])[:76 * 76].reshape(76, 76).astype(bool)
    V = np.unpackbits(raw[h:])[:76 * 76].reshape(76, 76).astype(np.uint8)
U = (np.arange(76) + .5) / 76; uu, vv = np.meshgrid(U, U); PU = np.stack([uu, vv, np.ones_like(uu)], -1)
SRC = np.float32([[0, 0], [1, 0], [0, 1], [1, 1]])
DIH = [lambda b: b, lambda b: b[:, ::-1], lambda b: b[::-1], lambda b: b[::-1, ::-1],
       lambda b: b.T, lambda b: b.T[:, ::-1], lambda b: b.T[::-1], lambda b: b.T[::-1, ::-1]]
TXT = np.frombuffer(b'DOLBY SR*D', np.uint8); COLS = [29, 28, 31, 30, 33, 32, 35, 34, 37, 36]

def pts(C):
    H = cv2.getPerspectiveTransform(SRC, np.float32(C)); p = PU @ H.T
    return (p[..., 0] / p[..., 2]).astype(np.float32), (p[..., 1] / p[..., 2]).astype(np.float32)
def dots(s, C, sh=0, sig=2, win=25):
    ss = s if sh == 0 else s + sh * (s - cv2.GaussianBlur(s, (0, 0), sig))
    x, y = pts(C); v = cv2.remap(ss, x, y, cv2.INTER_LINEAR)
    valid = (y > 1.5) & (y < s.shape[0] - 2.5) & (x > 0) & (x < s.shape[1] - 1)
    w = valid.astype(np.float32); mu = cv2.blur(np.where(valid, v, 0).astype(np.float32), (win, win)) / np.maximum(cv2.blur(w, (win, win)), 1e-3)
    return v - mu, valid
def score(d, valid, k):
    b = DIH[k]((d > 0).astype(np.uint8)); m = DIH[k](valid) & M
    return ((b == V) & m).sum() / max(m.sum(), 1)
def fit(s, C, k, steps):
    cur = score(*dots(s, C), k)
    for st in steps:
        for _ in range(4):
            imp = False
            for i in range(4):
                for j in range(2):
                    for dd in (-st, st):
                        C2 = C.copy(); C2[i, j] += dd; sc = score(*dots(s, C2), k)
                        if sc > cur: cur, C, imp = sc, C2, True
            if not imp: break
    return cur, C
def pack(fx):
    g = np.zeros((19, 38), np.uint8)
    for bit in range(8):
        col, row = bit // 4, bit % 4; g |= (fx[row::4, col::2].astype(np.uint8) << bit)
    return g
def conf(d):
    m = np.abs(d) / (d.std() + 1e-9); c = np.full((19, 38), np.inf)
    for bit in range(8):
        col, row = bit // 4, bit % 4; c = np.minimum(c, m[row::4, col::2])
    return c.astype(np.float32)
def best_read(s, C, k, valid_only=False):
    best = None
    for sh, sig, win in ((0, 1, 15), (2, 2, 25), (3, 2, 25), (5, 3, 25), (8, 3, 35)):
        d, valid = dots(s, C, sh, sig, win); G = pack(DIH[k]((d > 0).astype(np.uint8)))
        e = int((G[11, COLS] != TXT).sum())
        if best is None or e < best[0]: best = (e, G, DIH[k](d).copy(), DIH[k](valid).copy())
    return best

try:
    import os, sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from srd import blocks as _B
except Exception:
    _B=None; print('module srd introuvable : ni choix de la meilleure lecture par le code correcteur, ni recollage')
def _clean(G_, C_):
    """True si le bloc est entièrement cohérent (8 mots de code externes et 3 internes)."""
    gc, ok, iok, er = _B.correct(G_.copy(), C_)
    return ok == 8 and iok == 3
def pick_read(s, C, k):
    """Essaie plusieurs réglages de netteté ; garde le premier dont le bloc est entièrement cohérent, sinon celui
    qui lit le mieux le texte fixe « DOLBY SR*D »."""
    cands = []
    for sh, sig, win in ((0, 1, 15), (2, 2, 25), (3, 2, 25), (5, 3, 25), (8, 3, 35)):
        d, valid = dots(s, C, sh, sig, win); G_ = pack(DIH[k]((d > 0).astype(np.uint8))); dd = DIH[k](d).copy()
        cands.append((int((G_[11, COLS] != TXT).sum()), G_, dd))
    cands.sort(key=lambda t: t[0])
    if _B is not None and not a.fast:
        for e, G_, dd in cands:
            if _clean(G_, conf(dd)): return e, G_, dd
    return cands[0]
cap = cv2.VideoCapture(a.video); N = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
cap.set(cv2.CAP_PROP_POS_FRAMES, a.start); end = N if a.count == 0 else min(N, a.start + a.count)
K = None; PREV = []; out = []; partial = {}
for fi in range(a.start, end):
    ok, f = cap.read()
    if not ok: break
    g = cv2.cvtColor(f, cv2.COLOR_BGR2GRAY).astype(np.float32); Hh, Ww = g.shape; sc = Hh / 2160
    x0, x1 = int(a.x0 * Ww), int(a.x1 * Ww); s = g[:, x0:x1]
    th = (s > np.percentile(s, 97) * 0.9).astype(np.uint8)
    th = cv2.morphologyEx(th, cv2.MORPH_OPEN, np.ones((int(31 * sc), int(31 * sc)), np.uint8))
    n, _, st, cen = cv2.connectedComponentsWithStats(th)
    P = sorted([cen[i] for i in range(1, n) if 150 * sc < st[i][2] < 440 * sc and 120 * sc < st[i][3] < 340 * sc], key=lambda p: p[1])
    if len(P) < 2: continue
    pitch = np.median(np.diff([p[1] for p in P])); side = pitch / 4.75 * 2.45
    cs = [((P[i][0] + P[i + 1][0]) / 2, (P[i][1] + P[i + 1][1]) / 2) for i in range(len(P) - 1)]
    good = []
    for cx, cy in cs:
        C = np.array([[cx - side / 2, cy - side / 2], [cx + side / 2, cy - side / 2], [cx - side / 2, cy + side / 2], [cx + side / 2, cy + side / 2]])
        q = -1
        for Cp in PREV:
            if abs(Cp[:, 1].mean() - cy) < 40 * sc:
                q, C2 = fit(s, Cp.copy(), K, (2 * sc, sc, .5 * sc, .25 * sc))
                if q > 0.95: C = C2; break
        if q <= 0.95:
            best = None
            for k in (range(8) if K is None else (K,)):
                for ds in np.arange(-8, 9, 4) * sc:
                    for dy in np.arange(-20, 21, 4) * sc:
                        for dx in np.arange(-20, 21, 4) * sc:
                            Cc = C + [dx, dy] + np.array([[-1, -1], [1, -1], [-1, 1], [1, 1]]) * ds / 2
                            qq = score(*dots(s, Cc), k)
                            if best is None or qq > best[0]: best = (qq, k, Cc)
            if best[0] < 0.65: continue
            q, C = fit(s, best[2], best[1], (6 * sc, 3 * sc, 1.5 * sc, .75 * sc, .4 * sc, .2 * sc))
            if q > 0.93 and K is None: K = best[1]
        if q < 0.75 or K is None: continue   # les paquets du chargeur n'ont pas le texte DOLBY SR*D, donc un calage plus bas
        e, G, d = pick_read(s, C, K); out.append((fi, 0, q, G, conf(d))); good.append(C)
    if good: PREV = good
    good.sort(key=lambda C: C[:, 1].mean())
    if len(good) >= 2:   # blocs coupés en haut et en bas de l'image
        steps = np.round(np.diff([C[:, 1].mean() for C in good]) / pitch).clip(1)
        D = np.mean([(good[t + 1] - good[t]) / steps[t] for t in range(len(good) - 1)], 0)
        for side_name, Cp in (('top', good[0] - D), ('bot', good[-1] + D)):
            q, Cr = fit(s, Cp.copy(), K, (1.5 * sc, .6 * sc, .25 * sc))
            reads = []
            for sh, sig, win in ((3, 2, 25), (0, 1, 15), (2, 2, 25), (5, 3, 25), (8, 3, 35)):
                d, valid = dots(s, Cr, sh, sig, win); reads.append((DIH[K](d).astype(np.float32), DIH[K](valid).copy()))
            partial[(fi, side_name)] = (q, reads)
    if fi % 50 == 0: print(f'image {fi}/{end}  blocs entiers {sum(1 for o in out if o[1] == 0)}', flush=True)



full=[o for o in out if o[1]==0]
if not full: raise SystemExit("aucun bloc SR-D trouvé (bande mal cadrée ? vidéo AV1 illisible par OpenCV ? voir l'aide)")
stitched=[]
if _B is not None and partial and not a.no_stitch:
    def _num(g, c):
        gc, ok, iok, er = _B.correct(g.copy(), c)
        if ok < 8 or iok < 3: return None                                # bloc entièrement cohérent uniquement
        w = _B.dewhiten(_B.frame_words(gc, 235))
        if (int(w[1]) >> 8) not in (9, 11) or int(w[0]) not in (0x0000, 0x4000, 0x8000): return None
        raw = _B.frame_words(gc, 235)
        return ((int(raw[218]) >> 8 & 3) << 16) | ((int(raw[218]) & 0xFF) << 8) | (int(raw[219]) & 0xFF)
    nums={}
    for fi, kd, q, G_, C_ in full:
        n=_num(G_, C_)
        if n is not None: nums.setdefault(fi, []).append(n)
    have={n for v in nums.values() for n in v}
    for fi in sorted({f for f, sd in partial if sd == 'bot'}):
        if (fi+1,'top') not in partial: continue
        qb, rb = partial[(fi,'bot')]; qt, rt = partial[(fi+1,'top')]
        if (rb[0][1] | rt[0][1]).mean() < 0.85: continue
        near = [m for f in (fi, fi+1) for m in nums.get(f, [])]
        for (db, vb), (dt, vt) in zip(rb, rt):                            # mêmes réglages de netteté des deux côtés
            na = np.where(vb, db / (db[vb].std() + 1e-9), 0); nb = np.where(vt, dt / (dt[vt].std() + 1e-9), 0)
            d = np.where(np.abs(na) >= np.abs(nb), na, nb); G_ = pack((d > 0).astype(np.uint8)); C_ = conf(d)
            n = _num(G_, C_)
            if n is None: continue
            if n in have or (near and min(abs(n - m) for m in near) > 8): break   # déjà lu, ou numéro sans rapport
            stitched.append((fi+1, 1, min(qb, qt), G_, C_)); have.add(n); break
out = sorted(full + stitched, key=lambda o: (o[0], o[1]))
np.savez_compressed(a.out, grids=np.stack([o[3] for o in out]), conf=np.stack([o[4] for o in out]),
                    frame=np.array([o[0] for o in out]), kind=np.array([o[1] for o in out]), q=np.array([o[2] for o in out]))
print(f"{len(full)} blocs entiers, {len(stitched)} blocs recollés, écrits dans {a.out}")
