#!/usr/bin/env python3
"""
srd_prepare.py : prepare les photos DNG de piste Dolby SR-D pour l'analyse.

Pour chaque photo :
  - developpe le DNG en niveaux de gris,
  - trouve les perforations et redresse le film,
  - decoupe chaque bloc SR-D (image PNG legere, sans perte),
  - lit les 76 x 76 points de chaque bloc et mesure la qualite de lecture,
  - ecrit un rapport (qualite, bobine, numero de bloc) pour reperer les photos ratees.

Les resultats sont regroupes dans des fichiers zip d'environ 200 Mo, prets a envoyer.

Utilisation :
    python srd_prepare.py DOSSIER_DES_DNG -o DOSSIER_DE_SORTIE
Options :
    --zip-mo 200      taille maximale d'un zip (Mo)
    --workers N       nombre de processus en parallele (par defaut : nombre de coeurs)
    --no-crops        ne garder que les valeurs des points (beaucoup plus leger, moins complet)
"""
import argparse
import csv
import json
import os
import sys
import time
import zipfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import cv2 as cv
import numpy as np

try:
    import rawpy
except ImportError:
    sys.exit("Il manque rawpy : pip install rawpy")

PERF_PITCH_MM = 4.75          # pas des perforations 35 mm
BLOCK_MM = 2.45               # cote approximatif d'un bloc SR-D
GRID = 76                     # 76 x 76 points par bloc

# Motifs fixes d'un bloc (coins, logo, texte "DOLBY SR*D") : servent a caler la grille
# et a mesurer la qualite de lecture. (ligne, colonne, lignes de bits)
PATTERNS = [
    (0, 0, ['1110010', '1110010', '1110010', '0001101', '0001101', '1110010', '0001101']),
    (0, 69, ['1110010', '1110010', '1110010', '0001101', '0001101', '1110010', '0001101']),
    (69, 0, ['1110010', '1110010', '1110010', '0001101', '0001101', '1110010', '0001101']),
    (69, 69, ['1110010', '1110010', '1110010', '0001101', '0001101', '1110010', '0001101']),
    (32, 32, ['000000000000', '000000000000', '011000000110', '011100001110', '011110011110',
              '011110011110', '011110011110', '011110011110', '011100001110', '011000000110',
              '000000000000', '000000000000']),
    (44, 56, ['10000000001101110000', '10001000010011100011', '11110111000101011100',
              '10000010001000000010']),
]
KP = np.array([(t + y, l + x, int(v)) for t, l, rows in PATTERNS
               for y, row in enumerate(rows) for x, v in enumerate(row)])

DIH = [lambda b: b, lambda b: b[:, ::-1], lambda b: b[::-1, :], lambda b: b[::-1, ::-1],
       lambda b: b.T, lambda b: b.T[:, ::-1], lambda b: b.T[::-1, :], lambda b: b.T[::-1, ::-1]]

_U = (np.arange(GRID) + 0.5) / GRID


# ----------------------------------------------------------------- image
def load_gray(path):
    with rawpy.imread(str(path)) as raw:
        rgb = raw.postprocess(use_camera_wb=True, output_bps=8, gamma=(2.222, 4.5))
    return cv.cvtColor(rgb, cv.COLOR_RGB2GRAY)


def find_perfs(g):
    """Centres des perforations, alignes sur une droite."""
    hi = np.percentile(g, 99.5)
    best = None
    for frac in (0.80, 0.72, 0.65, 0.88):
        th = (g > hi * frac).astype(np.uint8)
        th = cv.morphologyEx(th, cv.MORPH_OPEN, np.ones((21, 21), np.uint8))
        n, _, st, cen = cv.connectedComponentsWithStats(th)
        cand = []
        for i in range(1, n):
            x, y, w, h, a = st[i]
            if 200 < w < 1100 and 150 < h < 800 and 0.45 < h / w < 0.95 and a > 0.7 * w * h:
                cand.append(cen[i])
        if len(cand) < 3:
            continue
        cand = np.array(cand)
        # garder la colonne la plus peuplee (meme x a peu pres)
        xs = cand[:, 0]
        groups = [cand[np.abs(xs - x) < 250] for x in xs]
        col = max(groups, key=len)
        col = col[np.argsort(col[:, 1])]
        if best is None or len(col) > len(best):
            best = col
    return best


def sample_h(img, C):
    src = np.float32([[0, 0], [1, 0], [0, 1], [1, 1]])
    H = cv.getPerspectiveTransform(src, np.float32(C))
    uu, vv = np.meshgrid(_U, _U)
    p = np.stack([uu, vv, np.ones_like(uu)], -1) @ H.T
    return cv.remap(img, (p[..., 0] / p[..., 2]).astype(np.float32),
                    (p[..., 1] / p[..., 2]).astype(np.float32), cv.INTER_LINEAR)


def binarize(v):
    return (v > cv.blur(v.astype(np.float32), (15, 15))).astype(np.uint8)


def score(bits):
    return float((bits[KP[:, 0], KP[:, 1]] == KP[:, 2]).mean())


def square(x0, y0, s):
    return np.array([[x0, y0], [x0 + s, y0], [x0, y0 + s], [x0 + s, y0 + s]], float)


def read_block(img, s0, dih=None):
    """Cale la grille 76 x 76 sur le bloc (recherche puis affinage des 4 coins)."""
    c = img.shape[0] / 2
    best = (-1, None, None)
    for ds in range(-12, 13, 3):
        for dy in range(-40, 41, 4):
            for dx in range(-40, 41, 4):
                C = square(c - s0 / 2 + dx, c - s0 / 2 + dy, s0 + ds)
                bits = binarize(sample_h(img, C))
                for k in ([dih] if dih is not None else range(8)):
                    sc = score(DIH[k](bits))
                    if sc > best[0]:
                        best = (sc, C, k)
    cur, C, k = best
    for step in (6, 3, 1.5, 0.75, 0.4):
        for _ in range(4):
            improved = False
            for i in range(4):
                for j in range(2):
                    for d in (-step, step):
                        C2 = C.copy()
                        C2[i, j] += d
                        sc = score(DIH[k](binarize(sample_h(img, C2))))
                        if sc > cur:
                            cur, C, improved = sc, C2, True
            if not improved:
                break
    vals = sample_h(img, C)
    return cur, k, C, vals


def byte(bits, c, r):
    v = n = 0
    for x in range(2):
        for y in range(4):
            if bits[r * 4 + y, c * 2 + x]:
                v |= 1 << n
            n += 1
    return v


def block_info(bits):
    asc = ''.join(chr(b) if 32 <= b < 127 else '?' for b in
                  (byte(bits, c, 11) for c in (29, 28, 31, 30, 33, 32, 35, 34, 37, 36)))
    lo, mid, up = byte(bits, 25, 11), byte(bits, 27, 11), byte(bits, 26, 11)
    # meme calcul que le lecteur de reference (dolby-srd) ; les trois octets bruts sont aussi gardes
    num = (((up << 16 & 3) | 65535) & (mid << 8 | 255 | 255 << 16) & (lo | 65535 << 8))
    return asc, num, up >> 2, f"{lo:02x}{mid:02x}{up:02x}"


# ----------------------------------------------------------------- une photo
def process(path, out_dir, keep_crops, dih_hint):
    t0 = time.time()
    g = load_gray(path)
    P = find_perfs(g)
    if P is None or len(P) < 2:
        return {"file": path.name, "error": "perforations introuvables"}, []
    d = np.diff(P, axis=0)
    ang = float(np.degrees(np.arctan2(d[:, 0].mean(), d[:, 1].mean())))
    pitch = float(np.median(np.hypot(d[:, 0], d[:, 1])))
    s0 = pitch / PERF_PITCH_MM * BLOCK_MM
    half = int(pitch * 0.31)
    gf = g.astype(np.float32)
    rows, dih = [], dih_hint
    for i in range(len(P) - 1):
        if abs(np.hypot(*d[i]) - pitch) > 0.15 * pitch:
            continue                       # perforation manquante entre les deux
        cx, cy = (P[i] + P[i + 1]) / 2
        y0, x0 = int(cy) - half, int(cx) - half
        if y0 < 0 or x0 < 0 or y0 + 2 * half > g.shape[0] or x0 + 2 * half > g.shape[1]:
            continue
        crop = gf[y0:y0 + 2 * half, x0:x0 + 2 * half]
        M = cv.getRotationMatrix2D((half, half), -ang, 1)
        crop = cv.warpAffine(crop, M, (2 * half, 2 * half), borderValue=0)
        work = cv.GaussianBlur(crop, (0, 0), max(0.5, s0 / GRID * 0.15))
        sc, k, C, vals = read_block(work, s0, dih)
        if dih is None and sc > 0.9:
            dih = k                        # meme orientation pour tout le film
        bits = DIH[k](binarize(vals))
        asc, num, reel, raw = block_info(bits)
        name = f"{path.stem}_b{i:02d}"
        if keep_crops:
            cv.imwrite(str(out_dir / "blocs" / f"{name}.png"), np.clip(crop, 0, 255).astype(np.uint8))
        rows.append({"file": path.name, "bloc": name, "position": i, "qualite": round(sc * 100, 2),
                     "orientation": k, "texte": asc, "bobine": reel, "numero": num, "octets_numero": raw,
                     "coins": json.dumps(np.round(C, 2).tolist()),
                     "_vals": DIH[k](np.clip(vals, 0, 255).astype(np.uint8))})
    info = {"file": path.name, "angle": round(ang, 2), "px_par_point": round(s0 / GRID, 2),
            "blocs": len(rows), "secondes": round(time.time() - t0, 1),
            "qualite_mediane": round(float(np.median([r["qualite"] for r in rows])), 2) if rows else 0}
    return info, rows


def _worker(args):
    path, out_dir, keep_crops, dih = args
    try:
        return process(Path(path), Path(out_dir), keep_crops, dih)
    except Exception as e:
        return {"file": Path(path).name, "error": repr(e)}, []


# ----------------------------------------------------------------- programme
def main():
    ap = argparse.ArgumentParser(description="Prepare les photos SR-D pour l'envoi.")
    ap.add_argument("dossier", help="dossier contenant les DNG")
    ap.add_argument("-o", "--sortie", required=True)
    ap.add_argument("--zip-mo", type=int, default=200)
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--no-crops", action="store_true")
    a = ap.parse_args()

    src = Path(a.dossier)
    files = sorted([p for p in src.iterdir() if p.suffix.lower() == ".dng"])
    if not files:
        sys.exit(f"Aucun DNG dans {src}")
    out = Path(a.sortie)
    (out / "blocs").mkdir(parents=True, exist_ok=True)

    print(f"{len(files)} photos a traiter avec {a.workers} processus...")
    # la premiere photo fixe l'orientation du film, les autres en profitent (plus rapide)
    info0, rows0 = _worker((files[0], out, not a.no_crops, None))
    dih = rows0[0]["orientation"] if rows0 else None
    infos, allrows = [info0], list(rows0)
    print(f"  {info0['file']} : {info0.get('blocs', 0)} blocs, qualite {info0.get('qualite_mediane', '-')} %")

    with ProcessPoolExecutor(a.workers) as ex:
        futs = [ex.submit(_worker, (f, out, not a.no_crops, dih)) for f in files[1:]]
        for n, fu in enumerate(as_completed(futs), 2):
            info, rows = fu.result()
            infos.append(info)
            allrows += rows
            msg = info.get("error") or f"{info['blocs']} blocs, qualite {info['qualite_mediane']} %"
            print(f"  [{n}/{len(files)}] {info['file']} : {msg}")

    allrows.sort(key=lambda r: r["bloc"])
    # valeurs des points : un seul fichier compact
    np.savez_compressed(out / "points.npz",
                        noms=np.array([r["bloc"] for r in allrows]),
                        valeurs=np.stack([r["_vals"] for r in allrows]) if allrows else np.zeros((0, 76, 76)))
    with open(out / "rapport.csv", "w", newline="", encoding="utf-8") as f:
        cols = ["file", "bloc", "position", "qualite", "orientation", "texte", "bobine", "numero",
                "octets_numero", "coins"]
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(allrows)
    with open(out / "photos.json", "w", encoding="utf-8") as f:
        json.dump(sorted(infos, key=lambda i: i["file"]), f, indent=1, ensure_ascii=False)

    # resume
    q = np.array([r["qualite"] for r in allrows]) if allrows else np.array([0])
    bad = sorted({i["file"] for i in infos if i.get("error") or i.get("qualite_mediane", 0) < 97})
    print(f"\n{len(allrows)} blocs extraits. Qualite : mediane {np.median(q):.1f} %, "
          f"{(q >= 98).mean() * 100:.0f} % des blocs au-dessus de 98 %.")
    if bad:
        print(f"Photos a verifier (qualite faible ou erreur) : {', '.join(bad[:20])}"
              + (" ..." if len(bad) > 20 else ""))

    # zips prets a envoyer
    items = [out / "rapport.csv", out / "photos.json", out / "points.npz"]
    items += sorted((out / "blocs").glob("*.png"))
    limit, part, size, z = a.zip_mo * 1_000_000, 1, 0, None
    for p in items:
        if z is None or (size + p.stat().st_size > limit and size > 0):
            if z:
                z.close()
            z = zipfile.ZipFile(out / f"envoi_{part:02d}.zip", "w", zipfile.ZIP_STORED)
            part, size = part + 1, 0
        z.write(p, p.relative_to(out))
        size += p.stat().st_size
    if z:
        z.close()
    print(f"Fichiers a m'envoyer : {', '.join(sorted(x.name for x in out.glob('envoi_*.zip')))}")


if __name__ == "__main__":
    import multiprocessing
    multiprocessing.set_start_method("spawn", force=True)   # evite les blocages de rawpy sous Linux
    main()
