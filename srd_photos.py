#!/usr/bin/env python3
"""Des photos de la piste SR-D au son : lit les sorties de srd_prepare.py (points.npz + rapport.csv).

Usage :
    python srd_photos.py dossier_envoi_01 [dossier_envoi_02 ...] -o film
    (chaque dossier est un zip « envoi_XX » décompressé ; il suffit qu'un seul contienne points.npz et rapport.csv)

Ce que fait le script, dans l'ordre :
  1. chaque bloc de chaque photo est lu de cinq façons (seuil local plus ou moins large, avec ou sans
     accentuation) ; on garde la première lecture que le code correcteur valide entièrement
     (8 mots de code externes et 3 internes cohérents) ;
  2. dans chaque photo, les blocs décodés donnent le numéro des blocs voisins qui ont échoué ;
  3. un bloc qui a échoué dans toutes les photos où il apparaît est retenté en ADDITIONNANT ses lectures
     (deux lectures médiocres en donnent souvent une bonne) ; il n'est gardé que s'il est cohérent et porte
     un numéro encore manquant ;
  4. écriture de  <sortie>_trames.npz        les trames audio vérifiées (EC9 ou EC11)
                  <sortie>_manquants.csv     les blocs encore manquants, avec les photos où ils apparaissent
                  <sortie>_5.1.wav, _stereo.wav, _centre.wav   le son (48 kHz), sauf avec --no-audio

Dépendances : numpy, opencv-python, reedsolo, et le dossier srd/ du dépôt à côté de ce script.
"""
import argparse, collections, csv, itertools, os, sys, time
import numpy as np, cv2

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from srd import blocks

READS = ((15, 0.0), (9, 0.0), (25, 0.0), (15, 1.5), (9, 1.5))          # (fenêtre du seuil local, accentuation)


def pack(bits):
    g = np.zeros((19, 38), np.uint8)
    for bit in range(8):
        col, row = bit // 4, bit % 4
        g |= (bits[row::4, col::2].astype(np.uint8) << bit)
    return g


def conf(d):
    m = np.abs(d) / (d.std() + 1e-9); c = np.full((19, 38), np.inf)
    for bit in range(8):
        col, row = bit // 4, bit % 4
        c = np.minimum(c, m[row::4, col::2])
    return c.astype(np.float32)


def soft(v, win, sh):
    """valeurs des points -> écart au seuil local, normalisé (positif = point clair)"""
    v = v.astype(np.float32)
    x = v if sh == 0 else v + sh * (v - cv2.GaussianBlur(v, (0, 0), 1.0))
    d = x - cv2.blur(x, (win, win))
    return d / (d.std() + 1e-9)


def decode_soft(d):
    """-> None, ou dict(kind='audio'|'autre', an, version, mode, reel, bits) pour un bloc entièrement cohérent"""
    gc, ok, iok, er = blocks.correct(pack(d > 0), conf(d))
    if ok < 8 or iok < 3:
        return None
    w = blocks.dewhiten(blocks.frame_words(gc, 235)); version, mode = int(w[1]) >> 8, int(w[0])
    raw = blocks.frame_words(gc, 235)
    an = ((int(raw[218]) >> 8 & 3) << 16) | ((int(raw[218]) & 0xFF) << 8) | (int(raw[219]) & 0xFF)
    bn, reel = blocks.block_number(gc)
    if version not in (9, 11) or mode not in (0x8000, 0x4000, 0x0000):
        return dict(kind='autre', an=an, version=version, mode=mode, reel=reel, bits=None)
    bits = ((w[2:217, None] >> np.arange(15, -1, -1)) & 1).reshape(-1).astype(np.int8)
    if version == 9:
        from srd import ec9
        if not ec9.validate(bits):
            return None
    return dict(kind='audio', an=an, version=version, mode=mode, reel=reel, bits=bits)


def _one(v):
    for win, sh in READS:
        r = decode_soft(soft(v, win, sh))
        if r is not None:
            return r
    return None


def _init(step):
    blocks.ERASURE_STEP = step


def load(dirs):
    """-> liste de dict(file, bloc, position, qualite, vals) sans doublon"""
    names, vals, rep = [], [], {}
    for d in dirs:
        p = os.path.join(d, 'points.npz')
        if os.path.exists(p):
            z = np.load(p); names += [str(n) for n in z['noms']]; vals.append(z['valeurs'])
        p = os.path.join(d, 'rapport.csv')
        if os.path.exists(p):
            for r in csv.DictReader(open(p, encoding='utf-8')):
                rep[r['bloc']] = r
    if not vals:
        raise SystemExit('points.npz introuvable dans les dossiers donnés')
    vals = np.concatenate(vals); seen, out = set(), []
    for i, n in enumerate(names):
        if n in seen:
            continue
        seen.add(n); r = rep.get(n, {})
        photo = r.get('file', n.rsplit('_b', 1)[0]); pos = int(r.get('position', n.rsplit('_b', 1)[1]))
        out.append(dict(file=photo, bloc=n, position=pos, qualite=float(r.get('qualite', 0) or 0), vals=vals[i]))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('dirs', nargs='+'); ap.add_argument('-o', '--out', default='srd')
    ap.add_argument('--workers', type=int, default=None, help='processus en parallèle (défaut : tous les coeurs)')
    ap.add_argument('--fast', action='store_true', help='moins d essais d effacements (2 fois plus rapide, quelques blocs en moins)')
    ap.add_argument('--no-audio', action='store_true', help='ne pas écrire les fichiers WAV')
    ap.add_argument('--no-ser', action='store_true', help='EC11 : sans le bruit de reconstruction')
    ap.add_argument('--conceal', choices=['repeat', 'predict'], default='repeat',
                    help="blocs manquants : 'repeat' rejoue le bloc précédent (défaut), 'predict' = prédiction linéaire")
    ap.add_argument('--level', choices=['loud', 'peak'], default='loud',
                    help="niveau de sortie : 'loud' = niveau moyen à -20 dBFS (défaut), 'peak' = crête à 0,89")
    ap.add_argument('--limit', type=int, default=0, help='ne traiter que les N premières photos (essai)')
    a = ap.parse_args()
    step = 4 if a.fast else 2; _init(step)
    B = load(a.dirs)
    photos = sorted({b['file'] for b in B})
    if a.limit:
        keep = set(photos[:a.limit]); B = [b for b in B if b['file'] in keep]; photos = sorted(keep)
    print(f'{len(B)} blocs dans {len(photos)} photos')

    # 1. chaque bloc seul
    t0 = time.time()
    import multiprocessing as mp
    with mp.get_context('spawn').Pool(a.workers, initializer=_init, initargs=(step,)) as pool:
        res = []
        for k, r in enumerate(pool.imap(_one, (b['vals'] for b in B), chunksize=16)):
            res.append(r)
            if (k + 1) % 500 == 0 or k + 1 == len(B):
                print(f'  {k + 1}/{len(B)} blocs lus, {sum(1 for x in res if x and x["kind"] == "audio")} trames audio '
                      f'({time.time() - t0:.0f} s)', flush=True)
    for b, r in zip(B, res):
        b['res'] = r
    audio = [b for b in B if b['res'] and b['res']['kind'] == 'audio']
    if not audio:
        raise SystemExit('aucune trame audio décodée')
    reel = collections.Counter(b['res']['reel'] for b in audio).most_common(1)[0][0]
    found = {}
    for b in audio:
        if b['res']['reel'] == reel:
            found.setdefault(b['res']['an'], b['res'])
    nums = sorted(found); runs, s = [], [nums[0]]
    for n in nums[1:]:
        if n - s[-1] > 64: runs.append(s); s = [n]
        else: s.append(n)
    runs.append(s); run = max(runs, key=len); lo, hi = run[0], run[-1]
    print(f'bobine {reel} : {len(run)} trames sur {hi - lo + 1} positions (blocs {lo} à {hi}) après la lecture bloc par bloc')

    # 2. numéro prévu des blocs non décodés, d'après les blocs décodés de la même photo
    byphoto = collections.defaultdict(list)
    for b in B:
        byphoto[b['file']].append(b)
    steps = collections.Counter()
    for bl in byphoto.values():
        ok = [b for b in bl if b['res'] and b['res']['kind'] == 'audio' and b['res']['reel'] == reel]
        for x, y in itertools.combinations(ok, 2):
            if abs(x['position'] - y['position']) == 1:
                steps[(y['res']['an'] - x['res']['an']) * (y['position'] - x['position'])] += 1
    direction = 1 if steps.get(1, 0) >= steps.get(-1, 0) else -1
    cand = collections.defaultdict(list)                       # numéro prévu -> blocs non décodés
    for bl in byphoto.values():
        ok = [b for b in bl if b['res'] and b['res']['kind'] == 'audio' and b['res']['reel'] == reel]
        for b in bl:
            if b['res'] is not None:
                continue
            for r in ok:                                        # chaque bloc décodé de la photo donne une prévision
                cand[r['res']['an'] + direction * (b['position'] - r['position'])].append(b)

    # 3. addition des lectures d'un même bloc
    missing = [n for n in range(lo, hi + 1) if n not in found]
    rec = 0; t0 = time.time()
    if missing: print(f'{len(missing)} blocs manquants : essai en additionnant les lectures (sur un seul coeur)...', flush=True)
    for im, m in enumerate(missing):
        if im and im % 10 == 0:
            print(f'  {im}/{len(missing)} blocs manquants traités, {rec} récupérés ({time.time() - t0:.0f} s)', flush=True)
        if m in found:
            continue
        pool_ = {id(b): b for k in (-1, 0, 1) for b in cand.get(m + k, [])}   # un bloc non audio peut décaler d'un cran
        idx = sorted(pool_.values(), key=lambda b: -b['qualite'])[:7]
        items = [p for p in itertools.combinations(idx, 2) if p[0]['file'] != p[1]['file']]
        items += [p for p in itertools.combinations(idx, 3) if len({x['file'] for x in p}) == 3]
        got = False
        for win, sh in READS[:4]:
            for it in items:
                d = sum(soft(b['vals'], win, sh) * max(0.05, (b['qualite'] - 50) / 50) ** 2 for b in it)
                r = decode_soft(d)
                if r is not None and r['kind'] == 'audio' and r['reel'] == reel and r['an'] not in found and lo <= r['an'] <= hi:
                    found[r['an']] = r; rec += 1; got = got or r['an'] == m
                if got: break
            if got: break
    left = [n for n in range(lo, hi + 1) if n not in found]
    print(f'{rec} blocs récupérés en additionnant les lectures ({time.time() - t0:.0f} s) ; '
          f'{hi - lo + 1 - len(left)} trames sur {hi - lo + 1} ({100 * (1 - len(left) / (hi - lo + 1)):.1f} %), {len(left)} manquantes')

    # 4. sorties
    ks = [n for n in range(lo, hi + 1) if n in found]
    np.savez_compressed(a.out + '_trames.npz', K=np.array(ks), B=np.array([found[n]['bits'] for n in ks], np.int8),
                        version=np.array([found[n]['version'] for n in ks]), mode=np.array([found[n]['mode'] for n in ks]), reel=reel)
    with open(a.out + '_manquants.csv', 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f); w.writerow(['bloc_audio', 'temps_s', 'photos (photo:position:qualite)'])
        for m in left:
            c = {b['bloc']: b for k in (-1, 0, 1) for b in cand.get(m + k, [])}
            w.writerow([m, round((m - lo) * 512 / 48000, 2),
                        ' ; '.join(f"{b['file']}:{b['position']}:{b['qualite']}" for b in sorted(c.values(), key=lambda b: b['file']))])
    print('écrit :', a.out + '_trames.npz,', a.out + '_manquants.csv')
    if not a.no_audio: print('décodage du son (quelques secondes par minute de film)...', flush=True)
    if a.no_audio:
        return
    from srd import decoder
    from srd.plc import conceal_tdac
    import srd_decode
    rng = None if a.no_ser else np.random.default_rng(0); specs, holes = [], []
    for t, n in enumerate(range(lo, hi + 1)):
        r = found.get(n); sp = None
        if r is not None:
            try:
                if r['version'] == 9:
                    from srd import ec9
                    sp = ec9.spectra(r['bits'], rng, not a.no_ser)
                else:
                    sp = decoder.spectra(r['bits'], rng, not a.no_ser)
            except Exception:
                sp = None
        if sp is None: holes.append(t)
        specs.append(sp)
    shape = next(s for s in specs if s is not None)
    if a.conceal == 'predict':
        specs = [np.zeros_like(shape) if s is None else s for s in specs]
        audio = conceal_tdac(specs, holes, decoder.synthesize, decoder.KC, decoder.KS, decoder.KBD, decoder.W0) if holes else decoder.synthesize(specs)
    else:
        # bloc manquant : on rejoue le spectre du bloc précédent (jusqu'à 8 blocs de suite, silence au-delà).
        # Mesuré sur des blocs retirés exprès : c'est ce qui conserve le mieux l'énergie et le timbre.
        last, gap, filled = np.zeros_like(shape), 0, []
        for s in specs:
            if s is not None: last, gap = s, 0; filled.append(s)
            else:
                gap += 1; filled.append(last if gap <= 8 else np.zeros_like(shape))
        audio = decoder.synthesize(filled)
    L, R, Ls, Rs, C, LFE = audio
    g = 0.707; st = np.vstack([L + g * C + g * Ls + 0.5 * LFE, R + g * C + g * Rs + 0.5 * LFE])
    amp, gain = srd_decode.level(audio, st, a.level)
    srd_decode.write_wav(a.out + '_5.1.wav', amp(np.vstack([L, R, C, LFE, Ls, Rs])))
    srd_decode.write_wav(a.out + '_stereo.wav', amp(st)); srd_decode.write_wav(a.out + '_centre.wav', amp(C))
    print(f'{audio.shape[1] / 48000:.1f} s de son écrites ({len(holes)} blocs masqués, gain {20 * np.log10(gain):+.1f} dB)')


if __name__ == '__main__':
    main()
