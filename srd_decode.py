#!/usr/bin/env python3
"""Decode Dolby SR-D (EC11 and EC9) audio from images of the blocks, or from block grids.

Inputs (any mix): folders of images, image files (JPEG, PNG, TIFF, or camera RAW such as DNG),
722-byte grid files (19 x 38, row-major, e.g. output.bin from davidferguson/dolby-srd srd-recognize),
or .npy arrays (n, 19, 38). Blocks seen in several overlapping photos are merged automatically.
Output: <out>_5.1.wav (L R C LFE Ls Rs), <out>_stereo.wav, <out>_centre.wav, 48 kHz 16-bit.
"""
import argparse, os, sys, wave, numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import Counter
from srd import blocks, decoder

def _read_image(path):
    from srd import images
    try: return path, images.grids_from_image(path), None
    except Exception as e: return path, [], str(e)

def load_grids(paths, workers=None, log=print):
    """-> list of (grid 19x38, confidence 19x38 or None)"""
    from srd import images
    files, grids = [], []
    for p in paths:
        if os.path.isdir(p):
            files += sorted(os.path.join(p, f) for f in os.listdir(p) if os.path.splitext(f)[1].lower() in images.IMG_EXT)
        elif os.path.splitext(p)[1].lower() in images.IMG_EXT: files.append(p)
        elif p.endswith('.npz'):                       # srd_extract_video.py output: grids + confidence maps
            z = np.load(p); grids += list(zip(z['grids'].astype(np.uint8), z['conf']))
        elif p.endswith('.npy'): grids += [(g, None) for g in np.load(p).astype(np.uint8).reshape(-1, 19, 38)]
        else: grids += [(g, None) for g in np.frombuffer(open(p, 'rb').read(), np.uint8).reshape(-1, 19, 38)]
    if files:
        log(f'reading {len(files)} image(s) on {workers or os.cpu_count()} core(s)...')
        done = 0; nblk = 0
        import multiprocessing as mp
        with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context('spawn')) as ex:
            futs = [ex.submit(_read_image, f) for f in files]
            for fu in as_completed(futs):
                path, gr, err = fu.result(); done += 1; nblk += len(gr)
                grids += gr
                if err: log(f'  {os.path.basename(path)}: {err}')
                if done % 10 == 0 or done == len(files): log(f'  {done}/{len(files)} images, {nblk} blocks read')
    if not grids: raise SystemExit('nothing to decode')
    return grids

def level(audio, st, mode='loud'):
    """-> function applying the output gain. 'peak': highest sample at 0.89 (at most +24 dB), as before.
    'loud': average level of the stereo mix at -20 dBFS, the few peaks above 0.7 being rounded off smoothly.
    With 'peak' one isolated loud sample sets the level of the whole file, so exports differ a lot in loudness."""
    if mode == 'peak':
        peak = max(np.max(np.abs(audio)), np.max(np.abs(st))) or 1.0; gain = min(0.89 / peak, 16.0)
        return (lambda x: x * gain), gain
    rms = float(np.sqrt(np.mean(st ** 2))) or 1.0; gain = min(0.1 / rms, 64.0)
    def f(x):
        y = x * gain; a = np.abs(y)
        return np.where(a > 0.7, np.sign(y) * (0.7 + 0.29 * np.tanh((a - 0.7) / 0.29)), y)
    return f, gain

def write_wav(path, x):
    x = np.atleast_2d(x)
    with wave.open(path, 'wb') as w:
        w.setnchannels(x.shape[0]); w.setsampwidth(2); w.setframerate(48000)
        w.writeframes((np.clip(x, -1, 1).T * 32767).astype('<i2').tobytes())

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('inputs', nargs='+'); ap.add_argument('-o', '--out', default='srd')
    ap.add_argument('--no-ser', action='store_true', help='disable noise reconstruction of quantization error')
    ap.add_argument('--workers', type=int, default=None, help='parallel processes for image reading')
    ap.add_argument('--max-gap', type=int, default=8, help='missing blocks concealed by repetition up to this gap')
    ap.add_argument('--level', choices=['loud', 'peak'], default='loud',
                    help="output level: 'loud' = average level at -20 dBFS (default), 'peak' = highest sample at 0.89")
    ap.add_argument('--conceal', choices=['auto', 'repeat', 'time'], default='auto',
                    help="missing blocks: 'repeat' the previous frame's spectrum, or 'time' = rebuild the waveform "
                         "from the audio before and after the gap (better when many blocks are missing, e.g. video "
                         "scans; experimental, usually worse). 'auto' = 'repeat'")
    ap.add_argument('--surround-delay', type=int, default=0, metavar='N',
                    help='delay the surround channels by N samples. Off by default: on film all channels are '
                         'aligned. 952 reproduces the surround delay found on the DTS disc of one reel')
    ap.add_argument('--restore', action='store_true', help='optional restoration of codec artifacts (see srd/restore.py)')
    a = ap.parse_args()
    run(a.inputs, a.out, ser=not a.no_ser, workers=a.workers, max_gap=a.max_gap, restore=a.restore,
        surround_delay=a.surround_delay, conceal=a.conceal, level_mode=a.level)

def run(inputs, out, ser=True, workers=None, log=print, max_gap=8, restore=False, surround_delay=0, conceal='auto', level_mode='loud'):
    decoder.SHRINK = 0.25 if restore else 0.0
    if restore: ser = False                          # restoration replaces the uniform noise fill
    class A: pass
    a = A(); a.out, a.no_ser = out, not ser
    grids = load_grids(inputs, workers, log)
    found, skipped, recovered = {}, 0, 0
    for g, conf in grids:
        b = blocks.audio_bits(g.copy(), conf)
        if b is None: skipped += 1; continue
        key = (b['reel'], b['audio_block'])
        if key in found and not found[key]['recovered']: continue      # prefer a clean copy
        found[key] = b; recovered += b['recovered']
    if not found: raise SystemExit('no decodable EC11 or EC9 audio block found')
    reel = Counter(r for r, _ in found).most_common(1)[0][0]
    allnums = sorted(n for r, n in found if r == reel)
    runs, s = [], [allnums[0]]                       # split into continuous runs (block counter wraps at 2^18)
    for n in allnums[1:]:
        if n - s[-1] > 64: runs.append(s); s = [n]
        else: s.append(n)
    runs.append(s); nums = max(runs, key=len)
    log(f'{len(grids)} blocks read, {len(found)} audio frames ({recovered} recovered with erasures), {skipped} non-audio/uncorrectable; reel {reel}; '
          f'{len(runs)} run(s), decoding the longest: blocks {nums[0]}..{nums[-1]} ({len(nums)} present); '
          f'modes {sorted(set(hex(v["mode"]) for v in found.values()))}, revisions {sorted(set(v["version"] for v in found.values()))}')
    rng = None if a.no_ser else np.random.default_rng(0)
    if conceal == 'auto':
        conceal = 'repeat'      # measured on deliberately removed blocks: repetition keeps energy and timbre best
    specs, last, gap, missing, holes = [], None, 0, 0, []
    for n in range(nums[0], nums[-1] + 1):
        b = found.get((reel, n))
        if b is None and conceal == 'time':
            holes.append(len(specs)); missing += 1; specs.append(None); continue
        if b is not None:
            if b['version'] == 9:
                from srd import ec9
                last = ec9.spectra(b['bits'], rng, ser)
            else:
                last = decoder.spectra(b['bits'], rng, ser)
            gap = 0
        else:
            gap += 1; missing += 1
        specs.append(last if (b is not None or gap <= max_gap) else np.zeros_like(last))
    if conceal == 'time':
        shape = next(s for s in specs if s is not None)
        specs = [np.zeros_like(shape) if s is None else s for s in specs]
        log(f'timeline {len(specs)} frames, {missing} missing (rebuilt from the audio around each gap)')
    else:
        log(f'timeline {len(specs)} frames, {missing} missing (previous frame repeated)')
    if restore:
        from srd.restore import restore as _restore
        specs = list(_restore(np.array(specs))); log('restoration applied')
    if conceal == 'time' and holes:
        # the missing block is predicted in the time domain, turned back into MDCT/MDST coefficients and
        # synthesized with its real neighbours, so that their overlap is cancelled properly
        from srd.plc import conceal_tdac
        audio = conceal_tdac(specs, holes, decoder.synthesize, decoder.KC, decoder.KS, decoder.KBD, decoder.W0)
    else:
        audio = decoder.synthesize(specs)
    if surround_delay:                       # stream channels 2 and 3 are the surrounds
        for ch in (2, 3):
            audio[ch, surround_delay:] = audio[ch, :-surround_delay]
            audio[ch, :surround_delay] = 0
        log(f'surround channels delayed by {surround_delay} samples')
    L, R, Ls, Rs, C, LFE = audio
    g = 0.707; st = np.vstack([L + g * C + g * Ls + 0.5 * LFE, R + g * C + g * Rs + 0.5 * LFE])
    amp, gain = level(audio, st, level_mode)
    write_wav(f'{a.out}_5.1.wav', amp(np.vstack([L, R, C, LFE, Ls, Rs])))
    write_wav(f'{a.out}_stereo.wav', amp(st)); write_wav(f'{a.out}_centre.wav', amp(C))
    log(f'{audio.shape[1] / 48000:.1f} s written to {a.out}_*.wav (gain {20 * np.log10(gain):+.1f} dB)')
    return [f'{a.out}_5.1.wav', f'{a.out}_stereo.wav', f'{a.out}_centre.wav']

if __name__ == '__main__':
    main()
