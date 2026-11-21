#!/usr/bin/env python3
"""Write codebooks.npz, the EC11 mantissa codebooks, from a dump of hardware you own.

The dump is the flash of a Dolby Cat. No. 673A card: a 128 KB 29C010 marked 84147, firmware 3.00.
An older tables.npz made by an earlier version of this tool is accepted as well, and converted.
Nothing leaves your machine, and codebooks.npz is in .gitignore.

Only the prefix codes are read: for each of 6 band groups, five Huffman codebooks (5, 7, 15, 31 and 63 levels)
and one codebook of three-level pairs. Everything else the decoder needs is rules written in this repository.
This file is the only place that knows where those codes sit in that firmware; the decoder itself reads the
plain format described in srd/codebooks.py.

Usage:  python tools/extract_tables.py your_dump.bin [codebooks.npz]
        python tools/extract_tables.py old_tables.npz [codebooks.npz]
"""
import hashlib
import sys

import numpy as np

SIZE = 131072
KNOWN = {'a86e58095957a4756521a7c71cbd5ab84a8f4e0f': 'Cat. No. 673A, 84147, firmware 3.00'}
HQ_FROM, HQ_TO = 0x2600, 0x2B80         # DSP memory words holding the quantizer data
BOOKS_AT = 0x2AFD                       # per band: pointer to that band's list of codebooks
SIZES = (5, 7, 15, 31, 63)
GROUP_FIRST_BAND = (0, 5, 9, 12, 15, 18)


def dsp_word(d, a):
    """one 32-bit word of DSP memory, from the records the flash stores it in"""
    off = (0x8CCA + (a - 0x2B11) * 4) if a >= 0x27DF else (0x15ABA + (a - 0x280D) * 4)
    return int.from_bytes(d[off:off + 4], 'big')


def s20(v):
    v &= 0xFFFFF; return (v - (1 << 20) if v & 0x80000 else v) / (1 << 19)


def s12(v):
    v &= 0xFFF; return (v - 4096 if v & 0x800 else v) / 2048


def main(src, dst='codebooks.npz'):
    if src.lower().endswith('.npz'):
        z = np.load(src)
        if 'hq' not in z.files: sys.exit('this .npz is not a tables.npz from an earlier version of this tool')
        hq, base = z['hq'], int(z['hq_base']); print(f'converting {src}')
    else:
        d = open(src, 'rb').read(); sha = hashlib.sha1(d).hexdigest()
        print(f'dump: {len(d)} bytes, sha1 {sha}')
        if len(d) != SIZE:
            sys.exit(f'this is not a 128 KB flash dump ({len(d)} bytes). Expected the 84147 of a 673A card.')
        print('match: ' + KNOWN[sha] if sha in KNOWN else
              'note: this dump is not one we have seen. The result is checked below; if the checks fail,\n'
              '      the addresses probably differ in your firmware version.')
        hq = np.array([dsp_word(d, a) for a in range(HQ_FROM, HQ_TO)], np.uint32); base = HQ_FROM

    def H(a):
        return int(hq[a - base])

    huff = np.zeros((6, 5, 63, 2), np.int32); pair3 = np.zeros((6, 16, 2), np.int32)
    try:
        for g, band in enumerate(GROUP_FIRST_BAND):
            lst = H(BOOKS_AT + band) & 0xFFFF
            for si, L in enumerate(SIZES):
                p = H(lst + si + 1) & 0xFFFF; P = H(p) & 0xFFFFF
                assert 1 <= P <= 12, 'codebook depth out of range'
                ent = [H(p + 2 + k) & 0xFFFFF for k in range(1 << P)]; k = 0; seen = set()
                while k < len(ent):                     # runs of equal entries are one code
                    run = 1
                    while run * 2 <= len(ent) - k and k % (run * 2) == 0 and all(e == ent[k] for e in ent[k:k + run * 2]):
                        run *= 2
                    ln = P - int(np.log2(run)); lev = (s20(ent[k]) * L + L - 1) / 2
                    assert abs(lev - round(lev)) < 1e-3 and 0 <= round(lev) < L, 'codebook value is not a quantizer level'
                    lev = int(round(lev)); assert lev not in seen, 'a level has two codes'
                    seen.add(lev); huff[g, si, lev] = (ln, k >> (P - ln)); k += run
                assert len(seen) == L, f'{L}-level codebook has {len(seen)} codes'
                assert abs(sum(2.0 ** -int(l) for l in huff[g, si, :L, 0]) - 1) < 1e-9, 'not a complete prefix code'
            p = H(lst) & 0xFFFF
            for k in range(16):                         # one entry per 4-bit window: code length and the pair it gives
                w = H(p + k); ln = w & 0xFF; x = round(s12(w >> 20) * 1.5) + 1; y = round(s12(w >> 8) * 1.5) + 1
                assert 1 <= ln <= 4 and 0 <= x <= 2 and 0 <= y <= 2, 'pair codebook does not look right'
                pair3[g, k] = (ln, 3 * x + y)
            k = 0
            while k < 16:                               # the windows must tile into whole codes
                ln = int(pair3[g, k, 0]); span = 1 << (4 - ln)
                assert k % span == 0 and (pair3[g, k:k + span] == pair3[g, k]).all(), 'pair codebook is not a prefix code'
                k += span
    except (AssertionError, IndexError) as e:
        sys.exit(f'extraction failed its checks: {e}\nThe dump is probably a different firmware version, or not a 673A flash.')
    np.savez(dst, huff=huff, pair3=pair3)
    print(f'ok -> {dst}  (30 Huffman codebooks, 6 pair codebooks; keep it local, it is in .gitignore)')


if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(*sys.argv[1:3])
