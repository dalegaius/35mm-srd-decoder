#!/usr/bin/env python3
"""Read the tables the decoder needs out of a dump of hardware you own.

The dump is the flash on a Dolby Cat. No. 673A card: a 128 KB 29C010 marked 84147, firmware 3.00.
This script never leaves your machine, and the file it writes (tables.npz) is in .gitignore.

What it reads:
  - the Reed-Solomon interleave tables, outer RS(87,67) x 8 and inner RS(177,157) x 3;
  - the ZR38000 quantizer data: dispatch table, step sizes and the 36 Huffman codebooks.

Usage:  python tools/extract_tables.py your_84147_dump.bin [tables.npz]
"""
import hashlib
import sys

import numpy as np

SIZE = 131072
KNOWN = {'a86e58095957a4756521a7c71cbd5ab84a8f4e0f': 'Cat. No. 673A, 84147, firmware 3.00'}
OUTER_AT, OUTER_N = 0x1343, 696         # 8 codewords x 87 grid addresses, 1-based
INNER_AT, INNER_N = 0x1B6B, 531         # 3 codewords x 177
HQ_FROM, HQ_TO = 0x2600, 0x2B80         # ZR38000 words holding the quantizer data


def be24(d, off, n):
    return np.array([int.from_bytes(d[off + 3 * i:off + 3 * i + 3], 'big') for i in range(n)], np.int32)


def zoran_word(d, a):
    """One 32-bit ZR38000 word, from the two Motorola S-records the firmware stores it in."""
    off = (0x8CCA + (a - 0x2B11) * 4) if a >= 0x27DF else (0x15ABA + (a - 0x280D) * 4)
    return int.from_bytes(d[off:off + 4], 'big')


def main(src, dst='tables.npz'):
    d = open(src, 'rb').read()
    sha = hashlib.sha1(d).hexdigest()
    print(f'dump: {len(d)} bytes, sha1 {sha}')
    if len(d) != SIZE:
        sys.exit(f'this is not a 128 KB flash dump ({len(d)} bytes). Expected the 84147 on a 673A card.')
    print('match: ' + KNOWN[sha] if sha in KNOWN else
          'note: this dump is not one we have seen. Extraction will be checked below; if the checks\n'
          '      fail, the addresses probably differ in your firmware version.')

    outer = be24(d, OUTER_AT, OUTER_N)
    inner = be24(d, INNER_AT, INNER_N)
    hq = np.array([zoran_word(d, a) for a in range(HQ_FROM, HQ_TO)], np.uint32)

    def H(a):
        return int(hq[a - HQ_FROM])

    # Checks: wrong addresses give nonsense, and it is better to find out here than in the audio.
    try:
        GRID = 19 * 38                       # zeros are padding and are allowed
        assert outer.min() >= 0 and outer.max() <= GRID, 'outer interleave out of grid range'
        assert inner.min() >= 0 and inner.max() <= GRID, 'inner interleave out of grid range'
        nz = outer[outer > 0]
        assert len(set(nz.tolist())) == len(nz), 'outer interleave has repeats'
        assert len(nz) > OUTER_N - 20, 'outer interleave is mostly empty'
        tri = [(H((H(0x2B2F + i) & 0xFFFF)) & 0xFFFFF) for i in range(37)]
        assert tri[:4] == [0, 0, 0, 0], 'quantizer dispatch: allocations 0-3 should be the zero handler'
        steps = [(H(0x2B54 + a) & 0xFFFFF) / 2.0 ** 18 for a in range(37)]
        assert abs(steps[0] - 1.0) < 1e-6 and abs(steps[9] - 0.2) < 1e-3, 'step table does not look right'
    except AssertionError as e:
        sys.exit(f'extraction failed its checks: {e}\n'
                 f'The dump is probably a different firmware version, or not a 673A flash.')

    np.savez(dst, outer=outer, inner=inner, hq=hq, hq_base=HQ_FROM)
    print(f'ok -> {dst}  (keep it local: it holds Dolby data and is in .gitignore)')


if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(*sys.argv[1:3])
