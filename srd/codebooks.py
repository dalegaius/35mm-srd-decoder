"""The EC11 mantissa codebooks, in a plain format.

Revision 11 codes most mantissas with prefix codes chosen by Dolby from signal statistics. They follow no rule,
so they cannot be written down as a formula, and this repository does not ship them. They are read from a local
file, codebooks.npz, holding two arrays:

  huff   [6 band groups, 5 sizes (5, 7, 15, 31, 63 levels), 63 levels, 2] = (code length, code value)
         level i of an L-level quantizer has the value (2 i - (L - 1)) / L; unused entries have length 0
  pair3  [6 band groups, 16 windows, 2] = (code length, pair) for each value of the next 4 bits; pair 3 a + b
         stands for the two values ((a - 1) * 2 / 3, (b - 1) * 2 / 3) of a coefficient (its two half-blocks).
         It is kept as a window table because in the lowest group two codes give the same pair

Band groups: bands 0-4, 5-8, 9-11, 12-14, 15-17, 18-19.
The file can be produced by tools/extract_tables.py from a dump of hardware you own. Any other way of obtaining
the same codes (for instance learning them from film and a reference soundtrack) gives the same file.
Revision 9 needs none of this.
"""
import os
import numpy as np

SIZES = (5, 7, 15, 31, 63)
GROUP = [0] * 5 + [1] * 4 + [2] * 3 + [3] * 3 + [4] * 3 + [5] * 2
_NEED = ('The EC11 mantissa codebooks are missing (codebooks.npz).\n\n'
         'Revision 9 (EC9) prints decode without them. Revision 11 (EC11) prints need them, and this repository\n'
         'does not ship them. Generate the file once from a dump of hardware you own:\n\n'
         '    python tools/extract_tables.py your_dump.bin\n\n'
         'See NOTICE.md for why it works this way.')
_C = None
_LUT = {}


def path():
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'codebooks.npz')


def available():
    return os.path.exists(path())


def _load():
    global _C
    if _C is None:
        if not available():
            raise SystemExit(_NEED)
        z = np.load(path()); _C = (z['huff'].astype(np.int64), z['pair3'].astype(np.int64))
    return _C


def _table(entries, values):
    """entries: (length, code) per symbol -> (P, value per P-bit window, code length per P-bit window)"""
    P = int(max(l for l, _ in entries)); vals = [None] * (1 << P); lens = [0] * (1 << P)
    for (l, c), v in zip(entries, values):
        if l == 0: continue
        for k in range(1 << (P - l)):
            vals[(c << (P - l)) + k] = v; lens[(c << (P - l)) + k] = int(l)
    if any(l == 0 for l in lens):
        raise SystemExit('codebooks.npz: a codebook is not a complete prefix code')
    return P, vals, lens


def huffman(band, levels):
    """-> (P, values[2^P], lengths[2^P]) for the Huffman quantizer with this many levels in this band"""
    key = (GROUP[band], levels)
    if key not in _LUT:
        e = _load()[0][GROUP[band], SIZES.index(levels), :levels]
        _LUT[key] = _table([tuple(x) for x in e], [(2 * i - (levels - 1)) / levels for i in range(levels)])
    return _LUT[key]


def pairs(band):
    """-> (4, (value 0, value 1)[16], lengths[16]) for the three-level pair quantizer in this band"""
    key = (GROUP[band], 3)
    if key not in _LUT:
        e = _load()[1][GROUP[band]]
        _LUT[key] = (4, [((int(k) // 3 - 1) * 2 / 3, (int(k) % 3 - 1) * 2 / 3) for _, k in e], [int(l) for l, _ in e])
    return _LUT[key]
