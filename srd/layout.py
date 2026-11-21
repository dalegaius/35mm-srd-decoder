"""Where the Reed-Solomon symbols sit in an SR-D block, generated from rules (no table).

A block is a grid of 19 rows x 38 columns of bytes. The four corners (4 bytes wide, 2 rows high) hold fixed
marks and carry no data, which leaves 690 cells. They are protected by a product code:

  outer code   8 codewords RS(87,67)    (codewords 2 to 7 are one symbol short: their first symbol is absent)
  inner code   3 codewords RS(177,157)  (codeword 2 is one symbol short in the same way)

The 157 message symbols of the three inner codewords are the bytes of the audio frame; the outer code covers
every cell.

The cells are visited in one fixed order, "the stream":

  1. the left 16 columns of rows 8, 9, 10, then the right 16 columns of rows 8, 9, 10. Each of these six
     segments is read cyclically, starting from the cell whose row-by-row rank in the block is a multiple of 8;
  2. rows 0 to 7, row by row;
  3. rows 12 to 18, row by row;
  4. row 11;
  5. the 6 centre columns of rows 8, 9, 10 (around the logo).

The stream is cut into triplets of cells. A triplet belongs to one inner codeword, the three codewords taking
turns. Its last two cells hold the next two message bytes of that codeword. Its first cell holds a parity byte:
  - an inner parity byte of the same codeword when the triplet number modulo 11 is 0, 3 or 7;
  - an outer parity byte otherwise;
  - and once all the parity bytes are placed (60 inner, 160 outer), one more message byte of the codeword.

Every cell that is not an outer parity byte is a message symbol of an outer codeword. Cell number s of the
stream goes to outer codeword s modulo 8, message symbols and parity symbols each taking the next free rank of
that codeword in stream order. Near the end the codewords are simply completed in ascending order: this applies
to the parity cells from stream cell 616 on (inner and outer), and to the message cells from stream cell 675 on.

`build()` returns the same two arrays the decoder used to read from a table: for each codeword, the grid
address (1 + 38 * row + column, 0 = absent) of each of its symbols. `selfcheck()` verifies the result on real
blocks: a layout that is right gives zero syndromes on every codeword of a cleanly read block.
"""
import numpy as np

ROWS, COLS = 19, 38
OUTER_N, OUTER_K, INNER_N, INNER_K = 87, 67, 177, 157
INNER_PARITY_SLOTS = (0, 3, 7)          # triplet number modulo 11
FILL_PARITY_FROM, FILL_MESSAGE_FROM = 616, 675
OUTER_FIRST = (0, 0, 1, 1, 1, 1, 1, 1)  # first rank used by each outer codeword
INNER_FIRST = (0, 0, 1)

_cache = None


def stream():
    """the 690 data cells (row, column) in stream order"""
    corner = lambda r, c: r in (0, 1, ROWS - 2, ROWS - 1) and (c < 4 or c >= COLS - 4)
    cells = [(r, c) for r in range(ROWS) for c in range(COLS) if not corner(r, c)]
    rank = {rc: n for n, rc in enumerate(cells)}
    out = []
    for c0 in (0, COLS - 16):
        for r in (8, 9, 10):
            seg = [(r, c) for c in range(c0, c0 + 16)]
            k = next(i for i, rc in enumerate(seg) if rank[rc] % 8 == 0)
            out += seg[k:] + seg[:k]
    out += [rc for rc in cells if rc[0] <= 7]
    out += [rc for rc in cells if 12 <= rc[0]]
    out += [(11, c) for c in range(COLS)]
    out += [(r, c) for r in (8, 9, 10) for c in range(16, COLS - 16)]
    assert len(out) == 690 and len(set(out)) == 690
    return out


def build():
    """-> (outer [8, 87], inner [3, 177]) grid addresses, 1-based, 0 = absent"""
    global _cache
    if _cache is not None:
        return _cache
    outer = np.zeros((8, OUTER_N), np.int32); inner = np.zeros((3, INNER_N), np.int32)
    addr = lambda rc: rc[0] * COLS + rc[1] + 1
    n_op = [0] * 8; n_om = [0] * 8          # outer parity / message symbols placed, per codeword
    n_ip = [0] * 3; n_im = list(INNER_FIRST)   # inner parity placed / next inner message rank
    turn = 0; cur = 0; left = 0             # inner codeword rotation, current codeword, cells left in its group
    for s, rc in enumerate(stream()):
        t, first = s // 3, s % 3 == 0
        if first and t % 11 in INNER_PARITY_SLOTS and sum(n_ip) < 60:
            j = turn if s < FILL_PARITY_FROM else next(x for x in range(3) if n_ip[x] < 20)
            inner[j][INNER_K + n_ip[j]] = addr(rc); n_ip[j] += 1
        elif first and sum(n_op) < 160:
            k = s % 8 if s < FILL_PARITY_FROM else next(x for x in range(8) if n_op[x] < 20)
            outer[k][OUTER_K + n_op[k]] = addr(rc); n_op[k] += 1
            continue                         # an outer parity cell is not an outer message symbol
        else:
            if first or left == 0:           # a new group of message bytes: 2 cells, or 3 when it starts a triplet
                while n_im[turn] >= INNER_K: turn = (turn + 1) % 3
                cur = turn; turn = (turn + 1) % 3
                left = min(3 if first else 2, INNER_K - n_im[cur])
            if left == 0:                    # the codeword filled up inside the triplet: move on to the next one
                while n_im[turn] >= INNER_K: turn = (turn + 1) % 3
                cur = turn; turn = (turn + 1) % 3; left = INNER_K - n_im[cur]
            inner[cur][n_im[cur]] = addr(rc); n_im[cur] += 1; left -= 1
        k = s % 8 if s < FILL_MESSAGE_FROM else next(x for x in range(8) if OUTER_FIRST[x] + n_om[x] < OUTER_K)
        outer[k][OUTER_FIRST[k] + n_om[k]] = addr(rc); n_om[k] += 1
    assert n_op == [20] * 8 and n_ip == [20] * 3 and n_im == [INNER_K] * 3
    assert [OUTER_FIRST[k] + n_om[k] for k in range(8)] == [OUTER_K] * 8
    _cache = (outer, inner)
    return _cache


def selfcheck(grids):
    """grids: iterable of 19x38 uint8 blocks as read off film, before any correction.
    -> (outer codewords with zero syndromes, outer codewords examined, inner clean, inner examined).
    With the right layout a good share of the codewords of decently read blocks is already clean; with a wrong
    one, none is (a random word passes with probability 256^-20)."""
    import reedsolo as rs
    rs.init_tables(prim=0x187, generator=2, c_exp=8)
    outer, inner = build(); res = [0, 0, 0, 0]
    for g in grids:
        g = np.asarray(g)
        for which, cws in ((0, outer), (2, inner)):
            for cw in cws:
                word = bytearray(int(g[divmod(int(a) - 1, COLS)]) if a else 0 for a in cw)
                res[which] += max(rs.rs_calc_syndromes(word, 20, fcr=120, generator=2)) == 0; res[which + 1] += 1
    return tuple(res)
