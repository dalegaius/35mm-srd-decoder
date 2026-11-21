"""EC11 exponent coding: a 4-bit first exponent, then one variable-length difference per following band.

The differences use a fixed prefix code, the same on every print:
  00 -> 0 | 010 -> +1 | 0110 -> -1 | 0111 -> +2 | 1000 -> -2 | 10010 -> +3 | 10011 -> -3
  then 6-bit codes 101000 to 111111 for +4, +8, +6, +5, +7, +10, +9, +15, +14, +13, +12, +11, -4 ... -15
It is stored here as a 64-entry table indexed by the next 6 bits (value, and number of bits used).
"""

T = [0] * 16 + [1] * 8 + [-1] * 4 + [2] * 4 + [-2] * 4 + [3, 3, -3, -3] + \
    [4, 8, 6, 5, 7, 10, 9, 15, 14, 13, 12, 11] + list(range(-4, -16, -1))
L = [2] * 16 + [3] * 8 + [4] * 12 + [5] * 4 + [6] * 24
assert len(T) == 64 and len(L) == 64


def peek6(bits, p):
    v = 0
    for k in range(6):
        v = v * 2 + (int(bits[p + k]) if p + k < len(bits) else 0)
    return v


def decode_channel(bits, p, n, first_bits=4):
    """first exponent on 4 bits, then n coded differences. Returns (exponents, end position)."""
    e = 0
    for k in range(first_bits):
        e = e * 2 + int(bits[p + k])
    p += first_bits
    out = [e]
    for _ in range(n):
        i = peek6(bits, p)
        e += T[i]
        p += L[i]
        out.append(e)
    return out, p
