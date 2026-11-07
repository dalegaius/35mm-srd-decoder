#!/usr/bin/env python3
"""
EC11 exponent VLC from the ZR38000 table H:$25EB (fredrun 3.1.12, flash 0xE07F record).

The 64-entry table is indexed by the next 6 bits; entries form a complete prefix code:
  00 -> 0 | 010 -> +1 | 0110 -> -1 | 0111 -> +2 | 1000 -> -2 | 10010 -> +3 | 10011 -> -3
  101000..111111 (6 bits) -> table values (4,8,6,5,7,10,9,15,14,13,12,11,-4..-15)
"""
import numpy as np

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
    """4-bit absolute exponent then n VLC deltas. Returns (exponents, end position)."""
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
