#!/usr/bin/env python3
"""
EC11 block parser, field order read from the ZR38000 decoder (fredrun 3.1.12, H:$31AE/$315F/$3248).

De-whitened payload bit 0 = frame word 2 (frame words 0,1 = 0x8000, 0x0B40; 0x0B40>>8 = 11 -> EC11).
Mode 11 config (H:$1A00 odd entries): flags 0xFAAD0, 20-band tables at $255E/$259B.

  bits    reader  field
  0-119   $3A05   30 x 4 bits, unsigned       -> $2CD  (coupling coord exponents, 5 ch x 6 composite bands)
  120-129 $3A05   5 x 2 bits                  -> $32C  (flag 0x2000)
  130-249 $39DE   30 x 4 bits, left-justified -> $331  (coupling coord mantissas)
  250-254 $3A1B   5 x 1-bit flags, each + 8-bit value if set (flag 0x10)
  next    $3A05   7 x 2 bits exponent strategies -> $2C5
  next    $3A05   1 x 6 bits                  -> $2CC -> [$3E0]
  next    $3264   7 exponent sets: 4-bit absolute + VLC deltas (table H:$25EB);
                  sets 1-5: 13 deltas (14 bands, discrete), set 6: 5 (6 LFE bands),
                  set 7: 5 (6 composite bands 14-19); first exponent += 3*strategy
  then            mantissas
"""
import numpy as np
from .ec11_vlc import decode_channel

NDELTA = [13] * 5 + [5, 5]   # [$44E]=13, [$3F1]=5, [$451]=[$450]-1=5


def bits_int(b, p, n):
    v = 0
    for k in range(n):
        v = v * 2 + int(b[p + k])
    return v


def parse(b):
    p = 0
    out = {}
    out['r1'] = [bits_int(b, p + 4 * i, 4) for i in range(30)]; p += 120
    out['c2'] = [bits_int(b, p + 2 * i, 2) for i in range(5)]; p += 10
    out['r2'] = [bits_int(b, p + 4 * i, 4) for i in range(30)]; p += 120
    steer = []
    for c in range(5):
        f = int(b[p]); p += 1
        v = None
        if f:
            v = bits_int(b, p, 8); p += 8
        steer.append(v)
    out['steer'] = steer
    out['strategy'] = [bits_int(b, p + 2 * i, 2) for i in range(7)]; p += 14
    out['f6'] = bits_int(b, p, 6); p += 6
    out['exp_start'] = p
    exps = []
    for k, n in enumerate(NDELTA):
        e, p = decode_channel(b, p, n)
        off = 3 * out['strategy'][k]          # $3264: first exponent = 4-bit + 3*strategy (op22 = x D2)
        exps.append([x + off for x in e])
    out['exps'] = exps
    out['mant_start'] = p
    return out


