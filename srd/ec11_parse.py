"""EC11 frame parser (SR-D revision 11).

Bit 0 is the first bit of frame word 2 (word 0 is the mode word, word 1 = 0x0B40: revision 11 in the high byte).

  bits      field
  0-119     30 x 4 bits   coupling level exponents, 5 channels x 6 coupled bands             ('r1')
  120-129    5 x 2 bits   one field per channel, added three times over to its exponents    ('c2')
  130-249   30 x 4 bits   coupling level mantissas                                          ('r2')
  250-...    5 flags      each followed by an 8-bit value when set (never seen set)         ('steer')
  next       7 x 2 bits   one field per exponent set, added three times over to it          ('strategy')
  next       6 bits       global allocation offset                                          ('f6')
  next       7 exponent sets: 5 channels x 14 bands, LFE x 6 bands, composite x 6 bands (bands 14-19);
             each is a 4-bit first exponent then coded differences (srd/ec11_vlc.py)
  then       mantissas (srd/decoder.py)
"""
import numpy as np
from .ec11_vlc import decode_channel

NDELTA = [13] * 5 + [5, 5]   # differences per exponent set


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
        off = 3 * out['strategy'][k]          # the whole set is raised by 3 x its 2-bit field
        exps.append([x + off for x in e])
    out['exps'] = exps
    out['mant_start'] = p
    return out


