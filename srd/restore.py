"""Optional restoration of codec artifacts (--restore). Works on the decoded coefficients of the whole
timeline before synthesis. Nothing here recovers lost information: it only makes the artifacts of the
coarse quantization less audible (spectral holes that flicker on and off, jumpy coupled highs).

1. Hole filling without "birdies": a coefficient that is zero in one frame but non-zero in the frames
   before and after gets the geometric mean of its neighbours' magnitudes (x0.7) with a random sign.
   Other zeros in coded bands get faint noise at a fraction of the band's level.
2. Coupled highs (above ~4.7 kHz): the level of each coupled band, per channel and half-block, is
   smoothed over 3 frames (median), with the correction limited to +-6 dB.
"""
import numpy as np
from .decoder import START, WIDTHS

def restore(spec, rng=None, hole_gain=0.7, floor_gain=0.15, smooth=3, hf_lift_db=0.0):
    """spec: (T, 6, 2, 256) array, modified copy returned."""
    rng = rng or np.random.default_rng(1)
    S = spec.copy(); A = np.abs(spec); T = len(S)
    # 1. holes between two non-zero neighbours
    Z = (spec == 0)
    prev = np.zeros_like(A); nxt = np.zeros_like(A)
    prev[1:] = A[:-1]; nxt[:-1] = A[1:]
    hole = Z & (prev > 0) & (nxt > 0)
    fill = hole_gain * np.sqrt(prev * nxt)
    S[hole] = fill[hole] * rng.choice([-1.0, 1.0], size=int(hole.sum()))
    # faint floor on remaining zeros inside bands that carry signal
    for b in range(2, 20):
        k0, k1 = START[b], START[b] + WIDTHS[b]
        band = A[..., k0:k1]
        lvl = np.sqrt((band ** 2).mean(axis=-1, keepdims=True))
        z = (S[..., k0:k1] == 0) & (lvl > 0)
        noise = rng.uniform(-1, 1, size=z.shape) * floor_gain * lvl
        S[..., k0:k1] = np.where(z, noise, S[..., k0:k1])
    # 2. smooth the level of coupled bands over time
    for b in range(14, 20):
        k0, k1 = START[b], START[b] + WIDTHS[b]
        e = np.sqrt((S[..., k0:k1] ** 2).mean(axis=-1)) + 1e-12          # (T, 6, 2)
        h = smooth // 2
        stack = [np.concatenate([e[:1].repeat(max(0, i), 0), e[:len(e) - max(0, i)]]) if i >= 0 else
                 np.concatenate([e[-i:], e[-1:].repeat(-i, 0)]) for i in range(-h, h + 1)]
        target = np.median(np.stack(stack), axis=0)
        g = np.clip(target / e, 0.5, 2.0)
        S[..., k0:k1] *= g[..., None]
    # 3. optional high-frequency lift (shelf from ~7.5 kHz, linear in dB up to hf_lift_db at 15 kHz)
    if hf_lift_db:
        f = (np.arange(256) + 0.5) * 93.75
        gdb = np.clip((f - 7500) / 7500, 0, 1) * hf_lift_db
        S *= 10 ** (gdb / 20)
    return S
