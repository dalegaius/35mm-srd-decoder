# What is solid, what is not

Work in progress. This file is the honest version of the README.

## Checked against an independent reference

**Revision 11**, on a 2005 feature reel for which the DTS disc of the same print was available
(43 seconds, 3,965 frames, sample-accurate alignment):

- every decoded block matches the DTS track: median correlation 0.998, none below 0.8;
- band energies agree within a few tenths of a dB on average from 0 to 20 kHz, on L, R, C and on both
  surrounds, for separate and coupled bands alike; 1.3 dB rms per block in the coupled bands;
- the LFE channel is found in the low end of the DTS surround tracks, where DTS stores it;
- channel order L, R, Ls, Rs, C, LFE.

One thing this comparison showed and that is not a decoding matter: on that DTS disc the surround
programme is 952 samples (19.8 ms) late relative to SR-D, while the bass stored in the same tracks is
on time. `--surround-delay 952` reproduces it; it is off by default.

## Checked on film, without a reference

- **Blocks and error correction**, both revisions: the symbol layout is generated from rules that
  reproduce the reference layout exactly, and that film alone confirms. Blocks are accepted only when
  all 11 codewords agree.
- **Revision 9 bitstream**: on clean frames the 9-bit value of band 0 of each channel is normalized,
  which a wrong allocation breaks at once. A Dolby
  leader with a different tone in each channel decodes to six pure tones at the expected frequencies.
- **Revision 9 on programme material**: a 1999 trailer (2 min 29 s, 98.8% of frames from JPEG
  photographs), a 2007 cinema advertisement, rating cards. The trailer sounds right and its spectrum is
  continuous across the 3 kHz coupling boundary.

## Not checked yet

- **Revision 9 levels above 3 kHz.** The coupling rule is the same as revision 11's, which is checked,
  but no revision 9 print has been compared to an independent soundtrack. A possible 1 to 2 dB step
  just above 3 kHz was seen on one trailer; it may be the mix.
- **The order of the 76 phase flags** of revision 9 is assumed.
- **Noise reconstruction.** In revision 11 it does not change band energies (within 0.2 dB); whether
  it improves the sound is a matter of taste. It is off for revision 9, which sounded better without.
- **The mode word** (0x0000, 0x4000, 0x8000): meaning unknown, no effect on decoding found.
- **Block switching** (short transforms), described in the literature: not found in either revision
  on the prints we read.
- **The window and the launchers** have not been retested on macOS with the current tools.

## Known limits

- **Missing blocks** are concealed by repeating the previous block. Measured on blocks removed on
  purpose, this keeps energy and timbre better than the more elaborate methods we tried (which remain
  in `srd/plc.py` as experiments). With about 1% of blocks missing, small stutters remain audible.
- **Compressed video scans.** Under about 98.5% of dots read correctly almost no block decodes.
  A learned demodulator (a linear filter trained on the blocks that already decode, with decision
  feedback from the neighbouring dots) lifted one such scan from 28% to 89% of blocks on sample
  frames. It is not part of this repository yet.

## Towards needing no dump

| Needed by the decoder | Where it comes from |
|---|---|
| Transform, window, de-whitening, Reed-Solomon parameters | formulas |
| Symbol layout of a block | **rules** (`srd/layout.py`), verifiable on film |
| Frame byte order | derived from the inner code |
| Band widths, base allocations, exponent code, quantizer selection and step sizes | rules and short constants that describe the format |
| Everything revision 9 needs | **all of the above: no dump** |
| The 36 mantissa codebooks of revision 11 | **a local `codebooks.npz`. The last obstacle.** |

The codebooks can be recovered without any dump, from film plus a reference soundtrack of the same
reel: the reference says roughly which level each code must give, and the prefix tree is rebuilt from
thousands of frames, one codebook at a time. On 43 seconds of film this gave 13 of the 36 codebooks,
all 13 with exactly the right codes. The rest need more film: the rarer quantizers do not occur often
enough in 43 seconds. This work is not in the repository yet.
