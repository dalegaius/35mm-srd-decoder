# What is solid, what is not

Work in progress. This file is the honest version of the README.

## Solid

- **Bitstream reading**: blocks, Reed-Solomon, frame assembly, de-whitening, EC11 header, exponents,
  allocation, quantizer. On a test reel, 99.3% of audio blocks are consumed to within 16 bits of their
  last data bit, median +1.
- **Synthesis**: evenly stacked TDAC, KBD alpha 4. The window matches the table Dolby broadcasts in
  module E3 on the film itself to within half an LSB over 124 samples.
- **Front channels and LFE**: correct to the ear and to measurement.

## Measured, not derived

- **Per-band gain of the coupled bands**: 2.4, 3.6, 6.2, 8.9, 9.9, 11.1 dB for bands 14 to 19. Fitted
  against a reference mix, not read out of the format. The real table most likely lives in module E3 on
  the film, in packets we never recovered. Anyone who captures those packets would replace this.

## Not understood

- **Surround channels.** Right level, roughly the right spectral shape, wrong waveform. Ruled out:
  - *steering*: the five flags are inactive on all 790 blocks examined;
  - *coupling master exponents*: 3 on every channel in 98% of blocks, so a constant;
  - *a fixed delay*: about 1100 samples improves the match but does not reproduce cleanly, which is why
    `--surround-delay` exists and defaults to off.

## Scope

- **EC11 only**, and only tested on the prints we had. The firmware also recognises an EC9 path, but in
  version 3.00 that branch is a bare return: EC9 prints are rejected rather than decoded. Other
  revisions may exist; the per-revision configuration table sits at H:$1A00, just past the end of the
  memory region we extracted.
- The **coupling mode** is in the header (mode 2, first coupled band 14, on the prints we tested) but
  the decoder assumes that value. Prints using another mode need that made dynamic.

## Towards needing no dump

| Needed by the decoder | Where it comes from now |
|---|---|
| KBD window, TDAC, de-whitening, quantizer step sizes, RS parameters | formulas, nothing to ship |
| RS interleave tables | the dump, but derivable from the film by finding which positions form valid codewords |
| Frame byte order | **derived from the inner code**: no table, no third-party file |
| Band widths, allocation table B, exponent VLC, dispatch | a few dozen constants that describe the format |
| **36 Huffman codebooks for the mantissas** | **the dump. This is the real obstacle.** |

The codebooks are trained tables with no construction rule, so recovering them from data alone would be
a research project of its own. Everything else can go.
