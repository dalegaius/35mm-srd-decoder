# The SR-D format, as this decoder reads it

This is a description of the format, not of any implementation. Where something is assumed rather than
established, it says so. `docs/STATUS.md` lists how each part was checked.

## 1. The block

- One block sits between each pair of perforations on the soundtrack side: 4 blocks per film frame,
  96 per second at 24 frames per second.
- A block is 76 x 76 dots. A byte is 2 dots wide and 4 dots high, which gives a grid of 19 rows of
  38 bytes, 722 bytes. A clear dot is a 1.
- The four corners (4 bytes wide, 2 rows high, 32 bytes in all) hold fixed alignment marks. A logo
  sits in the centre. That leaves 690 bytes of data.
- Row 11 ends with the text `DOLBY SR*D`, preceded by the block number (18 bits) and the reel number
  (6 bits). Trailers and advertisements usually carry reel 63; a 1999 trailer carries reel 1.
- The frame is whitened before printing (section 3), so a block looks like noise even when it codes silence.

## 2. Error correction

A product code over GF(256), field polynomial 0x187, first root alpha^120, 20 check symbols per
codeword:

- **outer code**: 8 codewords RS(87,67). Codewords 2 to 7 are one symbol short (their first symbol is
  absent and counts as zero);
- **inner code**: 3 codewords RS(177,157). Codeword 2 is one symbol short in the same way.

690 = 470 frame bytes + 160 outer check bytes + 60 inner check bytes. The outer code covers every
cell except its own check bytes; the 157 message symbols of each inner codeword are frame bytes.

### Where the symbols sit

The 690 cells are visited in one fixed order, the stream:

1. the left 16 columns of rows 8, 9, 10, then the right 16 columns of rows 8, 9, 10. Each of these six
   segments is read cyclically, starting from the cell whose row-by-row rank in the block (corners
   excluded) is a multiple of 8;
2. rows 0 to 7, row by row;
3. rows 12 to 18, row by row;
4. row 11;
5. the 6 centre columns of rows 8, 9, 10.

The stream is cut into triplets of cells. A triplet belongs to one inner codeword, the three
codewords taking turns. Its last two cells hold the next two frame bytes of that codeword. Its first
cell holds a check byte: an inner check byte of the same codeword when the triplet number modulo 11
is 0, 3 or 7, an outer check byte otherwise; once all 220 check bytes are placed, one more frame byte.

Cell number s of the stream belongs to outer codeword s modulo 8. Near the end, codewords are simply
completed in ascending order: from stream cell 616 for check bytes, from cell 675 for the others.
Those two thresholds are observations: they are right, but we know no deeper reason for them.

`srd/layout.py` implements exactly this. A check that needs nothing but film: with the right layout a
good share of the codewords of a decently read block already have zero syndromes; with a rule off by
one cell, none has.

### Decoding practice

A codeword may be corrected with errors and erasures, but a result is accepted only with a safety
margin (a 20-check word that "corrects" with 20 erasures proves nothing). A block is kept only if all
11 codewords are consistent on the final grid. About three blocks in every 96 are not audio.

## 3. The frame

- 16-bit words are formed from the message symbols of the three inner codewords in turn, two symbols
  per word: word w takes symbols (2p + o, 2p + 1 + o) of codeword w mod 3, with p = w div 3 and
  o = 0, 0, 1 for the three codewords.
- De-whitening: word w is XORed with (0xD944 + w x 0x2777) mod 65536.
- Word 0 is a mode word. Seen so far: 0x0000 (a 1999 trailer), 0x4000 (prints of 2004 to 2007),
  0x8000 (a 2005 feature). Its meaning is not known; the decoding is the same for all three.
- The high byte of word 1 is the codec revision: 9 or 11. We call them EC9 and EC11, after the names
  reported for them. No other value has been seen.
- Words 2 to 216 are the audio payload: 3,440 bits for six channels, about 322 kbit/s.

## 4. What both revisions share

- 48 kHz, 512 new samples per block, six channels in the order **L, R, Ls, Rs, C, LFE**.
- The transform is the evenly stacked TDAC of the patent: each block carries two sets of 256
  coefficients, a modified DCT for the first half-block and a modified DST for the second, N = 512,
  hop 256, Kaiser-Bessel-derived window with alpha = 4.
- A coefficient is a mantissa times 2^-exponent. Exponents are sent per band; the number of bits per
  mantissa is not sent but computed by the decoder from the exponents.
- Above a crossover frequency the five main channels share one composite signal (coupling). The
  composite only gives the **shape** of each band: it is brought to unit RMS over the band, both
  half-blocks together. Each channel then has its own **absolute level** for the band, sent as a
  small exponent and mantissa.

## 5. Revision 9 (EC9)

Everything has a fixed size and position. Bit 0 is the first bit of word 2.

| Bits | Field |
|---|---|
| 0-503 | 126 exponents of 4 bits: 39 for the first group (20 bands of L, then 19 composite bands), 20 each for R, Ls, Rs, C, 7 for LFE |
| 504-517 | 7 fields of 2 bits, one per exponent set (L, R, Ls, Rs, C, LFE, composite). An exponent is its 4 bits plus 3 times its set's field |
| 518-521 | 4-bit global allocation offset f |
| 522-996 | 95 coupling exponents of 5 bits: 5 channels x 19 composite bands |
| 997-1471 | 95 coupling mantissas of 5 bits |
| 1472-1547 | 76 phase-inversion flags: 4 channels (all but L) x 19 bands. The order within this field is assumed |
| 1548-3439 | mantissas |

**Bands.** 39 bands of widths 12 x 1, 5 x 2, 2 x 3, 6 x 4, 1 x 6, 9 x 10, 1 x 15, 3 x 18 coefficients
(Table I of the patent). Bands 0 to 19 (up to about 3 kHz) are sent per channel, bands 20 to 38 are
the composite.

**Allocation.** For each band, bits = base[band] + round((Xmax - X) / 4) + f - 8, limited to 0..9,
where Xmax is the largest exponent of the frame, and base is 5 for bands 0-11, 4 for 12-23, 3 for
24-33, 2, 2, 0, 0, 0 for the last five, plus 1 for composite bands. Bands 0 and 1 always get 9 bits.
A band whose 4-bit exponent is 15 is silent and starts with 0 bits. The 1,730 bits left after bands 0
and 1 are then handed out one bit per mantissa at a time, band 2 upward across the five channels and
then the composite bands, in repeated rounds, up to 9 bits, until a band no longer fits. What remains
goes to the last coefficients of that band, half in each pass.

**Mantissas.** Plain two's complement. Two passes: all the DCT values of all channels and bands, then
all the DST values (which have no band 0).

**Coupling.** Channel level for a composite band = (m / 32) x 2^-e, with the 5-bit e and m of that
channel and band, applied to the unit-RMS composite; negated when the channel's phase flag is set.

## 6. Revision 11 (EC11)

Same payload size, spent more cleverly: variable-length exponents and entropy-coded mantissas free
enough bits to keep the channels separate up to about 4.7 kHz.

| Bits | Field |
|---|---|
| 0-119 | 30 coupling exponents of 4 bits: 5 channels x 6 coupled bands |
| 120-129 | 5 fields of 2 bits, one per channel (c2) |
| 130-249 | 30 coupling mantissas of 4 bits |
| then | 5 flags, each followed by an 8-bit value when set (never seen set) |
| then | 7 fields of 2 bits, one per exponent set |
| then | 6-bit global allocation offset f |
| then | 7 exponent sets: 5 channels x 14 bands, LFE x 6 bands, composite x 6 bands |
| then | mantissas |

**Bands.** 20 bands of widths 1, 1, 1, 1, 1, 2, 2, 3, 3, 4, 6, 7, 8, 10, 13, 16, 20, 26, 36, 56.
Bands 0 to 13 per channel, bands 14 to 19 composite.

**Exponents.** Each set is a 4-bit first exponent, raised by 3 times the set's 2-bit field, followed
by one coded difference per following band: `00` 0, `010` +1, `0110` -1, `0111` +2, `1000` -2,
`10010` +3, `10011` -3, then 6-bit codes from `101000` for +4, +8, +6, +5, +7, +10, +9, +15, +14, +13,
+12, +11, -4 to -15.

**Allocation.** For each band an index a = Xmax - X + B[band] (+4 for composite bands) + (f - 32),
with B = 20 x 6, 17, 16 x 5, 13, 12 x 4, 9, 2, 0; raised to 8 if positive, limited to 0..36, and to
at most 4 (20 - X).

**Quantizers.** The index selects one:

| Index | Quantizer |
|---|---|
| 0-3 | nothing sent |
| 4-8 | three levels, the two values of a coefficient (DCT, DST) in one code of 1 to 4 bits |
| 9-11, 12-13, 14-18, 19-22, 23-26 | 5, 7, 15, 31, 63 levels, one Huffman code per value |
| 27-30, 31-34, 35-36 | plain two's complement on 7, 8, 9 bits |

Level i of an L-level quantizer has the value (2i - (L - 1)) / L. The Huffman codes depend on the band
group (bands 0-4, 5-8, 9-11, 12-14, 15-17, 18-19): 30 codebooks, plus 6 for the three-level pairs.
**These codebooks are the one part of the format this document cannot give**: they are statistical
codes with no rule. See `NOTICE.md`.

Bands 0 and 1 are always plain 9 bits. Mantissas are read in one pass, band by band, DCT and DST
values alternating; channel order L (with the composite after its band 13), R, Ls, Rs, C, LFE.

**Coupling.** Channel level for a coupled band = ((m + 16) / 32) x 2^-E with E = e + 3 x c2 of the
channel. When E reaches 24 the mantissa has no implicit leading one (level m / 32 x 2^-24), so that
a level of exactly zero can be sent.

## 7. Differences in one table

| | EC9 | EC11 |
|---|---|---|
| Bands per channel | 20, up to 3 kHz | 14, up to 4.7 kHz |
| Coupled bands | 19 | 6 |
| LFE bands | 7 | 6 |
| Exponents | 4 bits fixed | variable length |
| Mantissas | plain, 0 to 9 bits | none, 3-level pairs, Huffman, or plain |
| Mantissa order | all DCT, then all DST | one pass, alternating |
| Coupling level | 5 + 5 bits | 4 + 4 bits, exponent offset per channel |
| Layout | fixed | varies from frame to frame |
