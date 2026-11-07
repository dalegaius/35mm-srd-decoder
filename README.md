# 36mm-srd-decoder

An open decoder for **Dolby SR-D**, the 5.1 digital soundtrack printed between the perforations of
35 mm release prints. It goes from photographs of the film to 5.1 audio, with no Dolby hardware in the
playback path.

## Work in progress

This is research code, not a finished product. It decodes, and the result is good, but:

- it has only ever been tested on **EC11**, on prints we had to hand. Other revisions of the codec
  exist and are not handled;
- the **surround channels are not right yet** (see `docs/STATUS.md`);
- a few values are **measured rather than derived from the format**, in particular the per-band gain of
  the coupled bands;
- it still needs tables from a firmware dump, and **removing that dependency is active work**. Part of
  it is done already: the frame byte order, which used to require a third-party table, is now derived
  from the error-correcting code itself.

Expect rough edges, and please report what breaks.

## Where it stands

The bitstream itself is read exactly: on a test reel, 99.3% of audio blocks are consumed to within
16 bits of their last data bit, median +1 bit. Per-band normalisation checks out (the largest value in
each coded band is at or above 0.5, median 1.000).

Front channels and LFE sound right and measure right. Surrounds do not yet.

## You need a dump of your own hardware

This repository contains **no Dolby firmware and no tables extracted from it**. The decoder needs a few
(Reed-Solomon interleaves, the quantizer dispatch, step sizes, 36 Huffman codebooks), and you generate
them once from a dump of hardware you own:

```
python tools/extract_tables.py your_84147_dump.bin
```

The dump is the 128 KB flash marked **84147** on a Dolby **Cat. No. 673A** card, firmware 3.00. The
script checks what it reads and tells you if the addresses do not fit your version. The file it writes,
`tables.npz`, stays on your machine and is in `.gitignore`.

Read `NOTICE.md` for why the project is built this way. The long-term aim is to need no dump at all:
most of what the decoder uses is already formulas or derivable from the film, and the 36 Huffman
codebooks are the real remaining obstacle.

## Install

```
pip install -r requirements.txt
```

## Use

```
python srd_decode.py path/to/photos -o myreel
```

Or double-click `Decode SR-D (Windows).bat` or `Decode SR-D (Mac).command` for a window that asks for a
folder.

Writes `myreel_5.1.wav` (L R C LFE Ls Rs), `myreel_stereo.wav` and `myreel_centre.wav`, 48 kHz.
Inputs can be a folder of photographs (JPEG, PNG, TIFF or camera RAW), 722-byte grid files such as the
`output.bin` of davidferguson/dolby-srd's srd-recognize, or `.npy` arrays of shape (n, 19, 38).

Options:

| Option | Effect |
|---|---|
| `--no-ser` | disable the noise reconstruction of quantization error |
| `--max-gap N` | conceal up to N missing blocks by repeating the previous one |
| `--surround-delay N` | delay the surrounds by N samples; try 1100. Off by default, see below |
| `--restore` | soften codec artifacts; does not recover lost information |
| `--workers N` | read photographs in parallel |

### Photographing the blocks

The block finder locates and orients blocks on its own, including mirrored film. Aim for at least
5 pixels per dot, keep a third of overlap between consecutive frames, and focus on the blocks
themselves. Blocks seen in several photographs are merged.

## How it works

| Stage | Detail |
|---|---|
| Grid | 76 x 76 fixels, 2 x 4 fixels per byte, 722 bytes per block |
| Outer code | RS(87,67) x 8 over GF(256), poly 0x187, fcr 120 |
| Inner code | RS(177,157) x 3. **Frame word w is the message symbols of inner codeword w mod 3**, taken in pairs: bytes (sym[j][2p+o], sym[j][2p+1+o]) with p = w // 3 and o = (0, 0, 1). This is why no third-party offset table is needed. |
| De-whitening | word w XOR (0xD944 + w * 0x2777) mod 65536, from frame word 21 |
| Header | word 0 = mode, word 1 = 0x0B40: block type SRD, EC revision 11, coupling mode 2 (first coupled band 14) |
| Exponents | 7 sets, 4-bit start plus VLC deltas; 5 channels x 14 bands, 6 shared coupled bands, 6 LFE bands |
| Allocation | a = Xmax - X + B[band] (+4 coupled) + (f6 - 32); floor 8 if positive; clip 0-36; cap 4(20 - X) |
| Quantizer | by allocation: none, 3-level pairs, Huffman (5 to 63 levels, per band group), or plain 7 to 9 bits |
| Coupling | sqrt((m + 16) / 32 * 2^-e) per channel and coupled band, times a per-band gain of 2.4, 3.6, 6.2, 8.9, 9.9, 11.1 dB (measured, not derived) |
| Synthesis | **evenly stacked TDAC**: half-block 0 cos(2 pi / N (n + 128) k), half-block 1 sin(...), N = 512, hop 256, **KBD window alpha 4** |
| Channels | stream order L, R, Ls, Rs, C, LFE |

Two points cost the most to find and are worth stating plainly. The transform is **evenly stacked**,
not the odd stacking of AC-3: using the wrong one shifts every frequency by 47 Hz and makes speech
sound metallic. And the window is **KBD alpha 4**, which matches the window table Dolby broadcasts on
the film itself in module E3 to within half a least-significant bit.

The format is described in **US patent 5,583,962** (Dolby Laboratories, 1992), now expired.

## Surrounds

The surround channels have about the right level and roughly the right spectral shape, but their
waveform does not match what it should be. Shifting them by about 1100 samples improves things without
explaining anything, so `--surround-delay 1100` is provided and left off by default. Ruled out so far:
the steering flags (inactive on every block we examined) and the coupling master exponents (constant).
If you know this format, we would like to hear from you.

## Credits

**Mike Saur** and **Alexis Ghnassia**, equally. With thanks to **David Ferguson**, whose `dolby-srd`
documented the format publicly and provides a block reader, and to **furrtek** for the flash dump used
during development. See `NOTICE.md`.

## Licence

**GNU GPL v3 or later**, for the source code in this repository only. Any version you build on must stay
free software. See `LICENSE` and `NOTICE.md`.
