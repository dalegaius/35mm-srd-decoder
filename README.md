# 35mm-srd-decoder

An open decoder for **Dolby SR-D**, the 5.1 digital soundtrack printed between the perforations of
35 mm release prints. It goes from photographs or scans of the film to 5.1 audio, with no Dolby hardware
in the playback path.

Both codec revisions found on film are decoded:

| | Revision 9 ("EC9") | Revision 11 ("EC11") |
|---|---|---|
| Needs data from a hardware dump | **no** | yes, the mantissa codebooks only |
| Checked against an independent reference | not yet | yes, the DTS disc of the same reel |
| Tested on | a 1999 trailer, a 2007 cinema ad, rating cards, a Dolby leader | a 2005 feature reel, a 2006 trailer |

## Work in progress

This is research code, not a finished product. It decodes, and the result is good, but read
`docs/STATUS.md` before relying on it: it says what is established, what is measured, and what is still
assumed. Expect rough edges, and please report what breaks.

## Install

```
pip install -r requirements.txt
```

Python 3.9 or later.

## Use

### The window

Double-click `Decode SR-D (Windows).bat` or `Decode SR-D (Mac).command`, choose a folder. What happens
depends on what the folder holds:

| Folder contents | What runs |
|---|---|
| Camera RAW photos (`.dng`) | `srd_prepare.py`, then `srd_photos.py`: the path that recovers the most blocks |
| `points.npz` and `rapport.csv` from an earlier run | `srd_photos.py` only |
| Other images (JPEG, PNG, TIFF), or `.npz` block files from the video extractor | `srd_decode.py` |

The result goes to a `srd_decoded` folder inside the one you chose: `name_5.1.wav` (L R C LFE Ls Rs),
`name_stereo.wav`, `name_centre.wav`, 48 kHz 16 bit. The window shows the current step, the elapsed time
and an estimate of the time left.

### The command line

```
python srd_decode.py path/to/images -o myreel            # images, grid files, .npz block files
python srd_prepare.py photos_dng -o work --no-crops      # RAW photos: read every block once
python srd_photos.py work -o work/myreel                 # correct, combine the reads, decode
python srd_extract_video.py scan.mkv -o blocks.npz       # film scans delivered as video
```

`srd_decode.py` also reads 722-byte grid files such as the `output.bin` of davidferguson/dolby-srd's
srd-recognize, and `.npy` arrays of shape (n, 19, 38).

| Option of `srd_decode.py` | Effect |
|---|---|
| `--workers N` | read images in parallel |
| `--level loud` / `peak` | output level: average at -20 dBFS (default), or highest sample at 0.89 |
| `--max-gap N` | conceal up to N missing blocks in a row by repeating the previous one (default 8) |
| `--no-ser` | EC11: disable the noise reconstruction of quantization error |
| `--surround-delay N` | delay the surrounds by N samples; off by default |
| `--restore` | soften codec artifacts; does not recover lost information |

`srd_photos.py` is the better tool when a block appears in several photographs: a block that fails in
every photograph is retried by **adding its reads together** before thresholding, which often turns two
poor reads into a good one. It also writes `name_manquants.csv`, the blocks still missing and the
photographs they appear in, so that you know what to shoot again.

### Photographing the blocks

- At least 3.5 pixels per dot; 5 is comfortable.
- Each block should appear in **two** photographs. With six blocks per photograph, advance three blocks
  between shots.
- Blocks in the middle of the frame read at 94 to 99%, blocks at the edges often below 15%: what matters
  is that every block is near the middle of at least one photograph.
- Even, diffuse backlight, the camera parallel to the film, focus on the blocks.
- The block finder locates and orients blocks on its own, including mirrored film.

Film scans published as compressed video are a harder source: below about 98.5% of dots read correctly,
almost no block survives the error correction. `srd_extract_video.py` handles the good ones.

## The EC11 codebooks

Revision 11 codes its mantissas with 36 prefix codes that Dolby chose from statistics. They follow no
rule, so this repository cannot contain them as a formula, and it does not contain them as data. You
generate them once, from a dump of hardware you own:

```
python tools/extract_tables.py your_dump.bin
```

The dump is the 128 KB flash marked **84147** on a Dolby **Cat. No. 673A** card, firmware 3.00. The tool
checks what it reads and writes `codebooks.npz`, a plain file (code lengths and code values, format in
`srd/codebooks.py`) that stays on your machine and is in `.gitignore`. A `tables.npz` made by an earlier
version of this repository can be converted with the same command.

Everything else is rules written in this repository. **A revision 9 print needs no dump at all.**
Read `NOTICE.md` for why the project is built this way.

## How it works

`docs/FORMAT.md` describes the format in full. In short:

| Stage | Detail |
|---|---|
| Grid | 76 x 76 dots, 2 x 4 dots per byte, 722 bytes per block, 4 blocks per film frame |
| Error correction | outer RS(87,67) x 8 and inner RS(177,157) x 3 over GF(256), poly 0x187, first root 120. Symbol positions are generated from rules (`srd/layout.py`). A block is accepted only if all 11 codewords agree on the final grid |
| Frame | 16-bit words from the message symbols of the three inner codewords, de-whitened: word w XOR (0xD944 + w x 0x2777) mod 65536 |
| Header | word 0 = mode (0x0000, 0x4000 and 0x8000 seen), word 1 high byte = revision (9 or 11) |
| Audio | 3,440 bits per block for six channels, 512 samples at 48 kHz, stream order L, R, Ls, Rs, C, LFE |
| Synthesis | **evenly stacked TDAC**: half-block 0 cos(2 pi / N (n + 128) k), half-block 1 sin(...), N = 512, hop 256, **KBD window alpha 4** |
| Coupling | above 3 kHz (EC9) or 4.7 kHz (EC11) the channels share one composite signal, which only gives the shape of each band; every channel has its own absolute level per band |

Two points cost the most to find and are worth stating plainly. The transform is **evenly stacked**, not
the odd stacking of AC-3: using the wrong one shifts every frequency by 47 Hz and makes speech sound
metallic. And the window is **KBD alpha 4**.

The format is described in **US patent 5,583,962** (Dolby Laboratories, filed 1991), now expired. Revision
9 follows it closely.

## Credits

**Mike Saur** and **Alexis Ghnassia**, equally. With thanks to **David Ferguson**, whose `dolby-srd`
documented the format publicly and provides a block reader, and to **furrtek** for the flash dump used
during development. See `NOTICE.md`.

## Licence

**GNU GPL v3 or later**, for the source code in this repository only. Any version you build on must stay
free software. See `LICENSE` and `NOTICE.md`.
