# Changes

## November 2026

**Revision 9 (EC9) is decoded**, with no data from any hardware dump (`srd/ec9.py`).

**Less dependence on a dump for revision 11.**
- The positions of the Reed-Solomon symbols are generated from rules (`srd/layout.py`) instead of
  being read from a table.
- Quantizer selection and step sizes are rules.
- The only local data left is the mantissa codebooks, in a plain format: `codebooks.npz`
  (`srd/codebooks.py`). `tools/extract_tables.py` writes it from a dump, or converts a `tables.npz`
  made by an earlier version. **`tables.npz` is no longer read: run the tool once.**

**Coupling.** The measured per-band gains are gone. Each channel has an absolute level per coupled
band, applied to the unit-RMS composite; in revision 11 its exponent is raised by three times the
channel's 2-bit field, which the decoder used to ignore. Against a DTS soundtrack: 1.3 dB rms per
block, where the former rule gave 8 to 9 dB.

**Error correction.**
- A corrected codeword is accepted only with a safety margin. Before, a 20-check codeword could be
  "corrected" with 20 erasures, which accepts anything.
- A block is kept only if all 11 codewords are consistent on the final grid.
- Confidence maps from the readers are used for erasures.

**Photographs.** `srd_photos.py` works from the output of `srd_prepare.py`: several reads of each
block, block numbers predicted from decoded neighbours, and blocks that fail everywhere retried by
adding their reads together. On one reel: from 1,033 to 3,966 frames out of 4,014.

**Video scans.** `srd_extract_video.py` reads blocks from film scans delivered as video, stitching
blocks cut between two frames.

**Output.** Average level set to -20 dBFS by default (`--level peak` for the former behaviour).
Missing blocks are concealed by repeating the previous one.

**Mode word 0x0000** (seen on a 1999 trailer) is accepted alongside 0x4000 and 0x8000.

**Window.** `srd_gui.py` picks the right tools for the folder's contents and shows progress and an
estimate of the time left.
