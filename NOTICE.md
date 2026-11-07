# What this repository does and does not contain

## Does not contain

- **No Dolby firmware**, in whole or in part: no flash dumps, no code, no disassembly listings.
- **No tables extracted from that firmware.** The decoder needs a few of them (Reed-Solomon interleave
  tables, the quantizer dispatch table, step sizes and 36 Huffman codebooks). They are **not** in this
  repository and never will be. `tables.npz` is in `.gitignore`.
- **No film soundtrack material.** Test audio is a copyrighted work and stays out.

## What you provide

A dump of the flash on hardware **you own**: the 84147 on a Dolby Cat. No. 673A card, firmware 3.00
(128 KB). `tools/extract_tables.py` reads the tables out of it into a local `tables.npz`.

Nothing in this repository helps you obtain such a dump, and we do not distribute one.

## Why it is built this way

Understanding a format and writing your own decoder for it is one thing; redistributing a
manufacturer's data is another. The second is what this split avoids. The tables that carry real
creative weight are the 36 Huffman codebooks: several thousand values Dolby chose from statistics,
with no construction rule, which is the part that most resembles original work rather than a technical
fact. They stay on your machine.

The format itself is described in **US patent 5,583,962** (Dolby Laboratories, filed 1992), which has
expired.

## Credits

**Mike Saur** and **Alexis Ghnassia**, equally, for the work this decoder is built on:

- reading the blocks off the film and correcting them (outer and inner Reed-Solomon);
- the frame layout, and deriving the byte order from the inner code;
- the de-whitening key;
- the EC11 bitstream: header, coupling fields, exponent VLC, allocation rule;
- the quantizer: dispatch, step sizes, codebooks, three-level pairs;
- the synthesis: evenly stacked TDAC and the KBD alpha 4 window;
- the channel order, and the calibration of the coupled bands.

With thanks to:

- **David Ferguson** (`davidferguson/dolby-srd`), for documenting the format publicly and for the
  srd-recognize block reader; sample images by **Denis-Carl Robidoux**.
- **furrtek**, for the 84147 flash dump used during development (not redistributed here).

## Trademarks

Dolby, Dolby Digital and SR-D are trademarks of Dolby Laboratories. This project is independent and
has no connection with Dolby Laboratories. The names are used only to say what the software reads.
