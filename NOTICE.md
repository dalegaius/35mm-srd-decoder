# What this repository does and does not contain

## Does not contain

- **No Dolby firmware**, in whole or in part: no flash dump, no code, no disassembly listing.
- **No codebooks.** Revision 11 codes its mantissas with 30 Huffman codebooks and 6 codebooks of
  three-level pairs. They are not in this repository and never will be. `codebooks.npz` is in
  `.gitignore`.
- **No film soundtrack material.** Test audio is a copyrighted work and stays out.

## Does contain

Source code written by the authors, and the constants that describe the format: the block geometry, the
Reed-Solomon parameters and the rules that place the symbols in a block, the de-whitening sequence, the
band widths and base allocations of both revisions, the prefix code of the EC11 exponents, the order of
the fields in a frame.

It is fair to say where that knowledge comes from, because not all of it comes from the same place:

- **US patent 5,583,962** (Dolby Laboratories, filed 1991, expired). The transform, the band structure
  and the allocation of revision 9 are described there; the two small tables in `srd/ec9.py` are the ones
  printed in it.
- **Public documentation**, chiefly David Ferguson's `dolby-srd` wiki, for the block and its error
  correction.
- **The film itself.** Every rule here can be checked on prints: a wrong symbol layout leaves no codeword
  with zero syndromes (`srd/layout.py`, `selfcheck`), a wrong field order makes frames end at the wrong
  bit, a wrong level rule shows against a reference soundtrack.
- **The study of a decoder we own**, for interoperability. Several rules were first understood by
  studying how that decoder behaves, then written here in our own words and verified on film. No code
  and no data were copied from it. The one thing that cannot be re-expressed as a rule, the codebooks,
  is exactly what this repository leaves out.

We are not lawyers. If redistribution matters to you, take advice that fits your country.

## What you provide, for revision 11 only

A dump of the flash on hardware **you own**: the 84147 on a Dolby Cat. No. 673A card, firmware 3.00
(128 KB). `tools/extract_tables.py` reads the codebooks out of it into a local `codebooks.npz`. That
tool is the only file that knows where they sit in that firmware.

Nothing in this repository helps you obtain such a dump, and we do not distribute one.

**Revision 9 prints need nothing**: revision 9 uses plain binary mantissas.

## Why it is built this way

Understanding a format and writing your own decoder for it is one thing; redistributing a
manufacturer's data is another. The second is what this split avoids. The codebooks are several hundred
codes Dolby chose from statistics, with no construction rule: the part that most resembles original
work rather than a technical fact. They stay on your machine.

The aim is to need no dump at all. See `docs/STATUS.md` for how far that has gone.

## Credits

**Mike Saur** and **Alexis Ghnassia**, equally, for the work this decoder is built on:

- reading the blocks off the film and correcting them (outer and inner Reed-Solomon);
- the frame layout, the symbol layout rules, the de-whitening key;
- the EC11 bitstream: header, coupling fields, exponent code, allocation rule, quantizers;
- the EC9 bitstream: exponents, bit allocation, mantissa order, coupling;
- the synthesis: evenly stacked TDAC and the KBD alpha 4 window;
- the channel order, and the coupling rule checked against a DTS soundtrack.

With thanks to:

- **David Ferguson** (`davidferguson/dolby-srd`), for documenting the format publicly and for the
  srd-recognize block reader; sample images by **Denis-Carl Robidoux**.
- **furrtek**, for the 84147 flash dump used during development (not redistributed here).

## Trademarks

Dolby, Dolby Digital and SR-D are trademarks of Dolby Laboratories. DTS is a trademark of its owner.
This project is independent and has no connection with them. The names are used only to say what the
software reads.
