"""Stack single-page PDF figures top to bottom, preserving vector content.

Run with the base Python interpreter, supplying source PDFs in display order
and --output for the combined single-page PDF. Requires pypdf.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
from pathlib import Path

from pypdf import PdfReader, PdfWriter, Transformation


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sources", type=Path, nargs="+")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.resolve() in {path.resolve() for path in args.sources}:
        parser.error("The combined output must not overwrite a source PDF")
    source_hashes = [sha256(path.read_bytes()).hexdigest() for path in args.sources]
    readers = [PdfReader(path) for path in args.sources]
    if any(len(reader.pages) != 1 for reader in readers):
        parser.error("Every source must contain exactly one page")
    pages = [reader.pages[0] for reader in readers]
    for page in pages:
        if page.rotation:
            page.transfer_rotation_to_content()
    width = max(float(page.cropbox.width) for page in pages)
    scales = [width / float(page.cropbox.width) for page in pages]
    heights = [float(page.cropbox.height) * scale for page, scale in zip(pages, scales)]
    margin, gap = 12.0, 18.0
    height = sum(heights) + gap * (len(pages) - 1) + 2 * margin
    writer = PdfWriter()
    canvas = writer.add_blank_page(width=width + 2 * margin, height=height)
    cursor = height - margin
    for page, scale, panel_height in zip(pages, scales, heights):
        cursor -= panel_height
        transform = (Transformation()
                     .translate(-float(page.cropbox.left), -float(page.cropbox.bottom))
                     .scale(scale).translate(margin, cursor))
        canvas.merge_transformed_page(page, transform)
        cursor -= gap
    writer.add_metadata({"/Title": "Bottom-boundary topology and closed-field helicity"})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("wb") as stream:
        writer.write(stream)
    assert len(PdfReader(args.output).pages) == 1
    assert source_hashes == [sha256(path.read_bytes()).hexdigest() for path in args.sources]
    print(f"Saved one-page PDF with unchanged sources: {args.output}")


if __name__ == "__main__":
    main()
