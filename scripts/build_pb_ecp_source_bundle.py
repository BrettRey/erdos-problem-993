#!/usr/bin/env python3
"""Build the ECP source bundle for the Poisson--binomial paper.

The ECP sample (ejpecp/sample.tex, L439--442) asks for the bibliography inside
the document, so this copies the manuscript with the \\bibliographystyle and
\\bibliography lines replaced by its current .bbl, adds ejpecp.cls, and
zips the result. It refuses to build while the supplement DOI placeholder is
still in the source, unless --allow-placeholder is given.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys
import zipfile

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "paper" / "poisson_binomial"
PLACEHOLDER = "ZENODO-DOI-PENDING"
STEM = "variance-scaled-turan-first-descent"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--allow-placeholder", action="store_true")
    parser.add_argument("--out", type=Path, default=PAPER / "ecp_source_bundle.zip")
    args = parser.parse_args()

    tex = (PAPER / f"{STEM}.tex").read_text(encoding="utf-8")
    if PLACEHOLDER in tex and not args.allow_placeholder:
        print(f"refusing: {PLACEHOLDER} still in {STEM}.tex", file=sys.stderr)
        return 2
    bbl = (PAPER / f"{STEM}.bbl").read_text(encoding="utf-8")
    old = "\\bibliographystyle{amsplain}\n\\bibliography{references}\n"
    if tex.count(old) != 1:
        print("refusing: bibliography lines not found exactly once", file=sys.stderr)
        return 2
    inlined = tex.replace(old, bbl)
    with zipfile.ZipFile(args.out, "w", compression=zipfile.ZIP_DEFLATED) as z:
        z.writestr(f"{STEM}.tex", inlined)
        z.write(PAPER / "ejpecp.cls", "ejpecp.cls")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
