#!/usr/bin/env python3
"""Build the public preprint PDF of the Poisson--binomial paper.

The submitted manuscript uses \\documentclass[ECP]{ejpecp}, which prints the
ECP header and the template's placeholder DOI (10.1214/YY-TN). A public copy
must not carry those (DECISIONS 2026-10-04), so this copies the manuscript into
a temporary folder, switches on the class's preprint option, builds it with
pdflatex/bibtex, and refuses to write the PDF if the build has errors or the
text still shows journal furniture. The canonical .tex is not modified.
"""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "paper" / "poisson_binomial"
STEM = "variance-local-log-concavity-poisson-binomial"
JOURNAL_LINE = "\\documentclass[ECP]{ejpecp}\n"
PREPRINT_LINE = "\\documentclass[ECP,preprint]{ejpecp}\n"
FURNITURE = ("YY-TN", "Electron. Commun. Probab. 0", "Submitted to ECP", "ISSN")


def run(cmd: list[str], cwd: Path) -> None:
    subprocess.run(cmd, cwd=cwd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=PAPER / f"{STEM}-preprint.pdf")
    args = parser.parse_args()

    tex = (PAPER / f"{STEM}.tex").read_text(encoding="utf-8")
    if not tex.startswith(JOURNAL_LINE):
        print("refusing: the manuscript does not start with the journal class line", file=sys.stderr)
        return 2
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        (work / f"{STEM}.tex").write_text(PREPRINT_LINE + tex[len(JOURNAL_LINE):], encoding="utf-8")
        for name in ("references.bib", "ejpecp.cls"):
            shutil.copy(PAPER / name, work / name)
        run(["pdflatex", "-interaction=nonstopmode", f"{STEM}.tex"], work)
        run(["bibtex", STEM], work)
        for _ in range(2):
            run(["pdflatex", "-interaction=nonstopmode", f"{STEM}.tex"], work)
        log = (work / f"{STEM}.log").read_text(encoding="latin-1")
        errors = sum(line.startswith("!") for line in log.splitlines())
        undefined = log.lower().count("undefined")
        overfull = log.count("Overfull")
        if errors or undefined or overfull:
            print(f"refusing: {errors} errors, {undefined} undefined, {overfull} overfull", file=sys.stderr)
            return 1
        text = subprocess.run(["pdftotext", f"{STEM}.pdf", "-"], cwd=work,
                              capture_output=True, text=True, check=True).stdout
        found = [s for s in FURNITURE if s in text]
        if found:
            print(f"refusing: journal furniture in the text: {found}", file=sys.stderr)
            return 1
        pages = subprocess.run(["pdfinfo", f"{STEM}.pdf"], cwd=work,
                               capture_output=True, text=True, check=True).stdout
        npages = next(line.split()[-1] for line in pages.splitlines() if line.startswith("Pages:"))
        shutil.copy(work / f"{STEM}.pdf", args.out)
    digest = hashlib.sha256(args.out.read_bytes()).hexdigest()
    print(f"wrote {args.out} ({npages} pages, {args.out.stat().st_size} bytes)")
    print(f"sha256 {digest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
