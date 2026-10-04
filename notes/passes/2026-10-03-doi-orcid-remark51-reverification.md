# DOI, ORCID, Remark 5.1 cut, and a missed build error
<!-- SUMMARY: reserved Zenodo DOI inserted; wrong ORCID fixed; Remark 5.1 moved to the supplement; TeX errors in the supplement box found and fixed · status: clean build (errors checked) · updated: 2026-10-03 -->

## Changes

- **DOI.** Zenodo draft 23129480 has reserved DOI 10.5281/zenodo.23129480. It
  now replaces `ZENODO-DOI-PENDING` in the supplement description and appears
  in CERTIFICATE.md ("Scope and archive status"). The draft is not published.
- **ORCID.** The paper, the cover letter, the portal fields,
  `submission-info.md`, `.zenodo.json` and `paper/main.tex` (the Erdős
  manuscript) all carried 0000-0003-2407-9448. Brett's ORCID is
  0000-0003-0073-7195: it is in his house template, appears 360 times across
  the portfolio, and is what Zenodo's directory lists for "Reynolds, Brett".
  All six files are fixed.
  - Brett's two existing Erdős Zenodo records (19100781, 18745547) list no
    ORCID, so the wrong ID was not attached publicly there.
  - The E-JC submission of the Erdős manuscript did print the wrong ID.
- **Remark 5.1** (Pitman's ratio bound against (1.6)) was cut, as Brett
  approved. Its text is now `SUPPLEMENTARY-REMARKS.md`, which is in the archive
  and listed in CERTIFICATE.md.

## Missed build error (now fixed)

The restructure (commit ff4d15a, recommendation E) split the supplement
description into three paragraphs with blank lines. The class defines
`\sdescription` with plain `\def`, so its argument cannot contain a paragraph
break. That produced two TeX errors: "Paragraph ended before \sdescription was
complete" and "Extra }". In the PDF, **the first paragraph, with the archive
link and the computations the proof needs, was missing**.

From ff4d15a until now, every build-integrity record checked only
"undefined|Overfull" and so missed the errors. The pasted ChatGPT review saw
the symptom ("I couldn't find an embedded archive link").

Fix: the paragraphs are separated by `\endgraf\smallskip`, which a non-long
argument accepts. The log now has 0 errors, 0 undefined references and 0
overfull boxes, and the box renders all four parts. From here on, build checks
also count lines beginning with `!` in the log.

## Checks

13 pages. The PDF links resolve to the correct ORCID (checked after qpdf
decompression). The supplement archive was rebuilt (77 members) and verified.
