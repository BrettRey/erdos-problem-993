# EJMS template-compliance changes and re-check
<!-- SUMMARY: journal mode instead of preprint, template metadata block restored, Table 1 moved beside its first callout, to make two EJMS submission checkboxes literally true; no wording changed · status: clean build (errors checked); 13 pp · updated: 2026-10-04 -->

**Why.** EJMS asks the author to tick ten statements about the LaTeX source
before a submission can start. Each was checked against the source; eight
already held. Box 1 holds because the local `ejpecp.cls` (v1.12.3) is
byte-identical to the official copy at github.com/vtex-soft/texsupport.ims-ejpecp,
downloaded on 2026-10-04. Two did not hold as worded:

- **Box 2.** "Document class and option `\documentclass[ECP]{ejpecp}` was not
  changed." The source said `[ECP,preprint]`.
- **Box 9.** "All floating objects … appear in correct places next to their
  callouts." Table 1 was first cited on pages 2–3 but printed at the top of
  page 12.

Brett chose to comply literally (2026-10-04).

**Changes.**
1. **Class line.** `\documentclass[ECP]{ejpecp}`.
2. **Template metadata block.** The template sample's block "set by the
   Managing Editor" is restored: `\VOLUME{0}`, `\YEAR{2023}`, `\PAPERNUM{0}`,
   `\DOI{10.1214/YY-TN}`. Without it journal mode throws 14 undefined-control-sequence errors.
3. **Dates.** `\SUBMITTED{October 4, 2026}` and `\ACCEPTED{\ldots}`. Without
   them journal mode prints "FIXME!" on page 1.
4. **Table 1.** The float moved in the source from §5.3 to just after its
   first full callout in §1 ("Table 1 shows this for a family …"). It now
   prints at the top of page 4. Its content is unchanged.

**Effect on the text.** A word-multiset comparison of the old and new PDF
text shows the additions are journal furniture only:
- the header ("Electron. Commun. Probab. 0 (2023), article no. 0", the
  placeholder DOI, the ISSN, the logo text);
- the footers ("ECP 0 (2023), paper 0.");
- the line "Submitted to ECP on October 4, 2026, final version accepted on …".

Nothing was removed apart from math tokens that pdftotext splits differently.
Numbering is unchanged: Table 1 is still Table 1, and no equation or result
number moved.

**Passes.** The edit marked §1, §5.3 and the preamble as changed. No sentence
of the paper's prose changed, so each stale pass was judged against this
diff:

| Pass | Judgement |
|---|---|
| contribution-alignment | title and abstract unchanged; the metadata is journal furniture |
| level-category, terminological-hygiene, negative-claims, source-reread, external-review-triage | no words changed |
| coherence-cohesion, reader-pass | the table now sits beside the claim it supports, an improvement; read on the rendered page 4 |
| de-ai-prose, editorial-scar-tissue | no prose changed |
| adversarial-cold-read | page 1 gained only header furniture; the abstract and opening are unchanged |
| proofread | rendered pages 1 and 4 read |

**Checks.**
- Clean rebuild: 0 lines beginning `!`, 0 undefined, 0 overfull, 13 pages.
  There is one underfull vbox on page 4 (badness 7081, a page carrying a
  float); visually fine.
- Source bundle rebuilt (SHA-256 `926d432a…b2b9`); it compiles standalone with
  text identical to the manuscript PDF.
