# Re-verification after the genre restructure (recommendations A–H)
<!-- SUMMARY: checks after Section 1 was reorganized, Remarks added, Section 5 created, and the supplement and acknowledgement rewritten · status: no content lost; clean build · updated: 2026-10-03 -->

Brett accepted all of the recommendations in
`notes/ecp-genre-assessment-2026-10-03.md`.

## What was done

- **A. Section 1 reorganized.** It now holds the setting, gap and Theorem 1.1
  (restated so it can be read alone), statements of Prop 1.2, Cor 1.3,
  Prop 1.5 and Prop 1.6, a short Example 1.9, and the proof strategy. The new
  Section 5, "Proofs of the results of Section 1", has the four proofs, the
  balanced family with Table 1, Remark 5.1 (Pitman's ratio bound) and the
  Haar-unitary derivation. Section 1 now ends on p. 4 instead of p. 6.
- **B. Remarks.** Remarks 1.4 (the √V window and W − W′), 1.7 (reciprocal
  bound against ULC), 1.8 (real-rooted polynomials) and 5.1 (Pitman's ratio
  bound).
- **C. Balanced-family passage split** into four paragraphs. The bound
  derived from Pitman is displayed as (5.1), and θ is defined where it is used.
- **D. Roadmap trimmed.** The parameter ranges and the list of lost constants
  went; one sentence now says the argument does not determine κ⋆.
- **E. Supplement description** split into three parts: what the proof needs,
  the additional check, and formalization.
- **F. Abstract.** The curvature comparison is now its own clause, and the
  last sentence says "the proof of the main theorem".
- **G. Acknowledgement.**
  - Every specific credit is kept, and dates and model labels moved to
    CERTIFICATE.md ("Provenance of AI contributions").
  - The ChatGPT and Elicit outputs are now called "AI-generated reviews".
  - The second Codex review, which prompted this reorganization, is credited.
- **H. Example 1.9.** It states the result. The derivation and the
  N = 1998 program sit in §5.4 and the supplement.
- **I. Smaller fixes.** The §4.1 forward reference to S and T was replaced by
  a remark at the point in §4.2 where they are defined. Theorem 1.1 is now
  self-contained.

## Checks

| Check | Result |
|---|---|
| Nothing lost | A sentence-level diff of the old and new body (pre-restructure copy in the session scratchpad) shows every removed sentence either moved, reworded as listed above, or deliberately dropped: the roadmap's parameter ranges (still in §4) and the list of lost constants (replaced). The four proofs are verbatim apart from the intended edits (the forward reference "of Section 2" removed, the Darroch citation added, "(1.7) below" becoming "(1.7)"). |
| Numbering | Remarks share the theorem counter, so Prop 1.4 became 1.5, Prop 1.5 became 1.6, and Example 1.6 became 1.9. Every internal reference uses `\ref`. External current files were updated (CERTIFICATE.md, the CUE program's docstring, the portal AI disclosure, both replay records), and a final grep finds no stale numbers. Equation numbers in §§1–2 are unchanged, and the new bound is (5.1). |
| Build | 0 undefined references, 0 overfull boxes, 13 pages. Supplement rebuilt and verified. |
| Proofread | Read pp. 2–4 and pp. 10–11 in the PDF. |
| Cold read | Not rerun, by Brett's decision (2026-10-03). Pages 1–2 keep the setting, the gap and Theorem 1.1. The Prop 1.2 proof moved out. |
