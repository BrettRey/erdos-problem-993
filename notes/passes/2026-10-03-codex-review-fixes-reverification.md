# Re-verification after the Codex review fixes
<!-- SUMMARY: pass re-checks after the cross-family Codex review (no math error) and its four fixes · status: diff-scoped, no new findings · updated: 2026-10-03 -->

**Diff:**
- line 71: the zero extension now covers all k outside {0,…,n};
- outline: the formalization's scope is named exactly (HJ cubics, strict
  log-concavity, maximal-mass bound) instead of "results of … Newton …";
- abstract, last sentence: "including the cubic inequalities and the
  maximal-mass bound it uses".

Outside the manuscript, `scripts/verify_pb_cue_threshold.py` now rounds
outward.

**External review:** Codex (GPT-6.1-sol, OpenAI) found no mathematical error in
§§2–4 or in the formal statements. Record: `runs/pb-codex-review-20261003/`.
This is the cross-family check that the new §4.1 had lacked.

| Pass | Check | Result |
|---|---|---|
| negative-claims-audit, source-reread | The narrower scope claim matches `PBDeduction/Defs.lean` (`DeductionHyp.strict_lc`) and the theorem `pb_strict_lc`. Line 72 still credits Newton's inequalities for log-concavity, which is accurate. | Supported. |
| adversarial-cold-read, contribution-alignment | Only the method sentence of the abstract changed, and it now names its scope more precisely. Portal abstract re-synced. | Round-9 verdict stands. |
| numbers-audit | The CUE program's enclosures are now outward-rounded and contain the exact values. The paper's printed 0.215 and 0.0004 are unchanged and still supported. | Correct. |
| build-integrity, proofread, terminological-hygiene, level-category, reader, coherence, prose passes | 0 undefined references, 0 overfull boxes, 13 pages. Read each changed sentence. | No issue. |
