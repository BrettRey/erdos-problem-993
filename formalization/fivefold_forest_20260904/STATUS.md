# Fivefold forest verification

The full auxiliary finite-forest theorem is **Lean-verified and locally
rebuilt**, including the graph-to-count bridge. For every finite forest with
independence number at least five, `b₂ = 0` implies `5 b₃ ≤ e₃` and
`5 b₄ ≤ e₄`. This does not resolve Erdős #993 or the remaining `b₂ > 0` branch.

- Project: `a4a50054-798f-4d51-b4e7-87f8b9f029d0`
- Original task: `7b254d22-2403-4b14-a5d5-e11ddd91f8a8`
- Submitted: 2026-09-04, task creation `2026-09-04T22:29:04.063310` UTC
- Original task ended `OUT_OF_BUDGET` at `2026-09-05T01:00:17.064123` UTC.
- Brett requested the exact continuation prompt `resume`.
- Resumed task: `dd8a1e5b-085b-40b9-a8f7-cfac6e926293`, created
  `2026-09-05T01:23:50.970902` UTC; **COMPLETE** at
  `2026-09-05T01:36:31.237000` UTC (4 September, Toronto).
- Input: `aristotle_input.md`
- Input SHA-256: `50c29d7ec1a52f12832aac2bc394b6dad9857b046c82f2e4db8c1e3e1bcaef52`
- Compiled specification: `Statement.lean`
- Specification SHA-256: `e44ac947a12885cc7208cb00226898e96ede579fab083c4bf9f08623faa80cce`
- Local checker: Lean 4.28.0, existing repository Mathlib v4.28.0
- Exact compilation command: `lake env lean formalization/fivefold_forest_20260904/Statement.lean`
- Independent reviews: `../../reviews/fivefold-20260904-182512/`

The independent Codex review returned **VALID AS WRITTEN**. The separate
Opus review returned **VALID WITH EXPLICIT REPAIR**: cite the existing
upper bound in the independence-number-one star case and expose the use of
the zero blocked coefficients in the three-leaf recurrence. Those explanatory
changes are now incorporated in the live note. The submitted packet and the
review source remain frozen; neither theorem nor proof strategy changed.
Raw reviews, hashes, and the exact verdicts are in the review directory.
No backend task remains active. No further proving run was submitted.

Submission command:

```bash
python3 scripts/aristotle_cli.py formalize \
  formalization/fivefold_forest_20260904/aristotle_input.md
```

Final result and local verification:

```bash
python3 scripts/aristotle_cli.py result a4a50054-798f-4d51-b4e7-87f8b9f029d0 \
  --destination formalization/fivefold_forest_20260904/resumed-result.tar.gz
cd formalization/fivefold_forest_20260904/output-final_aristotle
lake build
lake env lean ../AxiomAudit.lean
```

- Final archive SHA-256:
  `07ce5a45aa4f8c0a0645f56f7619ebf49b100e85e8312138fd055f0022b405cd`.
- Main source: `output-final_aristotle/RequestProject/Main.lean:115`,
  `FivefoldForest.fivefold_of_forest`; SHA-256
  `8cfd458fe8fe0eb1126cd13136d2db432321abde83f492d86ad225269bdacdf1`.
- Build: exit 0, 3,144 jobs, all twelve returned Lean modules compiled from
  source; no warnings. The project pins Lean 4.28.0 and Mathlib commit
  `8f9d9cff6bd728b17a24e163c9402775d9e6a365`. Existing local dependency
  caches were reused after checking all nine pinned revisions and clean
  tracked package sources. Logs: `resumed-build.log` (local generated log).
- `AxiomAudit.lean` independently restates the counts and proves the requested
  statement from the returned theorem. It and the main theorem plus twenty
  supporting declarations depend only on `propext`, `Classical.choice`,
  `Quot.sound`. Exact successful output: `resumed-axioms.txt`; exit 0.
- Source scan found no `sorry`, `admit`, added `axiom`, `implemented_by`,
  `native_decide`, `unsafe`, `run_tac`, or `run_elab`. Definitions count actual
  independent sets contained in **maximum**, not merely maximal, sets; counts
  at negative sizes are zero. No hidden profile, order, or deficiency premise.
- The original out-of-budget artifact already built locally. The resumed
  source only removes unused hypotheses from three cone lemmas and updates
  their call sites; the final theorem and definitions are unchanged.
- Our first independent read-back harness hit an elaboration timeout, then
  needed an explicit finite-set equality argument. Only the audit harness
  was repaired; no returned proof was edited. Its final axiom audit is clean.
- Original archive and extracted source are retained separately. Aristotle's
  service summary mentions its own commits. At the verification checkpoint
  no local commit, push, or public upload had occurred. Brett authorized
  repository shipping on 5 September; Prove2Me publication remains separate.
  Novelty has not been established.
