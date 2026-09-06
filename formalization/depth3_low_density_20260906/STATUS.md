# Low-density depth-three formalization

Status: COMPLETE; both frozen graph targets independently rebuilt and audited.

Aristotle project: `be8de245-5246-41b8-aa5e-ca01f73d8008`.
Task: `2d0cdd1e-7763-4120-a7d1-e327f5f9a2c4`.
Submitted 6 September 2026 at approximately 01:58 UTC. Only `input/` was uploaded;
the pinned dependency manifest and source are included, not local build caches.
Brett explicitly instructed: if the proving run times out, resume it. Resume
the same project and frozen targets; a dropped monitoring connection alone is
not a stopped proving task and should only be reconnected.

User request: independently review and Lean-verify the new results.
Independent review run: `reviews/review-board-20260906-014743/`.

Frozen target: `input/RequestProject/DepthThreeSpec.lean`, SHA-256
`aaf283ca986d1119706e68169cff73fab3c050f99b66ce68736b353e6a6da090`.
The actual finite-forest blocked-shadow theorem AND the low-density graph
inequality are required for completion.

Pre-submission local checks (historical):

- Specification compiles under Lean 4.28.0.
- `DepthThree.low_density_algebra` compiles and depends only on
  `propext`, `Classical.choice`, and `Quot.sound`. This is the scalar assembly,
  not a proof of graph-level hypotheses.
- Existing code-reserve dependencies are supplied. The missing explicit
  graph-count/code bijection is a required new proof obligation.
- The unrelated `ExhaustiveSmallCodes.lean` source was omitted from the packet
  because it is not a dependency and uses `native_decide`; the original is intact.

Request SHA-256: `ef0334cb7ba45f55e095077e50d399af83b800c206a7869c90f838134edb72e5`.
## Local verification of the returned source

Completed 6 September 2026, approximately 02:42 UTC. The provider returned
COMPLETE_WITH_ERRORS while its report claimed COMPLETE. Local evidence resolves
that discrepancy: `lake build` from a fresh project build directory passed all
8,049 jobs, and `LocalAudit.lean` passed both exact target ascriptions and the
axiom audit. The full logs are retained as `local-build.log` and `local-audit.log`.

Proved declarations:

- `DepthThree.blocked_shadow_of_forest : DepthThree.BlockedShadowTarget`.
- `DepthThree.low_density_of_forest : DepthThree.LowDensityTarget`.

Both and every explicitly audited graph/count bridge use only `propext`,
`Classical.choice`, and `Quot.sound`. Source scan found no proof holes, new
axioms, native_decide, implemented_by, unsafe declarations, or externs. The
frozen specification and every supplied dependency are byte-for-byte unchanged.
The main agent read all six new modules, including the fibrewise bijection that
counts each partial independent set once, not once per maximum completion.

Returned source: `output-first/input_aristotle/`.
Archive: `result-20260906T0239.tar.gz`, SHA-256
`c2e5c7ab4626643aa087b806f709db492ba5ed0eee9b52bcc06e76686e9b9c76`.
`DECLARATION_MAP.md` in the returned source maps every supporting declaration.

This completes the two low-density packet targets, not the separate b1-zero
window theorem and not Erdős #993. The same project will now be continued for
the b1-zero target, reusing the locally verified bridge and low-density theorem.
