# b1-zero window formalization

Status: **graph-level Lean verification passed locally on 6 September 2026**.

The saved checkpoint proves the exact frozen theorem
`DepthThree.b1_zero_window : DepthThree.B1ZeroWindowTarget`:
for every finite tree with `33 <= n <= 38`, `17 <= alpha <= 19`,
`2 alpha <= n + 5`, and `b_1=0`, the actual graph counts satisfy
`s_2 s_4 < s_3^2`. No representation, enumeration, or density hypothesis
has been added to the target. This does not settle the b1-positive branch,
the full depth-three window, or Erdős #993.

## Local verification

- Verified source: `checkpoint-timeout-1/` (52 Lean modules).
- Fresh default `lake build`: **PASS, 8,078 jobs**, exit 0. Project modules
  were rebuilt from source; only the existing pinned dependency-package
  cache was reused. Warnings concern unused variables and section variables only.
  Complete output: [local-checkpoint-build.log](local-checkpoint-build.log).
- Exact target ascriptions for the b1-zero, low-density, and blocked-shadow
  graph theorems: **PASS**. Target definitions were printed and read back.
  The frozen specification, dependency manifest, and toolchain are unchanged.
- All three targets and the audited structural, count, representation,
  root-switch, and finite-certificate declarations use only `propext`,
  `Classical.choice`, and `Quot.sound`. Audit source: [LocalAudit.lean](LocalAudit.lean);
  complete output: [local-checkpoint-audit.log](local-checkpoint-audit.log).
- The full returned Lean source scan has no `sorry`, `admit`, extra `axiom`,
  `native_decide`, `implemented_by`, `unsafe`, `extern`, or elaboration/IO
  escape hatches. The finite check uses kernel reduction with proved
  enumeration coverage, not imported JSON rows.
- Source hashes: [returned-checkpoint-source-sha256.txt](returned-checkpoint-source-sha256.txt).
  Frozen specification SHA-256:
  `aaf283ca986d1119706e68169cff73fab3c050f99b66ce68736b353e6a6da090`.
  Final target module SHA-256:
  `aa7be0df446ceea33b4ca6f670d3f899bf2fd6ec8a4ef8dba1324820f75252fe`.
- Saved download: `result-at-timeout-20260906T0450.tar.gz`, SHA-256
  `5a27724e939002e7f1a7a052def5e2c3edfb69166676316dfcce5d2c3c4c2583`.
  The archive filename identifies the timed-out task; it was downloaded
  after the morning resumption. The provider status is not the proof warrant.

## What the formal proof establishes

The 21 new `RequestProject/B1Zero/` modules prove the allowed/forced/forbidden
decomposition, the pendant matching structure, actual finite-set convolutions,
the graph-to-rooted-corona representation, and the cardinality-preserving
root-switch domination. For at most three forbidden vertices, cone-ray density
bounds suffice. In the sole four-forbidden case, the graph has `alpha=19`,
nine forced vertices, and ten flexible pairs. The proof bounds its actual
blocked count by the union-bound numerator
`73 p0 + 24 p1 + 3 p2 + 15 q0 + 6 q1 + q2`. A kernel-checked inequality over
all 16,796 plane rooted forests, with proved coverage and evaluator semantics,
shows this numerator is at most one quarter of the extendable count. The
already verified low-density forest theorem completes the exact target.

This is a different sufficient route from the written proof's sharper
`52513/217404` graph bound. The formalization does **not** establish that
sharp graph bound or the special `b_3 >= 11 b_2` estimate, and does not need
either. The sharper encoded-forest certificate remains separately proved.
Independent mathematical review of the written route is preserved separately.

The [independent source-semantic follow-up](../../reviews/review-board-20260906-014743/lean-checkpoint-followup.md)
also passed after reading all 21 new modules and the common graph/count
layer. It performed no second build; the successful compilation and axiom
audit above are separate main-agent checks. Its full report and exact prompt
are preserved with hashes in the review manifest.

**Provider documentation qualification:** `checkpoint-timeout-1/DECLARATION_MAP.md`
lines 111--113 state that connectedness, vertex-count, deficiency, and b1
restrictions are unused. That is correct only for the general blocked-shadow
and low-density layer, not the B1Zero development. The final B1Zero declaration
explicitly uses the frozen tree-window hypotheses. The raw returned snapshot
is preserved byte-for-byte; this qualification corrects its interpretation,
not the Lean source.

## Submission history

The first submission attempt returned “too many projects in progress” and
created no project. After the low-density source passed local verification,
the same project `be8de245-5246-41b8-aa5e-ca01f73d8008` was continued at
approximately 02:43 UTC on 6 September 2026. No running task was canceled.

Initial task: `a5f2e62f-3d9e-49e2-a9e6-ce29a5673346`.
Accepted request SHA-256:
`e5426f9525ebab3ae130d40e183a095669833b5a338572a30708b21ac1e45107`.
Support archive SHA-256:
`6bd277ddac7e88cf6ed0957f5b0f28516b61938fe1d898253d7dbc986788b8fa`.
Only the request and the scoped support archive were added to the existing
project. The frozen specification and six verified graph modules are retained.
Brett explicitly asked for automatic resumption on timeout; continue this same
project with the target unchanged if its actual proving task times out.

## Resumption, 6 September 2026

The initial b1-zero task exhausted its time budget at 04:50:48 UTC, without a
final report. After Brett said “resume”, a fresh status check confirmed
OUT_OF_BUDGET and the same project was resumed at 09:59:49 UTC.
Resumed task: `9a3460e6-3088-4ed4-beed-1ee665e9684a`.
The continuation explicitly preserves the target, the completed graph/count
bridges, and all saved B1Zero work. No new target or extra hypothesis was added.
The saved checkpoint was downloaded separately for local inspection. It
already contained a complete target proof, which subsequently passed the
local build and axiom checks above; no further proof completion was needed.
After those checks and the independent source-semantic review passed, the
redundant resumed task was canceled at 10:09:36 UTC. Its status was confirmed
as CANCELED. Only this session's unnecessary continuation was stopped; the
verified local checkpoint, its archive, and all source/logs remain preserved.

User request: independently review and Lean-verify the new results.
Independent review run: `reviews/review-board-20260906-014743/`.
Final repository checks: all 55 core tests and all four proof-frontier tests
passed; the proof graph validates (16 claims, 30 edges, six bridges), and
`git diff --check` passes. No manuscript edit, commit, or push was made.

Subsequently, Brett requested commit and push. The ship gate, both packet
builds and exact-target audits, all 59 Python tests, and proof-graph validation
passed again. The source snapshots and verification logs are included in the
scoped repository change; download tarballs and generated caches remain local.
This does not alter the verified source or publish a revised manuscript.
The staged whitespace scan flags only inherited trailing spaces in the input
and returned copies of `indpoly.py`; both are byte-identical to the existing
repository helper and are preserved to keep the frozen source hashes valid.
All other staged files pass the whitespace check.

Frozen target: `input/RequestProject/DepthThreeSpec.lean`, SHA-256
`aaf283ca986d1119706e68169cff73fab3c050f99b66ce68736b353e6a6da090`.
The frozen completion requirement was the actual universal finite-tree window theorem, including
structural, graph-count/code, domination, and finite-enumeration completeness
bridges. These obligations are now discharged; a check of the 1,842 imported
coefficient rows alone would not have been sufficient.

The specification and `DepthThree.four_forbidden_algebra` compile under Lean
4.28.0; the latter has only the standard three axioms and is a scalar implication,
not the graph theorem. The independent review reproduced
the finite arithmetic; that is not a Lean proof. Existing code-reserve modules
are included for reuse; their actual graph-count interface has now been proved
and independently replayed in the previous low-density task.
The unrelated native-decide short-code enumeration was omitted from the packet;
the original source remains untouched.

Initial rejected request SHA-256:
`ca5d0300e597a9f6b1313a9c8ec18cf585356444c8fccc2ece900819fef909e5`.
The accepted continuation includes the additional supporting proofs below.
The following records describe supporting work before the final graph proof.

### Earlier local supporting results (historical)

- `DepthThree.ParameterCertificate.all_small_forbidden_rays` proves the exact
  rational bound for every ray in all twelve r<=3 parameter cases.
- `DepthThree.ParameterCertificate.five_ray_cone` proves the five-term monotone
  cone implication. Both compile with only the standard three axioms.
- `DepthThree.RootedForestCertificate.rooted_forest_bound` now proves the bound
  for every encoded plane rooted forest of size ten, using a proved recursive
  enumeration coverage theorem and kernel reduction. The enumeration has
  exactly 16,796 entries (also proved), encompassing the ordered versions of
  the 1,842 unordered types. No imported row list or native-decide oracle is
  used. The graph/encoding/evaluator bridge remains a separate obligation.
- Final finite-module replay: 22.68 seconds wall time, exit 0. The coverage
  and bound have only the standard three axioms; the exact enumeration count
  has no axioms. Source and complete output are retained in this packet.
- The assembled packet, including four reused forest-count/decomposition
  modules under `RequestProject.Fivefold`, passed `lake build` (8,049 jobs).
  Their only change from the verified fivefold project is the import prefix.
  Source scan found no proof holes, extra axioms, native_decide, or
  implemented_by. Full local output: `local-build.log`.
- Additional local replay passed for `FivefoldBridge.e_eq` and `b_eq`, which
  identify the internal finite-set helper counts with the frozen DepthThree
  graph counts, and for `small_forbidden_algebra`, the joint-density scalar
  assembly. All use only the standard three axioms. This count bridge does
  not supply the still-required matching-code bijection.
- `PolynomialCertificate.concentration_coeff_le` and `profile_eq_firstFive`
  passed locally with only the standard three axioms: the concentration step
  is coefficientwise, and the finite evaluator agrees with its full reversed
  polynomial recurrence. Graph interpretation remains separate.
- `all_small_forbidden_density_rays` passed in 18 seconds with the standard
  three axioms: every r<=3 cone ray is below the corresponding low-density
  threshold. The r=4 ratio 52513/217404 is below 1/4. The preferred final route
  therefore uses the proved low-density graph theorem after proving the
  structural and graph/domination/evaluator bridges; the separate b2 and
  special b3 bounds need not be formalized for completion.
- All support plus the six verified graph modules passed the combined default
  build (8,057 jobs); output is `local-extended-build.log`.
