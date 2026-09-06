# Summary of changes for run 2d0cdd1e-7763-4120-a7d1-e327f5f9a2c4
## Result: G0 COMPLETE — both frozen targets proved

`RequestProject/DepthThreeSpec.lean` is byte-for-byte unchanged; all original dependency modules are untouched. Five new modules were added, all covered by the existing `RequestProject.+` default Lake target. `lake build` succeeds from a cleaned project build directory with no `sorry`/`admit`, no added axioms, no `@[implemented_by]`, no `native_decide`, and no external certificates.

**Targets** (in `RequestProject/DepthThreeTargets.lean`):
- `DepthThree.blocked_shadow_of_forest : DepthThree.BlockedShadowTarget`
- `DepthThree.low_density_of_forest : DepthThree.LowDensityTarget`

`#print axioms` (emitted during the build) gives for both: `[propext, Classical.choice, Quot.sound]`.

**Supporting targets, all complete** (full table in `DECLARATION_MAP.md`):

- **G1 — the actual graph/count ↔ code bridge.** `RequestProject/GraphBridge.lean` proves that a finite graph has a maximum matching (`exists_maximum_matching`), that every finite forest underlies a `TreeMatching`, that bags partition the vertices with at most two vertices each, that `G.indepNum = Fintype.card D.Bag` (both directions, via König and the vertex-cover bound), and that `D.maxIndepSets` is exactly Mathlib's set of maximum independent sets (`isMaximumIndepSet_iff_mem_maxIndepSets`). `RequestProject/ExtendableCount.lean` then proves the cardinality-preserving fibrewise bijection `card_extSets_fiber` between extendable independent vertex sets with a given bag support and the projected partial words of the code of maximum independent sets (flips handled through `coverWord`), giving `card_extSets`, the guarded identities `e_eq_card_extSets` / `b_eq_card_blkSets`, and finally `e_eq_erasure : DepthThree.e D.G d = D.erasure d` for `d ≤ α`. This is what turns `TreeCodeBridge.erasure_depth_three_reserve` into the graph Pascal reserve; it was proved, not assumed.

- **G2 — unary-or-pair certificate and blocked shadow.** `RequestProject/BlockedCertificate.lean` reads extendability on the Boolean constraint system (`exists_maxIndep_superset_iff`) and, using `PosetCode.codeProj_idealCode` through `exists_idealIndicator_extending`, proves `exists_small_blocked_subset`: every blocked independent set contains a blocked subset of cardinality at most two. `RequestProject/BlockedShadow.lean` proves `card_extVerts_le_six` (an independent `(α−3)`-set leaves three bags empty, each of size ≤ 2) and `card_blocked_deletions_ge` (at least `α−4` deletions outside a chosen certificate stay blocked — deletions inside the certificate are not assumed blocked), and double counts the incidences to get `(α−4) b₂ ≤ 6 b₃`.

- **G3 — extendable incidence and positivity.** `RequestProject/LowDensity.lean` proves `4 e₄ ≤ (α−3) e₃` for every forest with `α ≥ 4` by double counting extendable one-vertex extensions, and `extSets_nonempty`, giving `e₂ > 0` and `e₄ > 0`.

- **G4 — numerical assembly.** `depthThreeStrict_of_lowDensity` combines the graph Pascal reserve, the blocked shadow bound, the incidence bound and positivity, casts to ℝ, and applies the existing scalar implication `DepthThree.low_density_algebra`.

Scope notes: disconnected forests and isolated vertices are included; there is no restriction on vertex count, matching deficiency, `b₁` or `b₂`; "extendable" is the specification's `Extendable` (containment in a *maximum* independent set); the `d ≤ indepNum` guard is carried explicitly; partial sets are counted once, not once per maximum completion. `DepthThree.B1ZeroWindowTarget` was not part of the request and is not attempted. Nothing here claims to resolve Erdős #993.
