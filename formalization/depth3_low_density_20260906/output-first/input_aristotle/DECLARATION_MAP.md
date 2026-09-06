# Declaration map and completion status

Toolchain: Lean 4.28.0, Mathlib as pinned by `lake-manifest.json`.
Build: `lake build` (all modules under `RequestProject/` are default Lake targets via the
`RequestProject.+` glob).

`RequestProject/DepthThreeSpec.lean` is **unchanged**.

## Completion status

| Target | Status |
| --- | --- |
| G0 `DepthThree.blocked_shadow_of_forest : DepthThree.BlockedShadowTarget` | **COMPLETE** |
| G0 `DepthThree.low_density_of_forest : DepthThree.LowDensityTarget` | **COMPLETE** |
| G1 graph/count ↔ code bridge and graph Pascal reserve | **COMPLETE** |
| G2 unary-or-pair certificate and blocked-shadow inequality | **COMPLETE** |
| G3 extendable incidence inequality and positivity | **COMPLETE** |
| G4 numerical assembly and final low-density theorem | **COMPLETE** |
| `DepthThree.B1ZeroWindowTarget` | not attempted (not requested) |

No `sorry`, no `admit`, no added `axiom`, no `@[implemented_by]`, no `native_decide`,
no external certificates.

`#print axioms` (emitted by `RequestProject/DepthThreeTargets.lean` during the build):

```
'DepthThree.blocked_shadow_of_forest' depends on axioms: [propext, Classical.choice, Quot.sound]
'DepthThree.low_density_of_forest'   depends on axioms: [propext, Classical.choice, Quot.sound]
```

## New modules

### `RequestProject/GraphBridge.lean` — G1, structural bridge

| Declaration | Statement |
| --- | --- |
| `MatchingBag.exists_maximum_matching` | a finite graph has a maximum matching (finite maximisation) |
| `MatchingBag.exists_treeMatching_of_isAcyclic` | every finite forest underlies a `TreeMatching` |
| `MatchingBag.isIndepSet_coe_iff` | Mathlib `IsIndepSet` on a `Finset`, unfolded |
| `MatchingBag.TreeMatching.bagSet` | the vertex set of a bag (matched pair, or singleton) |
| `MatchingBag.TreeMatching.bagOf` | the bag of a vertex; `bagOf_eq_iff` shows bags partition `V` |
| `MatchingBag.TreeMatching.card_bagSet_le` | every bag has at most two vertices |
| `MatchingBag.TreeMatching.bagOf_injOn_indep` | an independent set meets each bag at most once |
| `MatchingBag.TreeMatching.indepNum_eq_card_Bag` | `α(G) = #bags` (König, both directions) |
| `MatchingBag.TreeMatching.isMaximumIndepSet_iff_mem_maxIndepSets` | `D.maxIndepSets` = actual maximum independent sets |
| `MatchingBag.TreeMatching.extendable_iff` | `DepthThree.Extendable` in terms of `D.maxIndepSets` |
| `MatchingBag.TreeMatching.pick`, `pick_coverWord_mem`, `eq_pick_coverWord` | the unique vertex of a maximum independent set inside a bag |

### `RequestProject/ExtendableCount.lean` — G1, the counting bridge

| Declaration | Statement |
| --- | --- |
| `MatchingBag.TreeMatching.coverWord_agree` | on the bags met by `S`, the word of `S` equals that of any maximum independent set containing it (this is where the flips are handled, via `coverWord`) |
| `MatchingBag.TreeMatching.extSets` / `blkSets` | extendable / blocked independent sets of a given size, counted once each |
| `MatchingBag.TreeMatching.card_extSets_fiber` | **the cardinality-preserving bijection**: extendable sets with bag support `K` ↔ `codeProj K` of the code of maximum independent sets |
| `MatchingBag.TreeMatching.card_extSets` | `#extSets k = codeP D.treeCode k` |
| `MatchingBag.TreeMatching.e_eq_card_extSets`, `b_eq_card_blkSets` | the guarded specification counts `DepthThree.e`, `DepthThree.b` |
| `MatchingBag.TreeMatching.e_eq_erasure` | **`DepthThree.e D.G d = D.erasure d`** for `d ≤ α`; this is what turns `TreeCodeBridge.erasure_depth_three_reserve` into the graph Pascal reserve |

### `RequestProject/BlockedCertificate.lean` — G2, the certificate

| Declaration | Statement |
| --- | --- |
| `MatchingBag.exists_idealIndicator_extending` | a partial assignment violating no comparison inside its support extends to an order-ideal indicator (derived from `PosetCode.codeProj_idealCode`) |
| `MatchingBag.TreeMatching.subset_compl_coverOf_iff` | a vertex set sits in the independent set complementary to `coverOf x` iff `x` takes the prescribed value on each bag met |
| `MatchingBag.TreeMatching.exists_maxIndep_superset_iff` | extendability read on the Boolean constraint system |
| `MatchingBag.TreeMatching.exists_small_blocked_subset` | **every blocked independent set contains a blocked subset of cardinality ≤ 2** |

### `RequestProject/BlockedShadow.lean` — G2, the inequality

| Declaration | Statement |
| --- | --- |
| `MatchingBag.TreeMatching.blocked_mono` | supersets of blocked sets are blocked |
| `MatchingBag.TreeMatching.extVerts`, `card_extVerts_le_six` | an independent `(α-3)`-set has at most six one-vertex extensions |
| `MatchingBag.TreeMatching.card_blocked_deletions_ge` | a blocked `(α-2)`-set has at least `α-4` blocked one-vertex deletions (those outside a chosen certificate; deletions inside the certificate are *not* assumed blocked) |
| `MatchingBag.TreeMatching.blocked_shadow` | `(α-4) b₂ ≤ 6 b₃` |

### `RequestProject/LowDensity.lean` — G3 and G4

| Declaration | Statement |
| --- | --- |
| `MatchingBag.TreeMatching.extSets_nonempty` | `e_d > 0` for `d ≤ α` (positivity of `e₂`, `e₄`) |
| `MatchingBag.TreeMatching.extendable_incidence` | `4 e₄ ≤ (α-3) e₃` |
| `MatchingBag.TreeMatching.depthThreeStrict_of_lowDensity` | numerical assembly with the graph Pascal reserve, `blocked_shadow`, `extendable_incidence` and `DepthThree.low_density_algebra` |

### `RequestProject/DepthThreeTargets.lean` — G0

| Declaration | Statement |
| --- | --- |
| `DepthThree.blocked_shadow_of_forest` | `DepthThree.BlockedShadowTarget` |
| `DepthThree.low_density_of_forest` | `DepthThree.LowDensityTarget` |

## Reused, unchanged dependencies

`Codes.lean`, `CodeInvariance.lean`, `PosetCode.lean`, `KonigHall.lean`, `ForestLemmas.lean`,
`TreeMatching.lean`, `BagPoset.lean`, `PascalSmoothing.lean`, `PascalBridge.lean`,
`TreeCodeBridge.lean`, `UltraLogConcave.lean`, `DepthThreeAlgebra.lean`,
`DepthThreeSpec.lean` are untouched.

## Notes on the argument

* Disconnected forests and isolated vertices are covered: nothing in the development assumes
  connectedness or a positive vertex count.
* No restriction on the vertex count, the matching deficiency, `b₁` or `b₂` is used.
* "Extendable" always means "contained in a *maximum* independent set", never merely a
  maximal one; this is `DepthThree.Extendable`, unchanged.
* Counts at negative sizes are zero, as encoded by the `d ≤ indepNum` guard; the bridges
  `e_eq_card_extSets` and `b_eq_card_blkSets` carry that guard as a hypothesis.
* Partial sets are counted once, not once per maximum completion: `extSets`/`blkSets` are
  `Finset (Finset V)`, and the fibrewise bijection `card_extSets_fiber` maps each set to a
  single projected partial word.
