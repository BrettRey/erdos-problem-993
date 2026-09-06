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
| G0 `DepthThree.b1_zero_window : DepthThree.B1ZeroWindowTarget` | **COMPLETE** |
| B1 structural reduction (`b₁ = 0`: well covered, forced isolated, pendant bags, `r ≤ 4`) | **COMPLETE** |
| B2 `r ≤ 3` density closure | **COMPLETE** |
| B3 graph-to-rooted-forest representation and root-moving domination | **COMPLETE** |
| B4 `r = 4` density closure and final assembly | **COMPLETE** |

No `sorry`, no `admit`, no added `axiom`, no `@[implemented_by]`, no `native_decide`,
no external certificates.

`#print axioms` (emitted during the build by `RequestProject/DepthThreeTargets.lean` for the
first two targets, by `RequestProject/B1Zero/Target.lean` for `b1_zero_window`, and by
`RequestProject/B1Zero/SimpleCertificate.lean` for the enumeration certificate):

```
'DepthThree.blocked_shadow_of_forest' depends on axioms: [propext, Classical.choice, Quot.sound]
'DepthThree.low_density_of_forest'   depends on axioms: [propext, Classical.choice, Quot.sound]
'DepthThree.b1_zero_window'          depends on axioms: [propext, Classical.choice, Quot.sound]
'DepthThree.RootedForestCertificate.rooted_forest_simple_bound'
                                     depends on axioms: [propext, Classical.choice, Quot.sound]
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


## The `b₁ = 0` window: new modules under `RequestProject/B1Zero/`

All are default Lake targets (`RequestProject.+` glob) and contain no `sorry`.

### Structural reduction (B1)

| Declaration | Statement |
| --- | --- |
| `MatchingBag.TreeMatching.AllowedV`, `ForcedV`, `ForbiddenV`, `FlexV` (`Allowed.lean`) | vertices in some / every / no maximum independent set, and the flexible ones |
| `MatchingBag.TreeMatching.bagSet_flex_of_no_forced`, `bag_of_flex`, `bag_of_forbidden` (`Allowed.lean`) | the bag classification: a bag is a forced singleton, a forbidden vertex with its forced mate, or two flexible vertices |
| `MatchingBag.TreeMatching.card_verts_add_card_forced`, `card_flex` (`BagCount.lean`) | `\|V\| + \|Forced\| = 2α + \|Forbidden\|` and `\|Flex\| = 2(α - \|Forced\|)` |
| `MatchingBag.TreeMatching.B1Zero` (`Pendant.lean`) | the hypothesis `b₁ = 0`, in matching-bag form |
| `MatchingBag.TreeMatching.exists_leaf_in_flex_bag` (`Pendant.lean`) | **pendant lemma**: every flexible bag has a vertex with no allowed neighbour outside its bag |
| `MatchingBag.TreeMatching.exists_maxIndepSet_superset`, `extendable_iff_subset_allowed` (`WellCovered.lean`) | **well-coveredness**: with `b₁ = 0` an independent set is extendable iff it avoids the forbidden vertices |
| `SimpleGraph.IsAcyclic.ncard_edgeSet_add_one_le`, `IsAcyclic.card_adjPairs_le` (`ForestCount.lean`) | a finite forest has at most `\|V\| - 1` edges; edges between two disjoint sets are at most `\|A\| + \|B\| - 1` |
| `MatchingBag.TreeMatching.three_le_card_forcedNbrs`, `card_forcedNbrsSet_ge` (`ForcedDegree.lean`) | every forbidden vertex has `≥ 3` forced neighbours; `\|N_C(S)\| ≥ 2\|S\| + 1` |
| `MatchingBag.TreeMatching.card_flexBags`, `pflex_cone` (`FlexCone.lean`) | `#flexBags + \|Forced\| = α`; the matching-bag deletion inequality `(l+1) p_{l+1} ≤ 2(m-l) p_l` |
| `MatchingBag.TreeMatching.two_mul_card_forbidden_lt` (`DensitySmall.lean`) | `2r + 1 ≤ \|Forced\|` when `r ≥ 1`, hence `r ≤ 4` in the window |

### Counting (B2)

| Declaration | Statement |
| --- | --- |
| `MatchingBag.TreeMatching.card_indSubsets` (`Decompose.lean`) | the convolution `#{independent S ⊆ F ∪ Flex, \|S\| = k} = ∑_{i+j=k} C(\|F\|,i) p_j` |
| `MatchingBag.TreeMatching.card_extSets_eq_EE`, `card_blkSets_le` (`EBounds.lean`) | the extendable count, and the blocked bound `b_k ≤ ∑_t C(r,t) EE(c-2t-1, k-t)` |
| `DepthThree.B1Zero.density_numeric` (`Numeric.lean`) | the twelve-case finite parameter LP for `r ≤ 3` |
| `MatchingBag.TreeMatching.density_small` (`DensitySmall.lean`) | **`3(α-3) b₄ ≤ (α-7) e₄` for `r ≤ 3`** |

### The corona bridge (B3)

| Declaration | Statement |
| --- | --- |
| `MatchingBag.B1Zero.cntRev`, `cntRev_choose`, `cntRev_union` (`Convolution.lean`) | independent sets counted by codimension; binomial values on an independent set; the convolution across an edgeless splitting |
| `MatchingBag.B1Zero.ReachIn`, `compIn`, `IsTransversal` (`Reach.lean`) | connectivity inside a finite vertex set, components, and transversals |
| `MatchingBag.B1Zero.not_reachIn_of_adj` (`Reach.lean`) | acyclicity: two distinct neighbours of `v` are not joined avoiding `v` |
| `MatchingBag.B1Zero.transversal_erase`, `transversal_children`, `exists_transversal_subset` (`Reach.lean`) | the transversals produced by peeling a component and then its root; a set meeting every component contains a transversal |
| `MatchingBag.B1Zero.CoronaData` (`Corona.lean`) | a base set with private pendant neighbours; `cor X = X ∪ pend '' X` |
| `MatchingBag.B1Zero.CoronaData.card_le_of_indep`, `cor_no_cross` (`Corona.lean`) | independent subsets of a corona have at most `\|X\|` elements; no edges across coronas over non-adjacent base sets |
| `MatchingBag.B1Zero.CoronaData.vec`, `vec_union`, `vec_peel`, `vec_peel_root` (`CoronaCount.lean`) | the top-five codimension profile of a corona and **the two recursions matching `RootedForestCertificate.profile`** |
| `MatchingBag.B1Zero.exists_forest_repr` (`CoronaRepr.lean`) | **graph-to-forest representation**: every corona over a transversal-rooted base is the `profile` of a plane rooted `Forest` term of the same size |
| `MatchingBag.B1Zero.CoronaData.cntRev_move_roots` (`CoronaMove.lean`) | **root moving / domination**: avoiding base roots instead of pendant roots only increases the counts |
| `DepthThree.RootedForestCertificate.rooted_forest_simple_bound` (`SimpleCertificate.lean`) | **the finite certificate** `4(73p₀+24p₁+3p₂+15q₀+6q₁+q₂) ≤ 126p₀+84p₁+36p₂+9p₃+p₄` for all 16796 plane rooted forests of size ten, by `decide +kernel` |

### The four-forbidden case and the target (B4, G0)

| Declaration | Statement |
| --- | --- |
| `MatchingBag.TreeMatching.flexLeaf`, `flexBase`, `flexCorona`, `cor_flexBaseSet` (`FlexCorona.lean`) | the corona structure carried by the flexible vertices; `cor base = Flex` |
| `MatchingBag.TreeMatching.exists_attach_reach` (`FlexCorona.lean`) | in a tree, every flexible component contains a vertex adjacent to a forbidden vertex |
| `MatchingBag.TreeMatching.reachIn_base_of_reachIn_cor`, `exists_root_transversal` (`FlexCorona.lean`) | reachability projects to the base; a transversal of base components consisting of attachment vertices up to pendants |
| `MatchingBag.TreeMatching.flex_certificate` (`FlexCorona.lean`) | the certificate transported to the graph counts of the flexible part |
| `MatchingBag.TreeMatching.cntRev_split`, `card_extSets_eq_convolution` (`FourStructure.lean`) | `e₄ = 126p₀ + 84p₁ + 36p₂ + 9p₃ + p₄` |
| `MatchingBag.TreeMatching.card_fibre_le`, `sum_flexNbrs_le`, `card_blkSets_le_union_bound` (`FourStructure.lean`) | fibring the blocked sets over their forbidden part, the union bound over the four forbidden vertices, and `b₄ ≤ 73p₀+24p₁+3p₂+15q₀+6q₁+q₂` |
| `MatchingBag.TreeMatching.density_four_aux` (`FourStructure.lean`) | `4 b₄ ≤ e₄` when `(α, \|Forced\|, \|Forbidden\|, #flexBags) = (19, 9, 4, 10)` |
| `MatchingBag.TreeMatching.density_four`, `density_window` (`Target.lean`) | the density bound for `r = 4`, and in the whole window |
| `DepthThree.b1_zero_window` (`Target.lean`) | **`DepthThree.B1ZeroWindowTarget`** |

### Notes on the `b₁ = 0` route

* The route is the shorter one suggested in `request.md`: only the graph density
  `3(α-3) b₄ ≤ (α-7) e₄` is proved in the window, and the already proved
  `DepthThree.low_density_of_forest` then gives `s₂ s₄ < s₃²`.
* For `r = 4` the parameters are forced to `(|V|, α, |Forced|, r, #flexBags) = (33, 19, 9, 4, 10)`
  by `card_verts_add_card_forced` together with `|N_C(U)| ≥ 2|U| + 1`.
* Instead of the attachment-concentration step of the written note, the blocked count is
  bounded by a plain union bound over the four forbidden vertices.  This yields the weaker
  coefficient vector `(73, 24, 3 ; 15, 6, 1)` in place of `(58, 21, 3 ; 30, 9, 1)`, which is
  still below the quarter threshold for every plane rooted forest of size ten; the fresh
  kernel-checked enumeration `rooted_forest_simple_bound` records this.  The supplied
  `RootedForestCertificate.rooted_forest_bound` and `PolynomialCertificate` are therefore not
  used by the final proof, and neither are `ParameterCertificate`, `FivefoldBridge` and the
  `Fivefold` helpers; all remain in the build unchanged.
* No profile is assumed: `exists_forest_repr` builds the `Forest` encoding of an arbitrary
  graph corona by peeling components and roots, and `cntRev_move_roots` supplies the
  domination needed when a component attaches at a pendant.
