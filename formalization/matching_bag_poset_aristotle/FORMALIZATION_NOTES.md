# Formalization of `matching_bag_poset_reduction_2026-08-20.md`

All Lean sources live in `RequestProject/` and build without `sorry`.

## What is formalized

| Note | Statement | Lean |
| --- | --- | --- |
| §1 | a cover of size `q` meeting a matching of `q` disjoint edges contains exactly one endpoint of each edge and nothing else | `MatchingBag.cover_card_eq_matching_card` (`MatchingCover.lean`) |
| §2 | `p_k(Ω) = ∑_{|K|=k} |π_K(Ω)|` | `MatchingBag.codeP` (`Codes.lean`) |
| §2 (3) | `p_k(Ω(P)) = ∑_{|K|=k} |J(P[K])|` | `MatchingBag.codeP_idealCode` (`PosetCode.lean`) |
| §2 (4) | the bipartite graph `B(P)` | `MatchingBag.posetBip` (`PosetCode.lean`) |
| §2 (5) | `∑_k p_k(P) t^k = I(B(P);t)` | `MatchingBag.codeP_eq_indepSets_card`, `MatchingBag.indepPoly_eq_sum_codeP` |
| §2 (6) | a constant coordinate multiplies the profile polynomial by `1+t` | `MatchingBag.codeP_appendConst_zero`, `MatchingBag.codeP_appendConst_succ` (`ConstantCoordinate.lean`) |
| §3 (8) | `I(B(P);t) = ∑_{O ⊆ P} t^{|O|}(1+t)^{r-|↓O|}` | `MatchingBag.indepPoly_eq_sum_downClosure` (`Antichain.lean`) |
| §3 (9) | `I(B(P);t) = (1+t)^r A_P(t/(1+t))` | `MatchingBag.indepPoly_eq_sum_antichains`, `MatchingBag.indepPoly_eq_antichainCount` |
| §3 (10) | `E_d = ∑_j a_j C(M-j, d)` | `MatchingBag.coeff_erasure` |
| §5 | the seven-bit code has profile `(1,14,80,187,220,145,52,8)`, is log-concave but not ultra-log-concave | `MatchingBag.counterCode_profile`, `MatchingBag.counterCode_logConcave`, `MatchingBag.counterCode_not_ultraLogConcave` (`UltraLogConcave.lean`) |
| §5 | ordinary and normalized log-concavity for every nonempty binary code of length ≤ 4 | `MatchingBag.logConcave_length_one/two/three/four` (`ExhaustiveSmallCodes.lean`) |
| companion note | universal Pascal-smoothing lemma for erasure shadows of downward-closed families on `Fin M` | `PascalSmoothing.pascal_smoothing_shadow_lemma` and its corollaries (`PascalSmoothing.lean`) |
| §3 | the antichains of `P` form a downward-closed family | `MatchingBag.antichainsOf_downward_closed` (`PascalBridge.lean`) |
| §3 (10) | `E_d = ∑_j a_j C(M-j,d)` over `ℕ`, and `E_d = [t^{M-d}]((1+t)^c I(B(P);t))` | `MatchingBag.erasureProfile`, `MatchingBag.erasureProfile_eq_sum`, `MatchingBag.erasureProfile_eq_coeff` (`PascalBridge.lean`) |
| §3 | `E` is exactly the erasure shadow of a downward-closed family on `M` coordinates | `MatchingBag.erasureShadow_pascalFamily` (`PascalBridge.lean`) |
| §3 (7) | `E_3^2 ≥ E_2 E_4` for `M ≥ 4` (also in coefficient form, and strict when `E_2E_4 > 0`) | `MatchingBag.erasureProfile_depth_three_log_concave`, `MatchingBag.coeff_erasure_depth_three_log_concave`, `MatchingBag.erasureProfile_depth_three_strict` (`PascalBridge.lean`) |
| §3 (7a) | the reserve `27(M-3) E_3^2 ≥ 32(M-2) E_2 E_4` | `MatchingBag.erasureProfile_depth_three` (`PascalBridge.lean`) |
| §3 | Pascal-smoothing inequality for adjacent `E` values, log-concavity through defect depth eight, and at every interior defect for `M ≤ 33` | `MatchingBag.erasureProfile_pascal_smoothing`, `MatchingBag.erasureProfile_log_concave_depth_le_eight`, `MatchingBag.erasureProfile_log_concave_of_M_le_33` (`PascalBridge.lean`) |

## Conventions

* A binary code is a `Finset (ι → Bool)` on a finite index type; a projection is recorded as
  a partial word `ι → Option Bool` (`none` = erased coordinate).
* `Ω(P)` is the set of indicator functions of order *ideals* (down-closed subsets) of `P`.
  Accordingly, `B(P)` is oriented so that `i^1 j^0` is an edge exactly when `j ≤ i`, i.e.
  exactly when the partial assignment `x_i = 1, x_j = 0` violates down-closedness. This is the
  graph of formula (4) up to the direction of the order convention; independent sets of
  `B(P)` are precisely the extendable partial assignments, which is what the counts use.
* Formula (9) is stated without denominators as
  `I(B(P);t) = ∑_{A antichain} t^{|A|}(1+t)^{r-|A|}`, which is `(1+t)^r A_P(t/(1+t))`
  cleared of the `1+t` denominators.

## What is not formalized

* §4 (the blocked-path correction `b_d`) is descriptive and contains no self-contained claim
  to formalize; the note lists it as future work.

(The reduction of §1 from a tree with a maximum matching to a forest poset — equation (2) —
was sketched informally in the note; it is now fully formalized, see the traceability table
below.)

## The Pascal bridge (`PascalSmoothing.lean`, `PascalBridge.lean`)

`PascalSmoothing.lean` is the formalization of the companion note
`pascal_smoothing_shadow_lemma_2026-08-20.md`: for a downward-closed family `Δ` of subsets
of `Fin M` with face numbers `a_j`, the erasure shadow `E d = ∑_j a_j C(M-j,d)` satisfies
the denominator-cleared inequality
`8 (d+1)(m+1) E_{d-1} E_{d+1} ≤ 9 d m E_d^2` with `m = M - d`, together with its
consequences (log-concavity through defect depth eight, log-concavity for `M ≤ 33`, and the
depth-three reserve).  Its proof runs through the normalized face densities
`b_j = a_j / C(M,j)`, which are nonincreasing by local LYM, and the absorption identity
`E_d = C(M,d) · S_{M-d}` for the Pascal transform `S`.

`PascalBridge.lean` connects this to the poset side.  The Pascal development is phrased for
subsets of `Fin M`, whereas `P` is an arbitrary finite type, so the antichains of `P` are
*transported* along an embedding `posetEmb P c : P ↪ Fin (|P| + c)` (which exists because
`|P| ≤ M`).  Nothing is assumed about the two profiles: the transported family
`pascalFamily P c` is proved downward closed (`pascalFamily_downward_closed`), nonempty, and
face-number preserving (`faceCount_pascalFamily`), whence
`erasureShadow (pascalFamily P c) d = erasureProfile P c d` (`erasureShadow_pascalFamily`).
Transporting the Pascal inequalities across this identity gives equations (7) and (7a) for
every finite poset `P` and every `c`.

All theorems in these two files depend only on `propext`, `Classical.choice` and `Quot.sound`
(checked with `#print axioms`); in particular nothing here depends on the `native_decide`
computation in `ExhaustiveSmallCodes.lean`.

## The tree → forest-poset bridge (`FOLLOWUP_TREE_TO_FOREST_POSET_20260820.md`)

Source files, in dependency order:

* `KonigHall.lean` — matchings, vertex covers, König's theorem for bipartite graphs.
* `TreeMatching.lean` — the `TreeMatching` structure (forest + bipartition + maximum
  matching), the bags, the constraint system, and the cover ↔ solution bijection.
* `ForestLemmas.lean` — the graph-theoretic core: contracting a matching of a forest gives a
  forest, and the directed comparison relation has no closed chain.
* `BagPoset.lean` — forced/free variables, the poset `P` on the free variables, the
  solution ↔ order-ideal bijection, and acyclicity of the comparison and cover graphs.
* `TreeCodeBridge.lean` — equation (2), equation (6), and the inequality corollaries.

### Checkpoint traceability

| Checkpoint | Statement | Lean declaration (file) |
| --- | --- | --- |
| 1. Matching bags | a forest with a maximum matching yields an oriented `TreeMatching` | `MatchingBag.exists_treeMatching` (`TreeMatching.lean`) |
| 1 | the bag index type: matched edges plus unmatched singletons | `MatchingBag.TreeMatching.Idx`, `.Unm`, `.Bag` (`TreeMatching.lean`) |
| 1 | bags are pairwise disjoint; every vertex lies in exactly one bag | `MatchingBag.TreeMatching.bags_disjoint`, `.vertex_cases` (`TreeMatching.lean`) |
| 1 | an independent set contains at most one endpoint of a matched bag | `MatchingBag.TreeMatching.indep_not_both` (`TreeMatching.lean`) |
| 1 | cardinality bookkeeping: `#bags = |V| - |M|`, `#matched = 2|M|` | `MatchingBag.TreeMatching.card_Bag`, `.card_matchedVerts`, `.card_unmatchedVerts`, `.card_Idx` (`TreeMatching.lean`) |
| 1 | the feasible bag-assignment bijection: minimum covers ↔ solutions of the constraint system | `MatchingBag.TreeMatching.coverOf`, `.solOf`, `.coverOf_mem_minCovers`, `.solOf_mem_Sol`, `.solOf_coverOf`, `.coverOf_solOf`, `.minCovers_image_solOf` (`TreeMatching.lean`) |
| 2. Maximum layer (König) | König's theorem, proved from Mathlib's Hall theorem — not postulated | `MatchingBag.konig` (`KonigHall.lean`), with `MatchingBag.card_matching_le_cover`, `MatchingBag.hall_side` |
| 2 | a minimum cover has `|M|` vertices, hence exactly one endpoint of every matched edge and no unmatched vertex | `MatchingBag.TreeMatching.minCover_structure` (`TreeMatching.lean`), reusing `MatchingBag.cover_card_eq_matching_card` (`MatchingCover.lean`) |
| 2 | maximum independent sets are exactly the complements of minimum covers | `MatchingBag.TreeMatching.maxIndepSets_eq_image_compl` (`TreeMatching.lean`) |
| 3. Inequality system | equation (1): a non-matching edge `R_i — L_j` gives `x_i ≤ x_j` | `MatchingBag.TreeMatching.rel`, `.minCover_rel` (`TreeMatching.lean`) |
| 3 | edges to unmatched vertices give unary forced values | `MatchingBag.TreeMatching.ForcedZero`, `.ForcedOne`, `.minCover_forcedZero`, `.minCover_forcedOne` (`TreeMatching.lean`) |
| 3 | the constraint system and its solution set | `MatchingBag.TreeMatching.IsSol`, `.Sol` (`TreeMatching.lean`) |
| 4. Consistency | the unary values are consistent (`Sol ≠ ∅`), derived from maximality of `M` through König | `MatchingBag.TreeMatching.Sol_nonempty` (`TreeMatching.lean`) |
| 4. Forest structure | contracting pairwise disjoint edges of a forest gives a forest | `MatchingBag.bagGraph_isAcyclic` (`ForestLemmas.lean`) |
| 4 | the directed comparison relation has no closed chain (antisymmetry input) | `MatchingBag.no_closed_matching_chain` (`ForestLemmas.lean`) |
| 4 | the undirected comparison graph on matched bags is a forest | `MatchingBag.TreeMatching.compGraph_isAcyclic` (`BagPoset.lean`) |
| 4 | propagating and deleting forced variables; the poset `P` on the remaining free variables, with antisymmetry *proved* | `MatchingBag.TreeMatching.Forced`, `.forcedVal`, `.Free`, `instPartialOrderFree` (`BagPoset.lean`) |
| 4 | the cover graph of `P` is a forest | `MatchingBag.TreeMatching.coverGraph_isAcyclic` (`BagPoset.lean`), via `.rel_of_covBy` |
| 5. Code equivalence | mutually inverse maps between solutions and order-ideal indicators of `P` | `MatchingBag.TreeMatching.restrictFree`, `.extendFree`, `.isIdealIndicator_restrictFree`, `.extendFree_mem_Sol`, `.restrictFree_extendFree`, `.extendFree_restrictFree`, `.Sol_image_restrict`, `.restrictFree_injOn` (`BagPoset.lean`) |
| 5 | the bags split as free variables ⊕ constant coordinates | `MatchingBag.TreeMatching.Const`, `.numConst`, `.bagEquiv`, `.card_Bag_eq_card_Free_add_numConst` (`TreeCodeBridge.lean`) |
| 5 | **equation (2)**: `treeCode = codeRelabel bagEquiv bagFlip (appendConsts (idealCode P) Const)` | `MatchingBag.TreeMatching.treeCode_eq_codeRelabel` (`TreeCodeBridge.lean`) |
| 5 | the code of maximum independent sets is the bitwise complement of the cover code, with the same profile | `MatchingBag.TreeMatching.maxIndepCode_eq_flip_treeCode`, `.codeP_maxIndepCode` (`TreeCodeBridge.lean`) |
| 6. Polynomial consequence | **equation (6)**, coefficientwise and as polynomials: profile `= (1+t)^c I(B(P);t)` | `MatchingBag.TreeMatching.codeP_treeCode`, `.codePoly_treeCode` (`TreeCodeBridge.lean`) |
| 6 | the tree erasure profile equals the poset erasure profile of `PascalBridge.lean` | `MatchingBag.TreeMatching.erasure`, `.erasure_eq_erasureProfile` (`TreeCodeBridge.lean`) |
| 6 | depth-three log-concavity `e_2 e_4 ≤ e_3^2` | `MatchingBag.TreeMatching.erasure_depth_three` (`TreeCodeBridge.lean`) |
| 6 | depth-three reserve `32(M-2) e_2 e_4 ≤ 27(M-3) e_3^2` | `MatchingBag.TreeMatching.erasure_depth_three_reserve` (`TreeCodeBridge.lean`) |
| 6 | log-concavity through defect depth eight | `MatchingBag.TreeMatching.erasure_log_concave_depth_le_eight` (`TreeCodeBridge.lean`) |
| 6 | log-concavity at every interior defect when `M ≤ 33` | `MatchingBag.TreeMatching.erasure_log_concave_of_le_33` (`TreeCodeBridge.lean`) |

### Conventions used by the bridge

* A `TreeMatching V` packages a forest `G` on a finite vertex type, a proper Boolean
  colouring `col`, and a maximum matching `M` given as a `Finset (V × V)` whose pairs are
  oriented from the `true` side to the `false` side. `exists_treeMatching` shows this is no
  loss of generality: any forest with any maximum matching can be so oriented, because a
  forest is bipartite (`SimpleGraph.IsAcyclic.coloringTwo`) and reorienting pairs does not
  change the cardinality of the matching. Connectedness of `T` is never used, so all
  statements are proved for forests, which is strictly more general than trees.
* `M` in the erasure statements denotes `Fintype.card D.Bag = |V| - |M|`, the independence
  number of the forest, matching the note's use of `M` for the code length.
* `bagFlip` records the coordinate complementations needed to turn the *cover* code into an
  order-ideal code: a bag forced to `0` and every unmatched singleton bag is flipped.
  `codeRelabel` (`CodeInvariance.lean`) is proved to preserve the profile `codeP`, so the
  flips are harmless for every statement in §2, §3.
* The poset `P = D.Free` carries the order generated by the comparison relation `rel`,
  reflexive-transitively closed. Antisymmetry is *derived* from
  `no_closed_matching_chain`, i.e. from acyclicity of the forest, rather than assumed.

### Soundness

`lake build` succeeds for every module. `RequestProject/` contains no `sorry`, `admit`,
`axiom`, or `implemented_by`. `#print axioms` reports exactly
`[propext, Classical.choice, Quot.sound]` for every principal bridge theorem, namely
`treeCode_eq_codeRelabel`, `codeP_treeCode`, `codePoly_treeCode`,
`maxIndepCode_eq_flip_treeCode`, `codeP_maxIndepCode`, `erasure_eq_erasureProfile`,
`erasure_depth_three`, `erasure_depth_three_reserve`,
`erasure_log_concave_depth_le_eight`, `erasure_log_concave_of_le_33`,
`coverGraph_isAcyclic`, `compGraph_isAcyclic`, `Sol_image_restrict`, `konig`,
`bagGraph_isAcyclic`, `no_closed_matching_chain`, `exists_treeMatching`, `card_Bag`.
In particular no bridge theorem depends on `Lean.ofReduceBool`: the only `native_decide`
in the project is the pre-existing, self-contained enumeration in
`ExhaustiveSmallCodes.lean`, which nothing above imports.
