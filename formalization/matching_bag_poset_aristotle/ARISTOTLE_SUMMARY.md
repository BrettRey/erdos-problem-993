# Summary of changes for run fcbe1949-5ec5-4c84-8e73-a085b2352297
The task in `FOLLOWUP_TREE_TO_FOREST_POSET_20260820.md` is complete: the explicit tree/maximum-matching → forest-poset code bridge is formalized and consumes the existing `PascalBridge` results. The whole project builds with `lake build` and contains no `sorry`, `admit`, `axiom`, or `implemented_by`.

**New source files** (all under `RequestProject/`):
- `KonigHall.lean` — matchings, vertex covers, and **König's theorem derived from Mathlib's Hall theorem** (`konig`), not postulated.
- `TreeMatching.lean` — the `TreeMatching` structure (forest + bipartition + maximum matching), `exists_treeMatching` showing the orientation of matched edges is no loss of generality, the bags (`Idx`/`Unm`/`Bag`) with full cardinality bookkeeping (`card_Bag : #bags = |V| - |M|`), the constraint system (`rel`, `ForcedZero`, `ForcedOne`, `IsSol`, `Sol`), `minCover_structure`, and the explicit mutually inverse maps `coverOf`/`solOf` between minimum covers and solutions (`minCovers_image_solOf`).
- `ForestLemmas.lean` — the graph-theoretic core: contracting pairwise disjoint edges of a forest yields a forest (`bagGraph_isAcyclic`), and the directed comparison relation admits no closed chain (`no_closed_matching_chain`).
- `BagPoset.lean` — forced/free variables, the poset `P` on the free variables whose **antisymmetry is proved** from forest acyclicity, the solution ↔ order-ideal bijection (`Sol_image_restrict`), and forest structure of both the comparison graph (`compGraph_isAcyclic`) and the cover graph (`coverGraph_isAcyclic`).
- `TreeCodeBridge.lean` — the bag splitting `bagEquiv : Bag ≃ Free ⊕ Const`, **equation (2)** (`treeCode_eq_codeRelabel`), its independent-set form (`maxIndepCode_eq_flip_treeCode`), **equation (6)** coefficientwise and in polynomial form (`codeP_treeCode`, `codePoly_treeCode`: profile `= (1+t)^c · I(B(P);t)`), the identification of the tree erasure profile with the `PascalBridge` poset profile (`erasure_eq_erasureProfile`), and the four inequality consequences: depth-three `e₂e₄ ≤ e₃²`, the depth-three reserve `32(M-2)e₂e₄ ≤ 27(M-3)e₃²`, log-concavity through defect depth eight, and log-concavity at every interior defect for `M ≤ 33`.

All six required checkpoints are met. `#print axioms` reports exactly `[propext, Classical.choice, Quot.sound]` for every principal bridge theorem; nothing in the bridge depends on the pre-existing `native_decide` enumeration in `ExhaustiveSmallCodes.lean`.

Two honest notes on scope: connectedness of the tree is never used, so every statement is proved for arbitrary finite **forests**, which is strictly more general; and the "port-labelled quotient tree" of checkpoint 1 is realized as the contracted bag graph `SimpleGraph.fromRel` on matched bags, proved acyclic, together with the bag cardinality and cover/solution bijection results, rather than as a separately packaged labelled object.

`FORMALIZATION_NOTES.md` now carries a full checkpoint-by-checkpoint traceability table naming the exact Lean declaration and file for each requirement, plus the conventions and soundness record. All work is committed and pushed.

# Summary of changes for run 3dfa9612-a379-45a7-9af2-40b697b9aba4
I integrated the uploaded Pascal theorem into the library and closed the abstract poset extendable-profile inequalities end to end.

**What changed**

- `PascalSmoothing.lean` was moved into the library as `RequestProject/PascalSmoothing.lean` (unchanged content), so it is now a build target of the project.
- New file `RequestProject/PascalBridge.lean` connects it to the poset development of `RequestProject/Antichain.lean`.

**The bridge.** With `M = |P| + c`, the erasure profile `E d = Σ_j a_j · C(M-j, d)` is recorded as the natural-number sequence `MatchingBag.erasureProfile P c d`, and `erasureProfile_eq_coeff` proves `E d = [t^{M-d}]((1+t)^c · I(B(P);t))` for `d ≤ M` (so the definition of the note and the antichain sum agree). The Pascal development is phrased for subsets of `Fin M`, while `P` is an arbitrary finite type; rather than assuming `P = Fin r` or equality of the two profiles, the antichains of `P` are transported along an embedding `posetEmb P c : P ↪ Fin (|P| + c)`. The transported family `pascalFamily P c` is proved downward closed, nonempty, and face-number preserving (`faceCount_pascalFamily`), giving `erasureShadow (pascalFamily P c) d = erasureProfile P c d` — i.e. `E` is *precisely* the erasure shadow of a downward-closed family.

**Results proved (all six requested items), with no log-concavity assumed:**

1. `antichainsOf_downward_closed` — the antichains of `P` form a downward-closed family.
2. `erasureProfile_eq_sum` / `erasureProfile_eq_coeff` — `E d = Σ_j antichainCount P j · C(M-j, d)`, equal to the coefficient of `(1+X)^c · indepPoly P` at `M-d`.
3. `erasureProfile_pascal_smoothing` — the denominator-cleared inequality `8(d+1)(m+1)·E_{d-1}E_{d+1} ≤ 9dm·E_d²`, `m = M-d`, for `1 ≤ d ≤ M-1`.
4. `erasureProfile_log_concave_depth_le_eight` — log-concavity through defect depth eight, arbitrary `M`.
5. `erasureProfile_log_concave_of_M_le_33` — log-concavity at every interior defect when `M ≤ 33`.
6. `erasureProfile_depth_three` — `32(M-2)·E_2E_4 ≤ 27(M-3)·E_3²` for `M ≥ 4` (equation (7a)), plus `erasureProfile_depth_three_strict` (strict `E_2E_4 < E_3²` under the positivity hypothesis), `erasureProfile_depth_three_log_concave` (equation (7)) and its coefficient form `coeff_erasure_depth_three_log_concave`.

**Soundness.** `lake build` succeeds for the whole project; there is no `sorry`, `admit`, `axiom`, `implemented_by` or `native_decide` in the bridge or in the Pascal file, and no target was weakened. No definition mismatch arose — the displayed identities all hold as stated. `#print axioms` for `erasureShadow_pascalFamily`, `erasureProfile_eq_coeff`, `erasureProfile_pascal_smoothing`, `erasureProfile_log_concave_depth_le_eight`, `erasureProfile_log_concave_of_M_le_33`, `erasureProfile_depth_three`, `erasureProfile_depth_three_strict`, `erasureProfile_depth_three_log_concave`, `coeff_erasure_depth_three_log_concave` and `antichainsOf_downward_closed` each reports exactly `[propext, Classical.choice, Quot.sound]`; in particular nothing new depends on the documented `native_decide` computation in `ExhaustiveSmallCodes.lean`, which is untouched.

`FORMALIZATION_NOTES.md` was updated: the traceability table now lists the Pascal lemma, the bridge, and equations (7)/(7a) (previously listed as not formalized), and a new section explains the transport argument. All work is committed and pushed.

# Summary of changes for run 7353ac87-201d-4f7f-a697-817a6eff1e30
I formalized the note `matching_bag_poset_reduction_2026-08-20.md` in Lean 4 (Mathlib). Everything in `RequestProject/` builds with no `sorry` and no added axioms; a table of the correspondence with the note, plus the conventions used, is in `FORMALIZATION_NOTES.md`.

**Setup** (`RequestProject/Codes.lean`): a binary code is a `Finset (ι → Bool)`, projections `π_K` are recorded as partial words (`ι → Option Bool`), and the retained-coordinate profile is `codeP Ω k = ∑_{|K|=k} |π_K(Ω)|`.

**Section 2** (`PosetCode.lean`): the order-ideal code `Ω(P)` of a finite poset, the bipartite graph `B(P)`, and proofs of
- formula (3): `p_k(Ω(P)) = ∑_{|K|=k} |J(P[K])|` (`codeP_idealCode`), via the fact that a partial assignment extends iff it is an ideal of the induced subposet;
- formula (5): `p_k(Ω(P))` equals the number of size-`k` independent sets of `B(P)` (`codeP_eq_indepSets_card`), and its polynomial form `∑_k p_k t^k = I(B(P);t)` (`indepPoly_eq_sum_codeP`).

**Constant coordinate** (`ConstantCoordinate.lean`): adjoining a coordinate on which the code is constant multiplies the profile polynomial by `1+t` (`codeP_appendConst_zero`, `codeP_appendConst_succ`).

**Section 3** (`Antichain.lean`): formula (8) `I(B(P);t) = ∑_{O⊆P} t^{|O|}(1+t)^{r-|↓O|}`; formula (9) in the denominator-free form `I(B(P);t) = ∑_{A antichain} t^{|A|}(1+t)^{r-|A|}` and its collected version `∑_j a_j t^j (1+t)^{r-j}`; and formula (10) `E_d = [t^{M-d}]((1+t)^c I(B(P);t)) = ∑_j a_j C(M-j,d)` with `M = r + c`.

**Section 1** (`MatchingCover.lean`): the pigeonhole step used to define the variables `x_i` — a cover whose size equals that of a matching of pairwise disjoint edges contains exactly one endpoint of each matching edge and no other vertex.

**Section 5** (`UltraLogConcave.lean`, `ExhaustiveSmallCodes.lean`): the seven-bit code `{24,47,65,67,86,93,97,99}` has profile `(1,14,80,187,220,145,52,8)`; its profile is ordinary log-concave, but normalized (ultra) log-concavity fails at `k = 6` (`52²·C(7,5)·C(7,7) = 56784 < 56840 = 145·8·C(7,6)²`). I also verified exhaustively that ordinary and normalized log-concavity hold for every nonempty binary code of length at most four (lengths 1–3 by kernel evaluation; length 4 enumerates all 2^16−1 codes and is discharged by compiled evaluation, so it additionally uses the `Lean.ofReduceBool` axiom).

Not formalized, as recorded in the notes file: the informal §1 construction of the forest poset from a tree with a maximum matching (equation (2)); the log-concavity inequalities (7)/(7a), whose proof the note attributes to a companion document that is not part of this project; and §4, which is descriptive and lists future work.