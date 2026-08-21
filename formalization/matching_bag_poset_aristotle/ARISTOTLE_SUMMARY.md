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