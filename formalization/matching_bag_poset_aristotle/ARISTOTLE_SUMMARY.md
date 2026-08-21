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