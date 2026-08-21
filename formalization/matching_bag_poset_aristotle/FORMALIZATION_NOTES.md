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

* The reduction of §1 from a tree with a maximum matching to a forest poset (equation (2))
  is a construction sketched informally in the note; only the pigeonhole step above is
  formalized.
* The log-concavity statements (7) and (7a) of §3 are attributed there to the companion note
  `pascal_smoothing_shadow_lemma_2026-08-20.md`, which is not part of this project, and are
  not proved here.
* §4 (the blocked-path correction `b_d`) is descriptive and contains no self-contained claim
  to formalize; the note lists it as future work.
