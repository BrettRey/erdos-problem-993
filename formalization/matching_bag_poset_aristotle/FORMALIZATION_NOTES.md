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

* The reduction of §1 from a tree with a maximum matching to a forest poset (equation (2))
  is a construction sketched informally in the note; only the pigeonhole step above is
  formalized.
* §4 (the blocked-path correction `b_d`) is descriptive and contains no self-contained claim
  to formalize; the note lists it as future work.

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
