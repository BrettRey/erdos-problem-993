# Summary of changes for run dc1051b8-55b2-473f-9f3d-62a5f5e55e00
## Result

The requested uniform Laurent-polynomial theorem is **true and now fully proved in Lean 4 / Mathlib**, for all natural numbers `r, s ≥ 1` and all rationals `c, d ≥ 1` — no bounded replacement, no weakening, and no `sorry`, `admit`, `axiom`, `implemented_by` or `native_decide` anywhere in the project. `lake build` succeeds, and `#print axioms EvenConnector.Gpoly_centrallyUnimodal` reports only `[propext, Classical.choice, Quot.sound]`.

## What was built

Three source-labelled files plus a `README.md` containing the traceability table and a prose account of the proof.

* `RequestProject/Coeff.lean` — coefficient calculus for `ℚ[T;T⁻¹]` (Mathlib has no coefficient API for `LaurentPolynomial`, so `lcoeff` and its rules were built from scratch), binomial coefficients with an integer lower index (`czz`, with symmetry, Pascal and row-monotonicity), the definitions `XL = z + z⁻¹`, `HL n = zⁿ + z⁻ⁿ`, and the coefficients of `X ^ n`.
* `RequestProject/Main.lean` — the block `Gpoly r s c d`; the requested exact identity `c·d·A + d·B + c·C + D = (c-1)(d-1)A + (c-1)(A+C) + (d-1)(A+B) + (A+B+C+D)` (`Gpoly_decomposition`); the full coefficient computation (`lcoeff_Gpoly`, `lcoeff_Gpoly_odd`); the definition `CentrallyUnimodal` (symmetry under `m ↦ -m`, nonnegativity, and `coeff f (m+2) ≤ coeff f m` for `m ≥ 0`); the main theorem `Gpoly_centrallyUnimodal` and its corollary `Gpoly_coeff_sub_nonneg`.
* `RequestProject/Binomial.lean` — the combinatorial core.

## Audit of the suggested proof, and the repair

The suggested reduction is correct, and the unit-constant case is where the difficulty lies. Writing `N = r+s+1` and substituting the reflected index `k = N-1-j`, the whole coefficient inequality collapses to two statements about `k` with `2k+2 ≤ N`:

* `C(N,k) + C(r,k) ≤ C(N,k+1) + C(r,k+1)` (and its mirror in `s`), corresponding to the summands `A+C` and `A+B`;
* `C(N,k) + C(r,k) + C(s,k) + [k=0] ≤ C(N,k+1) + C(r,k+1) + C(s,k+1)`, corresponding to `A+B+C+D`.

The `[k=0]` term is the contribution of `H(r+s+1)`; this boundary case is exactly the one that the `j<s` / `j=s` / `s<j<N/2` split handles awkwardly, and it needs `r, s ≥ 1`. Rather than repairing that split case by case, the proof here replaces it by a uniform Catalan estimate valid in every range and every parity:

1. `C(n,k) - C(n,k+1) ≤ catalan k` for **all** `n, k` (the maximum over `n` is at `n = 2k`);
2. `catalan (k+1) ≤ C(n,k+1) - C(n,k)` whenever `2k+2 ≤ n` (the ascent is monotone in `n`);
3. `2·catalan k ≤ catalan (k+1)` for `k ≥ 1`.

Chaining these gives `C(N,k+1) - C(N,k) ≥ catalan(k+1) ≥ 2·catalan k ≥ (C(r,k)-C(r,k+1)) + (C(s,k)-C(s,k+1))`, which is the required inequality for `k ≥ 1`; `k = 0` is checked directly.

The formal definitions were cross-checked against independent hand computations of the coefficient sequences (e.g. `G(1,1,1,1)` gives `4,5,5,4` and `G(1,3,1,1)` gives `4,9,14,14,9,4`), and the coefficient extraction was cross-checked against a direct expansion of `X³`.
