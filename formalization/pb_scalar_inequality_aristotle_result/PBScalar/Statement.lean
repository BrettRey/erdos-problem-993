import PBScalar.Compact
import PBScalar.Large

/-!
# Proposition 3.1 of "Variance and local log-concavity of Poisson–binomial laws"

For `0 < δ < 1/4` let `K` be the largest integer `r ≥ 1` with `(r + 1) δ < 1`,
characterised here by `(K + 1) δ < 1 ≤ (K + 2) δ`. Put

* `a = (1 - 2δ) / (1 - δ)`,
* `R r = a ^ r * ∏_{j=1}^{r-1} (1 - j δ)`,
* `L r = (1 - δ) ^ (-r) * ∏_{j=2}^{r+1} (1 - j δ)`,
* weights on `ℤ`: `w 0 = 1`, `w r = R r`, `w (-r) = L r` for `1 ≤ r ≤ K`,
* `A δ K = (1/2) * ∑_{i,j = -K}^{K} w i * w j * (i - j) ^ 2`.

Target (Proposition 3.1): `(3 + δ) / (4 δ ^ 2) ≤ A δ K`.

The paper's proof (see `PROOF_CONTEXT.md`) substitutes `H = 1/δ - 1 > 3`,
symmetrises the weights, and proves the resulting one-variable inequality by
exact Bernstein expansions on `3 < H ≤ 16` (data in `data/`) and by explicit
polynomial identities for `H ≥ 16`. The Lean proof follows the paper for `H ≥ 16`
(`PBScalar/Large.lean`); on `3 < H ≤ 16` it covers the same 13 cells by a
different exact method, monotonicity in `δ` plus exact rational endpoint checks
(`PBScalar/Compact.lean`).

The definitions (`a`, `R`, `L`, `w`, `A`) live verbatim in `PBScalar/Defs.lean`
(moved there unchanged so that the proof files can import them); the theorem
statements below are unchanged.
-/

namespace PBScalar

open Finset

/-- **G1** (compact range): `1/17 ≤ δ < 1/4`, i.e. `3 < H ≤ 16` with `H = 1/δ - 1`. -/
theorem scalar_inequality_compact (δ : ℝ) (K : ℕ)
    (hδ17 : 1 / 17 ≤ δ) (hδ4 : δ < 1 / 4)
    (hK : ((K : ℝ) + 1) * δ < 1) (hK' : 1 ≤ ((K : ℝ) + 2) * δ) :
    (3 + δ) / (4 * δ ^ 2) ≤ A δ K :=
  compact_main δ K hδ17 hδ4 hK hK'

/-- **G2** (large range): `0 < δ ≤ 1/17`, i.e. `H ≥ 16`. -/
theorem scalar_inequality_large (δ : ℝ) (K : ℕ)
    (hδ0 : 0 < δ) (hδ17 : δ ≤ 1 / 17)
    (hK : ((K : ℝ) + 1) * δ < 1) (hK' : 1 ≤ ((K : ℝ) + 2) * δ) :
    (3 + δ) / (4 * δ ^ 2) ≤ A δ K :=
  large_main δ K hδ0 hδ17 hK hK'

/-- **G0** (Proposition 3.1): for `0 < δ < 1/4`. -/
theorem scalar_inequality (δ : ℝ) (K : ℕ)
    (hδ0 : 0 < δ) (hδ4 : δ < 1 / 4)
    (hK : ((K : ℝ) + 1) * δ < 1) (hK' : 1 ≤ ((K : ℝ) + 2) * δ) :
    (3 + δ) / (4 * δ ^ 2) ≤ A δ K := by
  rcases le_or_gt δ (1 / 17) with h | h
  · exact scalar_inequality_large δ K hδ0 h hK hK'
  · exact scalar_inequality_compact δ K h.le hδ4 hK hK'

end PBScalar
