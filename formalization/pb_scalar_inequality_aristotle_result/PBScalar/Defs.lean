import Mathlib

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
polynomial identities for `H ≥ 16`.

The definitions below may be adjusted for elaboration (coercions, `Finset`
idioms) provided the mathematics is unchanged; record any such change in the
README traceability table. The theorem statements must not be weakened.
-/

namespace PBScalar

open Finset

noncomputable section

/-- `a = (1 - 2δ) / (1 - δ)`. -/
def a (δ : ℝ) : ℝ := (1 - 2 * δ) / (1 - δ)

/-- Right lower bound `R r = a ^ r * ∏_{j=1}^{r-1} (1 - j δ)`. -/
def R (δ : ℝ) (r : ℕ) : ℝ :=
  a δ ^ r * ∏ j ∈ Finset.Ico 1 r, (1 - (j : ℝ) * δ)

/-- Left lower bound `L r = (1 - δ) ^ (-r) * ∏_{j=2}^{r+1} (1 - j δ)`. -/
def L (δ : ℝ) (r : ℕ) : ℝ :=
  ((1 - δ) ^ r)⁻¹ * ∏ j ∈ Finset.Icc 2 (r + 1), (1 - (j : ℝ) * δ)

/-- Weights on `ℤ`: `w 0 = 1`, `w r = R r` and `w (-r) = L r` for `r ≥ 1`. -/
def w (δ : ℝ) (i : ℤ) : ℝ :=
  if i = 0 then 1 else if 0 < i then R δ i.toNat else L δ (-i).toNat

/-- `A δ K = (1/2) * ∑_{i,j=-K}^{K} w i * w j * (i - j) ^ 2`. -/
def A (δ : ℝ) (K : ℕ) : ℝ :=
  (1 / 2 : ℝ) *
    ∑ i ∈ Finset.Icc (-(K : ℤ)) (K : ℤ), ∑ j ∈ Finset.Icc (-(K : ℤ)) (K : ℤ),
      w δ i * w δ j * ((i : ℝ) - (j : ℝ)) ^ 2

end

end PBScalar
