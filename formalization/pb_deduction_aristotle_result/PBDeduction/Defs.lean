import Mathlib

/-!
# Definitions for Sections 2–3 (moved verbatim from `PBDeduction/Statement.lean`)

These definitions are copied character-for-character from the original
`PBDeduction/Statement.lean` so that the proof files can import them. The
mathematics is unchanged.
-/

namespace PBDeduction

open Finset

noncomputable section

/-! ## Abstract mass functions -/

/-- Normalized Turán deficit `δ_k = 1 - f (k-1) f (k+1) / f k ^ 2`. -/
def deficit (f : ℤ → ℝ) (k : ℤ) : ℝ :=
  1 - f (k - 1) * f (k + 1) / f k ^ 2

/-- Pairwise form `(1/2) ∑_{i,j=0}^{n} f i f j (i - j)^2` of the variance of a
mass function on `{0, …, n}` (it equals the variance when `∑ f = 1`). -/
def pairVar (f : ℤ → ℝ) (n : ℕ) : ℝ :=
  (1 / 2 : ℝ) * ∑ i ∈ Finset.Icc (0 : ℤ) n, ∑ j ∈ Finset.Icc (0 : ℤ) n,
    f i * f j * ((i : ℝ) - (j : ℝ)) ^ 2

/-- Maximal mass `M = max_{0 ≤ k ≤ n} f k`. -/
def maxMass (f : ℤ → ℝ) (n : ℕ) : ℝ :=
  (Finset.Icc (0 : ℤ) n).sup' ⟨0, by simp⟩ f

/-- `D` is the first index in `{1, …, n}` at which `f` decreases. -/
def IsFirstDescent (f : ℤ → ℝ) (n : ℕ) (D : ℤ) : Prop :=
  1 ≤ D ∧ D ≤ n ∧ f D < f (D - 1) ∧ ∀ k : ℤ, 1 ≤ k → k < D → f (k - 1) ≤ f k

/-- The properties of a Poisson–binomial mass function that Sections 2–3 use:
positivity on the support, zero outside it, total mass one, strict
log-concavity (Newton's inequalities), and the Hillion–Johnson cubic
inequalities (their (78)–(79); the paper's (2.2)–(2.3)). -/
structure DeductionHyp (f : ℤ → ℝ) (n : ℕ) : Prop where
  pos : ∀ k : ℤ, 0 ≤ k → k ≤ n → 0 < f k
  zero_out : ∀ k : ℤ, k < 0 ∨ (n : ℤ) < k → f k = 0
  sum_one : ∑ k ∈ Finset.Icc (0 : ℤ) n, f k = 1
  strict_lc : ∀ k : ℤ, 0 ≤ k → k ≤ n → f (k - 1) * f (k + 1) < f k ^ 2
  hj_left : ∀ k : ℤ,
    (f (k - 1) ^ 2 - f (k - 2) * f k) * f (k + 1) ≤
      f (k - 1) * (f k ^ 2 - f (k - 1) * f (k + 1))
  hj_right : ∀ k : ℤ,
    (f (k + 1) ^ 2 - f k * f (k + 2)) * f (k - 1) ≤
      f (k + 1) * (f k ^ 2 - f (k - 1) * f (k + 1))

/-! ## Poisson–binomial laws -/

/-- Probability-generating polynomial `∏ i (1 - p i + p i X)`. -/
def pgf (n : ℕ) (p : Fin n → ℝ) : Polynomial ℝ :=
  ∏ i, (Polynomial.C (1 - p i) + Polynomial.C (p i) * Polynomial.X)

/-- The Poisson–binomial mass function `k ↦ P(W = k)`, as a function on `ℤ`. -/
def pbPmf (n : ℕ) (p : Fin n → ℝ) (k : ℤ) : ℝ :=
  if 0 ≤ k then (pgf n p).coeff k.toNat else 0

/-- The variance `V = ∑ i p i (1 - p i)`. -/
def pbVar (n : ℕ) (p : Fin n → ℝ) : ℝ :=
  ∑ i, p i * (1 - p i)

end

end PBDeduction
