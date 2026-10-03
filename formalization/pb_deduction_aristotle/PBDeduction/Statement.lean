import PBScalar.Statement

/-!
# Sections 2–3 of "Variance and local log-concavity of Poisson–binomial laws"

Theorem 1.1 of the paper: for a Poisson–binomial law with success
probabilities `0 < p i < 1`, variance `V = ∑ p i (1 - p i) ≥ 1`, mass function
`f` and first descent `D`, we have `V * δ_D ≥ 1/4`, where
`δ_k = 1 - f (k-1) f (k+1) / f k ^ 2`.

Sections 2–3 deduce this from Proposition 3.1 (the scalar inequality), which
is already formalized and kernel-checked in `PBScalar` as
`PBScalar.scalar_inequality`. **`PBScalar/` is a verified dependency: do not
modify it.**

Mass functions are functions `ℤ → ℝ`, zero outside `{0, …, n}`, so that the
paper's conventions `f (-1) = f (n+1) = 0` and "for every `k ∈ ℤ`" are literal.

Graded targets (see `PROMPT.md`):
* **G1** `deduction`: the deduction for an abstract mass function satisfying
  the hypotheses that the paper takes from the literature.
* **G2** `max_mass_bound`: the maximal-mass bound of Bobkov, Marsiglietti and
  Melbourne for any mass function on `{0, …, n}`.
* **G3** `pb_*` basics: positivity, support, total mass, variance, strict
  log-concavity, existence of the first descent.
* **G4** `pb_hj_left`, `pb_hj_right`: the Hillion–Johnson cubic inequalities.
* **G0** `theorem_1_1`: Theorem 1.1 for Poisson–binomial laws.
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

/-- **G1** (Sections 2–3, conditional): the deduction of Theorem 1.1 from
Proposition 3.1, the maximal-mass bound and the properties in `DeductionHyp`. -/
theorem deduction (f : ℤ → ℝ) (n : ℕ) (D : ℤ) (h : DeductionHyp f n)
    (hD : IsFirstDescent f n D) (hV : 1 ≤ pairVar f n)
    (hmax : (1 / maxMass f n ^ 2 - 1) / 12 ≤ pairVar f n) :
    1 / 4 ≤ pairVar f n * deficit f D := by
  sorry

/-- **G2** (Bobkov–Marsiglietti–Melbourne, Corollary 3.2; proof in Section 3):
`Var ≥ (M^{-2} - 1)/12` for every mass function on `{0, …, n}`. -/
theorem max_mass_bound (f : ℤ → ℝ) (n : ℕ) (hnn : ∀ k, 0 ≤ f k)
    (hsum : ∑ k ∈ Finset.Icc (0 : ℤ) n, f k = 1) :
    (1 / maxMass f n ^ 2 - 1) / 12 ≤ pairVar f n := by
  sorry

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

/-- **G3a** -/
theorem pb_pos (n : ℕ) (p : Fin n → ℝ) (hp : ∀ i, 0 < p i ∧ p i < 1) :
    ∀ k : ℤ, 0 ≤ k → k ≤ n → 0 < pbPmf n p k := by
  sorry

/-- **G3b** -/
theorem pb_zero_out (n : ℕ) (p : Fin n → ℝ) :
    ∀ k : ℤ, k < 0 ∨ (n : ℤ) < k → pbPmf n p k = 0 := by
  sorry

/-- **G3c** -/
theorem pb_sum_one (n : ℕ) (p : Fin n → ℝ) :
    ∑ k ∈ Finset.Icc (0 : ℤ) n, pbPmf n p k = 1 := by
  sorry

/-- **G3d** -/
theorem pb_pairVar (n : ℕ) (p : Fin n → ℝ) :
    pairVar (pbPmf n p) n = pbVar n p := by
  sorry

/-- **G3e** (Newton's inequalities, strict form). -/
theorem pb_strict_lc (n : ℕ) (p : Fin n → ℝ) (hp : ∀ i, 0 < p i ∧ p i < 1) :
    ∀ k : ℤ, 0 ≤ k → k ≤ n →
      pbPmf n p (k - 1) * pbPmf n p (k + 1) < pbPmf n p k ^ 2 := by
  sorry

/-- **G3f** -/
theorem pb_first_descent_exists (n : ℕ) (p : Fin n → ℝ)
    (hp : ∀ i, 0 < p i ∧ p i < 1) (hV : 1 ≤ pbVar n p) :
    ∃ D : ℤ, IsFirstDescent (pbPmf n p) n D := by
  sorry

/-- **G4a** (Hillion–Johnson, Theorem A.2, their (78)). -/
theorem pb_hj_left (n : ℕ) (p : Fin n → ℝ) (hp : ∀ i, 0 < p i ∧ p i < 1) :
    ∀ k : ℤ,
      (pbPmf n p (k - 1) ^ 2 - pbPmf n p (k - 2) * pbPmf n p k) * pbPmf n p (k + 1) ≤
        pbPmf n p (k - 1) *
          (pbPmf n p k ^ 2 - pbPmf n p (k - 1) * pbPmf n p (k + 1)) := by
  sorry

/-- **G4b** (Hillion–Johnson, Corollary A.3, their (79)). -/
theorem pb_hj_right (n : ℕ) (p : Fin n → ℝ) (hp : ∀ i, 0 < p i ∧ p i < 1) :
    ∀ k : ℤ,
      (pbPmf n p (k + 1) ^ 2 - pbPmf n p k * pbPmf n p (k + 2)) * pbPmf n p (k - 1) ≤
        pbPmf n p (k + 1) *
          (pbPmf n p k ^ 2 - pbPmf n p (k - 1) * pbPmf n p (k + 1)) := by
  sorry

/-- **G0** (Theorem 1.1). -/
theorem theorem_1_1 (n : ℕ) (p : Fin n → ℝ) (hp : ∀ i, 0 < p i ∧ p i < 1)
    (hV : 1 ≤ pbVar n p) (D : ℤ) (hD : IsFirstDescent (pbPmf n p) n D) :
    1 / 4 ≤ pbVar n p * deficit (pbPmf n p) D := by
  sorry

end

end PBDeduction
