import PBScalar.Statement
import PBDeduction.Defs
import PBDeduction.Deduction
import PBDeduction.MaxMass
import PBDeduction.PBHillionJohnson

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

/-! The definitions `deficit`, `pairVar`, `maxMass`, `IsFirstDescent`,
`DeductionHyp`, `pgf`, `pbPmf`, `pbVar` were moved verbatim (character for
character) to `PBDeduction/Defs.lean`, so that the proof files can import them.
See `README.md`. -/

/-- **G1** (Sections 2–3, conditional): the deduction of Theorem 1.1 from
Proposition 3.1, the maximal-mass bound and the properties in `DeductionHyp`. -/
theorem deduction (f : ℤ → ℝ) (n : ℕ) (D : ℤ) (h : DeductionHyp f n)
    (hD : IsFirstDescent f n D) (hV : 1 ≤ pairVar f n)
    (hmax : (1 / maxMass f n ^ 2 - 1) / 12 ≤ pairVar f n) :
    1 / 4 ≤ pairVar f n * deficit f D :=
  deduction_main h hD hV hmax

/-- **G2** (Bobkov–Marsiglietti–Melbourne, Corollary 3.2; proof in Section 3):
`Var ≥ (M^{-2} - 1)/12` for every mass function on `{0, …, n}`. -/
theorem max_mass_bound (f : ℤ → ℝ) (n : ℕ) (hnn : ∀ k, 0 ≤ f k)
    (hsum : ∑ k ∈ Finset.Icc (0 : ℤ) n, f k = 1) :
    (1 / maxMass f n ^ 2 - 1) / 12 ≤ pairVar f n :=
  max_mass_bound_proof f n hnn hsum

/-- **G3a** -/
theorem pb_pos (n : ℕ) (p : Fin n → ℝ) (hp : ∀ i, 0 < p i ∧ p i < 1) :
    ∀ k : ℤ, 0 ≤ k → k ≤ n → 0 < pbPmf n p k :=
  (pbPmf_nonneg_pos n p hp).2

/-- **G3b** -/
theorem pb_zero_out (n : ℕ) (p : Fin n → ℝ) :
    ∀ k : ℤ, k < 0 ∨ (n : ℤ) < k → pbPmf n p k = 0 :=
  pbPmf_zero_out n p

/-- **G3c** -/
theorem pb_sum_one (n : ℕ) (p : Fin n → ℝ) :
    ∑ k ∈ Finset.Icc (0 : ℤ) n, pbPmf n p k = 1 :=
  pbPmf_sum_one n p

/-- **G3d** -/
theorem pb_pairVar (n : ℕ) (p : Fin n → ℝ) :
    pairVar (pbPmf n p) n = pbVar n p :=
  pbPmf_pairVar n p

/-- **G3e** (Newton's inequalities, strict form). -/
theorem pb_strict_lc (n : ℕ) (p : Fin n → ℝ) (hp : ∀ i, 0 < p i ∧ p i < 1) :
    ∀ k : ℤ, 0 ≤ k → k ≤ n →
      pbPmf n p (k - 1) * pbPmf n p (k + 1) < pbPmf n p k ^ 2 :=
  (pbPmf_lc n p hp).2

/-- **G3f** -/
theorem pb_first_descent_exists (n : ℕ) (p : Fin n → ℝ)
    (hp : ∀ i, 0 < p i ∧ p i < 1) (hV : 1 ≤ pbVar n p) :
    ∃ D : ℤ, IsFirstDescent (pbPmf n p) n D :=
  pbPmf_first_descent_exists n p hp hV

/-- **G4a** (Hillion–Johnson, Theorem A.2, their (78)). -/
theorem pb_hj_left (n : ℕ) (p : Fin n → ℝ) (hp : ∀ i, 0 < p i ∧ p i < 1) :
    ∀ k : ℤ,
      (pbPmf n p (k - 1) ^ 2 - pbPmf n p (k - 2) * pbPmf n p k) * pbPmf n p (k + 1) ≤
        pbPmf n p (k - 1) *
          (pbPmf n p k ^ 2 - pbPmf n p (k - 1) * pbPmf n p (k + 1)) :=
  pbPmf_hj_left n p hp

/-- **G4b** (Hillion–Johnson, Corollary A.3, their (79)). -/
theorem pb_hj_right (n : ℕ) (p : Fin n → ℝ) (hp : ∀ i, 0 < p i ∧ p i < 1) :
    ∀ k : ℤ,
      (pbPmf n p (k + 1) ^ 2 - pbPmf n p k * pbPmf n p (k + 2)) * pbPmf n p (k - 1) ≤
        pbPmf n p (k + 1) *
          (pbPmf n p k ^ 2 - pbPmf n p (k - 1) * pbPmf n p (k + 1)) :=
  pbPmf_hj_right n p hp

/-- **G0** (Theorem 1.1). -/
theorem theorem_1_1 (n : ℕ) (p : Fin n → ℝ) (hp : ∀ i, 0 < p i ∧ p i < 1)
    (hV : 1 ≤ pbVar n p) (D : ℤ) (hD : IsFirstDescent (pbPmf n p) n D) :
    1 / 4 ≤ pbVar n p * deficit (pbPmf n p) D := by
  have h : DeductionHyp (pbPmf n p) n :=
    { pos := pb_pos n p hp
      zero_out := pb_zero_out n p
      sum_one := pb_sum_one n p
      strict_lc := pb_strict_lc n p hp
      hj_left := pb_hj_left n p hp
      hj_right := pb_hj_right n p hp }
  have hpv := pb_pairVar n p
  have hmax := max_mass_bound (pbPmf n p) n (pbPmf_nonneg_pos n p hp).1 (pb_sum_one n p)
  rw [← hpv] at hV ⊢
  exact deduction (pbPmf n p) n D h hD hV hmax

end

end PBDeduction
