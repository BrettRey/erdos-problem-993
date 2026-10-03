import PBDeduction.Defs

/-!
# Basic facts about normalized Turán deficits (Section 2 of the paper)

Elementary consequences of `DeductionHyp`: nonnegativity, the bounds
`0 < δ_k ≤ 1` on the support, the endpoint values `δ_0 = δ_n = 1`, the ratio
identity `q_{k+1} = q_k (1 - δ_k)`, the recurrence (2.4) and the reciprocal
bound (2.5).
-/

namespace PBDeduction

open Finset

variable {f : ℤ → ℝ} {n : ℕ}

lemma DeductionHyp.nonneg (h : DeductionHyp f n) (k : ℤ) : 0 ≤ f k := by
  by_cases h0 : 0 ≤ k
  · by_cases hn : k ≤ n
    · exact (h.pos k h0 hn).le
    · rw [h.zero_out k (Or.inr (by omega))]
  · rw [h.zero_out k (Or.inl (by omega))]

lemma DeductionHyp.deficit_pos (h : DeductionHyp f n) {k : ℤ} (h0 : 0 ≤ k) (hn : k ≤ n) :
    0 < deficit f k := by
  have hk := h.pos k h0 hn
  have := h.strict_lc k h0 hn
  unfold deficit
  have : f (k - 1) * f (k + 1) / f k ^ 2 < 1 := by
    rw [div_lt_one (by positivity)]; exact this
  linarith

lemma DeductionHyp.deficit_le_one (h : DeductionHyp f n) (k : ℤ) : deficit f k ≤ 1 := by
  unfold deficit
  have := h.nonneg (k - 1)
  have := h.nonneg (k + 1)
  have : 0 ≤ f (k - 1) * f (k + 1) / f k ^ 2 := by positivity
  linarith

lemma DeductionHyp.deficit_zero (h : DeductionHyp f n) : deficit f 0 = 1 := by
  unfold deficit
  rw [h.zero_out (0 - 1) (Or.inl (by omega))]; simp

lemma DeductionHyp.deficit_top (h : DeductionHyp f n) : deficit f n = 1 := by
  unfold deficit
  rw [h.zero_out ((n : ℤ) + 1) (Or.inr (by omega))]; simp

/-- The ratio identity `f (k+1) / f k = (f k / f (k-1)) (1 - δ_k)`, in
multiplied-out form: `f (k+1) f (k-1) = (1 - δ_k) f k ^ 2`. -/
lemma next_mul_prev (f : ℤ → ℝ) (k : ℤ) (hk : f k ≠ 0) :
    f (k + 1) * f (k - 1) = (1 - deficit f k) * f k ^ 2 := by
  unfold deficit
  field_simp
  ring

/-- **Paper (2.4)**, the recurrence for normalized Turán deficits: for
`1 ≤ k ≤ n - 1`, `δ_{k-1}(1-δ_k) ≤ δ_k` and `δ_{k+1}(1-δ_k) ≤ δ_k`. -/
theorem deficit_recurrence (h : DeductionHyp f n) {k : ℤ} (h1 : 1 ≤ k) (h2 : k ≤ n - 1) :
    deficit f (k - 1) * (1 - deficit f k) ≤ deficit f k ∧
      deficit f (k + 1) * (1 - deficit f k) ≤ deficit f k := by
  have hb := h.pos (k - 1) (by omega) (by omega)
  have hc := h.pos k (by omega) (by omega)
  have hd := h.pos (k + 1) (by omega) (by omega)
  have hl := h.hj_left k
  have hr := h.hj_right k
  have e1 : k - 1 - 1 = k - 2 := by ring
  have e2 : k - 1 + 1 = k := by ring
  have e3 : k + 1 - 1 = k := by ring
  have e4 : k + 1 + 1 = k + 2 := by ring
  unfold deficit
  rw [e1, e2, e3, e4]
  constructor
  · have key : (1 - f (k - 2) * f k / f (k - 1) ^ 2) * (1 - (1 - f (k - 1) * f (k + 1) / f k ^ 2))
        - (1 - f (k - 1) * f (k + 1) / f k ^ 2) =
        ((f (k - 1) ^ 2 - f (k - 2) * f k) * f (k + 1) -
          f (k - 1) * (f k ^ 2 - f (k - 1) * f (k + 1))) / (f (k - 1) * f k ^ 2) := by
      field_simp; ring
    have : ((f (k - 1) ^ 2 - f (k - 2) * f k) * f (k + 1) -
          f (k - 1) * (f k ^ 2 - f (k - 1) * f (k + 1))) / (f (k - 1) * f k ^ 2) ≤ 0 :=
      div_nonpos_of_nonpos_of_nonneg (by linarith) (by positivity)
    linarith
  · have key : (1 - f k * f (k + 2) / f (k + 1) ^ 2) * (1 - (1 - f (k - 1) * f (k + 1) / f k ^ 2))
        - (1 - f (k - 1) * f (k + 1) / f k ^ 2) =
        ((f (k + 1) ^ 2 - f k * f (k + 2)) * f (k - 1) -
          f (k + 1) * (f k ^ 2 - f (k - 1) * f (k + 1))) / (f (k + 1) * f k ^ 2) := by
      field_simp; ring
    have : ((f (k + 1) ^ 2 - f k * f (k + 2)) * f (k - 1) -
          f (k + 1) * (f k ^ 2 - f (k - 1) * f (k + 1))) / (f (k + 1) * f k ^ 2) ≤ 0 :=
      div_nonpos_of_nonpos_of_nonneg (by linarith) (by positivity)
    linarith

/-- From `y (1 - x) ≤ x` with `x, y > 0`: `1/x - 1 ≤ 1/y`. -/
lemma recip_step {x y : ℝ} (hx : 0 < x) (hy : 0 < y) (h : y * (1 - x) ≤ x) :
    1 / x - 1 ≤ 1 / y := by
  rw [div_sub_one hx.ne', div_le_div_iff₀ hx hy]
  nlinarith

/-- **Paper (2.5)**, the reciprocal bound: `|1/δ_{j+1} - 1/δ_j| ≤ 1` for
`0 ≤ j ≤ n - 1`. -/
theorem reciprocal_bound (h : DeductionHyp f n) {j : ℤ} (h0 : 0 ≤ j) (h1 : j ≤ n - 1) :
    |1 / deficit f (j + 1) - 1 / deficit f j| ≤ 1 := by
  have pj := h.deficit_pos h0 (by omega)
  have pj1 := h.deficit_pos (k := j + 1) (by omega) (by omega)
  have lj := h.deficit_le_one j
  have lj1 := h.deficit_le_one (j + 1)
  have r1 : 1 ≤ 1 / deficit f j := by rw [le_div_iff₀ pj]; linarith
  have r2 : 1 ≤ 1 / deficit f (j + 1) := by rw [le_div_iff₀ pj1]; linarith
  rw [abs_le]
  constructor
  · -- `1/δ_{j+1} ≥ 1/δ_j - 1`
    rcases eq_or_lt_of_le h0 with hj | hj
    · subst hj; rw [h.deficit_zero, div_one]; linarith
    · have := (deficit_recurrence h (k := j) (by omega) (by omega)).2
      have := recip_step pj pj1 this
      linarith
  · -- `1/δ_j ≥ 1/δ_{j+1} - 1`
    rcases eq_or_lt_of_le h1 with hj | hj
    · have : j + 1 = n := by omega
      rw [this, h.deficit_top, div_one]; linarith
    · have := (deficit_recurrence h (k := j + 1) (by omega) (by omega)).1
      rw [show j + 1 - 1 = j by ring] at this
      have := recip_step pj1 pj this
      linarith

end PBDeduction
