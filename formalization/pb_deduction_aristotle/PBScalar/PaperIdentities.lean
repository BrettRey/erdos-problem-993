import PBScalar.Quadratic

/-!
# Cross-check of a printed identity from the paper (not used by the main proof)

Paper, eq. `asymmetric-P`: on `3 < H ≤ 4` (`K = 3`, unsymmetrised weights),
`A(δ) - Q(H) = P(H) / (4 H^5 (H+1)^3)` with `δ = 1/(H+1)` and the displayed
degree-10 polynomial `P`. We also check the certificate data for the cell `[3, 4]`
(numerator in `t = H - 3`, its degree-10 Bernstein expansion, positivity of the 11
coefficients) and derive the cell by the paper's route (`cell3_paper`).
-/

namespace PBScalar

/-- The degree-10 polynomial `P(H)` printed in the paper. -/
def Pasym (H : ℝ) : ℝ :=
  -3 * H ^ 10 - 16 * H ^ 9 + 750 * H ^ 8 - 3676 * H ^ 7 + 6613 * H ^ 6 - 5460 * H ^ 5
    + 800 * H ^ 4 + 4696 * H ^ 3 - 7176 * H ^ 2 + 4272 * H - 912

/-- Paper, eq. `asymmetric-P`, verified as an exact rational-function identity. -/
theorem asymmetric_identity (H : ℝ) (hH : 0 < H) :
    A (1 / (H + 1)) 3 - (3 * H + 4) * (H + 1) / 4 =
      Pasym H / (4 * H ^ 5 * (H + 1) ^ 3) := by
  have h1 : H + 1 ≠ 0 := by positivity
  have h2 : 1 - 1 / (H + 1) = H / (H + 1) := by field_simp; ring
  rw [A_formula]
  simp only [s0, s1, s2, Finset.sum_range_succ, Finset.sum_range_zero, R_succ, L_succ,
    R_zero, L_zero, a, h2]
  unfold Pasym
  field_simp
  ring


/-- Cell `[3, 4]`: the certificate's numerator is `P(3 + t)`. -/
theorem cell3_numerator (t : ℝ) :
    Pasym (3 + t) =
        (22272 : ℝ) * t ^ 0 +
        (432960 : ℝ) * t ^ 1 +
        (1043808 : ℝ) * t ^ 2 +
        (1111360 : ℝ) * t ^ 3 +
        (641177 : ℝ) * t ^ 4 +
        (205806 : ℝ) * t ^ 5 +
        (31099 : ℝ) * t ^ 6 +
        (-580 : ℝ) * t ^ 7 +
        (-897 : ℝ) * t ^ 8 +
        (-106 : ℝ) * t ^ 9 +
        (-3 : ℝ) * t ^ 10 := by
  unfold Pasym
  ring

/-- Cell `[3, 4]`: Bernstein expansion of `P(3 + t)` (coefficients from the certificate). -/
theorem cell3_bernstein (t : ℝ) :
    (22272 : ℝ) * t ^ 0 +
        (432960 : ℝ) * t ^ 1 +
        (1043808 : ℝ) * t ^ 2 +
        (1111360 : ℝ) * t ^ 3 +
        (641177 : ℝ) * t ^ 4 +
        (205806 : ℝ) * t ^ 5 +
        (31099 : ℝ) * t ^ 6 +
        (-580 : ℝ) * t ^ 7 +
        (-897 : ℝ) * t ^ 8 +
        (-106 : ℝ) * t ^ 9 +
        (-3 : ℝ) * t ^ 10 =
      (22272 : ℝ) *
          1 * t ^ 0 * (1 - t) ^ 10 +
        (65568 : ℝ) *
          10 * t ^ 1 * (1 - t) ^ 9 +
        (1980896/15 : ℝ) *
          45 * t ^ 2 * (1 - t) ^ 8 +
        (3465128/15 : ℝ) *
          120 * t ^ 3 * (1 - t) ^ 7 +
        (26231027/70 : ℝ) *
          210 * t ^ 4 * (1 - t) ^ 6 +
        (12167515/21 : ℝ) *
          252 * t ^ 5 * (1 - t) ^ 5 +
        (30312004/35 : ℝ) *
          210 * t ^ 6 * (1 - t) ^ 4 +
        (6308231/5 : ℝ) *
          120 * t ^ 7 * (1 - t) ^ 3 +
        (27004552/15 : ℝ) *
          45 * t ^ 8 * (1 - t) ^ 2 +
        (12623096/5 : ℝ) *
          10 * t ^ 9 * (1 - t) ^ 1 +
        (3486896 : ℝ) *
          1 * t ^ 10 * (1 - t) ^ 0 := by
  ring

/-- Cell `[3, 4]`: all Bernstein coefficients are positive. -/
theorem cell3_coeffs_pos :
    0 < (22272 : ℝ) ∧
    0 < (65568 : ℝ) ∧
    0 < (1980896/15 : ℝ) ∧
    0 < (3465128/15 : ℝ) ∧
    0 < (26231027/70 : ℝ) ∧
    0 < (12167515/21 : ℝ) ∧
    0 < (30312004/35 : ℝ) ∧
    0 < (6308231/5 : ℝ) ∧
    0 < (27004552/15 : ℝ) ∧
    0 < (12623096/5 : ℝ) ∧
    0 < (3486896 : ℝ) := by
  norm_num

/-- Cell `[3, 4]` by the paper's route: `Q(H) ≤ A(δ)` with `K = 3`, `δ = 1/(H+1)`,
for `H = 3 + t ∈ [3, 4]`. -/
theorem cell3_paper (t : ℝ) (ht0 : 0 ≤ t) (ht1 : t ≤ 1) :
    (3 * (3 + t) + 4) * ((3 + t) + 1) / 4 ≤ A (1 / ((3 + t) + 1)) 3 := by
  have h := asymmetric_identity (3 + t) (by linarith)
  rw [cell3_numerator, cell3_bernstein] at h
  have hs : 0 ≤ 1 - t := by linarith
  have hpos : 0 ≤
      (22272 : ℝ) *
          1 * t ^ 0 * (1 - t) ^ 10 +
        (65568 : ℝ) *
          10 * t ^ 1 * (1 - t) ^ 9 +
        (1980896/15 : ℝ) *
          45 * t ^ 2 * (1 - t) ^ 8 +
        (3465128/15 : ℝ) *
          120 * t ^ 3 * (1 - t) ^ 7 +
        (26231027/70 : ℝ) *
          210 * t ^ 4 * (1 - t) ^ 6 +
        (12167515/21 : ℝ) *
          252 * t ^ 5 * (1 - t) ^ 5 +
        (30312004/35 : ℝ) *
          210 * t ^ 6 * (1 - t) ^ 4 +
        (6308231/5 : ℝ) *
          120 * t ^ 7 * (1 - t) ^ 3 +
        (27004552/15 : ℝ) *
          45 * t ^ 8 * (1 - t) ^ 2 +
        (12623096/5 : ℝ) *
          10 * t ^ 9 * (1 - t) ^ 1 +
        (3486896 : ℝ) *
          1 * t ^ 10 * (1 - t) ^ 0 := by
    generalize 1 - t = s at hs ⊢; positivity
  have hd : (0 : ℝ) < 4 * (3 + t) ^ 5 * ((3 + t) + 1) ^ 3 := by positivity
  have := div_nonneg hpos hd.le
  linarith

end PBScalar
