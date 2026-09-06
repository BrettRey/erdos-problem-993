import Mathlib.Tactic

namespace DepthThree

/-! A scalar implication only. This does not prove its graph hypotheses. -/
theorem low_density_algebra
    (a : ℕ) (ha : 17 ≤ a) (ha' : a ≤ 19)
    (e2 e3 e4 b2 b3 b4 : ℝ)
    (he2 : 0 < e2) (he3 : 0 ≤ e3) (he4 : 0 < e4)
    (hb2 : 0 ≤ b2) (hb3 : 0 ≤ b3) (hb4 : 0 ≤ b4)
    (hreserve : 32 * ((a : ℝ) - 2) * (e2 * e4) ≤
      27 * ((a : ℝ) - 3) * e3 ^ 2)
    (hblocked : ((a : ℝ) - 4) * b2 ≤ 6 * b3)
    (hextendable : 4 * e4 ≤ ((a : ℝ) - 3) * e3)
    (hdensity : 3 * ((a : ℝ) - 3) * b4 ≤ ((a : ℝ) - 7) * e4) :
    (e2 + b2) * (e4 + b4) < (e3 + b3) ^ 2 := by
  have hcross := mul_le_mul hblocked hextendable (by positivity : 0 ≤ 4 * e4)
    (by positivity : 0 ≤ 6 * b3)
  have hd2 := mul_le_mul_of_nonneg_left hdensity (le_of_lt he2)
  have hdb := mul_le_mul_of_nonneg_left hdensity hb2
  have hp : 0 < e2 * e4 := mul_pos he2 he4
  interval_cases a <;> norm_num at * <;> nlinarith [sq_nonneg b3]

#print axioms low_density_algebra

end DepthThree
