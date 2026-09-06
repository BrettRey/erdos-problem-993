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

theorem four_forbidden_algebra
    (e2 e3 e4 b2 b3 b4 : ℝ)
    (he2 : 0 < e2) (he4 : 0 < e4) (hb2 : 0 ≤ b2) (hb3 : 0 ≤ b3)
    (hreserve : 34 * (e2 * e4) ≤ 27 * e3 ^ 2)
    (hblocked : 11 * b2 ≤ b3) (hextendable : e4 ≤ 4 * e3)
    (hdensity : 217404 * b4 ≤ 52513 * e4) :
    (e2 + b2) * (e4 + b4) < (e3 + b3) ^ 2 := by
  have hcross := mul_le_mul hblocked hextendable (le_of_lt he4) hb3
  have hd2 := mul_le_mul_of_nonneg_left hdensity (le_of_lt he2)
  have hdb := mul_le_mul_of_nonneg_left hdensity hb2
  have hp : 0 < e2 * e4 := mul_pos he2 he4
  have hpb : 0 ≤ e4 * b2 := mul_nonneg (le_of_lt he4) hb2
  nlinarith [sq_nonneg b3]

#print axioms four_forbidden_algebra

/-! Assembly for the small-forbidden branch once the graph density bounds and
the cone certificate have been transported to these real quantities. -/
theorem small_forbidden_algebra
    (e2 e3 e4 b2 b3 b4 u v : ℝ)
    (he2 : 0 < e2) (he3 : 0 ≤ e3) (he4 : 0 < e4)
    (hb2 : 0 ≤ b2) (hb3 : 0 ≤ b3)
    (hu : 0 ≤ u) (hv : 0 ≤ v)
    (hreserve : 34 * (e2 * e4) ≤ 27 * e3 ^ 2)
    (h2 : b2 ≤ u * e2) (h4 : b4 ≤ v * e4)
    (hjoint : u + v + u * v < (7 : ℝ) / 27) :
    (e2 + b2) * (e4 + b4) < (e3 + b3) ^ 2 := by
  have h24 := (mul_le_mul_of_nonneg_left h4 hb2).trans
    (mul_le_mul_of_nonneg_right h2 (by positivity : 0 ≤ v * e4))
  have h2e := mul_le_mul_of_nonneg_right h2 (le_of_lt he4)
  have he4b := mul_le_mul_of_nonneg_left h4 (le_of_lt he2)
  have hj := mul_lt_mul_of_pos_right hjoint (mul_pos he2 he4)
  have hcross := mul_nonneg he3 hb3
  nlinarith [sq_nonneg b3]

#print axioms small_forbidden_algebra

end DepthThree
