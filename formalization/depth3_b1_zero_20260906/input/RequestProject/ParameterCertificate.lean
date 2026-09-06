import Mathlib.Tactic

namespace DepthThree.ParameterCertificate

theorem five_ray_cone
    (t0 t1 t2 t3 t4 q0 q1 q2 q3 q4 : ℝ)
    (ht0 : 0 ≤ t0) (ht1 : t0 ≤ t1) (ht2 : t1 ≤ t2)
    (ht3 : t2 ≤ t3) (ht4 : t3 ≤ t4)
    (hq0 : 0 ≤ q0 + q1 + q2 + q3 + q4)
    (hq1 : 0 ≤ q1 + q2 + q3 + q4) (hq2 : 0 ≤ q2 + q3 + q4)
    (hq3 : 0 ≤ q3 + q4) (hq4 : 0 ≤ q4) :
    0 ≤ q0 * t0 + q1 * t1 + q2 * t2 + q3 * t3 + q4 * t4 := by
  nlinarith [mul_nonneg ht0 hq0, mul_nonneg (sub_nonneg.mpr ht1) hq1,
    mul_nonneg (sub_nonneg.mpr ht2) hq2, mul_nonneg (sub_nonneg.mpr ht3) hq3,
    mul_nonneg (sub_nonneg.mpr ht4) hq4]

def beta (m j : ℕ) : ℚ := (m.choose j : ℚ) / 2 ^ j

def numeratorCoeff (r c j : ℕ) : ℚ :=
  match j with
  | 0 => r * ((c - 3).choose 2 : ℚ) + (r.choose 2 : ℚ) * (c - 5 : ℕ) + r.choose 3
  | 1 => r * (c - 3 : ℕ) + (r.choose 2 : ℚ)
  | 2 => r
  | _ => 0

def upper2 (a delta r : ℕ) : ℚ :=
  let c := r + delta
  let m := a - c
  r / (beta m 2 + c * beta m 1 + (c.choose 2 : ℚ))

def upper4Ray (a delta r k : ℕ) : ℚ :=
  let c := r + delta
  let m := a - c
  (∑ j ∈ Finset.Icc k 4, numeratorCoeff r c j * beta m j) /
    (∑ j ∈ Finset.Icc k 4, (c.choose (4 - j) : ℚ) * beta m j)

/-! Exact finite arithmetic only; graph coefficient and cone bridges are separate. -/
set_option maxHeartbeats 4000000 in
theorem all_small_forbidden_rays
    (a delta r k : ℕ) (ha : 17 ≤ a) (ha' : a ≤ 19)
    (hd : delta ≤ 5) (hn : 33 + delta ≤ 2 * a) (hn' : 2 * a ≤ 38 + delta)
    (hr : 1 ≤ r) (hr' : r ≤ 3) (hrd : r + 1 ≤ delta) (hk : k ≤ 4) :
    upper2 a delta r + upper4Ray a delta r k +
      upper2 a delta r * upper4Ray a delta r k < (7 : ℚ) / 27 := by
  interval_cases a <;> interval_cases delta <;> interval_cases r <;>
    norm_num at * <;>
    interval_cases k <;>
    norm_num [upper2, upper4Ray, beta, numeratorCoeff, Finset.sum_Icc_succ_top, Nat.choose]

#print axioms all_small_forbidden_rays
#print axioms five_ray_cone

/-! Once the low-density graph theorem is available, this sharper routing
avoids the separate b2 bound and the joint endpoint assembly. -/
set_option maxHeartbeats 4000000 in
theorem all_small_forbidden_density_rays
    (a delta r k : ℕ) (ha : 17 ≤ a) (ha' : a ≤ 19)
    (hd : delta ≤ 5) (hn : 33 + delta ≤ 2 * a) (hn' : 2 * a ≤ 38 + delta)
    (hr : 1 ≤ r) (hr' : r ≤ 3) (hrd : r + 1 ≤ delta) (hk : k ≤ 4) :
    upper4Ray a delta r k ≤ ((a : ℚ) - 7) / (3 * ((a : ℚ) - 3)) := by
  interval_cases a <;> interval_cases delta <;> interval_cases r <;>
    norm_num at * <;>
    interval_cases k <;>
    norm_num [upper4Ray, beta, numeratorCoeff, Finset.sum_Icc_succ_top, Nat.choose]

#print axioms all_small_forbidden_density_rays

end DepthThree.ParameterCertificate
