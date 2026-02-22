import Mathlib

namespace Formal
namespace Bridge

section CommRing

variable {α : Type*} [CommRing α]

def mismatch (p0 p1 b0 b1 : α) : α := p0 * b1 - p1 * b0

def lcSurplus (b0 b1 b2 : α) : α := b1 ^ 2 - b2 * b0

def lcP (p0 p1 pm : α) : α := p1 ^ 2 - pm * p0

def lcQ (q0 q1 qm : α) : α := q1 ^ 2 - qm * q0

def cross (p0 p1 pm q0 q1 qm : α) : α := 2 * p1 * q1 - pm * q0 - p0 * qm

def combined (p0 p1 b0 b1 b2 : α) : α :=
  lcSurplus b0 b1 b2 + mismatch p0 p1 b0 b1

theorem mismatch_eq_pq
    (p0 p1 q0 q1 b0 b1 : α)
    (hb0 : b0 = p0 + q0)
    (hb1 : b1 = p1 + q1) :
    mismatch p0 p1 b0 b1 = p0 * q1 - p1 * q0 := by
  subst hb0
  subst hb1
  simp [mismatch]
  ring

theorem lc_surplus_expansion
    (p0 p1 pm q0 q1 qm b0 b1 b2 : α)
    (hb0 : b0 = p0 + q0)
    (hb1 : b1 = p1 + q1)
    (hb2 : b2 = pm + qm) :
    lcSurplus b0 b1 b2 =
      lcP p0 p1 pm + lcQ q0 q1 qm + cross p0 p1 pm q0 q1 qm := by
  subst hb0
  subst hb1
  subst hb2
  simp [lcSurplus, lcP, lcQ, cross]
  ring

theorem combined_decomposition
    (p0 p1 pm q0 q1 qm b0 b1 b2 : α)
    (hb0 : b0 = p0 + q0)
    (hb1 : b1 = p1 + q1)
    (hb2 : b2 = pm + qm) :
    combined p0 p1 b0 b1 b2 =
      lcP p0 p1 pm + lcQ q0 q1 qm + cross p0 p1 pm q0 q1 qm + (p0 * q1 - p1 * q0) := by
  subst hb0
  subst hb1
  subst hb2
  simp [combined, lcSurplus, lcP, lcQ, cross, mismatch]
  ring

theorem rise_neg_identity
    (p0 p1 q0 q1 b0 b1 db dq : α)
    (hb0 : b0 = p0 + q0)
    (hb1 : b1 = p1 + q1)
    (hdb : db = b1 - b0)
    (hdq : dq = q1 - q0) :
    (b1 * db + mismatch p0 p1 b0 b1) = p1 * db + b1 * dq := by
  subst hb0
  subst hb1
  subst hdb
  subst hdq
  simp [mismatch]
  ring

theorem i_m1_eq_b1_add_b0_add_p0
    (i_m1 p0 p1 q0 q1 b0 b1 : α)
    (hi_m1 : i_m1 = p1 + 2 * p0 + q1 + q0)
    (hb0 : b0 = p0 + q0)
    (hb1 : b1 = p1 + q1) :
    i_m1 = b1 + b0 + p0 := by
  subst hi_m1
  subst hb0
  subst hb1
  ring

theorem i_m_eq_b2_add_b1_add_p1
    (i_m p1 pm q1 qm b1 b2 : α)
    (hi_m : i_m = pm + 2 * p1 + qm + q1)
    (hb1 : b1 = p1 + q1)
    (hb2 : b2 = pm + qm) :
    i_m = b2 + b1 + p1 := by
  subst hi_m
  subst hb1
  subst hb2
  ring

end CommRing

section Ordered

variable {α : Type*} [CommRing α] [LinearOrder α] [IsStrictOrderedRing α]

theorem cross_plus_mismatch_nonneg_of_abs_bound
    (p0 p1 pm q0 q1 qm : α)
    (habs : |p0 * q1 - p1 * q0| <= cross p0 p1 pm q0 q1 qm) :
    0 <= cross p0 p1 pm q0 q1 qm + (p0 * q1 - p1 * q0) := by
  have hleft : -(cross p0 p1 pm q0 q1 qm) <= p0 * q1 - p1 * q0 := (abs_le.mp habs).1
  nlinarith

theorem cross_minus_mismatch_nonneg_of_abs_bound
    (p0 p1 pm q0 q1 qm : α)
    (habs : |p0 * q1 - p1 * q0| <= cross p0 p1 pm q0 q1 qm) :
    0 <= cross p0 p1 pm q0 q1 qm - (p0 * q1 - p1 * q0) := by
  have hright : p0 * q1 - p1 * q0 <= cross p0 p1 pm q0 q1 qm := (abs_le.mp habs).2
  nlinarith

theorem combined_nonneg_of_component_bounds
    (p0 p1 pm q0 q1 qm b0 b1 b2 : α)
    (hb0 : b0 = p0 + q0)
    (hb1 : b1 = p1 + q1)
    (hb2 : b2 = pm + qm)
    (hlcP : 0 <= lcP p0 p1 pm)
    (hlcQ : 0 <= lcQ q0 q1 qm)
    (hcross_mismatch : 0 <= cross p0 p1 pm q0 q1 qm + (p0 * q1 - p1 * q0)) :
    0 <= combined p0 p1 b0 b1 b2 := by
  have hdecomp :
      combined p0 p1 b0 b1 b2 =
        lcP p0 p1 pm + lcQ q0 q1 qm + cross p0 p1 pm q0 q1 qm + (p0 * q1 - p1 * q0) :=
    combined_decomposition p0 p1 pm q0 q1 qm b0 b1 b2 hb0 hb1 hb2
  nlinarith [hlcP, hlcQ, hcross_mismatch, hdecomp]

theorem combined_nonneg_of_lc_and_abs
    (p0 p1 pm q0 q1 qm b0 b1 b2 : α)
    (hb0 : b0 = p0 + q0)
    (hb1 : b1 = p1 + q1)
    (hb2 : b2 = pm + qm)
    (hlcP : 0 <= lcP p0 p1 pm)
    (hlcQ : 0 <= lcQ q0 q1 qm)
    (habs : |p0 * q1 - p1 * q0| <= cross p0 p1 pm q0 q1 qm) :
    0 <= combined p0 p1 b0 b1 b2 := by
  have hcross_mismatch :
      0 <= cross p0 p1 pm q0 q1 qm + (p0 * q1 - p1 * q0) :=
    cross_plus_mismatch_nonneg_of_abs_bound p0 p1 pm q0 q1 qm habs
  exact
    combined_nonneg_of_component_bounds
      p0 p1 pm q0 q1 qm b0 b1 b2 hb0 hb1 hb2 hlcP hlcQ hcross_mismatch

theorem rise_neg_nonneg_of_combined_and_sign
    (p0 p1 b0 b1 b2 riseNeg : α)
    (hsplit : combined p0 p1 b0 b1 b2 = riseNeg + b0 * (b1 - b2))
    (hcombined : 0 <= combined p0 p1 b0 b1 b2)
    (hb0_nonneg : 0 <= b0)
    (hb1_sub_b2_nonpos : b1 - b2 <= 0) :
    0 <= riseNeg := by
  have hterm_nonpos : b0 * (b1 - b2) <= 0 :=
    mul_nonpos_of_nonneg_of_nonpos hb0_nonneg hb1_sub_b2_nonpos
  have hrise : riseNeg = combined p0 p1 b0 b1 b2 - b0 * (b1 - b2) := by
    nlinarith [hsplit]
  nlinarith [hcombined, hterm_nonpos, hrise]

end Ordered

section Field

variable {α : Type*} [Field α]

theorem lambda_eq_bridge_ratio
    (lam i_m1 i_m p0 p1 pm q0 q1 qm b0 b1 b2 : α)
    (hlam : lam = i_m1 / i_m)
    (hi_m1 : i_m1 = p1 + 2 * p0 + q1 + q0)
    (hi_m : i_m = pm + 2 * p1 + qm + q1)
    (hb0 : b0 = p0 + q0)
    (hb1 : b1 = p1 + q1)
    (hb2 : b2 = pm + qm) :
    lam = (b1 + b0 + p0) / (b2 + b1 + p1) := by
  have hi_m1' : i_m1 = b1 + b0 + p0 :=
    i_m1_eq_b1_add_b0_add_p0 i_m1 p0 p1 q0 q1 b0 b1 hi_m1 hb0 hb1
  have hi_m' : i_m = b2 + b1 + p1 :=
    i_m_eq_b2_add_b1_add_p1 i_m p1 pm q1 qm b1 b2 hi_m hb1 hb2
  calc
    lam = i_m1 / i_m := hlam
    _ = (b1 + b0 + p0) / i_m := by simp [hi_m1']
    _ = (b1 + b0 + p0) / (b2 + b1 + p1) := by simp [hi_m']

end Field

end Bridge
end Formal
