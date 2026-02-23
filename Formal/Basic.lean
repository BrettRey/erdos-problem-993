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

section Route1

variable {α : Type*} [CommRing α] [LinearOrder α] [IsStrictOrderedRing α]

theorem route1_of_E
    (muP m exact_slack_B exact_excess_D : α)
    (hchain : muP - (m - 2) = exact_slack_B - exact_excess_D)
    (hE : exact_excess_D <= exact_slack_B) :
    m - 2 <= muP := by
  nlinarith [hchain, hE]

theorem E_of_route1
    (muP m exact_slack_B exact_excess_D : α)
    (hchain : muP - (m - 2) = exact_slack_B - exact_excess_D)
    (hroute1 : m - 2 <= muP) :
    exact_excess_D <= exact_slack_B := by
  nlinarith [hchain, hroute1]

theorem route1_iff_E
    (muP m exact_slack_B exact_excess_D : α)
    (hchain : muP - (m - 2) = exact_slack_B - exact_excess_D) :
    (m - 2 <= muP) ↔ (exact_excess_D <= exact_slack_B) := by
  constructor
  · intro hroute1
    exact E_of_route1 muP m exact_slack_B exact_excess_D hchain hroute1
  · intro hE
    exact route1_of_E muP m exact_slack_B exact_excess_D hchain hE

end Route1

section Route1Tail

variable {α : Type*} [Field α] [LinearOrder α] [IsStrictOrderedRing α]

theorem mu_ge_of_tail_certificate
    (lam Pval dP m headGain tailDef : α)
    (hP : 0 < Pval)
    (hbalance : lam * dP - (m - 2) * Pval = headGain - tailDef)
    (hcert : tailDef <= headGain) :
    m - 2 <= (lam * dP) / Pval := by
  have hnum_nonneg : 0 <= lam * dP - (m - 2) * Pval := by
    nlinarith [hbalance, hcert]
  have hmul : (m - 2) * Pval <= lam * dP := by
    nlinarith [hnum_nonneg]
  exact (le_div_iff₀ hP).2 hmul

theorem route1_of_tail_certificate
    (muP m exact_slack_B exact_excess_D : α)
    (lam Pval dP headGain tailDef : α)
    (hchain : muP - (m - 2) = exact_slack_B - exact_excess_D)
    (hmu : muP = (lam * dP) / Pval)
    (hP : 0 < Pval)
    (hbalance : lam * dP - (m - 2) * Pval = headGain - tailDef)
    (hcert : tailDef <= headGain) :
    exact_excess_D <= exact_slack_B := by
  have hroute1 : m - 2 <= muP := by
    calc
      m - 2 <= (lam * dP) / Pval :=
        mu_ge_of_tail_certificate lam Pval dP m headGain tailDef hP hbalance hcert
      _ = muP := by simp [hmu]
  exact E_of_route1 muP m exact_slack_B exact_excess_D hchain hroute1

theorem mu_ge_of_tail2_certificate
    (lam Pval dP m p1Term pmTerm tailPlus tailDef : α)
    (hP : 0 < Pval)
    (hdecomp :
      lam * dP - (m - 2) * Pval = -tailDef + p1Term + pmTerm + tailPlus)
    (htailPlus_nonneg : 0 <= tailPlus)
    (hcert : tailDef <= p1Term + pmTerm) :
    m - 2 <= (lam * dP) / Pval := by
  have hnum_nonneg : 0 <= lam * dP - (m - 2) * Pval := by
    nlinarith [hdecomp, htailPlus_nonneg, hcert]
  have hmul : (m - 2) * Pval <= lam * dP := by
    nlinarith [hnum_nonneg]
  exact (le_div_iff₀ hP).2 hmul

theorem route1_of_tail2_certificate
    (muP m exact_slack_B exact_excess_D : α)
    (lam Pval dP p1Term pmTerm tailPlus tailDef : α)
    (hchain : muP - (m - 2) = exact_slack_B - exact_excess_D)
    (hmu : muP = (lam * dP) / Pval)
    (hP : 0 < Pval)
    (hdecomp :
      lam * dP - (m - 2) * Pval = -tailDef + p1Term + pmTerm + tailPlus)
    (htailPlus_nonneg : 0 <= tailPlus)
    (hcert : tailDef <= p1Term + pmTerm) :
    exact_excess_D <= exact_slack_B := by
  have hroute1 : m - 2 <= muP := by
    calc
      m - 2 <= (lam * dP) / Pval :=
        mu_ge_of_tail2_certificate lam Pval dP m p1Term pmTerm tailPlus tailDef
          hP hdecomp htailPlus_nonneg hcert
      _ = muP := by simp [hmu]
  exact E_of_route1 muP m exact_slack_B exact_excess_D hchain hroute1

end Route1Tail

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

section OrderedField

variable {α : Type*} [Field α] [LinearOrder α] [IsStrictOrderedRing α]

theorem hard_ratio_of_rise_nonneg
    (p1 b1 db dq : α)
    (hb1 : 0 < b1)
    (hdb : 0 < db)
    (hrise : 0 <= p1 * db + b1 * dq) :
    (-dq / db) <= p1 / b1 := by
  have hmul : (-dq) * b1 <= p1 * db := by
    nlinarith [hrise]
  have hbound : -dq <= (p1 / b1) * db := by
    have htmp : -dq <= (p1 * db) / b1 :=
      (le_div_iff₀ hb1).2 (by simpa [mul_comm, mul_left_comm, mul_assoc] using hmul)
    calc
      -dq <= (p1 * db) / b1 := htmp
      _ = (p1 / b1) * db := by
          field_simp [ne_of_gt hb1]
  exact (div_le_iff₀ hdb).2 hbound

theorem rise_nonneg_of_star_conditions
    (p1 q1 b1 db dq : α)
    (hb1 : b1 = p1 + q1)
    (hp1_ge_q1 : q1 <= p1)
    (hq1_nonneg : 0 <= q1)
    (hdq_nonpos : dq <= 0)
    (hdb_ge_neg2dq : -2 * dq <= db) :
    0 <= p1 * db + b1 * dq := by
  have hp1_nonneg : 0 <= p1 := le_trans hq1_nonneg hp1_ge_q1
  have hb1_le_two_p1 : b1 <= 2 * p1 := by
    nlinarith [hb1, hp1_ge_q1]
  have hmul : b1 * dq >= (2 * p1) * dq := by
    nlinarith [hb1_le_two_p1, hdq_nonpos]
  have hdb2 : 0 <= db + 2 * dq := by
    nlinarith [hdb_ge_neg2dq]
  have hcore_nonneg : 0 <= p1 * (db + 2 * dq) := mul_nonneg hp1_nonneg hdb2
  have hrewrite : p1 * db + (2 * p1) * dq = p1 * (db + 2 * dq) := by
    ring
  nlinarith [hmul, hcore_nonneg, hrewrite]

theorem p1_ge_q1_of_h1_and_hard
    (p0 p1 q0 q1 db dq : α)
    (hdb : db = (p1 - p0) + dq)
    (hdq : dq = q1 - q0)
    (hdq_neg : dq < 0)
    (hdb_pos : 0 < db)
    (h1_local : q0 <= p0) :
    q1 <= p1 := by
  nlinarith [hdb, hdq, hdq_neg, hdb_pos, h1_local]

theorem db_ge_neg2dq_of_h2
    (p0 p1 q0 q1 db dq : α)
    (hdb : db = (p1 - p0) + dq)
    (hdq : dq = q1 - q0)
    (h2_local : 3 * (q0 - q1) <= p1 - p0) :
    -2 * dq <= db := by
  nlinarith [hdb, hdq, h2_local]

theorem rise_nonneg_of_h1_h2_hard
    (p0 p1 q0 q1 b1 db dq : α)
    (hb1 : b1 = p1 + q1)
    (hdb : db = (p1 - p0) + dq)
    (hdq : dq = q1 - q0)
    (hdq_neg : dq < 0)
    (hdb_pos : 0 < db)
    (hq1_nonneg : 0 <= q1)
    (h1_local : q0 <= p0)
    (h2_local : 3 * (q0 - q1) <= p1 - p0) :
    0 <= p1 * db + b1 * dq := by
  have hp1_ge_q1 : q1 <= p1 :=
    p1_ge_q1_of_h1_and_hard p0 p1 q0 q1 db dq hdb hdq hdq_neg hdb_pos h1_local
  have hdb_ge_neg2dq : -2 * dq <= db :=
    db_ge_neg2dq_of_h2 p0 p1 q0 q1 db dq hdb hdq h2_local
  exact
    rise_nonneg_of_star_conditions
      p1 q1 b1 db dq hb1 hp1_ge_q1 hq1_nonneg (le_of_lt hdq_neg) hdb_ge_neg2dq

theorem hard_ratio_of_h1_h2_hard
    (p0 p1 q0 q1 b1 db dq : α)
    (hb1 : b1 = p1 + q1)
    (hdb : db = (p1 - p0) + dq)
    (hdq : dq = q1 - q0)
    (hdq_neg : dq < 0)
    (hdb_pos : 0 < db)
    (hb1_pos : 0 < b1)
    (hq1_nonneg : 0 <= q1)
    (h1_local : q0 <= p0)
    (h2_local : 3 * (q0 - q1) <= p1 - p0) :
    (-dq / db) <= p1 / b1 := by
  have hrise_nonneg : 0 <= p1 * db + b1 * dq :=
    rise_nonneg_of_h1_h2_hard
      p0 p1 q0 q1 b1 db dq
      hb1 hdb hdq hdq_neg hdb_pos hq1_nonneg h1_local h2_local
  exact hard_ratio_of_rise_nonneg p1 b1 db dq hb1_pos hdb_pos hrise_nonneg

theorem rise_nonneg_of_hard_ratio
    (p1 b1 db dq : α)
    (hb1 : 0 < b1)
  (hdb : 0 < db)
  (hratio : (-dq / db) <= p1 / b1) :
  0 <= p1 * db + b1 * dq := by
  have hbound : -dq <= (p1 / b1) * db :=
    (div_le_iff₀ hdb).1 hratio
  have hmul : (-dq) * b1 <= ((p1 / b1) * db) * b1 :=
    mul_le_mul_of_nonneg_right hbound (le_of_lt hb1)
  have hcore : (-dq) * b1 <= p1 * db := by
    calc
      (-dq) * b1 <= ((p1 / b1) * db) * b1 := hmul
      _ = p1 * db := by
          field_simp [ne_of_gt hb1]
  nlinarith [hcore]

theorem rise_nonneg_iff_hard_ratio
    (p1 b1 db dq : α)
    (hb1 : 0 < b1)
    (hdb : 0 < db) :
    0 <= p1 * db + b1 * dq ↔ (-dq / db) <= p1 / b1 := by
  constructor
  · intro hrise
    exact hard_ratio_of_rise_nonneg p1 b1 db dq hb1 hdb hrise
  · intro hratio
    exact rise_nonneg_of_hard_ratio p1 b1 db dq hb1 hdb hratio

theorem hard_ratio_of_combined_and_sign
    (p0 p1 b0 b1 b2 db dq riseNeg : α)
    (hsplit : combined p0 p1 b0 b1 b2 = riseNeg + b0 * (b1 - b2))
    (hcombined : 0 <= combined p0 p1 b0 b1 b2)
    (hb0_nonneg : 0 <= b0)
    (hb1_sub_b2_nonpos : b1 - b2 <= 0)
  (hb1 : 0 < b1)
  (hdb : 0 < db)
  (hrise_def : riseNeg = p1 * db + b1 * dq) :
  (-dq / db) <= p1 / b1 := by
  have hrise_nonneg : 0 <= riseNeg :=
    rise_neg_nonneg_of_combined_and_sign
      p0 p1 b0 b1 b2 riseNeg hsplit hcombined hb0_nonneg hb1_sub_b2_nonpos
  have hrise_nonneg' : 0 <= p1 * db + b1 * dq := by
    simpa [hrise_def] using hrise_nonneg
  exact hard_ratio_of_rise_nonneg p1 b1 db dq hb1 hdb hrise_nonneg'

theorem hard_ratio_of_combined_and_Esign
    (p0 p1 b0 b1 b2 db dq riseNeg : α)
    (hsplit : combined p0 p1 b0 b1 b2 = riseNeg + b0 * (b1 - b2))
    (hcombined : 0 <= combined p0 p1 b0 b1 b2)
    (hb0_nonneg : 0 <= b0)
    (hE_sign : b1 - b2 <= 0)
    (hb1 : 0 < b1)
    (hdb : 0 < db)
    (hrise_def : riseNeg = p1 * db + b1 * dq) :
    (-dq / db) <= p1 / b1 := by
  exact
    hard_ratio_of_combined_and_sign
      p0 p1 b0 b1 b2 db dq riseNeg
      hsplit hcombined hb0_nonneg hE_sign hb1 hdb hrise_def

theorem rise_nonneg_of_det_factorization
    (p1 b1 db dq dAttach dBoundary : α)
    (hfactor : p1 * db + b1 * dq = dAttach * dBoundary)
    (hAttach_nonneg : 0 <= dAttach)
    (hBoundary_nonneg : 0 <= dBoundary) :
    0 <= p1 * db + b1 * dq := by
  calc
    0 <= dAttach * dBoundary := mul_nonneg hAttach_nonneg hBoundary_nonneg
    _ = p1 * db + b1 * dq := by simp [hfactor]

theorem hard_ratio_of_det_factorization
    (p1 b1 db dq dAttach dBoundary : α)
    (hfactor : p1 * db + b1 * dq = dAttach * dBoundary)
    (hAttach_nonneg : 0 <= dAttach)
    (hBoundary_nonneg : 0 <= dBoundary)
    (hb1 : 0 < b1)
    (hdb : 0 < db) :
    (-dq / db) <= p1 / b1 := by
  have hrise_nonneg : 0 <= p1 * db + b1 * dq :=
    rise_nonneg_of_det_factorization
      p1 b1 db dq dAttach dBoundary hfactor hAttach_nonneg hBoundary_nonneg
  exact hard_ratio_of_rise_nonneg p1 b1 db dq hb1 hdb hrise_nonneg

theorem rise_nonneg_of_matrix_factorization
    (p1 b1 db dq : α)
    (m11 m12 m21 m22 a11 a12 a21 a22 j11 j12 j21 j22 : α)
    (hentry11 : m11 = a11 * j11 + a12 * j21)
    (hentry12 : m12 = a11 * j12 + a12 * j22)
    (hentry21 : m21 = a21 * j11 + a22 * j21)
    (hentry22 : m22 = a21 * j12 + a22 * j22)
    (hdetM : m11 * m22 - m12 * m21 = p1 * db + b1 * dq)
    (hdetA_nonneg : 0 <= a11 * a22 - a12 * a21)
    (hdetJ_nonneg : 0 <= j11 * j22 - j12 * j21) :
    0 <= p1 * db + b1 * dq := by
  have hdet_factor :
      p1 * db + b1 * dq =
        (a11 * a22 - a12 * a21) * (j11 * j22 - j12 * j21) := by
    calc
      p1 * db + b1 * dq = m11 * m22 - m12 * m21 := by nlinarith [hdetM]
      _ = (a11 * j11 + a12 * j21) * (a21 * j12 + a22 * j22) -
            (a11 * j12 + a12 * j22) * (a21 * j11 + a22 * j21) := by
            simp [hentry11, hentry12, hentry21, hentry22]
      _ = (a11 * a22 - a12 * a21) * (j11 * j22 - j12 * j21) := by ring
  exact
    rise_nonneg_of_det_factorization
      p1 b1 db dq
      (a11 * a22 - a12 * a21) (j11 * j22 - j12 * j21)
      hdet_factor hdetA_nonneg hdetJ_nonneg

end OrderedField

section StrongC2

variable {α : Type*} [CommRing α] [LinearOrder α] [IsStrictOrderedRing α]

theorem strong_c2_of_lc_and_abs
    (lhs p0 p1 pm q0 q1 qm b0 b1 b2 : α)
    (hsplit : lhs = combined p0 p1 b0 b1 b2)
    (hb0 : b0 = p0 + q0)
    (hb1 : b1 = p1 + q1)
    (hb2 : b2 = pm + qm)
    (hlcP : 0 <= lcP p0 p1 pm)
    (hlcQ : 0 <= lcQ q0 q1 qm)
    (habs : |p0 * q1 - p1 * q0| <= cross p0 p1 pm q0 q1 qm) :
    0 <= lhs := by
  have hcombined : 0 <= combined p0 p1 b0 b1 b2 :=
    combined_nonneg_of_lc_and_abs p0 p1 pm q0 q1 qm b0 b1 b2
      hb0 hb1 hb2 hlcP hlcQ habs
  nlinarith [hsplit, hcombined]

end StrongC2

section Counterexample

open scoped BigOperators

theorem hard_regime_not_imply_b1_sub_b2_nonpos_example1 :
    let dq : ℚ := -14
    let db : ℚ := 658
    let b1 : ℚ := 3640
    let b2 : ℚ := 2844
    dq < 0 ∧ 0 < db ∧ ¬ (b1 - b2 <= 0) := by
  norm_num

theorem hard_regime_not_imply_b1_sub_b2_nonpos_example2 :
    let dq : ℚ := -48
    let db : ℚ := 1793
    let b1 : ℚ := 14297
    let b2 : ℚ := 11314
    dq < 0 ∧ 0 < db ∧ ¬ (b1 - b2 <= 0) := by
  norm_num

end Counterexample

section ConstructorSchema

set_option linter.unusedSectionVars false

variable {α : Type*} [CommRing α] [LinearOrder α] [IsStrictOrderedRing α]
variable {σ : Type*}

def HardBad (p1 b1 db dq : α) : Prop :=
  dq < 0 ∧ 0 < db ∧ p1 * db + b1 * dq < 0

def RouteBad (exact_slack_B exact_excess_D : α) : Prop :=
  exact_slack_B < exact_excess_D

theorem hard_bad_impossible_of_nonneg
    (p1 b1 db dq : α)
    (hE_hard : 0 <= p1 * db + b1 * dq) :
    ¬ HardBad p1 b1 db dq := by
  intro hbad
  exact not_lt_of_ge hE_hard hbad.2.2

theorem route_bad_impossible_of_order
    (exact_slack_B exact_excess_D : α)
    (hE_route1 : exact_excess_D <= exact_slack_B) :
    ¬ RouteBad exact_slack_B exact_excess_D := by
  intro hbad
  exact (not_lt_of_ge hE_route1) (by simpa [RouteBad] using hbad)

theorem no_bad_task_if_reachable_invariant
    (Reachable : σ → Prop)
    (M : σ → α)
    (hInv : ∀ s, Reachable s → 0 <= M s) :
    ∀ s, Reachable s → ¬ (M s < 0) := by
  intro s hs hbad
  exact not_lt_of_ge (hInv s hs) hbad

end ConstructorSchema

section FinalClosure

variable {α : Type*} [Field α] [LinearOrder α] [IsStrictOrderedRing α]

theorem route1_and_hard_ratio_of_E_lemmas
    (muP m exact_slack_B exact_excess_D : α)
    (p0 p1 b0 b1 b2 db dq riseNeg : α)
    (hchain : muP - (m - 2) = exact_slack_B - exact_excess_D)
    (hE_route1 : exact_excess_D <= exact_slack_B)
    (hsplit : combined p0 p1 b0 b1 b2 = riseNeg + b0 * (b1 - b2))
    (hcombined : 0 <= combined p0 p1 b0 b1 b2)
    (hb0_nonneg : 0 <= b0)
    (hE_sign : b1 - b2 <= 0)
    (hb1 : 0 < b1)
    (hdb : 0 < db)
    (hrise_def : riseNeg = p1 * db + b1 * dq) :
    (m - 2 <= muP) ∧ ((-dq / db) <= p1 / b1) := by
  constructor
  · exact route1_of_E muP m exact_slack_B exact_excess_D hchain hE_route1
  · exact
      hard_ratio_of_combined_and_Esign
        p0 p1 b0 b1 b2 db dq riseNeg
        hsplit hcombined hb0_nonneg hE_sign hb1 hdb hrise_def

theorem route1_and_hard_ratio_of_minimal_E
    (muP m exact_slack_B exact_excess_D : α)
    (p1 b1 db dq : α)
    (hchain : muP - (m - 2) = exact_slack_B - exact_excess_D)
    (hE_route1 : exact_excess_D <= exact_slack_B)
    (hE_hard : 0 <= p1 * db + b1 * dq)
    (hb1 : 0 < b1)
    (hdb : 0 < db) :
    (m - 2 <= muP) ∧ ((-dq / db) <= p1 / b1) := by
  constructor
  · exact route1_of_E muP m exact_slack_B exact_excess_D hchain hE_route1
  · exact hard_ratio_of_rise_nonneg p1 b1 db dq hb1 hdb hE_hard

end FinalClosure

end Bridge
end Formal
