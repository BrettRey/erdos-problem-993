import Mathlib.Tactic

/-!
  Axiom-style fixed-r certificate bridge targets.

  This file is intentionally outside the `Formal` Lean library.  It is a
  narrow proof-agent target.  The theorem statements use rational arithmetic so
  the first pass avoids extra polynomial infrastructure.
-/

namespace AxiomFixedR

/--
Adjacent coefficient margins at a known unimodal mode imply a global margin.
The hypotheses `hleft_side` and `hright_side` encode the only part of
unimodality needed here.
-/
theorem global_margin_of_adjacent_margins
    (F : Nat -> Rat) (m : Nat) (eta : Rat)
    (hleft_margin : eta * F m <= F m - F (m - 1))
    (hright_margin : eta * F m <= F m - F (m + 1))
    (hleft_side : forall k, k < m -> F k <= F (m - 1))
    (hright_side : forall k, m < k -> F k <= F (m + 1)) :
    forall k, k ≠ m -> eta * F m <= F m - F k := by
  intro k hk
  by_cases hkm : k < m
  · have hk_le : F k <= F (m - 1) := hleft_side k hkm
    linarith
  · have hmk : m < k := Nat.lt_of_le_of_ne (Nat.le_of_not_gt hkm) (Ne.symm hk)
    have hk_le : F k <= F (m + 1) := hright_side k hmk
    linarith

/--
Conservative tie-fugacity perturbation bound.

This is the formal version of the hand proof in
`notes/fixed_r_fugacity_shift_sublemma_2026-05-21.md`.
-/
theorem fugacity_shift_bound
    (Fm Fprev Gm Gprev lam0 lam T : Rat)
    (hFm : 0 < Fm)
    (hGm_nonneg : 0 <= Gm)
    (hGprev_nonneg : 0 <= Gprev)
    (hT_nonneg : 0 <= T)
    (hT_half : T <= (1 / 2 : Rat))
    (hlam0_def : lam0 = Fprev / Fm)
    (hlam_def : lam = (Fprev + Gprev) / (Fm + Gm))
    (hlam0_low : (3 / 4 : Rat) <= lam0)
    (hlam0_high : lam0 <= 2)
    (hGm_bound : Gm <= T * Fm)
    (hGprev_bound : Gprev <= T * Fm) :
    |lam - lam0| <= 4 * T /\ (1 / 2 : Rat) <= lam := by
  have hFm_ne : Fm ≠ 0 := ne_of_gt hFm
  have hden_pos : 0 < Fm + Gm := by linarith
  have hden_ne : Fm + Gm ≠ 0 := ne_of_gt hden_pos
  have hlam0_nonneg : 0 <= lam0 := by linarith
  have hTF_nonneg : 0 <= T * Fm := by nlinarith
  have hFprev_eq : Fprev = lam0 * Fm := by
    rw [hlam0_def]
    field_simp [hFm_ne]
  have hFprev_lower : (3 / 4 : Rat) * Fm <= Fprev := by
    rw [hFprev_eq]
    nlinarith
  have hGm_half : Gm <= (1 / 2 : Rat) * Fm := by
    nlinarith
  have hden_upper : Fm + Gm <= (3 / 2 : Rat) * Fm := by
    nlinarith
  have hhalf_numer : (1 / 2 : Rat) * (Fm + Gm) <= Fprev + Gprev := by
    nlinarith
  have hlam_lower : (1 / 2 : Rat) <= lam := by
    rw [hlam_def]
    exact (le_div_iff₀ hden_pos).2 hhalf_numer

  have hdiff : lam - lam0 = (Gprev - lam0 * Gm) / (Fm + Gm) := by
    rw [hlam_def, hFprev_eq]
    field_simp [hden_ne]
    ring
  have h_lamG_bound : lam0 * Gm <= 2 * T * Fm := by
    have hmul := mul_le_mul hlam0_high hGm_bound hGm_nonneg (by norm_num : (0 : Rat) <= 2)
    nlinarith
  have hnum_upper : Gprev - lam0 * Gm <= 3 * T * Fm := by
    nlinarith
  have hnum_lower : -(3 * T * Fm) <= Gprev - lam0 * Gm := by
    nlinarith
  have hnum_abs : |Gprev - lam0 * Gm| <= 3 * T * Fm := by
    exact abs_le.mpr ⟨hnum_lower, hnum_upper⟩
  have habs_den : |Gprev - lam0 * Gm| <= 4 * T * (Fm + Gm) := by
    nlinarith
  have habs_shift : |(Gprev - lam0 * Gm) / (Fm + Gm)| <= 4 * T := by
    rw [abs_div, abs_of_pos hden_pos]
    exact (div_le_iff₀ hden_pos).2 habs_den
  constructor
  · rw [hdiff]
    exact habs_shift
  · exact hlam_lower

/--
The mean-shift target after the calculus/variance argument has been packaged as
a Lipschitz certificate.  The harder follow-up is to prove the Lipschitz
certificate from the finite-support Gibbs distribution.
-/
theorem mean_shift_bound_from_lipschitz_certificate
    (N : Nat) (lam0 lam : Rat)
    (mu : Rat -> Rat)
    (hcert :
      |mu lam - mu lam0| <= (((N : Rat) ^ 2) / 2) * |lam - lam0|) :
    |mu lam - mu lam0| <= (((N : Rat) ^ 2) / 2) * |lam - lam0| := by
  exact hcert

end AxiomFixedR
