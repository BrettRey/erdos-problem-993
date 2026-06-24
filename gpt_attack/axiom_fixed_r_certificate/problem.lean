import Mathlib.Tactic

/-!
  Axiom-style fixed-r certificate bridge targets.

  This file is intentionally outside the `Formal` Lean library.  It is a
  narrow proof-agent target.  Coefficient-perturbation lemmas use rational
  arithmetic; the finite-support mean-shift bridge is stated over `ℝ`, where
  the calculus API lives.
-/

namespace AxiomFixedR

open BigOperators

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
If the hub-off coefficients have a global margin at `m`, and the hub-on
perturbation is everywhere smaller than that margin relative to its value at
`m`, then the sum has no larger coefficient away from `m`.
-/
theorem sum_mode_of_global_margin_and_perturbation
    (F G : Nat -> Rat) (m : Nat) (eta : Rat)
    (hglobal : ∀ k, k ≠ m -> eta * F m <= F m - F k)
    (hG_perturb : ∀ k, k ≠ m -> G k - G m <= eta * F m) :
    ∀ k, k ≠ m -> F k + G k <= F m + G m := by
  intro k hk
  have hF := hglobal k hk
  have hG := hG_perturb k hk
  linarith

/--
Certificate-shaped mode-preservation statement: adjacent hub-off margins plus
one-sided hub-on perturbation bounds imply that `m` remains a mode of `F + G`.
-/
theorem sum_mode_of_adjacent_margins_and_perturbation
    (F G : Nat -> Rat) (m : Nat) (eta : Rat)
    (hleft_margin : eta * F m <= F m - F (m - 1))
    (hright_margin : eta * F m <= F m - F (m + 1))
    (hleft_side : forall k, k < m -> F k <= F (m - 1))
    (hright_side : forall k, m < k -> F k <= F (m + 1))
    (hG_perturb : ∀ k, k ≠ m -> G k - G m <= eta * F m) :
    ∀ k, k ≠ m -> F k + G k <= F m + G m := by
  exact sum_mode_of_global_margin_and_perturbation F G m eta
    (global_margin_of_adjacent_margins F m eta hleft_margin hright_margin
      hleft_side hright_side)
    hG_perturb

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
Budgeted form of `fugacity_shift_bound`: if the certificate has already chosen
a permissible shift budget `delta`, it is enough to prove `4T <= delta`.
-/
theorem fugacity_shift_bound_with_budget
    (Fm Fprev Gm Gprev lam0 lam T delta : Rat)
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
    (hGprev_bound : Gprev <= T * Fm)
    (hbudget : 4 * T <= delta) :
    |lam - lam0| <= delta /\ (1 / 2 : Rat) <= lam := by
  have h := fugacity_shift_bound Fm Fprev Gm Gprev lam0 lam T
    hFm hGm_nonneg hGprev_nonneg hT_nonneg hT_half hlam0_def hlam_def
    hlam0_low hlam0_high hGm_bound hGprev_bound
  constructor
  · linarith
  · exact h.2

/-- Mean of a probability distribution supported on `{0, ..., N}`. -/
def finiteMean (N : Nat) (p : Fin (N + 1) -> ℝ) : ℝ :=
  ∑ k, p k * (k.val : ℝ)

/-- Second moment of a probability distribution supported on `{0, ..., N}`. -/
def finiteSecondMoment (N : Nat) (p : Fin (N + 1) -> ℝ) : ℝ :=
  ∑ k, p k * (k.val : ℝ) ^ 2

/-- Variance, written as a central second moment. -/
def finiteVariance (N : Nat) (p : Fin (N + 1) -> ℝ) : ℝ :=
  ∑ k, p k * ((k.val : ℝ) - finiteMean N p) ^ 2

theorem finiteVariance_nonneg
    (N : Nat) (p : Fin (N + 1) -> ℝ)
    (hp_nonneg : ∀ k, 0 <= p k) :
    0 <= finiteVariance N p := by
  unfold finiteVariance
  exact Finset.sum_nonneg (by
    intro k _
    exact mul_nonneg (hp_nonneg k) (sq_nonneg _))

theorem finiteVariance_eq_second_sub_mean_sq
    (N : Nat) (p : Fin (N + 1) -> ℝ)
    (hp_sum : (∑ k, p k) = 1) :
    finiteVariance N p = finiteSecondMoment N p - (finiteMean N p) ^ 2 := by
  set m : ℝ := finiteMean N p
  have hm : m = ∑ k, p k * (k.val : ℝ) := by
    simp [m, finiteMean]
  have hlinear : (∑ x, m * p x * (x.val : ℝ) * 2) = m ^ 2 * 2 := by
    calc
      (∑ x, m * p x * (x.val : ℝ) * 2)
          = m * (∑ x, p x * (x.val : ℝ)) * 2 := by
            rw [Finset.mul_sum, Finset.sum_mul]
            ring_nf
      _ = m ^ 2 * 2 := by
        rw [← hm]
        ring
  have hconst : (∑ x : Fin (N + 1), m ^ 2 * p x) = m ^ 2 := by
    calc
      (∑ x : Fin (N + 1), m ^ 2 * p x)
          = m ^ 2 * (∑ x : Fin (N + 1), p x) := by
            rw [Finset.mul_sum]
      _ = m ^ 2 := by
        rw [hp_sum]
        ring
  unfold finiteVariance finiteSecondMoment
  simp only [finiteMean]
  rw [← hm]
  ring_nf
  simp [Finset.sum_add_distrib, Finset.sum_neg_distrib]
  rw [show (∑ x : Fin (N + 1), p x * (x.val : ℝ) * m * 2) = m ^ 2 * 2 by
    calc
      (∑ x : Fin (N + 1), p x * (x.val : ℝ) * m * 2)
          = (∑ x : Fin (N + 1), m * p x * (x.val : ℝ) * 2) := by
            apply Finset.sum_congr rfl
            intro x _
            ring
      _ = m ^ 2 * 2 := hlinear]
  rw [show (∑ x : Fin (N + 1), p x * m ^ 2) = m ^ 2 by
    calc
      (∑ x : Fin (N + 1), p x * m ^ 2)
          = (∑ x : Fin (N + 1), m ^ 2 * p x) := by
            apply Finset.sum_congr rfl
            intro x _
            ring
      _ = m ^ 2 := hconst]
  ring

theorem finiteMean_nonneg
    (N : Nat) (p : Fin (N + 1) -> ℝ)
    (hp_nonneg : ∀ k, 0 <= p k) :
    0 <= finiteMean N p := by
  unfold finiteMean
  exact Finset.sum_nonneg (by
    intro k _
    exact mul_nonneg (hp_nonneg k) (by exact_mod_cast k.val.zero_le))

theorem finiteMean_le_degree
    (N : Nat) (p : Fin (N + 1) -> ℝ)
    (hp_nonneg : ∀ k, 0 <= p k)
    (hp_sum : (∑ k, p k) = 1) :
    finiteMean N p <= (N : ℝ) := by
  unfold finiteMean
  calc
    (∑ k, p k * (k.val : ℝ))
        <= ∑ k : Fin (N + 1), p k * (N : ℝ) := by
          apply Finset.sum_le_sum
          intro k _
          have hkN_nat : k.val <= N := Nat.le_of_lt_succ k.isLt
          have hkN : (k.val : ℝ) <= (N : ℝ) := by
            exact_mod_cast hkN_nat
          exact mul_le_mul_of_nonneg_left hkN (hp_nonneg k)
    _ = (N : ℝ) := by
      rw [← Finset.sum_mul, hp_sum]
      ring

theorem finiteSecondMoment_le_degree_mul_mean
    (N : Nat) (p : Fin (N + 1) -> ℝ)
    (hp_nonneg : ∀ k, 0 <= p k) :
    finiteSecondMoment N p <= (N : ℝ) * finiteMean N p := by
  unfold finiteSecondMoment finiteMean
  calc
    (∑ k, p k * (k.val : ℝ) ^ 2)
        <= ∑ k : Fin (N + 1), p k * ((N : ℝ) * (k.val : ℝ)) := by
          apply Finset.sum_le_sum
          intro k _
          have hk_nonneg : 0 <= (k.val : ℝ) := by
            exact_mod_cast k.val.zero_le
          have hkN_nat : k.val <= N := Nat.le_of_lt_succ k.isLt
          have hkN : (k.val : ℝ) <= (N : ℝ) := by
            exact_mod_cast hkN_nat
          have hk_sq : (k.val : ℝ) ^ 2 <= (N : ℝ) * (k.val : ℝ) := by
            nlinarith
          exact mul_le_mul_of_nonneg_left hk_sq (hp_nonneg k)
    _ = (N : ℝ) * (∑ k, p k * (k.val : ℝ)) := by
      rw [Finset.mul_sum]
      apply Finset.sum_congr rfl
      intro k _
      ring

/-- Popoviciu's variance bound for a distribution supported on `{0, ..., N}`. -/
theorem finiteVariance_le_degree_sq_div_four
    (N : Nat) (p : Fin (N + 1) -> ℝ)
    (hp_nonneg : ∀ k, 0 <= p k)
    (hp_sum : (∑ k, p k) = 1) :
    finiteVariance N p <= ((N : ℝ) ^ 2) / 4 := by
  have hvar_eq := finiteVariance_eq_second_sub_mean_sq N p hp_sum
  have hsecond := finiteSecondMoment_le_degree_mul_mean N p hp_nonneg
  have hm_nonneg := finiteMean_nonneg N p hp_nonneg
  have hm_le := finiteMean_le_degree N p hp_nonneg hp_sum
  rw [hvar_eq]
  calc
    finiteSecondMoment N p - (finiteMean N p) ^ 2
        <= (N : ℝ) * finiteMean N p - (finiteMean N p) ^ 2 := by
          linarith
    _ <= ((N : ℝ) ^ 2) / 4 := by
      nlinarith [sq_nonneg (finiteMean N p - (N : ℝ) / 2)]

theorem derivative_bound_from_gibbs_variance
    (N : Nat) (mu variance : ℝ -> ℝ) (t : ℝ)
    (ht : (1 / 2 : ℝ) <= t)
    (hvar_nonneg : 0 <= variance t)
    (hvar_le : variance t <= ((N : ℝ) ^ 2) / 4)
    (hderiv : deriv mu t = variance t / t) :
    |deriv mu t| <= ((N : ℝ) ^ 2) / 2 := by
  have ht_pos : 0 < t := by
    linarith
  rw [hderiv, abs_div, abs_of_nonneg hvar_nonneg, abs_of_pos ht_pos]
  exact (div_le_iff₀ ht_pos).2 (by nlinarith [sq_nonneg (N : ℝ)])

theorem abs_sub_le_of_deriv_bound_ordered
    (f : ℝ -> ℝ) (a b C : ℝ)
    (hab : a <= b)
    (hcont : ContinuousOn f (Set.Icc a b))
    (hdiff : DifferentiableOn ℝ f (Set.Ioo a b))
    (hderiv : ∀ x, x ∈ Set.Ioo a b -> |deriv f x| <= C) :
    |f b - f a| <= C * |b - a| := by
  by_cases h_eq : a = b
  · subst b
    simp
  · have hab_lt : a < b := lt_of_le_of_ne hab h_eq
    obtain ⟨c, hc_mem, hc_deriv⟩ := exists_deriv_eq_slope f hab_lt hcont hdiff
    have hb_sub_pos : 0 < b - a := by
      linarith
    have hdiff_eq : f b - f a = deriv f c * (b - a) := by
      rw [hc_deriv]
      field_simp [ne_of_gt hb_sub_pos]
    rw [hdiff_eq, abs_mul, abs_of_pos hb_sub_pos]
    exact mul_le_mul_of_nonneg_right (hderiv c hc_mem) (le_of_lt hb_sub_pos)

/--
Mean-shift bound from the Gibbs derivative identity and the finite-support
variance bound.  This replaces the old opaque Lipschitz-certificate target by
explicit proof obligations: continuity/differentiability of `mu`, the Gibbs
identity `deriv mu t = variance t / t`, and the variance bound on the interval.
-/
theorem mean_shift_bound_from_derivative_identity
    (N : Nat) (lam0 lam : ℝ) (mu variance : ℝ -> ℝ)
    (hlam0_low : (1 / 2 : ℝ) <= lam0)
    (hlam_low : (1 / 2 : ℝ) <= lam)
    (hcont : ContinuousOn mu (Set.uIcc lam0 lam))
    (hdiff : DifferentiableOn ℝ mu (Set.uIoo lam0 lam))
    (hvar_nonneg : ∀ t, t ∈ Set.uIoo lam0 lam -> 0 <= variance t)
    (hvar_le : ∀ t, t ∈ Set.uIoo lam0 lam -> variance t <= ((N : ℝ) ^ 2) / 4)
    (hderiv_id : ∀ t, t ∈ Set.uIoo lam0 lam -> deriv mu t = variance t / t) :
    |mu lam - mu lam0| <= (((N : ℝ) ^ 2) / 2) * |lam - lam0| := by
  by_cases horder : lam0 <= lam
  · have hcont' : ContinuousOn mu (Set.Icc lam0 lam) := by
      simpa [Set.uIcc, horder] using hcont
    have hdiff' : DifferentiableOn ℝ mu (Set.Ioo lam0 lam) := by
      simpa [Set.uIoo, horder] using hdiff
    exact abs_sub_le_of_deriv_bound_ordered mu lam0 lam (((N : ℝ) ^ 2) / 2) horder
      hcont' hdiff'
      (by
        intro t ht
        have ht_low : (1 / 2 : ℝ) <= t := by
          have hlt : lam0 < t := ht.1
          linarith
        exact derivative_bound_from_gibbs_variance N mu variance t ht_low
          (hvar_nonneg t (by simpa [Set.uIoo, horder] using ht))
          (hvar_le t (by simpa [Set.uIoo, horder] using ht))
          (hderiv_id t (by simpa [Set.uIoo, horder] using ht)))
  · have hrev : lam <= lam0 := le_of_not_ge horder
    have hcont' : ContinuousOn mu (Set.Icc lam lam0) := by
      simpa [Set.uIcc, hrev] using hcont
    have hdiff' : DifferentiableOn ℝ mu (Set.Ioo lam lam0) := by
      simpa [Set.uIoo, hrev] using hdiff
    have hbound :
        |mu lam0 - mu lam| <= (((N : ℝ) ^ 2) / 2) * |lam0 - lam| :=
      abs_sub_le_of_deriv_bound_ordered mu lam lam0 (((N : ℝ) ^ 2) / 2) hrev
        hcont' hdiff'
        (by
          intro t ht
          have ht_low : (1 / 2 : ℝ) <= t := by
            have hlt : lam < t := ht.1
            linarith
          exact derivative_bound_from_gibbs_variance N mu variance t ht_low
            (hvar_nonneg t (by simpa [Set.uIoo, hrev] using ht))
            (hvar_le t (by simpa [Set.uIoo, hrev] using ht))
            (hderiv_id t (by simpa [Set.uIoo, hrev] using ht)))
    simpa [abs_sub_comm] using hbound

/--
Mean-shift bound for a finite Gibbs distribution family.  The finite-support
probability hypotheses supply nonnegativity and Popoviciu's variance bound;
the only remaining analytic obligation is the standard Gibbs derivative
identity for the concrete mean function.
-/
theorem mean_shift_bound_from_finite_gibbs_distribution
    (N : Nat) (lam0 lam : ℝ) (mu : ℝ -> ℝ)
    (p : ℝ -> Fin (N + 1) -> ℝ)
    (hlam0_low : (1 / 2 : ℝ) <= lam0)
    (hlam_low : (1 / 2 : ℝ) <= lam)
    (hcont : ContinuousOn mu (Set.uIcc lam0 lam))
    (hdiff : DifferentiableOn ℝ mu (Set.uIoo lam0 lam))
    (hp_nonneg : ∀ t, t ∈ Set.uIoo lam0 lam -> ∀ k, 0 <= p t k)
    (hp_sum : ∀ t, t ∈ Set.uIoo lam0 lam -> (∑ k, p t k) = 1)
    (hderiv_id :
      ∀ t, t ∈ Set.uIoo lam0 lam -> deriv mu t = finiteVariance N (p t) / t) :
    |mu lam - mu lam0| <= (((N : ℝ) ^ 2) / 2) * |lam - lam0| := by
  exact mean_shift_bound_from_derivative_identity N lam0 lam mu
    (fun t => finiteVariance N (p t))
    hlam0_low hlam_low hcont hdiff
    (by
      intro t ht
      exact finiteVariance_nonneg N (p t) (hp_nonneg t ht))
    (by
      intro t ht
      exact finiteVariance_le_degree_sq_div_four N (p t) (hp_nonneg t ht) (hp_sum t ht))
    hderiv_id

/-! ### Gibbs derivative identity for the polynomial-weight family

The block below discharges the `hderiv_id` assumption of
`mean_shift_bound_from_finite_gibbs_distribution` for the concrete
polynomial-weight Gibbs family.

Key structural fact: because Lean's field division satisfies
`(∑ i, f i) / a = ∑ i, f i / a` unconditionally (both sides are `0` when
`a = 0`), `gibbsMean N w` is *globally* equal, as a function, to the rational
function `gibbsA N w / gibbsZ N w`.  Hence `HasDerivAt` at a point `t` needs
only the pointwise hypotheses `t ≠ 0` and `gibbsZ N w t ≠ 0`; no local
nonvanishing or interval hypothesis is required for the identity itself.
-/

/-- Gibbs partition function with polynomial weights. -/
def gibbsZ (N : Nat) (w : Fin (N + 1) -> ℝ) (t : ℝ) : ℝ :=
  ∑ k, w k * t ^ k.val

/-- First-moment numerator `∑ k, k * w k * t ^ k`. -/
def gibbsA (N : Nat) (w : Fin (N + 1) -> ℝ) (t : ℝ) : ℝ :=
  ∑ k, w k * (k.val : ℝ) * t ^ k.val

/-- Second-moment numerator `∑ k, k ^ 2 * w k * t ^ k`. -/
def gibbsB (N : Nat) (w : Fin (N + 1) -> ℝ) (t : ℝ) : ℝ :=
  ∑ k, w k * (k.val : ℝ) ^ 2 * t ^ k.val

/-- Gibbs probability weights `w k * t ^ k / Z(t)`. -/
noncomputable def gibbsProb (N : Nat) (w : Fin (N + 1) -> ℝ) (t : ℝ)
    (k : Fin (N + 1)) : ℝ :=
  (w k * t ^ k.val) / gibbsZ N w t

/-- Mean of the Gibbs distribution. -/
noncomputable def gibbsMean (N : Nat) (w : Fin (N + 1) -> ℝ) (t : ℝ) : ℝ :=
  finiteMean N (gibbsProb N w t)

theorem gibbsProb_sum_eq_one
    (N : Nat) (w : Fin (N + 1) -> ℝ) (t : ℝ)
    (hZ : gibbsZ N w t ≠ 0) :
    (∑ k, gibbsProb N w t k) = 1 := by
  unfold gibbsProb
  rw [← Finset.sum_div]
  show gibbsZ N w t / gibbsZ N w t = 1
  exact div_self hZ

theorem gibbsProb_nonneg
    (N : Nat) (w : Fin (N + 1) -> ℝ) (t : ℝ) (k : Fin (N + 1))
    (hw : 0 <= w k) (ht : 0 <= t) (hZ : 0 <= gibbsZ N w t) :
    0 <= gibbsProb N w t k := by
  unfold gibbsProb
  exact div_nonneg (mul_nonneg hw (pow_nonneg ht _)) hZ

/-- The Gibbs mean is unconditionally the rational function `A / Z`. -/
theorem finiteMean_gibbsProb
    (N : Nat) (w : Fin (N + 1) -> ℝ) (t : ℝ) :
    finiteMean N (gibbsProb N w t) = gibbsA N w t / gibbsZ N w t := by
  unfold finiteMean gibbsProb gibbsA
  rw [Finset.sum_div]
  exact Finset.sum_congr rfl fun k _ => by ring

/-- The Gibbs second moment is unconditionally the rational function `B / Z`. -/
theorem finiteSecondMoment_gibbsProb
    (N : Nat) (w : Fin (N + 1) -> ℝ) (t : ℝ) :
    finiteSecondMoment N (gibbsProb N w t) = gibbsB N w t / gibbsZ N w t := by
  unfold finiteSecondMoment gibbsProb gibbsB
  rw [Finset.sum_div]
  exact Finset.sum_congr rfl fun k _ => by ring

/-- `Z` is a polynomial in `t`; its derivative is the term-by-term one. -/
theorem hasDerivAt_gibbsZ
    (N : Nat) (w : Fin (N + 1) -> ℝ) (t : ℝ) :
    HasDerivAt (fun s => gibbsZ N w s)
      (∑ k : Fin (N + 1), w k * ((k.val : ℝ) * t ^ (k.val - 1))) t := by
  unfold gibbsZ
  exact HasDerivAt.fun_sum fun k _ => (hasDerivAt_pow k.val t).const_mul (w k)

/-- `A` is a polynomial in `t`; its derivative is the term-by-term one. -/
theorem hasDerivAt_gibbsA
    (N : Nat) (w : Fin (N + 1) -> ℝ) (t : ℝ) :
    HasDerivAt (fun s => gibbsA N w s)
      (∑ k : Fin (N + 1), w k * (k.val : ℝ) * ((k.val : ℝ) * t ^ (k.val - 1))) t := by
  unfold gibbsA
  exact HasDerivAt.fun_sum fun k _ =>
    (hasDerivAt_pow k.val t).const_mul (w k * (k.val : ℝ))

/-- `t * Z'(t) = A(t)`, with the `k = 0` term handled by the vanishing
coefficient (no division by `t` and no `t ≠ 0` hypothesis needed). -/
theorem gibbsZ_polyDeriv_mul_self
    (N : Nat) (w : Fin (N + 1) -> ℝ) (t : ℝ) :
    (∑ k : Fin (N + 1), w k * ((k.val : ℝ) * t ^ (k.val - 1))) * t
      = gibbsA N w t := by
  unfold gibbsA
  rw [Finset.sum_mul]
  refine Finset.sum_congr rfl fun k _ => ?_
  rcases Nat.eq_zero_or_pos k.val with h0 | hpos
  · simp [h0]
  · have hp : t ^ (k.val - 1) * t = t ^ k.val := by
      rw [← pow_succ, Nat.sub_add_cancel hpos]
    calc w k * ((k.val : ℝ) * t ^ (k.val - 1)) * t
        = w k * (k.val : ℝ) * (t ^ (k.val - 1) * t) := by ring
      _ = w k * (k.val : ℝ) * t ^ k.val := by rw [hp]

/-- `t * A'(t) = B(t)`, same boundary handling as above. -/
theorem gibbsA_polyDeriv_mul_self
    (N : Nat) (w : Fin (N + 1) -> ℝ) (t : ℝ) :
    (∑ k : Fin (N + 1), w k * (k.val : ℝ) * ((k.val : ℝ) * t ^ (k.val - 1))) * t
      = gibbsB N w t := by
  unfold gibbsB
  rw [Finset.sum_mul]
  refine Finset.sum_congr rfl fun k _ => ?_
  rcases Nat.eq_zero_or_pos k.val with h0 | hpos
  · simp [h0]
  · have hp : t ^ (k.val - 1) * t = t ^ k.val := by
      rw [← pow_succ, Nat.sub_add_cancel hpos]
    calc w k * (k.val : ℝ) * ((k.val : ℝ) * t ^ (k.val - 1)) * t
        = w k * (k.val : ℝ) ^ 2 * (t ^ (k.val - 1) * t) := by ring
      _ = w k * (k.val : ℝ) ^ 2 * t ^ k.val := by rw [hp]

/--
Gibbs derivative identity, `HasDerivAt` form.  Only the pointwise hypotheses
`t ≠ 0` and `gibbsZ N w t ≠ 0` are needed: `gibbsMean N w` equals `A / Z` as a
function on all of `ℝ`, so differentiability of the quotient at `t` transfers
directly, with no neighbourhood or sign hypotheses.
-/
theorem hasDerivAt_gibbsMean
    (N : Nat) (w : Fin (N + 1) -> ℝ) (t : ℝ)
    (ht : t ≠ 0) (hZ : gibbsZ N w t ≠ 0) :
    HasDerivAt (fun s => gibbsMean N w s)
      (finiteVariance N (gibbsProb N w t) / t) t := by
  have hdiv :=
    (hasDerivAt_gibbsA N w t).fun_div (hasDerivAt_gibbsZ N w t) hZ
  have hfun :
      (fun s => gibbsMean N w s) = fun y => gibbsA N w y / gibbsZ N w y := by
    funext s
    unfold gibbsMean
    exact finiteMean_gibbsProb N w s
  have hsum := gibbsProb_sum_eq_one N w t hZ
  have hval :
      finiteVariance N (gibbsProb N w t) / t
        = ((∑ k : Fin (N + 1), w k * (k.val : ℝ) * ((k.val : ℝ) * t ^ (k.val - 1)))
              * gibbsZ N w t
            - gibbsA N w t
              * (∑ k : Fin (N + 1), w k * ((k.val : ℝ) * t ^ (k.val - 1))))
            / gibbsZ N w t ^ 2 := by
    rw [finiteVariance_eq_second_sub_mean_sq N (gibbsProb N w t) hsum,
      finiteSecondMoment_gibbsProb N w t, finiteMean_gibbsProb N w t,
      ← gibbsA_polyDeriv_mul_self N w t, ← gibbsZ_polyDeriv_mul_self N w t]
    field_simp
  rw [hfun, hval]
  exact hdiv

/--
Gibbs derivative identity in `deriv` form: this is the statement assumed as
`hderiv_id` in `mean_shift_bound_from_finite_gibbs_distribution`, now proved
for the concrete polynomial-weight family.
-/
theorem deriv_gibbsMean_eq_variance_div
    (N : Nat) (w : Fin (N + 1) -> ℝ) (t : ℝ)
    (ht : 0 < t) (hZ : gibbsZ N w t ≠ 0) :
    deriv (fun s => gibbsMean N w s) t =
      finiteVariance N (gibbsProb N w t) / t :=
  (hasDerivAt_gibbsMean N w t ht.ne' hZ).deriv

/--
Mean-shift bound for the concrete polynomial-weight Gibbs family, with the
`hderiv_id`, continuity, differentiability, and probability hypotheses of
`mean_shift_bound_from_finite_gibbs_distribution` all discharged.  The only
remaining inputs are nonnegative weights and strict positivity of the
partition function on the closed fugacity interval; the certificate scripts
supply both (weights are independence-polynomial coefficients, and `Z(t) > 0`
follows from `w 0 = 1 > 0` with `t > 0`, or from any positive coefficient).
-/
theorem mean_shift_bound_for_gibbs_family
    (N : Nat) (w : Fin (N + 1) -> ℝ) (lam0 lam : ℝ)
    (hlam0_low : (1 / 2 : ℝ) <= lam0)
    (hlam_low : (1 / 2 : ℝ) <= lam)
    (hw_nonneg : ∀ k, 0 <= w k)
    (hZ_pos : ∀ t, t ∈ Set.uIcc lam0 lam -> 0 < gibbsZ N w t) :
    |gibbsMean N w lam - gibbsMean N w lam0|
      <= (((N : ℝ) ^ 2) / 2) * |lam - lam0| := by
  have hsub : Set.uIoo lam0 lam ⊆ Set.uIcc lam0 lam :=
    Set.uIoo_subset_uIcc_self
  have ht_pos : ∀ t ∈ Set.uIoo lam0 lam, 0 < t := by
    intro t ht
    simp only [Set.uIoo, Set.mem_Ioo] at ht
    have hmin : (1 / 2 : ℝ) <= lam0 ⊓ lam := le_inf hlam0_low hlam_low
    linarith [ht.1]
  have hfun :
      (fun s => gibbsMean N w s) = fun y => gibbsA N w y / gibbsZ N w y := by
    funext s
    unfold gibbsMean
    exact finiteMean_gibbsProb N w s
  refine mean_shift_bound_from_finite_gibbs_distribution N lam0 lam
    (fun s => gibbsMean N w s) (fun t => gibbsProb N w t)
    hlam0_low hlam_low ?_ ?_ ?_ ?_ ?_
  · -- continuity of `A / Z` on the closed interval, where `Z` is positive
    rw [hfun]
    intro s hs
    exact (((hasDerivAt_gibbsA N w s).differentiableAt.continuousAt).div
      ((hasDerivAt_gibbsZ N w s).differentiableAt.continuousAt)
      (hZ_pos s hs).ne').continuousWithinAt
  · -- differentiability on the open interval
    intro t ht
    exact ((hasDerivAt_gibbsMean N w t (ht_pos t ht).ne'
      (hZ_pos t (hsub ht)).ne').differentiableAt).differentiableWithinAt
  · -- nonnegativity of the Gibbs weights
    intro t ht k
    exact gibbsProb_nonneg N w t k (hw_nonneg k) (ht_pos t ht).le
      (hZ_pos t (hsub ht)).le
  · -- normalization
    intro t ht
    exact gibbsProb_sum_eq_one N w t (hZ_pos t (hsub ht)).ne'
  · -- the Gibbs derivative identity, no longer an assumption
    intro t ht
    exact deriv_gibbsMean_eq_variance_div N w t (ht_pos t ht)
      (hZ_pos t (hsub ht)).ne'

/--
Positivity of the partition function from a single positive coefficient.
For independence-polynomial coefficients the constant term is `w 0 = 1 > 0`,
so with nonnegative coefficients and `0 <= t` every term of `Z(t)` is
nonnegative and the `k = 0` term `w 0 * t ^ 0 = w 0` is strictly positive.
-/
theorem gibbsZ_pos
    (N : Nat) (w : Fin (N + 1) -> ℝ) (t : ℝ)
    (hw_nonneg : ∀ k, 0 <= w k)
    (hw0 : 0 < w ⟨0, Nat.succ_pos N⟩)
    (ht : 0 <= t) :
    0 < gibbsZ N w t := by
  unfold gibbsZ
  refine Finset.sum_pos'
    (fun k _ => mul_nonneg (hw_nonneg k) (pow_nonneg ht _)) ?_
  refine ⟨⟨0, Nat.succ_pos N⟩, Finset.mem_univ _, ?_⟩
  show 0 < w ⟨0, Nat.succ_pos N⟩ * t ^ (0 : ℕ)
  rw [pow_zero, mul_one]
  exact hw0

/--
Mean-shift bound for the independence-polynomial situation: the external
`hZ_pos` hypothesis of `mean_shift_bound_for_gibbs_family` is discharged by
`gibbsZ_pos`, using only nonnegative coefficients, a positive constant term,
and the fugacity lower bounds `1/2 <= lam0`, `1/2 <= lam` (which force
`0 <= t` on the closed interval).
-/
theorem mean_shift_bound_for_independence_polynomial
    (N : Nat) (w : Fin (N + 1) -> ℝ) (lam0 lam : ℝ)
    (hlam0_low : (1 / 2 : ℝ) <= lam0)
    (hlam_low : (1 / 2 : ℝ) <= lam)
    (hw_nonneg : ∀ k, 0 <= w k)
    (hw0 : 0 < w ⟨0, Nat.succ_pos N⟩) :
    |gibbsMean N w lam - gibbsMean N w lam0|
      <= (((N : ℝ) ^ 2) / 2) * |lam - lam0| := by
  refine mean_shift_bound_for_gibbs_family N w lam0 lam
    hlam0_low hlam_low hw_nonneg ?_
  intro t ht
  simp only [Set.uIcc, Set.mem_Icc] at ht
  have hmin : (1 / 2 : ℝ) <= lam0 ⊓ lam := le_inf hlam0_low hlam_low
  exact gibbsZ_pos N w t hw_nonneg hw0 (by linarith [ht.1])

/--
Route-2 arithmetic bridge from a reserve at `lam0` and a mean-shift budget.
The checked Gibbs mean-shift theorem supplies the perturbation bound; the
remaining inequality is the certificate reserve comparison.
-/
theorem route2_from_mean_shift_reserve
    (N : Nat) (w : Fin (N + 1) -> ℝ) (lam0 lam m eps : ℝ)
    (hlam0_low : (1 / 2 : ℝ) <= lam0)
    (hlam_low : (1 / 2 : ℝ) <= lam)
    (hw_nonneg : ∀ k, 0 <= w k)
    (hw0 : 0 < w ⟨0, Nat.succ_pos N⟩)
    (hreserve : m - (4 / 3 : ℝ) + eps <= gibbsMean N w lam0)
    (hbudget : (((N : ℝ) ^ 2) / 2) * |lam - lam0| <= eps) :
    m - (3 / 2 : ℝ) <= gibbsMean N w lam := by
  have hshift := mean_shift_bound_for_independence_polynomial N w lam0 lam
    hlam0_low hlam_low hw_nonneg hw0
  have hdist : |gibbsMean N w lam - gibbsMean N w lam0| <= eps := by
    linarith
  have hdiff_lower :
      -eps <= gibbsMean N w lam - gibbsMean N w lam0 :=
    (abs_le.mp hdist).1
  linarith

/--
Abstract fixed-`r` certificate composition (Route-2 shape).

This composes the checked bridge lemmas without modelling the tree family:

* `global_margin_of_adjacent_margins` turns the adjacent hub-off margin
  certificate facts into the global margin;
* `fugacity_shift_bound` turns exact rational coefficient bounds into
  `|lam - lam0| <= 4T` and `1/2 <= lam`;
* `mean_shift_bound_for_independence_polynomial` turns the fugacity shift
  into a mean shift of at most `(N^2 / 2) * 4T = 2N^2T`.

The remaining certificate facts are exactly the quantities the scripts check
with rational arithmetic: adjacent margins and side conditions for `F`,
tie-fugacity coefficient bounds for `G`, reserve at the rational tie fugacity
`lam0`, and the budget inequality `2N^2T <= eps + 1/6`.  The `1/6` slack is
the gap between the reserve constant `4/3` and the Route-2 constant `3/2`.

`w` is the abstract coefficient array of the polynomial whose Gibbs mean is
being controlled.  Linking `w`, `F`, and `G` to the spider family remains a
separate modelling step.
-/
theorem fixed_r_certificate_composition
    (F G : Nat -> Rat) (m N : Nat) (eta T eps lam0 lam : Rat)
    (w : Fin (N + 1) -> ℝ)
    (hleft_margin : eta * F m <= F m - F (m - 1))
    (hright_margin : eta * F m <= F m - F (m + 1))
    (hleft_side : forall k, k < m -> F k <= F (m - 1))
    (hright_side : forall k, m < k -> F k <= F (m + 1))
    (hFm : 0 < F m)
    (hGm_nonneg : 0 <= G m)
    (hGprev_nonneg : 0 <= G (m - 1))
    (hT_nonneg : 0 <= T)
    (hT_half : T <= (1 / 2 : Rat))
    (hlam0_def : lam0 = F (m - 1) / F m)
    (hlam_def : lam = (F (m - 1) + G (m - 1)) / (F m + G m))
    (hlam0_low : (3 / 4 : Rat) <= lam0)
    (hlam0_high : lam0 <= 2)
    (hGm_bound : G m <= T * F m)
    (hGprev_bound : G (m - 1) <= T * F m)
    (hw_nonneg : forall k, 0 <= w k)
    (hw0 : 0 < w ⟨0, Nat.succ_pos N⟩)
    (hreserve : (m : ℝ) - 4 / 3 + (eps : ℝ) <= gibbsMean N w (lam0 : ℝ))
    (hbudget : 2 * (N : Rat) ^ 2 * T <= eps + 1 / 6) :
    (forall k, k ≠ m -> eta * F m <= F m - F k)
      /\ (m : ℝ) - 3 / 2 <= gibbsMean N w (lam : ℝ) := by
  constructor
  · exact global_margin_of_adjacent_margins F m eta
      hleft_margin hright_margin hleft_side hright_side
  · have hshift :=
      fugacity_shift_bound (F m) (F (m - 1)) (G m) (G (m - 1)) lam0 lam T
        hFm hGm_nonneg hGprev_nonneg hT_nonneg hT_half
        hlam0_def hlam_def hlam0_low hlam0_high hGm_bound hGprev_bound
    have habs : |(lam : ℝ) - (lam0 : ℝ)| <= 4 * (T : ℝ) := by
      exact_mod_cast hshift.1
    have hlam_half : (1 / 2 : ℝ) <= (lam : ℝ) := by
      have h : ((1 / 2 : Rat) : ℝ) <= (lam : ℝ) := by
        exact_mod_cast hshift.2
      norm_num at h
      exact h
    have hlam0_half : (1 / 2 : ℝ) <= (lam0 : ℝ) := by
      have hrat : (1 / 2 : Rat) <= lam0 := by linarith
      have h : ((1 / 2 : Rat) : ℝ) <= (lam0 : ℝ) := by
        exact_mod_cast hrat
      norm_num at h
      exact h
    have hmean :=
      mean_shift_bound_for_independence_polynomial N w
        (lam0 : ℝ) (lam : ℝ) hlam0_half hlam_half hw_nonneg hw0
    have hmean_lower := (abs_le.mp hmean).1
    have hmul :
        (((N : ℝ) ^ 2) / 2) * |(lam : ℝ) - (lam0 : ℝ)|
          <= (((N : ℝ) ^ 2) / 2) * (4 * (T : ℝ)) :=
      mul_le_mul_of_nonneg_left habs (by positivity)
    have hmulT :
        (((N : ℝ) ^ 2) / 2) * (4 * (T : ℝ)) = 2 * (N : ℝ) ^ 2 * (T : ℝ) := by
      ring
    have hbudgetR : 2 * (N : ℝ) ^ 2 * (T : ℝ) <= (eps : ℝ) + 1 / 6 := by
      have h : ((2 * (N : Rat) ^ 2 * T : Rat) : ℝ)
          <= ((eps + 1 / 6 : Rat) : ℝ) := by
        exact_mod_cast hbudget
      norm_num at h
      norm_num
      exact h
    linarith

/--
Abstract fixed-`r` certificate composition with the `B` versus `F^-` split
from `fixed_r_certificate_target.tex`.

Here `wB` is the coefficient array for the target polynomial `B`, while
`wMinus` is the coefficient array for the hub-off comparison polynomial
`F^-`.  The reserve is assumed for `wMinus` at `lam0`, and the static
comparison hypothesis transfers it to `wB` at `lam0`.  The checked mean-shift
bound for `wB` then transfers from `lam0` to `lam`.
-/
theorem fixed_r_certificate_composition_split
    (F G : Nat -> Rat) (m N : Nat) (eta T eps delta lam0 lam : Rat)
    (wB wMinus : Fin (N + 1) -> ℝ)
    (hleft_margin : eta * F m <= F m - F (m - 1))
    (hright_margin : eta * F m <= F m - F (m + 1))
    (hleft_side : forall k, k < m -> F k <= F (m - 1))
    (hright_side : forall k, m < k -> F k <= F (m + 1))
    (hG_perturb : ∀ k, k ≠ m -> G k - G m <= eta * F m)
    (hFm : 0 < F m)
    (hGm_nonneg : 0 <= G m)
    (hGprev_nonneg : 0 <= G (m - 1))
    (hT_nonneg : 0 <= T)
    (hT_half : T <= (1 / 2 : Rat))
    (hlam0_def : lam0 = F (m - 1) / F m)
    (hlam_def : lam = (F (m - 1) + G (m - 1)) / (F m + G m))
    (hlam0_low : (3 / 4 : Rat) <= lam0)
    (hlam0_high : lam0 <= 2)
    (hGm_bound : G m <= T * F m)
    (hGprev_bound : G (m - 1) <= T * F m)
    (hwB_nonneg : forall k, 0 <= wB k)
    (hwB0 : 0 < wB ⟨0, Nat.succ_pos N⟩)
    (hreserve_minus :
      (m : ℝ) - 4 / 3 + (eps : ℝ) <= gibbsMean N wMinus (lam0 : ℝ))
    (hstatic :
      |gibbsMean N wB (lam0 : ℝ) - gibbsMean N wMinus (lam0 : ℝ)| <= (delta : ℝ))
    (hbudget : 2 * (N : Rat) ^ 2 * T + delta <= eps + 1 / 6) :
    (forall k, k ≠ m -> eta * F m <= F m - F k)
      /\ (forall k, k ≠ m -> F k + G k <= F m + G m)
      /\ (m : ℝ) - 3 / 2 <= gibbsMean N wB (lam : ℝ) := by
  constructor
  · exact global_margin_of_adjacent_margins F m eta
      hleft_margin hright_margin hleft_side hright_side
  · constructor
    · exact sum_mode_of_adjacent_margins_and_perturbation F G m eta
        hleft_margin hright_margin hleft_side hright_side hG_perturb
    · have hshift :=
        fugacity_shift_bound (F m) (F (m - 1)) (G m) (G (m - 1)) lam0 lam T
          hFm hGm_nonneg hGprev_nonneg hT_nonneg hT_half
          hlam0_def hlam_def hlam0_low hlam0_high hGm_bound hGprev_bound
      have habs : |(lam : ℝ) - (lam0 : ℝ)| <= 4 * (T : ℝ) := by
        exact_mod_cast hshift.1
      have hlam_half : (1 / 2 : ℝ) <= (lam : ℝ) := by
        have h : ((1 / 2 : Rat) : ℝ) <= (lam : ℝ) := by
          exact_mod_cast hshift.2
        norm_num at h
        exact h
      have hlam0_half : (1 / 2 : ℝ) <= (lam0 : ℝ) := by
        have hrat : (1 / 2 : Rat) <= lam0 := by linarith
        have h : ((1 / 2 : Rat) : ℝ) <= (lam0 : ℝ) := by
          exact_mod_cast hrat
        norm_num at h
        exact h
      have hmean :=
        mean_shift_bound_for_independence_polynomial N wB
          (lam0 : ℝ) (lam : ℝ) hlam0_half hlam_half hwB_nonneg hwB0
      have hmean_lower := (abs_le.mp hmean).1
      have hstatic_lower := (abs_le.mp hstatic).1
      have hmul :
          (((N : ℝ) ^ 2) / 2) * |(lam : ℝ) - (lam0 : ℝ)|
            <= (((N : ℝ) ^ 2) / 2) * (4 * (T : ℝ)) :=
        mul_le_mul_of_nonneg_left habs (by positivity)
      have hmulT :
          (((N : ℝ) ^ 2) / 2) * (4 * (T : ℝ)) = 2 * (N : ℝ) ^ 2 * (T : ℝ) := by
        ring
      have hbudgetR :
          2 * (N : ℝ) ^ 2 * (T : ℝ) + (delta : ℝ) <= (eps : ℝ) + 1 / 6 := by
        have h : ((2 * (N : Rat) ^ 2 * T + delta : Rat) : ℝ)
            <= ((eps + 1 / 6 : Rat) : ℝ) := by
          exact_mod_cast hbudget
        norm_num at h
        norm_num
        exact h
      linarith

/--
Exact certificate fields consumed by `fixed_r_certificate_composition`,
bundled as a record.  Constructing a term of this structure from exact
rational data is the certificate check for one instance; the unpacking
theorem `route2_of_certificate` then yields the global-margin and Route-2
conclusions with no further hypotheses.  Field names match the hypothesis
names of `fixed_r_certificate_composition` one-for-one.
-/
structure Route2Certificate where
  N : Nat
  F : Nat -> Rat
  G : Nat -> Rat
  m : Nat
  eta : Rat
  T : Rat
  eps : Rat
  lam0 : Rat
  lam : Rat
  w : Fin (N + 1) -> ℝ
  hleft_margin : eta * F m <= F m - F (m - 1)
  hright_margin : eta * F m <= F m - F (m + 1)
  hleft_side : ∀ k, k < m -> F k <= F (m - 1)
  hright_side : ∀ k, m < k -> F k <= F (m + 1)
  hFm : 0 < F m
  hGm_nonneg : 0 <= G m
  hGprev_nonneg : 0 <= G (m - 1)
  hT_nonneg : 0 <= T
  hT_half : T <= (1 / 2 : Rat)
  hlam0_def : lam0 = F (m - 1) / F m
  hlam_def : lam = (F (m - 1) + G (m - 1)) / (F m + G m)
  hlam0_low : (3 / 4 : Rat) <= lam0
  hlam0_high : lam0 <= 2
  hGm_bound : G m <= T * F m
  hGprev_bound : G (m - 1) <= T * F m
  hw_nonneg : ∀ k, 0 <= w k
  hw0 : 0 < w ⟨0, Nat.succ_pos N⟩
  hreserve : (m : ℝ) - 4 / 3 + (eps : ℝ) <= gibbsMean N w (lam0 : ℝ)
  hbudget : 2 * (N : Rat) ^ 2 * T <= eps + 1 / 6

/--
Unpacking theorem: a `Route2Certificate` instance yields the conclusion of
`fixed_r_certificate_composition` directly.
-/
theorem route2_of_certificate (c : Route2Certificate) :
    (∀ k, k ≠ c.m -> c.eta * c.F c.m <= c.F c.m - c.F k)
      /\ (c.m : ℝ) - 3 / 2 <= gibbsMean c.N c.w (c.lam : ℝ) :=
  fixed_r_certificate_composition c.F c.G c.m c.N c.eta c.T c.eps c.lam0
    c.lam c.w c.hleft_margin c.hright_margin c.hleft_side c.hright_side
    c.hFm c.hGm_nonneg c.hGprev_nonneg c.hT_nonneg c.hT_half c.hlam0_def
    c.hlam_def c.hlam0_low c.hlam0_high c.hGm_bound c.hGprev_bound
    c.hw_nonneg c.hw0 c.hreserve c.hbudget

/--
Certificate fields for the split version of the Route-2 criterion.  This is
closer to `fixed_r_certificate_target.tex`: it records hub-on mode preservation,
the `F^-` reserve, the static comparison between `B` and `F^-` at `lam0`, and
the widened perturbation budget.
-/
structure Route2SplitCertificate where
  N : Nat
  F : Nat -> Rat
  G : Nat -> Rat
  m : Nat
  eta : Rat
  T : Rat
  eps : Rat
  delta : Rat
  lam0 : Rat
  lam : Rat
  wB : Fin (N + 1) -> ℝ
  wMinus : Fin (N + 1) -> ℝ
  hleft_margin : eta * F m <= F m - F (m - 1)
  hright_margin : eta * F m <= F m - F (m + 1)
  hleft_side : ∀ k, k < m -> F k <= F (m - 1)
  hright_side : ∀ k, m < k -> F k <= F (m + 1)
  hG_perturb : ∀ k, k ≠ m -> G k - G m <= eta * F m
  hFm : 0 < F m
  hGm_nonneg : 0 <= G m
  hGprev_nonneg : 0 <= G (m - 1)
  hT_nonneg : 0 <= T
  hT_half : T <= (1 / 2 : Rat)
  hlam0_def : lam0 = F (m - 1) / F m
  hlam_def : lam = (F (m - 1) + G (m - 1)) / (F m + G m)
  hlam0_low : (3 / 4 : Rat) <= lam0
  hlam0_high : lam0 <= 2
  hGm_bound : G m <= T * F m
  hGprev_bound : G (m - 1) <= T * F m
  hwB_nonneg : ∀ k, 0 <= wB k
  hwB0 : 0 < wB ⟨0, Nat.succ_pos N⟩
  hreserve_minus : (m : ℝ) - 4 / 3 + (eps : ℝ) <= gibbsMean N wMinus (lam0 : ℝ)
  hstatic :
    |gibbsMean N wB (lam0 : ℝ) - gibbsMean N wMinus (lam0 : ℝ)| <= (delta : ℝ)
  hbudget : 2 * (N : Rat) ^ 2 * T + delta <= eps + 1 / 6

/--
Unpacking theorem for the split Route-2 certificate.
-/
theorem route2_of_split_certificate (c : Route2SplitCertificate) :
    (∀ k, k ≠ c.m -> c.eta * c.F c.m <= c.F c.m - c.F k)
      /\ (∀ k, k ≠ c.m -> c.F k + c.G k <= c.F c.m + c.G c.m)
      /\ (c.m : ℝ) - 3 / 2 <= gibbsMean c.N c.wB (c.lam : ℝ) :=
  fixed_r_certificate_composition_split c.F c.G c.m c.N c.eta c.T c.eps c.delta
    c.lam0 c.lam c.wB c.wMinus c.hleft_margin c.hright_margin
    c.hleft_side c.hright_side c.hG_perturb c.hFm c.hGm_nonneg
    c.hGprev_nonneg c.hT_nonneg c.hT_half c.hlam0_def c.hlam_def
    c.hlam0_low c.hlam0_high c.hGm_bound c.hGprev_bound c.hwB_nonneg
    c.hwB0 c.hreserve_minus c.hstatic c.hbudget

/-- The Route-2 mean inequality for one indexed instance. -/
def Route2Target (N m : Nat) (wB : Fin (N + 1) -> ℝ) (lam : Rat) : Prop :=
  (m : ℝ) - 3 / 2 <= gibbsMean N wB (lam : ℝ)

/-- The mean-inequality part of a split certificate. -/
theorem route2_target_of_split_certificate (c : Route2SplitCertificate) :
    Route2Target c.N c.m c.wB c.lam :=
  (route2_of_split_certificate c).2.2

/--
Parametrized split certificate for a family whose target data have already
been fixed externally.  This avoids equality casts in the family wrapper:
the certificate scripts can emit a value directly at
`(N a, m a, lam a, wB a)`.
-/
structure Route2SplitCertificateFor
    (N m : Nat) (lam : Rat) (wB : Fin (N + 1) -> ℝ) where
  F : Nat -> Rat
  G : Nat -> Rat
  eta : Rat
  T : Rat
  eps : Rat
  delta : Rat
  lam0 : Rat
  wMinus : Fin (N + 1) -> ℝ
  hleft_margin : eta * F m <= F m - F (m - 1)
  hright_margin : eta * F m <= F m - F (m + 1)
  hleft_side : ∀ k, k < m -> F k <= F (m - 1)
  hright_side : ∀ k, m < k -> F k <= F (m + 1)
  hG_perturb : ∀ k, k ≠ m -> G k - G m <= eta * F m
  hFm : 0 < F m
  hGm_nonneg : 0 <= G m
  hGprev_nonneg : 0 <= G (m - 1)
  hT_nonneg : 0 <= T
  hT_half : T <= (1 / 2 : Rat)
  hlam0_def : lam0 = F (m - 1) / F m
  hlam_def : lam = (F (m - 1) + G (m - 1)) / (F m + G m)
  hlam0_low : (3 / 4 : Rat) <= lam0
  hlam0_high : lam0 <= 2
  hGm_bound : G m <= T * F m
  hGprev_bound : G (m - 1) <= T * F m
  hwB_nonneg : ∀ k, 0 <= wB k
  hwB0 : 0 < wB ⟨0, Nat.succ_pos N⟩
  hreserve_minus : (m : ℝ) - 4 / 3 + (eps : ℝ) <= gibbsMean N wMinus (lam0 : ℝ)
  hstatic :
    |gibbsMean N wB (lam0 : ℝ) - gibbsMean N wMinus (lam0 : ℝ)| <= (delta : ℝ)
  hbudget : 2 * (N : Rat) ^ 2 * T + delta <= eps + 1 / 6

/-- Unpacking theorem for a parametrized split certificate. -/
theorem route2_of_split_certificate_for
    {N m : Nat} {lam : Rat} {wB : Fin (N + 1) -> ℝ}
    (c : Route2SplitCertificateFor N m lam wB) :
    (∀ k, k ≠ m -> c.eta * c.F m <= c.F m - c.F k)
      /\ (∀ k, k ≠ m -> c.F k + c.G k <= c.F m + c.G m)
      /\ Route2Target N m wB lam :=
  fixed_r_certificate_composition_split c.F c.G m N c.eta c.T c.eps c.delta
    c.lam0 lam wB c.wMinus c.hleft_margin c.hright_margin
    c.hleft_side c.hright_side c.hG_perturb c.hFm c.hGm_nonneg
    c.hGprev_nonneg c.hT_nonneg c.hT_half c.hlam0_def c.hlam_def
    c.hlam0_low c.hlam0_high c.hGm_bound c.hGprev_bound c.hwB_nonneg
    c.hwB0 c.hreserve_minus c.hstatic c.hbudget

/-- The mean-inequality part of a parametrized split certificate. -/
theorem route2_target_of_split_certificate_for
    {N m : Nat} {lam : Rat} {wB : Fin (N + 1) -> ℝ}
    (c : Route2SplitCertificateFor N m lam wB) :
    Route2Target N m wB lam :=
  (route2_of_split_certificate_for c).2.2

/-- The Route-2 target for a whole indexed family. -/
def Route2FamilyTarget
    (N m : Nat -> Nat) (lam : Nat -> Rat)
    (wB : (a : Nat) -> Fin (N a + 1) -> ℝ) (a : Nat) : Prop :=
  Route2Target (N a) (m a) (wB a) (lam a)

/--
Family-level threshold split.  The finite prefix is discharged by exact
instance checks; the tail is discharged by a certificate-producing function
`a ↦ Route2SplitCertificateFor`.
-/
theorem route2_family_from_finite_and_tail
    (A : Nat)
    (N m : Nat -> Nat) (lam : Nat -> Rat)
    (wB : (a : Nat) -> Fin (N a + 1) -> ℝ)
    (hfinite :
      ∀ a, 1 <= a -> a < A -> Route2FamilyTarget N m lam wB a)
    (htail :
      ∀ a, A <= a -> Route2SplitCertificateFor (N a) (m a) (lam a) (wB a)) :
    ∀ a, 1 <= a -> Route2FamilyTarget N m lam wB a := by
  intro a ha
  by_cases hlt : a < A
  · exact hfinite a ha hlt
  · exact route2_target_of_split_certificate_for (htail a (Nat.le_of_not_gt hlt))

/--
Bundled fixed-`r` family certificate.  A constructor for this structure is the
abstract Lean shape of the tex criterion: a finite exact prefix plus a
certificate-producing tail.
-/
structure Route2FamilyCertificate where
  A : Nat
  N : Nat -> Nat
  m : Nat -> Nat
  lam : Nat -> Rat
  wB : (a : Nat) -> Fin (N a + 1) -> ℝ
  hfinite :
    ∀ a, 1 <= a -> a < A -> Route2FamilyTarget N m lam wB a
  htail :
    ∀ a, A <= a -> Route2SplitCertificateFor (N a) (m a) (lam a) (wB a)

/-- Unpacking theorem for a bundled fixed-`r` family certificate. -/
theorem route2_of_family_certificate (c : Route2FamilyCertificate) :
    ∀ a, 1 <= a -> Route2FamilyTarget c.N c.m c.lam c.wB a :=
  route2_family_from_finite_and_tail c.A c.N c.m c.lam c.wB c.hfinite c.htail

end AxiomFixedR
