import PBScalar.Defs

/-!
# Elementary facts about the weights `R r` and `L r`

* one-step recursions for `R` and `L`;
* `L r = ∏_{s=1}^r (1 - s/H)` with `H = 1/δ - 1` (paper, eq. `b-def`);
* `R r ≥ L r ≥ 0` for `(r+1) δ ≤ 1` (paper, eq. `C-ratio`);
* monotonicity of `R r`, `L r` in `δ`;
* the Weierstrass product inequality.
-/

namespace PBScalar

open Finset

lemma R_zero (δ : ℝ) : R δ 0 = 1 := by simp [R]

lemma R_succ (δ : ℝ) (r : ℕ) : R δ (r + 1) = R δ r * (a δ * (1 - (r : ℝ) * δ)) := by
  unfold R
  rcases Nat.eq_zero_or_pos r with rfl | hr
  · simp
  · rw [Finset.prod_Ico_succ_top (by omega : 1 ≤ r)]
    ring

lemma L_zero (δ : ℝ) : L δ 0 = 1 := by simp [L]

lemma L_succ (δ : ℝ) (r : ℕ) :
    L δ (r + 1) = L δ r * ((1 - ((r : ℝ) + 2) * δ) / (1 - δ)) := by
  unfold L
  rw [Finset.prod_Icc_succ_top (by omega : 2 ≤ r + 1 + 1)]
  push_cast
  rw [pow_succ, mul_inv]
  ring

/-- **Shared lemma** (paper, eq. `b-def`): `L r = ∏_{s=1}^r (1 - s/H)` where `H = 1/δ - 1`. -/
theorem L_eq_prod (δ : ℝ) (hδ0 : 0 < δ) (hδ1 : δ < 1) (r : ℕ) :
    L δ r = ∏ s ∈ Finset.Icc 1 r, (1 - (s : ℝ) / (1 / δ - 1)) := by
  induction r with
  | zero => simp [L_zero]
  | succ r ih =>
    rw [L_succ, ih, Finset.prod_Icc_succ_top (by omega : 1 ≤ r + 1)]
    congr 1
    have h1 : (1 - δ) ≠ 0 := by linarith
    have h2 : 1 / δ - 1 = (1 - δ) / δ := by field_simp
    rw [h2]
    field_simp
    push_cast
    ring

/-- `R r ≥ L r ≥ 0` whenever `(r+1) δ ≤ 1`. -/
theorem L_nonneg_le_R (δ : ℝ) (hδ0 : 0 < δ) (r : ℕ) (hr : ((r : ℝ) + 1) * δ ≤ 1) :
    0 ≤ L δ r ∧ L δ r ≤ R δ r := by
  induction r with
  | zero => simp [L_zero, R_zero]
  | succ r ih =>
    push_cast at hr
    have hr' : ((r : ℝ) + 1) * δ ≤ 1 := by nlinarith
    obtain ⟨h0, hle⟩ := ih hr'
    have hδ1 : 0 < 1 - δ := by nlinarith
    have hy : 0 ≤ (1 - ((r : ℝ) + 2) * δ) / (1 - δ) := by
      apply div_nonneg _ hδ1.le; nlinarith
    have hyx : (1 - ((r : ℝ) + 2) * δ) / (1 - δ) ≤ a δ * (1 - (r : ℝ) * δ) := by
      rw [a, div_mul_eq_mul_div, div_le_div_iff_of_pos_right hδ1]
      have : (0:ℝ) ≤ r := r.cast_nonneg
      nlinarith [mul_nonneg this (mul_nonneg hδ0.le hδ0.le)]
    rw [L_succ, R_succ]
    refine ⟨mul_nonneg h0 hy, ?_⟩
    exact mul_le_mul hle hyx hy (h0.trans hle)

/-- **Shared lemma** (paper, Section 4): `R r ≥ L r` for `1 ≤ r ≤ K`. -/
theorem R_ge_L (δ : ℝ) (K : ℕ) (hδ0 : 0 < δ) (hK : ((K : ℝ) + 1) * δ < 1) (r : ℕ)
    (hr : r ≤ K) : L δ r ≤ R δ r := by
  refine (L_nonneg_le_R δ hδ0 r ?_).2
  have : (r : ℝ) ≤ K := by exact_mod_cast hr
  nlinarith

/-- Monotonicity in `δ`: if `0 < δ ≤ δ₁` and `(r+1) δ₁ ≤ 1` then
`0 ≤ L δ₁ r ≤ L δ r` and `0 ≤ R δ₁ r ≤ R δ r`. -/
theorem LR_antitone (δ δ₁ : ℝ) (hδ0 : 0 < δ) (hδδ₁ : δ ≤ δ₁) (r : ℕ)
    (hr : ((r : ℝ) + 1) * δ₁ ≤ 1) :
    (0 ≤ L δ₁ r ∧ L δ₁ r ≤ L δ r) ∧ (0 ≤ R δ₁ r ∧ R δ₁ r ≤ R δ r) := by
  induction r with
  | zero => simp [L_zero, R_zero]
  | succ r ih =>
    push_cast at hr
    have hr' : ((r : ℝ) + 1) * δ₁ ≤ 1 := by nlinarith
    obtain ⟨⟨hL0, hL⟩, ⟨hR0, hR⟩⟩ := ih hr'
    have hrn : (0:ℝ) ≤ r := r.cast_nonneg
    have hδ1 : 0 < 1 - δ₁ := by nlinarith
    have hδ1' : 0 < 1 - δ := by linarith
    -- left factors
    have hy0 : 0 ≤ (1 - ((r : ℝ) + 2) * δ₁) / (1 - δ₁) := by
      apply div_nonneg _ hδ1.le; nlinarith
    have hy : (1 - ((r : ℝ) + 2) * δ₁) / (1 - δ₁) ≤ (1 - ((r : ℝ) + 2) * δ) / (1 - δ) := by
      rw [div_le_div_iff₀ hδ1 hδ1']
      nlinarith
    -- right factors
    have ha0 : 0 ≤ a δ₁ := by
      unfold a; apply div_nonneg _ hδ1.le; nlinarith
    have ha : a δ₁ ≤ a δ := by
      unfold a; rw [div_le_div_iff₀ hδ1 hδ1']; nlinarith
    have hx0 : 0 ≤ 1 - (r : ℝ) * δ₁ := by nlinarith
    have hx : 1 - (r : ℝ) * δ₁ ≤ 1 - (r : ℝ) * δ := by nlinarith
    rw [L_succ, L_succ, R_succ, R_succ]
    refine ⟨⟨mul_nonneg hL0 hy0, mul_le_mul hL hy hy0 (hL0.trans hL)⟩,
      ⟨mul_nonneg hR0 (mul_nonneg ha0 hx0), ?_⟩⟩
    exact mul_le_mul hR (mul_le_mul ha hx hx0 (ha0.trans ha)) (mul_nonneg ha0 hx0)
      (hR0.trans hR)

/-- **Shared lemma**: Weierstrass product inequality `∏ (1 - x_s) ≥ 1 - ∑ x_s`
for `x_s ∈ [0, 1]`. -/
theorem weierstrass_prod {ι : Type*} (s : Finset ι) (x : ι → ℝ)
    (h0 : ∀ i ∈ s, 0 ≤ x i) (h1 : ∀ i ∈ s, x i ≤ 1) :
    1 - ∑ i ∈ s, x i ≤ ∏ i ∈ s, (1 - x i) := by
  classical
  induction s using Finset.induction_on with
  | empty => simp
  | insert j s hj ih =>
    rw [Finset.sum_insert hj, Finset.prod_insert hj]
    have ih' := ih (fun i hi => h0 i (Finset.mem_insert_of_mem hi))
      (fun i hi => h1 i (Finset.mem_insert_of_mem hi))
    have hxj0 := h0 j (Finset.mem_insert_self j s)
    have hxj1 := h1 j (Finset.mem_insert_self j s)
    have hS : 0 ≤ ∑ i ∈ s, x i := Finset.sum_nonneg (fun i hi => h0 i (Finset.mem_insert_of_mem hi))
    have hP : 0 ≤ ∏ i ∈ s, (1 - x i) := Finset.prod_nonneg
      (fun i hi => by linarith [h1 i (Finset.mem_insert_of_mem hi)])
    nlinarith [mul_le_mul_of_nonneg_left ih' (by linarith : (0:ℝ) ≤ 1 - x j)]

lemma sum_Icc_cast_div (H : ℝ) (r : ℕ) :
    ∑ s ∈ Finset.Icc 1 r, (s : ℝ) / H = (r : ℝ) * (r + 1) / (2 * H) := by
  induction r with
  | zero => simp
  | succ r ih =>
    rw [Finset.sum_Icc_succ_top (by omega), ih]
    push_cast
    ring

/-- Paper, eq. `bonferroni`: `L r ≥ 1 - r(r+1)/(2H)` when `r ≤ H = 1/δ - 1`. -/
theorem L_ge_lambda (δ : ℝ) (hδ0 : 0 < δ) (hδ1 : δ < 1) (r : ℕ)
    (hr : (r : ℝ) ≤ 1 / δ - 1) :
    1 - (r : ℝ) * (r + 1) / (2 * (1 / δ - 1)) ≤ L δ r := by
  rw [L_eq_prod δ hδ0 hδ1, ← sum_Icc_cast_div]
  have hH : 0 < 1 / δ - 1 := by
    rw [sub_pos, lt_div_iff₀ hδ0]; linarith
  apply weierstrass_prod
  · intro i _; positivity
  · intro i hi
    rw [div_le_one hH]
    have : i ≤ r := (Finset.mem_Icc.mp hi).2
    have : (i : ℝ) ≤ r := by exact_mod_cast this
    linarith

end PBScalar
