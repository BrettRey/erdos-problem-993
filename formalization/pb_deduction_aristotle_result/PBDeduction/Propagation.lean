import PBDeduction.Basic
import PBScalar.Weights

/-!
# Lemma 2.1 and the mass bounds (2.8)–(2.10)
-/

namespace PBDeduction

open Finset

variable {f : ℤ → ℝ} {n : ℕ}

/-- One step of the propagation: `y (1 - x) ≤ x`, `x ≤ b < 1` give `y ≤ b/(1-b)`. -/
lemma step_bound {x y b : ℝ} (hx : x ≤ b) (hb : b < 1) (hy : y * (1 - x) ≤ x) :
    y ≤ b / (1 - b) := by
  have hx1 : x < 1 := lt_of_le_of_lt hx hb
  have h1 : y ≤ x / (1 - x) := by rw [le_div_iff₀ (by linarith)]; exact hy
  have h2 : x / (1 - x) ≤ b / (1 - b) := by
    rw [div_le_div_iff₀ (by linarith) (by linarith)]; nlinarith
  linarith

lemma bound_succ {δ : ℝ} {r : ℕ} (h1 : 0 < 1 - r * δ) :
    δ / (1 - r * δ) / (1 - δ / (1 - r * δ)) = δ / (1 - (r + 1) * δ) := by
  have : 1 - δ / (1 - r * δ) = (1 - (r + 1) * δ) / (1 - r * δ) := by
    rw [eq_div_iff h1.ne', sub_mul, div_mul_cancel₀ _ h1.ne']; ring
  rw [this, div_div_div_cancel_right₀ h1.ne']

lemma bound_lt_one {δ : ℝ} {r : ℕ} (hδ0 : 0 < δ) (hr : ((r : ℝ) + 1) * δ < 1) :
    δ / (1 - r * δ) < 1 := by
  rw [div_lt_one (by nlinarith)]; nlinarith

/-- **Lemma 2.1** of the paper. Let `D ∈ {1,…,n}` with `δ = δ_D > 0`. For every
`r` with `(r+1) δ < 1` and each sign, `D ± r ∈ {1, …, n-1}` and
`δ_{D±r} ≤ δ/(1 - rδ) < 1`. -/
theorem propagation (h : DeductionHyp f n) {D : ℤ} (hD1 : 1 ≤ D) (hDn : D ≤ n) {δ : ℝ}
    (hδ : deficit f D = δ) (hδ0 : 0 < δ) (r : ℕ) (hr : ((r : ℝ) + 1) * δ < 1) :
    (1 ≤ D + r ∧ D + r ≤ n - 1 ∧ deficit f (D + r) ≤ δ / (1 - r * δ)) ∧
    (1 ≤ D - r ∧ D - r ≤ n - 1 ∧ deficit f (D - r) ≤ δ / (1 - r * δ)) ∧
    δ / (1 - r * δ) < 1 := by
  have hlt := bound_lt_one hδ0 hr
  induction r with
  | zero =>
    simp only [Nat.cast_zero, add_zero, sub_zero, zero_mul, div_one] at hr hlt ⊢
    have hDn' : D ≤ n - 1 := by
      rcases lt_or_eq_of_le hDn with h' | h'
      · omega
      · exfalso; rw [h', h.deficit_top] at hδ; linarith
    exact ⟨⟨hD1, hDn', hδ.le⟩, ⟨hD1, hDn', hδ.le⟩, hlt⟩
  | succ r ih =>
    push_cast at hr hlt ⊢
    have hr' : ((r : ℝ) + 1) * δ < 1 := by nlinarith
    obtain ⟨⟨a1, a2, a3⟩, ⟨b1, b2, b3⟩, c⟩ := ih hr' (bound_lt_one hδ0 hr')
    have hbnd := bound_succ (δ := δ) (r := r) (by nlinarith)
    refine ⟨?_, ?_, hlt⟩
    · have hrec := (deficit_recurrence h a1 a2).2
      have hb := step_bound a3 c hrec
      rw [hbnd] at hb
      have e : D + (r + 1 : ℤ) = D + r + 1 := by ring
      rw [e]
      refine ⟨by omega, ?_, hb⟩
      rcases lt_or_eq_of_le (show D + r + 1 ≤ n by omega) with h' | h'
      · omega
      · exfalso; rw [h', h.deficit_top] at hb; linarith
    · have hrec := (deficit_recurrence h (k := D - r) b1 b2).1
      have hb := step_bound b3 c hrec
      rw [hbnd] at hb
      have e : D - (r + 1 : ℤ) = D - r - 1 := by ring
      rw [e]
      have hD0 : 0 ≤ D - r - 1 := by omega
      refine ⟨?_, by omega, hb⟩
      rcases lt_or_eq_of_le hD0 with h' | h'
      · omega
      · exfalso; rw [← h', h.deficit_zero] at hb; linarith

/-! ## Mass bounds -/

/-- Adjacent mass ratio `q_k = f k / f (k-1)`. -/
noncomputable def ratio (f : ℤ → ℝ) (k : ℤ) : ℝ := f k / f (k - 1)

lemma ratio_succ (f : ℤ → ℝ) (k : ℤ) (h0 : 0 < f (k - 1)) (h1 : 0 < f k) :
    ratio f (k + 1) = ratio f k * (1 - deficit f k) := by
  unfold ratio deficit
  rw [show k + 1 - 1 = k by ring]
  field_simp
  ring

lemma one_sub_bound {δ : ℝ} {s : ℕ} (h1 : 0 < 1 - s * δ) :
    1 - δ / (1 - s * δ) = (1 - (s + 1) * δ) / (1 - s * δ) := by
  rw [eq_div_iff h1.ne', sub_mul, div_mul_cancel₀ _ h1.ne']; ring

lemma a_nonneg_of {δ : ℝ} {K : ℕ} (hδ0 : 0 < δ) (hK : ((K : ℝ) + 1) * δ < 1) (hK1 : 1 ≤ K) :
    0 ≤ PBScalar.a δ := by
  unfold PBScalar.a
  have : (1 : ℝ) ≤ K := by exact_mod_cast hK1
  exact div_nonneg (by nlinarith) (by nlinarith)

section Window

variable (h : DeductionHyp f n) {D : ℤ} (hD : IsFirstDescent f n D) {δ : ℝ}
  (hδ : deficit f D = δ) (hδ0 : 0 < δ) {K : ℕ} (hK : ((K : ℝ) + 1) * δ < 1)

include h hD hδ hδ0 hK

lemma window_bounds : 1 ≤ D - K ∧ D + K ≤ n - 1 :=
  ⟨((propagation h hD.1 hD.2.1 hδ hδ0 K hK).2.1).1,
    ((propagation h hD.1 hD.2.1 hδ hδ0 K hK).1).2.1⟩

lemma window_pos (i : ℤ) (h1 : D - K - 1 ≤ i) (h2 : i ≤ D + K + 1) : 0 < f i := by
  have := window_bounds h hD hδ hδ0 hK
  exact h.pos i (by omega) (by omega)

lemma deficit_window_bound (s : ℕ) (hs : s ≤ K) :
    deficit f (D + s) ≤ δ / (1 - s * δ) ∧ deficit f (D - s) ≤ δ / (1 - s * δ) := by
  have hsK : ((s : ℝ) + 1) * δ < 1 := by
    have : (s : ℝ) ≤ K := by exact_mod_cast hs
    nlinarith
  have := propagation h hD.1 hD.2.1 hδ hδ0 s hsK
  exact ⟨this.1.2.2, this.2.1.2.2⟩

/-- `1 - δ_{D±s} ≥ (1-(s+1)δ)/(1-sδ) ≥ 0` for `s ≤ K`. -/
lemma one_sub_deficit_window (s : ℕ) (hs : s ≤ K) :
    (1 - (s + 1) * δ) / (1 - s * δ) ≤ 1 - deficit f (D + s) ∧
    (1 - (s + 1) * δ) / (1 - s * δ) ≤ 1 - deficit f (D - s) ∧
    0 ≤ (1 - (s + 1) * δ) / (1 - s * δ) := by
  have hsK : ((s : ℝ) + 1) * δ < 1 := by
    have : (s : ℝ) ≤ K := by exact_mod_cast hs
    nlinarith
  have h1 : 0 < 1 - s * δ := by nlinarith
  have := deficit_window_bound h hD hδ hδ0 hK s hs
  rw [← one_sub_bound h1]
  refine ⟨by linarith [this.1], by linarith [this.2], ?_⟩
  rw [one_sub_bound h1]
  exact div_nonneg (by linarith) h1.le

/-- **Paper (2.8)**: `q_D ≥ a` (with `q_{D-1} ≥ 1` from the first descent). -/
theorem ratio_D_ge (hK1 : 1 ≤ K) : PBScalar.a δ ≤ ratio f D := by
  have hwb := window_bounds h hD hδ hδ0 hK
  have p2 := window_pos h hD hδ hδ0 hK (D - 1 - 1) (by omega) (by omega)
  have p1 := window_pos h hD hδ hδ0 hK (D - 1) (by omega) (by omega)
  have hrs := ratio_succ f (D - 1) p2 p1
  rw [show D - 1 + 1 = D by ring] at hrs
  rw [hrs]
  have hq : 1 ≤ ratio f (D - 1) := by
    unfold ratio
    rw [le_div_iff₀ p2, one_mul]
    exact hD.2.2.2 (D - 1) (by omega) (by omega)
  have hb := (one_sub_deficit_window h hD hδ hδ0 hK 1 hK1).2
  simp only [Nat.cast_one] at hb
  have ha : PBScalar.a δ = (1 - (1 + 1) * δ) / (1 - 1 * δ) := by
    unfold PBScalar.a; ring_nf
  rw [ha]
  calc (1 - (1 + 1) * δ) / (1 - 1 * δ) = 1 * ((1 - (1 + 1) * δ) / (1 - 1 * δ)) := by ring
    _ ≤ ratio f (D - 1) * (1 - deficit f (D - 1)) :=
        mul_le_mul hq (by simpa using hb.1) hb.2 (by linarith)

/-- Right ratio bound `q_{D+j} ≥ a (1 - jδ)` for `0 ≤ j ≤ K - 1`. -/
theorem ratio_right_ge (hK1 : 1 ≤ K) (j : ℕ) (hj : j + 1 ≤ K) :
    PBScalar.a δ * (1 - j * δ) ≤ ratio f (D + j) := by
  have ha0 : 0 ≤ PBScalar.a δ := a_nonneg_of hδ0 hK hK1
  have hwb := window_bounds h hD hδ hδ0 hK
  induction j with
  | zero => simpa using ratio_D_ge h hD hδ hδ0 hK hK1
  | succ j ih =>
    have ih := ih (by omega)
    push_cast
    have p0 := window_pos h hD hδ hδ0 hK (D + j - 1) (by omega) (by omega)
    have p1 := window_pos h hD hδ hδ0 hK (D + j) (by omega) (by omega)
    have hrs := ratio_succ f (D + j) p0 p1
    rw [show D + ((j : ℤ) + 1) = D + j + 1 by ring, hrs]
    have hb := one_sub_deficit_window h hD hδ hδ0 hK j (by omega)
    have hjK : ((j : ℝ) + 1 + 1) * δ < 1 := by
      have : (j : ℝ) + 1 ≤ K := by exact_mod_cast (show j + 1 ≤ K by omega)
      nlinarith
    have h1 : 0 < 1 - j * δ := by nlinarith
    calc PBScalar.a δ * (1 - (j + 1) * δ)
        = PBScalar.a δ * (1 - j * δ) * ((1 - (j + 1) * δ) / (1 - j * δ)) := by
          field_simp
      _ ≤ ratio f (D + j) * (1 - deficit f (D + j)) :=
          mul_le_mul ih hb.1 hb.2.2 (le_trans (mul_nonneg ha0 h1.le) ih)

/-- Left ratio bound: `q_{D-j} (1-(j+1)δ)/(1-δ) ≤ q_D` for `0 ≤ j ≤ K`. -/
theorem ratio_left_le (j : ℕ) (hj : j ≤ K) :
    ratio f (D - j) * ((1 - (j + 1) * δ) / (1 - δ)) ≤ ratio f D := by
  have hwb := window_bounds h hD hδ hδ0 hK
  have hδ1 : 0 < 1 - δ := by
    have : (0 : ℝ) ≤ K := by positivity
    nlinarith
  induction j with
  | zero => simp [div_self hδ1.ne']
  | succ j ih =>
    have ih := ih (by omega)
    push_cast
    have p0 := window_pos h hD hδ hδ0 hK (D - j - 1 - 1) (by omega) (by omega)
    have p1 := window_pos h hD hδ hδ0 hK (D - j - 1) (by omega) (by omega)
    have hrs := ratio_succ f (D - j - 1) p0 p1
    rw [show D - j - 1 + 1 = D - j by ring] at hrs
    have hb := (one_sub_deficit_window h hD hδ hδ0 hK (j + 1) hj)
    push_cast at hb
    rw [show D - ((j : ℤ) + 1) = D - j - 1 by ring] at hb ⊢
    have hjK : ((j : ℝ) + 1 + 1) * δ < 1 := by
      have : (j : ℝ) + 1 ≤ K := by exact_mod_cast hj
      nlinarith
    have h1 : 0 < 1 - (j + 1) * δ := by nlinarith
    have hq0 : 0 ≤ ratio f (D - j - 1) := by unfold ratio; positivity
    have hc0 : 0 ≤ (1 - (j + 1) * δ) / (1 - δ) := div_nonneg h1.le hδ1.le
    calc ratio f (D - j - 1) * ((1 - (j + 1 + 1) * δ) / (1 - δ))
        = ratio f (D - j - 1) * ((1 - (j + 1 + 1) * δ) / (1 - (j + 1) * δ)) *
            ((1 - (j + 1) * δ) / (1 - δ)) := by
          field_simp
      _ ≤ ratio f (D - j - 1) * (1 - deficit f (D - j - 1)) *
            ((1 - (j + 1) * δ) / (1 - δ)) := by
          apply mul_le_mul_of_nonneg_right _ hc0
          exact mul_le_mul_of_nonneg_left hb.2.1 hq0
      _ = ratio f (D - j) * ((1 - (j + 1) * δ) / (1 - δ)) := by rw [hrs]
      _ ≤ ratio f D := ih

/-- **Paper (2.9)** (right mass bound): `f (c + r) ≥ M R_r` for `0 ≤ r ≤ K`,
where `c = D - 1` and `M = f c`. -/
theorem right_mass_bound (hK1 : 1 ≤ K) (r : ℕ) (hr : r ≤ K) :
    f (D - 1) * PBScalar.R δ r ≤ f (D - 1 + r) := by
  have hwb := window_bounds h hD hδ hδ0 hK
  induction r with
  | zero => simp [PBScalar.R_zero]
  | succ r ih =>
    have ih := ih (by omega)
    have hq := ratio_right_ge h hD hδ hδ0 hK hK1 r (by omega)
    have p0 := window_pos h hD hδ hδ0 hK (D + r - 1) (by omega) (by omega)
    have hR0 : 0 ≤ PBScalar.R δ r := by
      have : ((r : ℝ) + 1) * δ ≤ 1 := by
        have : (r : ℝ) ≤ K := by exact_mod_cast (show r ≤ K by omega)
        nlinarith
      exact le_trans (PBScalar.L_nonneg_le_R δ hδ0 r this).1 (PBScalar.L_nonneg_le_R δ hδ0 r this).2
    have hM0 : 0 ≤ f (D - 1) := h.nonneg _
    have e1 : D - 1 + ((r + 1 : ℕ) : ℤ) = D + r := by push_cast; ring
    have e2 : D - 1 + (r : ℤ) = D + r - 1 := by ring
    rw [e1, PBScalar.R_succ]
    rw [e2] at ih
    have : f (D + r) = ratio f (D + r) * f (D + r - 1) := by
      unfold ratio; field_simp
    rw [this]
    calc f (D - 1) * (PBScalar.R δ r * (PBScalar.a δ * (1 - r * δ)))
        = (f (D - 1) * PBScalar.R δ r) * (PBScalar.a δ * (1 - r * δ)) := by ring
      _ ≤ f (D + r - 1) * ratio f (D + r) := by
          apply mul_le_mul ih hq
          · have ha0 : 0 ≤ PBScalar.a δ := a_nonneg_of hδ0 hK hK1
            have : (r : ℝ) + 1 ≤ K := by exact_mod_cast (show r + 1 ≤ K by omega)
            exact mul_nonneg ha0 (by nlinarith)
          · exact p0.le
      _ = ratio f (D + r) * f (D + r - 1) := by ring

/-- **Paper (2.10)** (left mass bound): `f (c - r) ≥ M L_r` for `0 ≤ r ≤ K`,
where `c = D - 1` and `M = f c`. -/
theorem left_mass_bound (r : ℕ) (hr : r ≤ K) :
    f (D - 1) * PBScalar.L δ r ≤ f (D - 1 - r) := by
  have hwb := window_bounds h hD hδ hδ0 hK
  have hδ1 : 0 < 1 - δ := by
    have : (0 : ℝ) ≤ K := by positivity
    nlinarith
  have hDD : f D < f (D - 1) := hD.2.2.1
  induction r with
  | zero => simp [PBScalar.L_zero]
  | succ r ih =>
    have ih := ih (by omega)
    have hq := ratio_left_le h hD hδ hδ0 hK (r + 1) hr
    push_cast at hq
    have pD := window_pos h hD hδ hδ0 hK (D - 1) (by omega) (by omega)
    have p0 := window_pos h hD hδ hδ0 hK (D - r - 1) (by omega) (by omega)
    have p1 := window_pos h hD hδ hδ0 hK (D - r - 1 - 1) (by omega) (by omega)
    have hqD : ratio f D < 1 := by
      unfold ratio; rw [div_lt_one pD]; exact hDD
    have hrK : ((r : ℝ) + 1 + 1) * δ < 1 := by
      have : (r : ℝ) + 1 ≤ K := by exact_mod_cast hr
      nlinarith
    have hc0 : 0 ≤ (1 - (r + 1 + 1) * δ) / (1 - δ) := div_nonneg (by nlinarith) hδ1.le
    -- `f (D-r-2) ≥ f (D-r-1) * c_{r+1}`
    have key : f (D - r - 1) * ((1 - (r + 1 + 1) * δ) / (1 - δ)) ≤ f (D - r - 1 - 1) := by
      rw [show D - ((r : ℤ) + 1) = D - r - 1 by ring] at hq
      unfold ratio at hq
      have : f (D - r - 1) / f (D - r - 1 - 1) * ((1 - (r + 1 + 1) * δ) / (1 - δ)) < 1 :=
        lt_of_le_of_lt hq hqD
      rw [div_mul_eq_mul_div, div_lt_one p1] at this
      exact this.le
    have hL0 : 0 ≤ PBScalar.L δ r := by
      have : ((r : ℝ) + 1) * δ ≤ 1 := by nlinarith
      exact (PBScalar.L_nonneg_le_R δ hδ0 r this).1
    have hM0 : 0 ≤ f (D - 1) := h.nonneg _
    rw [PBScalar.L_succ]
    have e1 : D - 1 - ((r + 1 : ℕ) : ℤ) = D - r - 1 - 1 := by push_cast; ring
    have e2 : D - 1 - (r : ℤ) = D - r - 1 := by ring
    rw [e1]
    rw [e2] at ih
    calc f (D - 1) * (PBScalar.L δ r * ((1 - ((r : ℝ) + 2) * δ) / (1 - δ)))
        = (f (D - 1) * PBScalar.L δ r) * ((1 - (r + 1 + 1) * δ) / (1 - δ)) := by ring
      _ ≤ f (D - r - 1) * ((1 - (r + 1 + 1) * δ) / (1 - δ)) :=
          mul_le_mul_of_nonneg_right ih hc0
      _ ≤ f (D - r - 1 - 1) := key

end Window

end PBDeduction
