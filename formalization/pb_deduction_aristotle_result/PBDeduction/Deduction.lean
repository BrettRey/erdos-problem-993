import PBDeduction.Propagation
import PBScalar.Statement

/-!
# Section 3: maximal mass, pairwise variance, and the deduction of Theorem 1.1
-/

namespace PBDeduction

open Finset

variable {f : ℤ → ℝ} {n : ℕ}

/-- Unimodality: `c = D - 1` is a mode, so `M = f c = maxMass f n`. -/
theorem maxMass_eq (h : DeductionHyp f n) {D : ℤ} (hD : IsFirstDescent f n D) :
    maxMass f n = f (D - 1) := by
  obtain ⟨hD1, hDn, hdesc, hinc⟩ := hD
  have left : ∀ m : ℕ, 0 ≤ D - 1 - m → f (D - 1 - m) ≤ f (D - 1) := by
    intro m
    induction m with
    | zero => intro _; simp
    | succ m ih =>
      intro hm
      push_cast at hm ⊢
      have := hinc (D - 1 - m) (by omega) (by omega)
      rw [show D - 1 - (m : ℤ) - 1 = D - 1 - (m + 1) by ring] at this
      exact this.trans (ih (by omega))
  have right : ∀ m : ℕ, D + m ≤ n → f (D + m) < f (D + m - 1) ∧ f (D + m) ≤ f (D - 1) := by
    intro m
    induction m with
    | zero => intro _; simp only [Nat.cast_zero, add_zero]; exact ⟨hdesc, hdesc.le⟩
    | succ m ih =>
      intro hm
      push_cast at hm ⊢
      obtain ⟨i1, i2⟩ := ih (by omega)
      have hlc := h.strict_lc (D + m) (by omega) (by omega)
      have p0 := h.pos (D + m) (by omega) (by omega)
      have p1 := h.pos (D + m - 1) (by omega) (by omega)
      rw [show D + ((m : ℤ) + 1) = D + m + 1 by ring, show D + m + 1 - 1 = D + m by ring]
      have : f (D + m - 1) * f (D + m + 1) < f (D + m - 1) * f (D + m) := by
        calc f (D + m - 1) * f (D + m + 1) < f (D + m) ^ 2 := hlc
          _ = f (D + m) * f (D + m) := by ring
          _ < f (D + m - 1) * f (D + m) := by
            rw [mul_comm (f (D + m - 1))]; exact mul_lt_mul_of_pos_left i1 p0
      have h3 := lt_of_mul_lt_mul_left this p1.le
      exact ⟨h3, h3.le.trans i2⟩
  have hle : ∀ k ∈ Icc (0 : ℤ) n, f k ≤ f (D - 1) := by
    intro k hk
    rw [mem_Icc] at hk
    by_cases hkD : k ≤ D - 1
    · obtain ⟨m, hm⟩ : ∃ m : ℕ, k = D - 1 - m := ⟨(D - 1 - k).toNat, by omega⟩
      subst hm; exact left m hk.1
    · obtain ⟨m, hm⟩ : ∃ m : ℕ, k = D + m := ⟨(k - D).toNat, by omega⟩
      subst hm; exact (right m hk.2).2
  exact le_antisymm (Finset.sup'_le _ _ hle) (Finset.le_sup' f (by rw [mem_Icc]; omega))

/-- Restriction of the pairwise variance to a window `{c-K, …, c+K}` inside the
support, with masses bounded below by `M w_i`. -/
theorem pairVar_ge_window (hnn : ∀ k, 0 ≤ f k) (c : ℤ) (K : ℕ) (M : ℝ) (w : ℤ → ℝ)
    (hc0 : 0 ≤ c - K) (hcn : c + K ≤ n) (hM : 0 ≤ M)
    (hw0 : ∀ i ∈ Icc (-(K : ℤ)) K, 0 ≤ w i)
    (hfw : ∀ i ∈ Icc (-(K : ℤ)) K, M * w i ≤ f (c + i)) :
    M ^ 2 * ((1 / 2 : ℝ) * ∑ i ∈ Icc (-(K : ℤ)) K, ∑ j ∈ Icc (-(K : ℤ)) K,
      w i * w j * ((i : ℝ) - (j : ℝ)) ^ 2) ≤ pairVar f n := by
  unfold pairVar
  have hsub : Icc (c - K) (c + K) ⊆ Icc (0 : ℤ) n := by
    intro x hx; rw [mem_Icc] at hx ⊢; omega
  have step1 : ∑ i ∈ Icc (c - K) (c + K), ∑ j ∈ Icc (c - K) (c + K),
      f i * f j * ((i : ℝ) - (j : ℝ)) ^ 2 ≤
      ∑ i ∈ Icc (0 : ℤ) n, ∑ j ∈ Icc (0 : ℤ) n, f i * f j * ((i : ℝ) - (j : ℝ)) ^ 2 := by
    calc _ ≤ ∑ i ∈ Icc (c - K) (c + K), ∑ j ∈ Icc (0 : ℤ) n,
          f i * f j * ((i : ℝ) - (j : ℝ)) ^ 2 := by
          apply sum_le_sum; intro i _
          apply sum_le_sum_of_subset_of_nonneg hsub
          intro j _ _; have := hnn i; have := hnn j; positivity
      _ ≤ _ := by
          apply sum_le_sum_of_subset_of_nonneg hsub
          intro i _ _; apply sum_nonneg; intro j _
          have := hnn i; have := hnn j; positivity
  have hmap : Icc (c - K) (c + K) = (Icc (-(K : ℤ)) K).map (addLeftEmbedding c) := by
    rw [map_add_left_Icc]; congr 1
  have step2 : ∑ i ∈ Icc (c - K) (c + K), ∑ j ∈ Icc (c - K) (c + K),
      f i * f j * ((i : ℝ) - (j : ℝ)) ^ 2 =
      ∑ i ∈ Icc (-(K : ℤ)) K, ∑ j ∈ Icc (-(K : ℤ)) K,
        f (c + i) * f (c + j) * ((i : ℝ) - (j : ℝ)) ^ 2 := by
    rw [hmap, sum_map]
    apply sum_congr rfl; intro i _
    rw [sum_map]
    apply sum_congr rfl; intro j _
    simp only [addLeftEmbedding_apply]
    push_cast; ring
  have step3 : M ^ 2 * ∑ i ∈ Icc (-(K : ℤ)) K, ∑ j ∈ Icc (-(K : ℤ)) K,
      w i * w j * ((i : ℝ) - (j : ℝ)) ^ 2 ≤
      ∑ i ∈ Icc (-(K : ℤ)) K, ∑ j ∈ Icc (-(K : ℤ)) K,
        f (c + i) * f (c + j) * ((i : ℝ) - (j : ℝ)) ^ 2 := by
    rw [mul_sum]
    apply sum_le_sum; intro i hi
    rw [mul_sum]
    apply sum_le_sum; intro j hj
    have h1 := hfw i hi
    have h2 := hfw j hj
    have h3 : 0 ≤ M * w i := mul_nonneg hM (hw0 i hi)
    have h4 : 0 ≤ M * w j := mul_nonneg hM (hw0 j hj)
    have : (M * w i) * (M * w j) ≤ f (c + i) * f (c + j) := mul_le_mul h1 h2 h4 (h3.trans h1)
    calc M ^ 2 * (w i * w j * ((i : ℝ) - (j : ℝ)) ^ 2)
        = (M * w i) * (M * w j) * ((i : ℝ) - (j : ℝ)) ^ 2 := by ring
      _ ≤ _ := mul_le_mul_of_nonneg_right this (sq_nonneg _)
  have := step2 ▸ step1
  nlinarith

/-- Existence of `K`: `(K+1) δ < 1 ≤ (K+2) δ` and `K ≥ 3` for `0 < δ < 1/4`. -/
lemma exists_K {δ : ℝ} (hδ0 : 0 < δ) (hδ4 : δ < 1 / 4) :
    ∃ K : ℕ, ((K : ℝ) + 1) * δ < 1 ∧ 1 ≤ ((K : ℝ) + 2) * δ ∧ 3 ≤ K := by
  set N := ⌈1 / δ⌉₊ with hN
  have h4 : (4 : ℝ) < 1 / δ := by rw [lt_div_iff₀ hδ0]; linarith
  have hN5 : 4 < N := Nat.lt_ceil.mpr (by exact_mod_cast h4)
  have hNle : 1 / δ ≤ N := Nat.le_ceil _
  have hNlt : (N : ℝ) < 1 / δ + 1 := Nat.ceil_lt_add_one (by positivity)
  refine ⟨N - 2, ?_, ?_, by omega⟩
  · rw [Nat.cast_sub (by omega)]; push_cast
    have : (N : ℝ) - 2 + 1 < 1 / δ := by linarith
    rw [lt_div_iff₀ hδ0] at this; linarith
  · rw [Nat.cast_sub (by omega)]; push_cast
    rw [div_le_iff₀ hδ0] at hNle; linarith

/-- **Paper (3.2)**: `pairVar f n ≥ M² A(δ)` with `M = f (D-1)`, under
`0 < δ = δ_D` and `(K+1) δ < 1`, `K ≥ 1`. -/
theorem pairVar_ge_A (h : DeductionHyp f n) {D : ℤ} (hD : IsFirstDescent f n D) {δ : ℝ}
    (hδ : deficit f D = δ) (hδ0 : 0 < δ) {K : ℕ} (hK : ((K : ℝ) + 1) * δ < 1)
    (hK1 : 1 ≤ K) :
    f (D - 1) ^ 2 * PBScalar.A δ K ≤ pairVar f n := by
  have hwb := window_bounds h hD hδ hδ0 hK
  have hRL : ∀ r : ℕ, r ≤ K → 0 ≤ PBScalar.L δ r ∧ PBScalar.L δ r ≤ PBScalar.R δ r := by
    intro r hr
    apply PBScalar.L_nonneg_le_R δ hδ0 r
    have : (r : ℝ) ≤ K := by exact_mod_cast hr
    nlinarith
  unfold PBScalar.A
  apply pairVar_ge_window h.nonneg (D - 1) K (f (D - 1)) (PBScalar.w δ) (by omega) (by omega)
    (h.nonneg _)
  · intro i hi
    rw [mem_Icc] at hi
    unfold PBScalar.w
    split_ifs with h0 h1
    · norm_num
    · have := hRL i.toNat (by omega); linarith
    · exact (hRL (-i).toNat (by omega)).1
  · intro i hi
    rw [mem_Icc] at hi
    unfold PBScalar.w
    split_ifs with h0 h1
    · subst h0; simp
    · have := right_mass_bound h hD hδ hδ0 hK hK1 i.toNat (by omega)
      rwa [show D - 1 + (i.toNat : ℤ) = D - 1 + i by omega] at this
    · have := left_mass_bound h hD hδ hδ0 hK (-i).toNat (by omega)
      rwa [show D - 1 - ((-i).toNat : ℤ) = D - 1 + i by omega] at this

/-- Final scalar step: from `V ≥ M² A`, `A ≥ (3+δ)/(4δ²)` and
`V ≥ (M⁻² - 1)/12`, conclude `V δ ≥ 1/4`. -/
lemma final_algebra {V M A δ : ℝ} (hδ0 : 0 < δ) (hV : 0 ≤ V) (hM : 0 < M)
    (hpv : M ^ 2 * A ≤ V) (hA : (3 + δ) / (4 * δ ^ 2) ≤ A)
    (hmax : (1 / M ^ 2 - 1) / 12 ≤ V) : 1 / 4 ≤ V * δ := by
  by_contra hcon
  push_neg at hcon
  set X := M ^ 2 with hX
  have hX0 : 0 < X := by positivity
  set B := (3 + δ) / (4 * δ ^ 2) with hB
  have hB0 : 0 < B := by positivity
  have hB4 : B * (4 * δ ^ 2) = 3 + δ := by rw [hB]; field_simp
  have h1 : 1 ≤ X * (1 + 12 * V) := by
    have : 1 / X ≤ 1 + 12 * V := by linarith
    rw [div_le_iff₀ hX0] at this; linarith
  have h2 : X * B ≤ V := le_trans (mul_le_mul_of_nonneg_left hA hX0.le) hpv
  have h3 : B ≤ V * (1 + 12 * V) := by
    calc B = B * 1 := by ring
      _ ≤ B * (X * (1 + 12 * V)) := mul_le_mul_of_nonneg_left h1 hB0.le
      _ = (X * B) * (1 + 12 * V) := by ring
      _ ≤ V * (1 + 12 * V) := mul_le_mul_of_nonneg_right h2 (by linarith)
  have h4 : 3 + δ ≤ V * (1 + 12 * V) * (4 * δ ^ 2) := by
    rw [← hB4]; exact mul_le_mul_of_nonneg_right h3 (by positivity)
  have h5 : 0 ≤ V * δ := by positivity
  nlinarith

/-- **G1** (Sections 2–3): the deduction of Theorem 1.1 from Proposition 3.1. -/
theorem deduction_main (h : DeductionHyp f n) {D : ℤ} (hD : IsFirstDescent f n D)
    (hV : 1 ≤ pairVar f n)
    (hmax : (1 / maxMass f n ^ 2 - 1) / 12 ≤ pairVar f n) :
    1 / 4 ≤ pairVar f n * deficit f D := by
  have hδ0 : 0 < deficit f D := h.deficit_pos (by have := hD.1; omega) hD.2.1
  by_cases hδ4 : 1 / 4 ≤ deficit f D
  · nlinarith
  push_neg at hδ4
  obtain ⟨K, hK, hK', hK3⟩ := exists_K hδ0 hδ4
  have hA := PBScalar.scalar_inequality _ K hδ0 hδ4 hK hK'
  have hpv := pairVar_ge_A h hD rfl hδ0 hK (by omega)
  rw [maxMass_eq h hD] at hmax
  have hwb := window_bounds h hD rfl hδ0 hK
  have hM : 0 < f (D - 1) := h.pos _ (by omega) (by omega)
  exact final_algebra hδ0 (by linarith) hM hpv hA hmax

end PBDeduction
