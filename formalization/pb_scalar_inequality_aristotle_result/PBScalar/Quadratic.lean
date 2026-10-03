import PBScalar.Weights

/-!
# The quadratic form `A`

* `Aw v K = (1/2) ∑_{i,j=-K}^{K} v i v j (i-j)^2` for an arbitrary weight function `v`;
* monotonicity of `Aw` in nonnegative weights;
* the pairwise identity `Aw v K = (∑ v)(∑ i² v) - (∑ i v)²`;
* consequences for `A δ K`: an explicit formula, the symmetrisation bound
  `A δ K ≥ S·T`, and antitonicity in `δ` for fixed `K`.
-/

namespace PBScalar

open Finset

/-- The quadratic form for an arbitrary weight function `v : ℤ → ℝ`. -/
noncomputable def Aw (v : ℤ → ℝ) (K : ℕ) : ℝ :=
  (1 / 2 : ℝ) *
    ∑ i ∈ Finset.Icc (-(K : ℤ)) (K : ℤ), ∑ j ∈ Finset.Icc (-(K : ℤ)) (K : ℤ),
      v i * v j * ((i : ℝ) - (j : ℝ)) ^ 2

lemma A_eq_Aw (δ : ℝ) (K : ℕ) : A δ K = Aw (w δ) K := rfl

/-- **Shared lemma**: `A` is nondecreasing in nonnegative weights. -/
theorem Aw_mono (v u : ℤ → ℝ) (K : ℕ)
    (hv : ∀ i ∈ Finset.Icc (-(K : ℤ)) (K : ℤ), 0 ≤ v i)
    (hvu : ∀ i ∈ Finset.Icc (-(K : ℤ)) (K : ℤ), v i ≤ u i) :
    Aw v K ≤ Aw u K := by
  unfold Aw
  apply mul_le_mul_of_nonneg_left _ (by norm_num)
  apply Finset.sum_le_sum
  intro i hi
  apply Finset.sum_le_sum
  intro j hj
  apply mul_le_mul_of_nonneg_right _ (sq_nonneg _)
  exact mul_le_mul (hvu i hi) (hvu j hj) (hv j hj) ((hv i hi).trans (hvu i hi))

/-- **Shared lemma** (pairwise identity):
`(1/2)∑∑ v_i v_j (i-j)^2 = (∑ v)(∑ i² v) - (∑ i v)²`. -/
theorem Aw_pairwise (v : ℤ → ℝ) (K : ℕ) :
    Aw v K = (∑ i ∈ Finset.Icc (-(K : ℤ)) (K : ℤ), v i) *
        (∑ i ∈ Finset.Icc (-(K : ℤ)) (K : ℤ), (i : ℝ) ^ 2 * v i) -
      (∑ i ∈ Finset.Icc (-(K : ℤ)) (K : ℤ), (i : ℝ) * v i) ^ 2 := by
  unfold Aw
  set s := Finset.Icc (-(K : ℤ)) (K : ℤ)
  have h : ∀ i j : ℤ, v i * v j * ((i : ℝ) - (j : ℝ)) ^ 2 =
      ((i : ℝ) ^ 2 * v i) * v j + v i * ((j : ℝ) ^ 2 * v j) -
        2 * (((i : ℝ) * v i) * ((j : ℝ) * v j)) := by intro i j; ring
  simp_rw [h, Finset.sum_sub_distrib, Finset.sum_add_distrib, ← Finset.mul_sum,
    ← Finset.sum_mul]
  ring

/-- Splitting a symmetric sum over `[-K, K]`. -/
lemma sum_Icc_symm (f : ℤ → ℝ) (K : ℕ) :
    ∑ i ∈ Finset.Icc (-(K : ℤ)) (K : ℤ), f i =
      f 0 + ∑ r ∈ Finset.range K, (f ((r : ℤ) + 1) + f (-((r : ℤ) + 1))) := by
  induction K with
  | zero => simp
  | succ K ih =>
    have hset : Finset.Icc (-((K + 1 : ℕ) : ℤ)) ((K + 1 : ℕ) : ℤ) =
        insert ((K : ℤ) + 1) (insert (-((K : ℤ) + 1)) (Finset.Icc (-(K : ℤ)) (K : ℤ))) := by
      ext x; simp only [Finset.mem_Icc, Finset.mem_insert]; push_cast; omega
    rw [hset, Finset.sum_insert (by simp only [Finset.mem_insert, Finset.mem_Icc]; omega),
      Finset.sum_insert (by simp only [Finset.mem_Icc]; omega), ih,
      Finset.sum_range_succ]
    ring

lemma w_zero (δ : ℝ) : w δ 0 = 1 := by simp [w]

lemma w_pos (δ : ℝ) (r : ℕ) : w δ ((r : ℤ) + 1) = R δ (r + 1) := by
  have h : ((r : ℤ) + 1).toNat = r + 1 := by omega
  simp only [w, h]
  rw [if_neg (by omega), if_pos (by omega)]

lemma w_neg (δ : ℝ) (r : ℕ) : w δ (-((r : ℤ) + 1)) = L δ (r + 1) := by
  have h : (-(-((r : ℤ) + 1))).toNat = r + 1 := by omega
  simp only [w, h]
  rw [if_neg (by omega), if_neg (by omega)]

/-- Right-hand sums used in the explicit formula for `A`. -/
noncomputable def s0 (δ : ℝ) (K : ℕ) : ℝ :=
  1 + ∑ r ∈ Finset.range K, (R δ (r + 1) + L δ (r + 1))

noncomputable def s1 (δ : ℝ) (K : ℕ) : ℝ :=
  ∑ r ∈ Finset.range K, ((r : ℝ) + 1) * (R δ (r + 1) - L δ (r + 1))

noncomputable def s2 (δ : ℝ) (K : ℕ) : ℝ :=
  ∑ r ∈ Finset.range K, ((r : ℝ) + 1) ^ 2 * (R δ (r + 1) + L δ (r + 1))

/-- Explicit formula: `A δ K = s0 · s2 - s1²`. -/
theorem A_formula (δ : ℝ) (K : ℕ) : A δ K = s0 δ K * s2 δ K - s1 δ K ^ 2 := by
  rw [A_eq_Aw, Aw_pairwise, sum_Icc_symm, sum_Icc_symm, sum_Icc_symm]
  have e1 : ∀ r : ℕ, (((r : ℤ) + 1 : ℤ) : ℝ) = (r : ℝ) + 1 := by intro r; push_cast; ring
  have e2 : ∀ r : ℕ, ((-((r : ℤ) + 1) : ℤ) : ℝ) = -((r : ℝ) + 1) := by intro r; push_cast; ring
  simp only [w_zero, w_pos, w_neg, e1, e2, s0, s1, s2, Int.cast_zero, neg_sq, neg_mul,
    ← sub_eq_add_neg, ← mul_add, ← mul_sub, zero_pow two_ne_zero, zero_mul, zero_add]

/-- Symmetric weights `w_{±r} = L r` (and `w_0 = L 0 = 1`). -/
noncomputable def wsym (δ : ℝ) (i : ℤ) : ℝ := L δ i.natAbs

/-- `S = 1 + 2 ∑_{r=1}^K L r` (paper, eq. `S-T`). -/
noncomputable def Ssym (δ : ℝ) (K : ℕ) : ℝ := 1 + 2 * ∑ r ∈ Finset.range K, L δ (r + 1)

/-- `T = 2 ∑_{r=1}^K r² L r` (paper, eq. `S-T`). -/
noncomputable def Tsym (δ : ℝ) (K : ℕ) : ℝ :=
  2 * ∑ r ∈ Finset.range K, ((r : ℝ) + 1) ^ 2 * L δ (r + 1)

/-- **Shared lemma**: `A_sym = S · T`. -/
theorem Aw_wsym (δ : ℝ) (K : ℕ) : Aw (wsym δ) K = Ssym δ K * Tsym δ K := by
  rw [Aw_pairwise, sum_Icc_symm, sum_Icc_symm, sum_Icc_symm]
  have h1 : ∀ r : ℕ, ((r : ℤ) + 1).natAbs = r + 1 := by intro r; omega
  have h2 : ∀ r : ℕ, (-((r : ℤ) + 1)).natAbs = r + 1 := by intro r; omega
  simp only [wsym, h1, h2, Int.natAbs_zero, L_zero, Ssym, Tsym]
  push_cast
  have h3 : ∑ r ∈ Finset.range K, (((r : ℝ) + 1) * L δ (r + 1) +
      (-((r : ℝ) + 1)) * L δ (r + 1)) = 0 := by
    apply Finset.sum_eq_zero; intro r _; ring
  rw [h3, Finset.mul_sum, Finset.mul_sum]
  simp only [zero_pow (two_ne_zero), zero_mul, zero_add, sub_zero]
  congr 1
  · congr 1; apply Finset.sum_congr rfl; intro r _; ring
  · apply Finset.sum_congr rfl; intro r _; ring

/-- Symmetrisation (paper, Section 4): `A δ K ≥ A_sym = S · T`. -/
theorem A_ge_ST (δ : ℝ) (K : ℕ) (hδ0 : 0 < δ) (hK : ((K : ℝ) + 1) * δ < 1) :
    Ssym δ K * Tsym δ K ≤ A δ K := by
  rw [← Aw_wsym, A_eq_Aw]
  have hle : ∀ i ∈ Finset.Icc (-(K : ℤ)) (K : ℤ), ((i.natAbs : ℕ) : ℝ) ≤ K := by
    intro i hi
    have := Finset.mem_Icc.mp hi
    exact_mod_cast (by omega : i.natAbs ≤ K)
  apply Aw_mono
  · intro i hi
    exact (L_nonneg_le_R δ hδ0 _ (by nlinarith [hle i hi])).1
  · intro i hi
    unfold wsym w
    split_ifs with h0 hpos
    · subst h0; simp [L_zero]
    · have : i.natAbs = i.toNat := by omega
      rw [this]
      exact R_ge_L δ K hδ0 hK _ (by have := Finset.mem_Icc.mp hi; omega)
    · have : i.natAbs = (-i).toNat := by omega
      rw [this]

/-- For fixed `K`, `A δ K` is antitone in `δ` as long as `(K+1) δ₁ ≤ 1`. -/
theorem A_antitone (δ δ₁ : ℝ) (K : ℕ) (hδ0 : 0 < δ) (hδδ₁ : δ ≤ δ₁)
    (hK : ((K : ℝ) + 1) * δ₁ ≤ 1) : A δ₁ K ≤ A δ K := by
  rw [A_eq_Aw, A_eq_Aw]
  have key : ∀ i ∈ Finset.Icc (-(K : ℤ)) (K : ℤ), 0 ≤ w δ₁ i ∧ w δ₁ i ≤ w δ i := by
    intro i hi
    have hi' := Finset.mem_Icc.mp hi
    unfold w
    split_ifs with h0 hpos
    · simp
    · have hr : ((i.toNat : ℕ) : ℝ) ≤ K := by exact_mod_cast (by omega : i.toNat ≤ K)
      exact (LR_antitone δ δ₁ hδ0 hδδ₁ _ (by nlinarith)).2
    · have hr : (((-i).toNat : ℕ) : ℝ) ≤ K := by exact_mod_cast (by omega : (-i).toNat ≤ K)
      exact (LR_antitone δ δ₁ hδ0 hδδ₁ _ (by nlinarith)).1
  exact Aw_mono _ _ K (fun i hi => (key i hi).1) (fun i hi => (key i hi).2)

end PBScalar
