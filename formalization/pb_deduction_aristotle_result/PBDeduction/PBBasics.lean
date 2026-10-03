import PBDeduction.MaxMass

/-!
# Poisson–binomial basics (G3)

Everything is reduced to the one-step recursion
`f^{[n+1]}_k = (1 - q) f^{[n]}_k + q f^{[n]}_{k-1}` (adding one Bernoulli factor),
`pbPmf_succ`, and proved by induction on the number of factors.
-/

namespace PBDeduction

open Finset Polynomial

noncomputable section

/-- Adding one Bernoulli(`q`) summand: `(conv q g) k = (1-q) g k + q g (k-1)`. -/
def conv (q : ℝ) (g : ℤ → ℝ) (k : ℤ) : ℝ := (1 - q) * g k + q * g (k - 1)

lemma pgf_succ (n : ℕ) (p : Fin (n + 1) → ℝ) :
    pgf (n + 1) p = (C (1 - p 0) + C (p 0) * X) * pgf n (Fin.tail p) := by
  unfold pgf; rw [Fin.prod_univ_succ]; rfl

lemma pbPmf_zero (p : Fin 0 → ℝ) (k : ℤ) : pbPmf 0 p k = if k = 0 then 1 else 0 := by
  unfold pbPmf pgf
  simp only [Finset.univ_eq_empty, Finset.prod_empty, Polynomial.coeff_one]
  by_cases h : k = 0
  · subst h; simp
  · rw [if_neg h]
    split_ifs with h1 h2
    · omega
    · rfl
    · rfl

lemma pbPmf_succ (n : ℕ) (p : Fin (n + 1) → ℝ) :
    pbPmf (n + 1) p = conv (p 0) (pbPmf n (Fin.tail p)) := by
  funext k
  unfold conv pbPmf
  rw [pgf_succ]
  set Q := pgf n (Fin.tail p)
  by_cases hk : 0 ≤ k
  · rw [if_pos hk, if_pos hk]
    rcases eq_or_lt_of_le hk with h0 | h0
    · subst h0
      rw [if_neg (by omega)]
      simp only [Int.toNat_zero]
      rw [add_mul, coeff_add, coeff_C_mul, mul_assoc, coeff_C_mul, coeff_X_mul_zero, mul_zero,
        add_zero]
    · rw [if_pos (by omega)]
      obtain ⟨m, hm⟩ : ∃ m : ℕ, k = m + 1 := ⟨(k - 1).toNat, by omega⟩
      subst hm
      have e1 : ((m : ℤ) + 1).toNat = m + 1 := by omega
      have e2 : ((m : ℤ) + 1 - 1).toNat = m := by omega
      rw [e1, e2, add_mul, coeff_add, coeff_C_mul, mul_assoc, coeff_C_mul, coeff_X_mul]
  · rw [if_neg hk, if_neg hk, if_neg (by omega)]; ring

/-! ## Support and total mass -/

theorem pbPmf_zero_out (n : ℕ) (p : Fin n → ℝ) :
    ∀ k : ℤ, k < 0 ∨ (n : ℤ) < k → pbPmf n p k = 0 := by
  induction n with
  | zero =>
    intro k hk; rw [pbPmf_zero, if_neg (by omega)]
  | succ n ih =>
    intro k hk
    rw [pbPmf_succ]
    unfold conv
    rw [ih (Fin.tail p) k (by omega), ih (Fin.tail p) (k - 1) (by omega)]; ring

lemma Icc_succ_top (m : ℕ) :
    Icc (0 : ℤ) ((m + 1 : ℕ) : ℤ) = insert ((m : ℤ) + 1) (Icc (0 : ℤ) m) := by
  ext x; simp only [mem_Icc, mem_insert]; push_cast; omega

lemma sum_Icc_succ_top' (F : ℤ → ℝ) (m : ℕ) :
    ∑ k ∈ Icc (0 : ℤ) ((m + 1 : ℕ) : ℤ), F k = ∑ k ∈ Icc (0 : ℤ) m, F k + F ((m : ℤ) + 1) := by
  rw [Icc_succ_top, sum_insert (by simp), add_comm]

lemma sum_Icc_succ_shift (F : ℤ → ℝ) (m : ℕ) :
    ∑ k ∈ Icc (0 : ℤ) ((m + 1 : ℕ) : ℤ), F k = F 0 + ∑ k ∈ Icc (0 : ℤ) m, F (k + 1) := by
  induction m with
  | zero => rw [sum_Icc_succ_top' F 0]; simp
  | succ m ih =>
    rw [sum_Icc_succ_top', ih, sum_Icc_succ_top' (fun k => F (k + 1)) m]
    push_cast; ring_nf

/-- The moment recursion for `conv`. -/
lemma sum_conv (q : ℝ) (g h : ℤ → ℝ) (n : ℕ) (hg1 : g (-1) = 0) (hgn : g ((n : ℤ) + 1) = 0) :
    ∑ k ∈ Icc (0 : ℤ) ((n + 1 : ℕ) : ℤ), h k * conv q g k =
      (1 - q) * ∑ k ∈ Icc (0 : ℤ) n, h k * g k + q * ∑ k ∈ Icc (0 : ℤ) n, h (k + 1) * g k := by
  unfold conv
  have e : ∀ k, h k * ((1 - q) * g k + q * g (k - 1)) =
      (1 - q) * (h k * g k) + q * (h k * g (k - 1)) := fun k => by ring
  simp only [e, sum_add_distrib, ← mul_sum]
  rw [sum_Icc_succ_top' (fun k => h k * g k), sum_Icc_succ_shift (fun k => h k * g (k - 1))]
  simp only [add_sub_cancel_right, zero_sub, hg1, hgn, mul_zero, zero_add, add_zero]

lemma pbPmf_neg_one (n : ℕ) (p : Fin n → ℝ) : pbPmf n p (-1) = 0 :=
  pbPmf_zero_out n p (-1) (Or.inl (by omega))

lemma pbPmf_top (n : ℕ) (p : Fin n → ℝ) : pbPmf n p ((n : ℤ) + 1) = 0 :=
  pbPmf_zero_out n p _ (Or.inr (by omega))

theorem pbPmf_sum_one (n : ℕ) (p : Fin n → ℝ) :
    ∑ k ∈ Finset.Icc (0 : ℤ) n, pbPmf n p k = 1 := by
  induction n with
  | zero => simp [pbPmf_zero]
  | succ n ih =>
    rw [pbPmf_succ]
    have := sum_conv (p 0) (pbPmf n (Fin.tail p)) (fun _ => 1) n (pbPmf_neg_one _ _)
      (pbPmf_top _ _)
    simp only [one_mul] at this
    rw [this, ih]; ring

/-! ## Mean and variance -/

lemma pbPmf_mean (n : ℕ) (p : Fin n → ℝ) :
    ∑ k ∈ Finset.Icc (0 : ℤ) n, (k : ℝ) * pbPmf n p k = ∑ i, p i := by
  induction n with
  | zero => simp [pbPmf_zero]
  | succ n ih =>
    rw [pbPmf_succ]
    have := sum_conv (p 0) (pbPmf n (Fin.tail p)) (fun k => (k : ℝ)) n (pbPmf_neg_one _ _)
      (pbPmf_top _ _)
    simp only at this
    rw [this, Fin.sum_univ_succ]
    have e : ∀ k : ℤ, ((k + 1 : ℤ) : ℝ) * pbPmf n (Fin.tail p) k =
        (k : ℝ) * pbPmf n (Fin.tail p) k + pbPmf n (Fin.tail p) k := fun k => by push_cast; ring
    simp only [e, sum_add_distrib, ih, pbPmf_sum_one]
    simp only [Fin.tail]; ring

lemma pbPmf_second (n : ℕ) (p : Fin n → ℝ) :
    ∑ k ∈ Finset.Icc (0 : ℤ) n, (k : ℝ) ^ 2 * pbPmf n p k =
      (∑ i, p i) ^ 2 + ∑ i, p i * (1 - p i) := by
  induction n with
  | zero => simp [pbPmf_zero]
  | succ n ih =>
    rw [pbPmf_succ]
    have := sum_conv (p 0) (pbPmf n (Fin.tail p)) (fun k => (k : ℝ) ^ 2) n (pbPmf_neg_one _ _)
      (pbPmf_top _ _)
    simp only at this
    rw [this, Fin.sum_univ_succ, Fin.sum_univ_succ]
    have e : ∀ k : ℤ, ((k + 1 : ℤ) : ℝ) ^ 2 * pbPmf n (Fin.tail p) k =
        (k : ℝ) ^ 2 * pbPmf n (Fin.tail p) k + 2 * ((k : ℝ) * pbPmf n (Fin.tail p) k)
          + pbPmf n (Fin.tail p) k := fun k => by push_cast; ring
    simp only [e, sum_add_distrib, ← mul_sum, ih, pbPmf_sum_one, pbPmf_mean]
    simp only [Fin.tail]; ring

theorem pbPmf_pairVar (n : ℕ) (p : Fin n → ℝ) :
    pairVar (pbPmf n p) n = pbVar n p := by
  rw [pairVar_eq_central _ _ (pbPmf_sum_one n p)]
  have hm : ∑ j ∈ Icc (0 : ℤ) n, pbPmf n p j * (j : ℝ) = ∑ i, p i := by
    rw [← pbPmf_mean]; apply sum_congr rfl; intro j _; ring
  rw [hm]
  have e : ∀ k : ℤ, pbPmf n p k * ((k : ℝ) - ∑ i, p i) ^ 2 =
      (k : ℝ) ^ 2 * pbPmf n p k - 2 * (∑ i, p i) * ((k : ℝ) * pbPmf n p k)
        + (∑ i, p i) ^ 2 * pbPmf n p k := fun k => by ring
  simp only [e, sum_add_distrib, sum_sub_distrib, ← mul_sum, pbPmf_second, pbPmf_mean,
    pbPmf_sum_one]
  unfold pbVar; ring

end

end PBDeduction
