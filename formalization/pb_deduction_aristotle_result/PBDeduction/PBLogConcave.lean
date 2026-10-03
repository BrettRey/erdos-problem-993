import PBDeduction.PBBasics

/-!
# Positivity, strict log-concavity (Newton) and the first descent (G3)
-/

namespace PBDeduction

open Finset

noncomputable section

/-! ## Positivity -/

theorem pbPmf_nonneg_pos (n : ℕ) (p : Fin n → ℝ) (hp : ∀ i, 0 < p i ∧ p i < 1) :
    (∀ k, 0 ≤ pbPmf n p k) ∧ (∀ k : ℤ, 0 ≤ k → k ≤ n → 0 < pbPmf n p k) := by
  induction n with
  | zero =>
    refine ⟨fun k => ?_, fun k h0 h1 => ?_⟩
    · rw [pbPmf_zero]; split_ifs <;> norm_num
    · rw [pbPmf_zero, if_pos (by omega)]; norm_num
  | succ n ih =>
    obtain ⟨ih1, ih2⟩ := ih (Fin.tail p) (fun i => hp i.succ)
    obtain ⟨q0, q1⟩ := hp 0
    rw [pbPmf_succ]
    unfold conv
    refine ⟨fun k => ?_, fun k h0 h1 => ?_⟩
    · have := ih1 k; have := ih1 (k - 1)
      have : 0 ≤ 1 - p 0 := by linarith
      positivity
    · have a1 := ih1 k; have a2 := ih1 (k - 1)
      have : 0 < 1 - p 0 := by linarith
      by_cases hk : k ≤ n
      · have := ih2 k h0 hk
        have : 0 < (1 - p 0) * pbPmf n (Fin.tail p) k := by positivity
        have : 0 ≤ p 0 * pbPmf n (Fin.tail p) (k - 1) := by positivity
        linarith
      · have := ih2 (k - 1) (by omega) (by push_cast at h1; omega)
        have : 0 ≤ (1 - p 0) * pbPmf n (Fin.tail p) k := by positivity
        have : 0 < p 0 * pbPmf n (Fin.tail p) (k - 1) := by positivity
        linarith

/-! ## Log-concavity -/

/-- The "cross" inequality `g_{k-2} g_{k+1} ≤ g_{k-1} g_k` for a nonnegative
log-concave sequence that is positive exactly on `{0, …, n}`. -/
lemma cross_ineq (g : ℤ → ℝ) (n : ℕ) (hnn : ∀ k, 0 ≤ g k)
    (hpos : ∀ k : ℤ, 0 ≤ k → k ≤ n → 0 < g k)
    (hzero : ∀ k : ℤ, k < 0 ∨ (n : ℤ) < k → g k = 0)
    (hlc : ∀ k, g (k - 1) * g (k + 1) ≤ g k ^ 2) (k : ℤ) :
    g (k - 2) * g (k + 1) ≤ g (k - 1) * g k := by
  by_cases hk : 1 ≤ k ∧ k ≤ n
  · have p1 := hpos (k - 1) (by omega) (by omega)
    have p2 := hpos k (by omega) hk.2
    have l1 := hlc (k - 1)
    have l2 := hlc k
    rw [show k - 1 - 1 = k - 2 by ring, show k - 1 + 1 = k by ring] at l1
    have hprod : g (k - 2) * g (k + 1) * (g (k - 1) * g k) ≤
        (g (k - 1) * g k) * (g (k - 1) * g k) := by
      have := mul_le_mul l1 l2 (mul_nonneg (hnn _) (hnn _)) (sq_nonneg _)
      nlinarith
    exact le_of_mul_le_mul_right (by linarith) (mul_pos p1 p2)
  · have : g (k - 2) * g (k + 1) = 0 := by
      rcases not_and_or.mp hk with h | h
      · rw [hzero (k - 2) (Or.inl (by omega))]; ring
      · rw [hzero (k + 1) (Or.inr (by omega))]; ring
    rw [this]; exact mul_nonneg (hnn _) (hnn _)

/-- The Turán expansion for `conv`. -/
lemma conv_turan (q : ℝ) (g : ℤ → ℝ) (k : ℤ) :
    conv q g k ^ 2 - conv q g (k - 1) * conv q g (k + 1) =
      (1 - q) ^ 2 * (g k ^ 2 - g (k - 1) * g (k + 1)) +
      q ^ 2 * (g (k - 1) ^ 2 - g (k - 2) * g k) +
      q * (1 - q) * (g (k - 1) * g k - g (k - 2) * g (k + 1)) := by
  unfold conv
  rw [show k - 1 - 1 = k - 2 by ring, show k + 1 - 1 = k by ring]
  ring

theorem pbPmf_lc (n : ℕ) (p : Fin n → ℝ) (hp : ∀ i, 0 < p i ∧ p i < 1) :
    (∀ k, pbPmf n p (k - 1) * pbPmf n p (k + 1) ≤ pbPmf n p k ^ 2) ∧
    (∀ k : ℤ, 0 ≤ k → k ≤ n → pbPmf n p (k - 1) * pbPmf n p (k + 1) < pbPmf n p k ^ 2) := by
  induction n with
  | zero =>
    have h : ∀ k, pbPmf 0 p (k - 1) * pbPmf 0 p (k + 1) = 0 := by
      intro k
      rw [pbPmf_zero, pbPmf_zero]
      split_ifs <;> first | omega | ring
    refine ⟨fun k => by rw [h]; positivity, fun k h0 h1 => ?_⟩
    rw [h, show k = 0 by omega, pbPmf_zero, if_pos rfl]; norm_num
  | succ n ih =>
    set g := pbPmf n (Fin.tail p)
    have hp' : ∀ i, 0 < Fin.tail p i ∧ Fin.tail p i < 1 := fun i => hp i.succ
    obtain ⟨ih1, ih2⟩ := ih (Fin.tail p) hp'
    obtain ⟨hnn, hpos⟩ := pbPmf_nonneg_pos n (Fin.tail p) hp'
    have hcross := cross_ineq g n hnn hpos (pbPmf_zero_out n _) ih1
    obtain ⟨q0, q1⟩ := hp 0
    rw [pbPmf_succ]
    have hq : 0 < (1 - p 0) := by linarith
    have key : ∀ k, conv (p 0) g k ^ 2 - conv (p 0) g (k - 1) * conv (p 0) g (k + 1) =
        (1 - p 0) ^ 2 * (g k ^ 2 - g (k - 1) * g (k + 1)) +
        p 0 ^ 2 * (g (k - 1) ^ 2 - g (k - 2) * g k) +
        p 0 * (1 - p 0) * (g (k - 1) * g k - g (k - 2) * g (k + 1)) := conv_turan (p 0) g
    have t1 : ∀ k, 0 ≤ (1 - p 0) ^ 2 * (g k ^ 2 - g (k - 1) * g (k + 1)) := fun k =>
      mul_nonneg (by positivity) (by linarith [ih1 k])
    have t2 : ∀ k, 0 ≤ p 0 ^ 2 * (g (k - 1) ^ 2 - g (k - 2) * g k) := by
      intro k
      have := ih1 (k - 1)
      rw [show k - 1 - 1 = k - 2 by ring, show k - 1 + 1 = k by ring] at this
      exact mul_nonneg (by positivity) (by linarith)
    have t3 : ∀ k, 0 ≤ p 0 * (1 - p 0) * (g (k - 1) * g k - g (k - 2) * g (k + 1)) := fun k =>
      mul_nonneg (by positivity) (by linarith [hcross k])
    refine ⟨fun k => ?_, fun k h0 h1 => ?_⟩
    · have := key k; linarith [t1 k, t2 k, t3 k]
    · have := key k
      by_cases hk : k ≤ n
      · have s := ih2 k h0 hk
        have : 0 < (1 - p 0) ^ 2 * (g k ^ 2 - g (k - 1) * g (k + 1)) :=
          mul_pos (by positivity) (by linarith)
        linarith [t2 k, t3 k]
      · have s := ih2 (k - 1) (by omega) (by push_cast at h1; omega)
        rw [show k - 1 - 1 = k - 2 by ring, show k - 1 + 1 = k by ring] at s
        have : 0 < p 0 ^ 2 * (g (k - 1) ^ 2 - g (k - 2) * g k) :=
          mul_pos (by positivity) (by linarith)
        linarith [t1 k, t3 k]

/-! ## First descent -/

/-- `f_{n-1} = f_n ∑ (1 - p i)/p i`. -/
lemma pbPmf_penultimate (n : ℕ) (p : Fin n → ℝ) (hp : ∀ i, 0 < p i) :
    pbPmf n p ((n : ℤ) - 1) = pbPmf n p n * ∑ i, (1 - p i) / p i := by
  induction n with
  | zero => simp [pbPmf_zero]
  | succ n ih =>
    have ih := ih (Fin.tail p) (fun i => hp i.succ)
    rw [pbPmf_succ, Fin.sum_univ_succ]
    unfold conv
    push_cast
    rw [show (n : ℤ) + 1 - 1 = n by ring, ih, pbPmf_top]
    have := hp 0
    simp only [Fin.tail]
    field_simp
    ring

theorem pbPmf_first_descent_exists (n : ℕ) (p : Fin n → ℝ)
    (hp : ∀ i, 0 < p i ∧ p i < 1) (hV : 1 ≤ pbVar n p) :
    ∃ D : ℤ, IsFirstDescent (pbPmf n p) n D := by
  have hn : 1 ≤ n := by
    rcases Nat.eq_zero_or_pos n with h | h
    · subst h; simp [pbVar] at hV; linarith
    · exact h
  have hpos := (pbPmf_nonneg_pos n p hp).2
  have hT : pbVar n p < ∑ i, (1 - p i) / p i := by
    unfold pbVar
    have : (Finset.univ : Finset (Fin n)).Nonempty := ⟨⟨0, hn⟩, by simp⟩
    apply Finset.sum_lt_sum_of_nonempty this
    intro i _
    obtain ⟨a, b⟩ := hp i
    rw [lt_div_iff₀ a]
    have : 0 < (1 - p i) * (1 - p i ^ 2) := mul_pos (by linarith) (by nlinarith)
    nlinarith
  have hdesc : pbPmf n p n < pbPmf n p ((n : ℤ) - 1) := by
    rw [pbPmf_penultimate n p (fun i => (hp i).1)]
    have := hpos n (by omega) le_rfl
    have : 1 < ∑ i, (1 - p i) / p i := by linarith
    nlinarith
  have hex : ∃ m : ℕ, 1 ≤ m ∧ m ≤ n ∧ pbPmf n p m < pbPmf n p ((m : ℤ) - 1) :=
    ⟨n, hn, le_rfl, hdesc⟩
  classical
  let m := Nat.find hex
  have hm := Nat.find_spec hex
  refine ⟨m, by exact_mod_cast hm.1, by exact_mod_cast hm.2.1, hm.2.2, ?_⟩
  intro k hk1 hkm
  obtain ⟨j, rfl⟩ : ∃ j : ℕ, k = j := ⟨k.toNat, by omega⟩
  have hjm : j < m := by exact_mod_cast hkm
  have := Nat.find_min hex hjm
  push_neg at this
  exact this (by exact_mod_cast hk1) (by have := hm.2.1; omega)

end

end PBDeduction
