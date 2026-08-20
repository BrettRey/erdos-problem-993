import Mathlib

/-!
# Laurent coefficient sequences and central unimodality

This file develops the Laurent-polynomial theory that the clan audit needs:
coefficient extraction, the notion of a *centrally unimodal* (symmetric,
nonnegative, step-2 nonincreasing on `m ≥ 0`) Laurent polynomial, and closure
of that notion under products.

All Laurent polynomials here have rational coefficients; `LaurentPolynomial ℚ`
is definitionally the finitely supported maps `ℤ → ℚ`, which is the
representation the request allows ("finite maps from `Int` to `Rat`").
-/

namespace ClanAudit

open LaurentPolynomial Finset

/-- Coefficient of a Laurent polynomial. -/
def coeffL (f : LaurentPolynomial ℚ) (m : ℤ) : ℚ := f m

@[simp] theorem coeffL_add (f g : LaurentPolynomial ℚ) (m : ℤ) :
    coeffL (f + g) m = coeffL f m + coeffL g m := rfl

@[simp] theorem coeffL_sub (f g : LaurentPolynomial ℚ) (m : ℤ) :
    coeffL (f - g) m = coeffL f m - coeffL g m := rfl

@[simp] theorem coeffL_smul (c : ℚ) (f : LaurentPolynomial ℚ) (m : ℤ) :
    coeffL (c • f) m = c * coeffL f m := rfl

@[simp] theorem coeffL_zero (m : ℤ) : coeffL (0 : LaurentPolynomial ℚ) m = 0 := rfl

@[simp] theorem coeffL_T (t m : ℤ) :
    coeffL (T t : LaurentPolynomial ℚ) m = if m = t then 1 else 0 := by
  unfold coeffL; rw [LaurentPolynomial.T]
  simp [AddMonoidAlgebra.single, eq_comm]

@[simp] theorem coeffL_one (m : ℤ) :
    coeffL (1 : LaurentPolynomial ℚ) m = if m = 0 then 1 else 0 := by
  have h : (1 : LaurentPolynomial ℚ) = T 0 := by simp
  rw [h, coeffL_T]

theorem coeffL_mul_T (f : LaurentPolynomial ℚ) (t m : ℤ) :
    coeffL (f * T t) m = coeffL f (m - t) := by
  unfold coeffL
  rw [LaurentPolynomial.T, AddMonoidAlgebra.mul_single_apply, mul_one, sub_eq_add_neg]

theorem coeffL_sum {ι : Type*} (s : Finset ι) (f : ι → LaurentPolynomial ℚ) (m : ℤ) :
    coeffL (∑ i ∈ s, f i) m = ∑ i ∈ s, coeffL (f i) m := by
  classical
  induction s using Finset.induction with
  | empty => simp
  | insert a s ha ih => rw [Finset.sum_insert ha, coeffL_add, ih, Finset.sum_insert ha]

theorem coeffL_ext {f g : LaurentPolynomial ℚ} (h : ∀ m, coeffL f m = coeffL g m) : f = g :=
  Finsupp.ext h

/-- Every Laurent polynomial has a bound `N` outside of which all coefficients vanish. -/
theorem exists_bound (f : LaurentPolynomial ℚ) :
    ∃ N : ℕ, ∀ m : ℤ, N < m.natAbs → coeffL f m = 0 := by
  classical
  refine ⟨(f.support.image Int.natAbs).sup id, fun m hm => ?_⟩
  by_contra h
  have hmem : m ∈ f.support := Finsupp.mem_support_iff.mpr h
  have hmem2 : m.natAbs ∈ f.support.image Int.natAbs := Finset.mem_image_of_mem _ hmem
  have hle : m.natAbs ≤ (f.support.image Int.natAbs).sup id := Finset.le_sup (f := id) hmem2
  omega

/-! ### Centrally unimodal Laurent polynomials -/

/-- A Laurent polynomial is *centrally unimodal* when it is symmetric under `z ↦ z⁻¹`,
has nonnegative coefficients, and its coefficients are nonincreasing in steps of two
away from the centre. -/
structure IsCU (f : LaurentPolynomial ℚ) : Prop where
  symm : ∀ m : ℤ, coeffL f (-m) = coeffL f m
  nonneg : ∀ m : ℤ, 0 ≤ coeffL f m
  decr : ∀ m : ℤ, 0 ≤ m → coeffL f (m + 2) ≤ coeffL f m

theorem IsCU.mono {f : LaurentPolynomial ℚ} (hf : IsCU f) :
    ∀ (a b : ℤ), 0 ≤ a → a ≤ b → (b - a) % 2 = 0 → coeffL f b ≤ coeffL f a := by
  intro a b ha hab hpar
  obtain ⟨k, hk⟩ : ∃ k : ℕ, b = a + 2 * k := ⟨((b - a) / 2).toNat, by omega⟩
  subst hk
  clear hab hpar
  induction k with
  | zero => simp
  | succ n ih =>
      have h1 : coeffL f (a + 2 * n + 2) ≤ coeffL f (a + 2 * n) := hf.decr _ (by positivity)
      calc coeffL f (a + 2 * ((n : ℤ) + 1)) = coeffL f (a + 2 * n + 2) := by
            congr 1; ring
        _ ≤ coeffL f (a + 2 * n) := h1
        _ ≤ coeffL f a := ih

theorem IsCU.coeff_abs {f : LaurentPolynomial ℚ} (hf : IsCU f) (a : ℤ) :
    coeffL f a = coeffL f (a.natAbs : ℤ) := by
  rcases Int.natAbs_eq a with h | h
  · exact congrArg _ h
  · calc coeffL f a = coeffL f (-(a.natAbs : ℤ)) := by rw [← h]
      _ = coeffL f (a.natAbs : ℤ) := hf.symm _

theorem IsCU.mono_abs {f : LaurentPolynomial ℚ} (hf : IsCU f) {a b : ℤ}
    (hab : a.natAbs ≤ b.natAbs) (hpar : (b - a) % 2 = 0) : coeffL f b ≤ coeffL f a := by
  have h1 : coeffL f a = coeffL f (a.natAbs : ℤ) := hf.coeff_abs a
  have h2 : coeffL f b = coeffL f (b.natAbs : ℤ) := hf.coeff_abs b
  rw [h1, h2]
  exact hf.mono _ _ (by positivity) (by exact_mod_cast hab) (by omega)

theorem IsCU.add {f g : LaurentPolynomial ℚ} (hf : IsCU f) (hg : IsCU g) : IsCU (f + g) where
  symm m := by simp [hf.symm, hg.symm]
  nonneg m := by simpa using add_nonneg (hf.nonneg m) (hg.nonneg m)
  decr m hm := by simpa using add_le_add (hf.decr m hm) (hg.decr m hm)

theorem IsCU.smul {f : LaurentPolynomial ℚ} {c : ℚ} (hc : 0 ≤ c) (hf : IsCU f) : IsCU (c • f) where
  symm m := by simp [hf.symm]
  nonneg m := by simpa using mul_nonneg hc (hf.nonneg m)
  decr m hm := by simpa using mul_le_mul_of_nonneg_left (hf.decr m hm) hc

/-! ### Boxes -/

/-- The "box" `z^(-j) + z^(-j+2) + ⋯ + z^j`. -/
noncomputable def box (j : ℕ) : LaurentPolynomial ℚ :=
  ∑ k ∈ Finset.range (j + 1), T ((j : ℤ) - 2 * k)

theorem coeffL_box (j : ℕ) (m : ℤ) :
    coeffL (box j) m = if m.natAbs ≤ j ∧ ((j : ℤ) - m) % 2 = 0 then 1 else 0 := by
  classical
  rw [box, coeffL_sum]
  by_cases h : m.natAbs ≤ j ∧ ((j : ℤ) - m) % 2 = 0
  · obtain ⟨h1, h2⟩ := h
    have hk : (((j : ℤ) - m) / 2).toNat ∈ Finset.range (j + 1) := by
      simp only [Finset.mem_range]; omega
    rw [Finset.sum_eq_single (((j : ℤ) - m) / 2).toNat]
    · rw [coeffL_T, if_pos (by omega), if_pos ⟨h1, h2⟩]
    · intro b _ hb
      rw [coeffL_T, if_neg (by omega)]
    · intro hc; exact absurd hk hc
  · rw [if_neg h, Finset.sum_eq_zero]
    intro k hk
    simp only [Finset.mem_range] at hk
    rw [coeffL_T, if_neg]
    intro hcon
    exact h ⟨by omega, by omega⟩

theorem coeffL_box_self (j : ℕ) : coeffL (box j) (j : ℤ) = 1 := by
  rw [coeffL_box, if_pos ⟨by simp, by simp⟩]

theorem coeffL_box_step (j : ℕ) {m : ℤ} (hm : 0 ≤ m) :
    coeffL (box j) m - coeffL (box j) (m + 2) = if m = (j : ℤ) then 1 else 0 := by
  rw [coeffL_box, coeffL_box]
  by_cases hpar : ((j : ℤ) - m) % 2 = 0
  · have hpar2 : ((j : ℤ) - (m + 2)) % 2 = 0 := by omega
    by_cases h1 : m.natAbs ≤ j
    · by_cases h2 : (m + 2).natAbs ≤ j
      · rw [if_pos ⟨h1, hpar⟩, if_pos ⟨h2, hpar2⟩, if_neg (by omega)]; ring
      · rw [if_pos ⟨h1, hpar⟩, if_neg (fun hc => h2 hc.1), if_pos (by omega)]; ring
    · have h2 : ¬ (m + 2).natAbs ≤ j := by omega
      rw [if_neg (fun hc => h1 hc.1), if_neg (fun hc => h2 hc.1), if_neg (by omega)]; ring
  · have hpar2 : ¬ ((j : ℤ) - (m + 2)) % 2 = 0 := by omega
    rw [if_neg (fun hc => hpar hc.2), if_neg (fun hc => hpar2 hc.2), if_neg (by omega)]; ring

theorem isCU_box (j : ℕ) : IsCU (box j) where
  symm m := by rw [coeffL_box, coeffL_box]; congr 1; simp only [eq_iff_iff]; omega
  nonneg m := by rw [coeffL_box]; split <;> norm_num
  decr m hm := by
    have h1 := coeffL_box_step j hm
    have h2 : (0:ℚ) ≤ if m = (j:ℤ) then 1 else 0 := by split <;> norm_num
    linarith

theorem coeffL_mul_box (f : LaurentPolynomial ℚ) (j : ℕ) (m : ℤ) :
    coeffL (f * box j) m = ∑ k ∈ Finset.range (j + 1), coeffL f (m - (j : ℤ) + 2 * k) := by
  rw [box, Finset.mul_sum, coeffL_sum]
  refine Finset.sum_congr rfl fun k _ => ?_
  rw [coeffL_mul_T]
  congr 1
  ring

/-- Convolution with a box preserves central unimodality. -/
theorem IsCU.mul_box {f : LaurentPolynomial ℚ} (hf : IsCU f) (j : ℕ) : IsCU (f * box j) where
  symm m := by
    rw [coeffL_mul_box, coeffL_mul_box,
      ← Finset.sum_range_reflect (fun k => coeffL f (m - (j:ℤ) + 2 * k)) (j + 1)]
    refine Finset.sum_congr rfl fun k hk => ?_
    simp only [Finset.mem_range] at hk
    have hcast : ((j + 1 - 1 - k : ℕ) : ℤ) = (j : ℤ) - k := by omega
    have h : -m - (j : ℤ) + 2 * (k : ℤ) = -(m - (j:ℤ) + 2 * ((j + 1 - 1 - k : ℕ) : ℤ)) := by
      rw [hcast]; ring
    rw [h, hf.symm]
  nonneg m := by
    rw [coeffL_mul_box]
    exact Finset.sum_nonneg fun k _ => hf.nonneg _
  decr m hm := by
    have key : coeffL (f * box j) m - coeffL (f * box j) (m + 2)
        = coeffL f (m - (j:ℤ)) - coeffL f (m + (j:ℤ) + 2) := by
      rw [coeffL_mul_box, coeffL_mul_box,
        Finset.sum_range_succ' (fun k => coeffL f (m - (j:ℤ) + 2 * k)) j,
        Finset.sum_range_succ (fun k => coeffL f (m + 2 - (j:ℤ) + 2 * k)) j]
      have hterm : ∀ k ∈ Finset.range j,
          coeffL f (m - (j:ℤ) + 2 * (((k : ℕ) + 1 : ℕ) : ℤ)) = coeffL f (m + 2 - (j:ℤ) + 2 * k) := by
        intro k _; congr 1; push_cast; ring
      rw [Finset.sum_congr rfl hterm]
      have h2 : m + 2 - (j:ℤ) + 2 * j = m + (j:ℤ) + 2 := by ring
      rw [h2]
      simp only [Nat.cast_zero, mul_zero, add_zero]
      ring
    have hle : coeffL f (m + (j:ℤ) + 2) ≤ coeffL f (m - (j:ℤ)) :=
      hf.mono_abs (by omega) (by omega)
    linarith

/-! ### Products of centrally unimodal polynomials -/

private theorem IsCU.mul_aux : ∀ (N : ℕ) (g : LaurentPolynomial ℚ), IsCU g →
    (∀ m : ℤ, N < m.natAbs → coeffL g m = 0) →
    ∀ f : LaurentPolynomial ℚ, IsCU f → IsCU (f * g) := by
  intro N
  induction N using Nat.strong_induction_on with
  | _ N ih =>
    match N, ih with
    | 0, _ =>
      intro g hg hsupp f hf
      have hgeq : g = (coeffL g 0) • box 0 := by
        refine coeffL_ext fun m => ?_
        rw [coeffL_smul, coeffL_box]
        rcases eq_or_ne m 0 with rfl | hm
        · norm_num
        · rw [if_neg (by omega), hsupp m (by omega)]; ring
      rw [hgeq, mul_smul_comm]
      exact (hf.mul_box 0).smul (hg.nonneg 0)
    | 1, _ =>
      intro g hg hsupp f hf
      have hgeq : g = (coeffL g 0) • box 0 + (coeffL g 1) • box 1 := by
        refine coeffL_ext fun m => ?_
        rw [coeffL_add, coeffL_smul, coeffL_smul, coeffL_box, coeffL_box]
        rcases eq_or_ne m 0 with rfl | hm0
        · norm_num
        rcases eq_or_ne m 1 with rfl | hm1
        · norm_num
        rcases eq_or_ne m (-1) with rfl | hmm1
        · have hs : coeffL g (-1) = coeffL g 1 := hg.symm 1
          rw [hs]; norm_num
        · rw [if_neg (by omega), if_neg (by omega), hsupp m (by omega)]; ring
      rw [hgeq, mul_add, mul_smul_comm, mul_smul_comm]
      exact ((hf.mul_box 0).smul (hg.nonneg 0)).add ((hf.mul_box 1).smul (hg.nonneg 1))
    | (n + 2), ih =>
      intro g hg hsupp f hf
      set a : ℚ := coeffL g ((n : ℤ) + 2) with ha
      set b : ℚ := coeffL g ((n : ℤ) + 1) with hb
      set g' : LaurentPolynomial ℚ := g - a • box (n + 2) - b • box (n + 1) with hg'
      have hcoeff : ∀ m : ℤ, coeffL g' m
          = coeffL g m - a * coeffL (box (n+2)) m - b * coeffL (box (n+1)) m := by
        intro m; rw [hg']; simp
      have ha0 : 0 ≤ a := hg.nonneg _
      have hb0 : 0 ≤ b := hg.nonneg _
      have hcast2 : (((n + 2 : ℕ) : ℤ)) = (n : ℤ) + 2 := by push_cast; ring
      have hcast1 : (((n + 1 : ℕ) : ℤ)) = (n : ℤ) + 1 := by push_cast; ring
      have hg'CU : IsCU g' := by
        constructor
        · intro m
          rw [hcoeff, hcoeff, hg.symm, (isCU_box (n+2)).symm, (isCU_box (n+1)).symm]
        · intro m
          rw [hcoeff, coeffL_box, coeffL_box]
          by_cases h1 : m.natAbs ≤ n + 2 ∧ (((n + 2 : ℕ) : ℤ) - m) % 2 = 0
          · by_cases h2 : m.natAbs ≤ n + 1 ∧ (((n + 1 : ℕ) : ℤ) - m) % 2 = 0
            · exfalso; omega
            · rw [if_pos h1, if_neg h2]
              have hle : a ≤ coeffL g m := by
                rw [ha, ← hcast2]
                exact hg.mono_abs (by simpa using h1.1) (by omega)
              linarith
          · by_cases h2 : m.natAbs ≤ n + 1 ∧ (((n + 1 : ℕ) : ℤ) - m) % 2 = 0
            · rw [if_neg h1, if_pos h2]
              have hle : b ≤ coeffL g m := by
                rw [hb, ← hcast1]
                exact hg.mono_abs (by simpa using h2.1) (by omega)
              linarith
            · rw [if_neg h1, if_neg h2]
              have := hg.nonneg m; linarith
        · intro m hm
          rw [hcoeff, hcoeff]
          have e1 := coeffL_box_step (n+2) hm
          have e2 := coeffL_box_step (n+1) hm
          have e3 := hg.decr m hm
          rw [hcast2] at e1
          rw [hcast1] at e2
          have f1 : a * (coeffL (box (n+2)) m - coeffL (box (n+2)) (m + 2))
              = a * (if m = (n : ℤ) + 2 then 1 else 0) := by rw [e1]
          have f2 : b * (coeffL (box (n+1)) m - coeffL (box (n+1)) (m + 2))
              = b * (if m = (n : ℤ) + 1 then 1 else 0) := by rw [e2]
          rcases eq_or_ne m ((n : ℤ) + 2) with rfl | h1
          · have hz : coeffL g ((n : ℤ) + 2 + 2) = 0 := hsupp _ (by omega)
            rw [if_pos rfl] at f1
            rw [if_neg (by omega)] at f2
            linarith [ha]
          · rcases eq_or_ne m ((n : ℤ) + 1) with rfl | h2
            · have hz : coeffL g ((n : ℤ) + 1 + 2) = 0 := hsupp _ (by omega)
              rw [if_neg (by omega)] at f1
              rw [if_pos rfl] at f2
              linarith [hb]
            · rw [if_neg h1] at f1
              rw [if_neg h2] at f2
              linarith
      have hsupp' : ∀ m : ℤ, n < m.natAbs → coeffL g' m = 0 := by
        intro m hm
        rw [hcoeff, coeffL_box, coeffL_box]
        rcases lt_trichotomy m.natAbs (n + 2) with h | h | h
        · have hmn : m.natAbs = n + 1 := by omega
          have hval : coeffL g m = b := by
            rw [hg.coeff_abs m, hmn, hb]; norm_cast
          rw [hval, if_neg (by omega), if_pos ⟨by omega, by omega⟩]
          ring
        · have hval : coeffL g m = a := by
            rw [hg.coeff_abs m, h, ha]; norm_cast
          rw [hval, if_pos ⟨by omega, by omega⟩, if_neg (by omega)]
          ring
        · rw [hsupp m (by omega), if_neg (by omega), if_neg (by omega)]
          ring
      have hmain : IsCU (f * g') := ih n (by omega) g' hg'CU hsupp' f hf
      have hsplit : f * g = f * g' + a • (f * box (n+2)) + b • (f * box (n+1)) := by
        rw [hg', mul_sub, mul_sub, mul_smul_comm, mul_smul_comm]
        ring
      rw [hsplit]
      exact (hmain.add ((hf.mul_box (n+2)).smul ha0)).add ((hf.mul_box (n+1)).smul hb0)

theorem IsCU.mul {f g : LaurentPolynomial ℚ} (hf : IsCU f) (hg : IsCU g) : IsCU (f * g) := by
  obtain ⟨N, hN⟩ := exists_bound g
  exact IsCU.mul_aux N g hg hN f hf

end ClanAudit
