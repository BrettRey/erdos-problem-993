import RequestProject.Defs

/-!
# Log-concavity for `ℤ`-indexed sequences

It is convenient to run the convolution argument on `ℤ`-indexed sequences,
where the index shifts appearing in the Cauchy–Binet identity cause no
truncated-subtraction trouble.  This file sets up the zero-extension
`LC.ZExt` of a sequence `ℕ → ℝ` to `ℤ → ℝ`, and proves the key
"inner product dominates outer product" lemma
`LC.ZLogConcave.inner_mul_ge`:

  for `a ≤ b ≤ c ≤ d` with `a + d = b + c`, one has `F b * F c ≥ F a * F d`.
-/

open scoped BigOperators

namespace LC

/-- Zero-extension of a sequence `ℕ → ℝ` to `ℤ → ℝ`. -/
noncomputable def ZExt (f : ℕ → ℝ) : ℤ → ℝ := fun k => if 0 ≤ k then f k.toNat else 0

@[simp] lemma ZExt_natCast (f : ℕ → ℝ) (n : ℕ) : ZExt f (n : ℤ) = f n := by
  simp [ZExt]

lemma ZExt_of_neg (f : ℕ → ℝ) {k : ℤ} (hk : k < 0) : ZExt f k = 0 := by
  simp [ZExt, not_le.mpr hk]

lemma ZExt_nonneg {f : ℕ → ℝ} (hf : Nonneg f) : ∀ k, 0 ≤ ZExt f k := by
  intro k
  by_cases h : 0 ≤ k
  · simpa [ZExt, h] using hf k.toNat
  · simp [ZExt, h]

lemma ZExt_logConcave {f : ℕ → ℝ} (hf0 : Nonneg f) (hf : LogConcave f) :
    ∀ k : ℤ, ZExt f (k + 1) * ZExt f (k + 1) ≥ ZExt f k * ZExt f (k + 2) := by
  intro k
  rcases lt_trichotomy k 0 with hk | hk | hk
  · -- `k < 0`
    rw [ZExt_of_neg f hk, zero_mul]
    exact mul_nonneg (ZExt_nonneg hf0 _) (ZExt_nonneg hf0 _)
  · subst hk
    have h := hf 0
    simpa using (by simpa using h : ZExt f 1 * ZExt f 1 ≥ ZExt f 0 * ZExt f 2)
  · -- `0 < k`, so `k = m` for a natural number `m`
    obtain ⟨m, rfl⟩ : ∃ m : ℕ, k = (m : ℤ) := ⟨k.toNat, by omega⟩
    have h := hf m
    have e1 : ((m : ℤ) + 1) = ((m + 1 : ℕ) : ℤ) := by push_cast; ring
    have e2 : ((m : ℤ) + 2) = ((m + 2 : ℕ) : ℤ) := by push_cast; ring
    rw [e1, e2, ZExt_natCast, ZExt_natCast, ZExt_natCast]
    exact h

lemma ZExt_noInternalZeros {f : ℕ → ℝ} (hf : NoInternalZeros f) :
    ∀ i j k : ℤ, i ≤ j → j ≤ k → ZExt f i ≠ 0 → ZExt f k ≠ 0 → ZExt f j ≠ 0 := by
  intro i j k hij hjk hi hk
  have hi0 : 0 ≤ i := by
    by_contra h
    exact hi (ZExt_of_neg f (by omega))
  have hj0 : 0 ≤ j := le_trans hi0 hij
  have hk0 : 0 ≤ k := le_trans hj0 hjk
  simp only [ZExt, if_pos hi0, if_pos hj0, if_pos hk0] at hi hk ⊢
  exact hf i.toNat j.toNat k.toNat (by omega) (by omega) hi hk

section

variable {F : ℤ → ℝ}
variable (hpos : ∀ k, 0 ≤ F k)
variable (hlc : ∀ k : ℤ, F (k + 1) * F (k + 1) ≥ F k * F (k + 2))
variable (hnz : ∀ i j k : ℤ, i ≤ j → j ≤ k → F i ≠ 0 → F k ≠ 0 → F j ≠ 0)

include hpos hlc hnz

/-- The one-step "shift" inequality: for `a ≤ b`, `F (a+1) * F b ≥ F a * F (b+1)`. -/
theorem step_mul_ge : ∀ {a b : ℤ}, a ≤ b → F (a + 1) * F b ≥ F a * F (b + 1) := by
  intro a b hab
  induction b, hab using Int.le_induction with
  | base => exact le_of_eq (mul_comm _ _)
  | succ b hab ih =>
    by_cases hb : F (b + 1) = 0
    · have : F a * F (b + 1 + 1) = 0 := by
        by_cases ha : F a = 0
        · simp [ha]
        by_cases hd : F (b + 1 + 1) = 0
        · simp [hd]
        exact absurd (hnz a (b + 1) (b + 1 + 1) (by omega) (by omega) ha hd) (by simpa using hb)
      rw [this]
      exact mul_nonneg (hpos _) (hpos _)
    · have hbpos : 0 < F (b + 1) := lt_of_le_of_ne (hpos _) (Ne.symm hb)
      have h1 : F (a + 1) * (F (b + 1) * F (b + 1)) ≥ F (a + 1) * (F b * F (b + 2)) :=
        mul_le_mul_of_nonneg_left (hlc b) (hpos _)
      have h2 : (F (a + 1) * F b) * F (b + 2) ≥ (F a * F (b + 1)) * F (b + 2) :=
        mul_le_mul_of_nonneg_right ih (hpos _)
      have h3 : F (b + 1) * (F (a + 1) * F (b + 1)) ≥ F (b + 1) * (F a * F (b + 1 + 1)) := by
        have e : b + 1 + 1 = b + 2 := by ring
        rw [e]; nlinarith [h1, h2]
      exact le_of_mul_le_mul_left h3 hbpos

/-- For `a ≤ b ≤ c ≤ d` with `a + d = b + c`, the "inner" product dominates the
"outer" product: `F b * F c ≥ F a * F d`. -/
theorem inner_mul_ge {a b c d : ℤ} (hab : a ≤ b) (hbc : b ≤ c) (hcd : c ≤ d)
    (hsum : a + d = b + c) : F b * F c ≥ F a * F d := by
  obtain ⟨s, hs⟩ : ∃ s : ℕ, b = a + (s : ℤ) := ⟨(b - a).toNat, by omega⟩
  have hd : d = c + (s : ℤ) := by omega
  subst hs; subst hd
  clear hab hcd hsum
  induction s generalizing a with
  | zero => simp
  | succ s ih =>
    have hac : a + 1 + (s : ℤ) ≤ c := by push_cast at hbc ⊢; omega
    have h1 : F (a + 1 + (s : ℤ)) * F c ≥ F (a + 1) * F (c + (s : ℤ)) := ih hac
    have h2 : F (a + 1) * F (c + (s : ℤ)) ≥ F a * F (c + (s : ℤ) + 1) :=
      step_mul_ge hpos hlc hnz (by omega)
    have e1 : a + ((s : ℤ) + 1) = a + 1 + (s : ℤ) := by ring
    have e2 : c + ((s : ℤ) + 1) = c + (s : ℤ) + 1 := by ring
    push_cast
    rw [e1, e2]
    linarith

end

end LC
