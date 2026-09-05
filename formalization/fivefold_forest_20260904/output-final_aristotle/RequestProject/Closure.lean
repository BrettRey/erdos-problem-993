import RequestProject.Union

/-!
# Product closure

The fivefold margins are preserved under disjoint unions, once the exceptional
`K₁,₄` profile is allowed as an alternative for the factors.
-/

open Finset

namespace FivefoldForest

open scoped Classical

variable {V : Type*} {G : SimpleGraph V} {S S₁ S₂ : Finset V}

lemma cast_iD (G : SimpleGraph V) (T : Finset V) (d : ℕ) :
    (iD G T d : ℤ) = (eD G T d : ℤ) + (bD G T d : ℤ) := by
  exact_mod_cast congrArg (Nat.cast : ℕ → ℤ) (eD_add_bD G T d).symm

lemma eD_split_three (h : Split G S S₁ S₂) :
    (eD G S 3 : ℤ) = (eD G S₁ 0 : ℤ) * (eD G S₂ 3 : ℤ) + (eD G S₁ 1 : ℤ) * (eD G S₂ 2 : ℤ)
      + (eD G S₁ 2 : ℤ) * (eD G S₂ 1 : ℤ) + (eD G S₁ 3 : ℤ) * (eD G S₂ 0 : ℤ) := by
  have := eD_split h 3
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, zero_add] at this
  exact_mod_cast congrArg (Nat.cast : ℕ → ℤ) this

lemma eD_split_four (h : Split G S S₁ S₂) :
    (eD G S 4 : ℤ) = (eD G S₁ 0 : ℤ) * (eD G S₂ 4 : ℤ) + (eD G S₁ 1 : ℤ) * (eD G S₂ 3 : ℤ)
      + (eD G S₁ 2 : ℤ) * (eD G S₂ 2 : ℤ) + (eD G S₁ 3 : ℤ) * (eD G S₂ 1 : ℤ)
      + (eD G S₁ 4 : ℤ) * (eD G S₂ 0 : ℤ) := by
  have := eD_split h 4
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, zero_add] at this
  exact_mod_cast congrArg (Nat.cast : ℕ → ℤ) this

lemma iD_split_three (h : Split G S S₁ S₂) :
    (iD G S 3 : ℤ) = (iD G S₁ 0 : ℤ) * (iD G S₂ 3 : ℤ) + (iD G S₁ 1 : ℤ) * (iD G S₂ 2 : ℤ)
      + (iD G S₁ 2 : ℤ) * (iD G S₂ 1 : ℤ) + (iD G S₁ 3 : ℤ) * (iD G S₂ 0 : ℤ) := by
  have := iD_split h 3
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, zero_add] at this
  exact_mod_cast congrArg (Nat.cast : ℕ → ℤ) this

lemma iD_split_four (h : Split G S S₁ S₂) :
    (iD G S 4 : ℤ) = (iD G S₁ 0 : ℤ) * (iD G S₂ 4 : ℤ) + (iD G S₁ 1 : ℤ) * (iD G S₂ 3 : ℤ)
      + (iD G S₁ 2 : ℤ) * (iD G S₂ 2 : ℤ) + (iD G S₁ 3 : ℤ) * (iD G S₂ 1 : ℤ)
      + (iD G S₁ 4 : ℤ) * (iD G S₂ 0 : ℤ) := by
  have := iD_split h 4
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, zero_add] at this
  exact_mod_cast congrArg (Nat.cast : ℕ → ℤ) this

/-- The margin of a disjoint union at defect three. -/
lemma margin3_split (h : Split G S S₁ S₂) (h1 : Eligible G S₁) (h2 : Eligible G S₂) :
    margin G S 3 =
      (eD G S₁ 0 : ℤ) * margin G S₂ 3 + (eD G S₂ 0 : ℤ) * margin G S₁ 3
        + (eD G S₁ 1 : ℤ) * (eD G S₂ 2 : ℤ) + (eD G S₁ 2 : ℤ) * (eD G S₂ 1 : ℤ) := by
  have e3 := eD_split_three h
  have i3 := iD_split_three h
  rw [cast_iD, cast_iD, cast_iD, cast_iD, cast_iD, cast_iD, cast_iD, cast_iD, cast_iD] at i3
  rw [bD_zero, h1.2.1, h1.2.2, bD_zero, h2.2.1, h2.2.2] at i3
  simp only [Nat.cast_zero, add_zero] at i3
  simp only [margin]
  linarith [i3, e3]

/-- The margin of a disjoint union at defect four. -/
lemma margin4_split (h : Split G S S₁ S₂) (h1 : Eligible G S₁) (h2 : Eligible G S₂) :
    margin G S 4 =
      (eD G S₁ 0 : ℤ) * margin G S₂ 4 + (eD G S₂ 0 : ℤ) * margin G S₁ 4
        + (eD G S₁ 1 : ℤ) * margin G S₂ 3 + (eD G S₂ 1 : ℤ) * margin G S₁ 3
        + (eD G S₁ 2 : ℤ) * (eD G S₂ 2 : ℤ) := by
  have e4 := eD_split_four h
  have i4 := iD_split_four h
  rw [cast_iD, cast_iD, cast_iD, cast_iD, cast_iD, cast_iD, cast_iD, cast_iD, cast_iD, cast_iD,
    cast_iD] at i4
  rw [bD_zero, h1.2.1, h1.2.2, bD_zero, h2.2.1, h2.2.2] at i4
  simp only [Nat.cast_zero, add_zero] at i4
  simp only [margin]
  linarith [i4, e4]

/-- **Product closure.**  A disjoint union of two eligible parts, each of which either
satisfies both fivefold bounds or has the exceptional star profile, satisfies both
fivefold bounds. -/
lemma split_good (h : Split G S S₁ S₂) (h1 : Eligible G S₁) (h2 : Eligible G S₂)
    (hf1 : S₁.card ≤ 2 * alpha G S₁) (hf2 : S₂.card ≤ 2 * alpha G S₂)
    (g1 : Good G S₁ ∨ StarProf G S₁) (g2 : Good G S₂ ∨ StarProf G S₂) : Good G S := by
  have m3 := margin3_split h h1 h2
  have m4 := margin4_split h h1 h2
  -- nonnegativity and the incidence bounds for the two parts
  have p10 : (1 : ℤ) ≤ (eD G S₁ 0 : ℤ) := by exact_mod_cast eD_pos (G := G) (S := S₁) (Nat.zero_le _)
  have p20 : (1 : ℤ) ≤ (eD G S₂ 0 : ℤ) := by exact_mod_cast eD_pos (G := G) (S := S₂) (Nat.zero_le _)
  have q1 : (eD G S₁ 0 : ℤ) ≤ 2 * (eD G S₁ 1 : ℤ) := by
    exact_mod_cast e0_le_two_e1 (G := G) (S := S₁) h1.1 hf1
  have q2 : (eD G S₂ 0 : ℤ) ≤ 2 * (eD G S₂ 1 : ℤ) := by
    exact_mod_cast e0_le_two_e1 (G := G) (S := S₂) h2.1 hf2
  have r1 : (eD G S₁ 1 : ℤ) ≤ (eD G S₁ 0 : ℤ) + 6 * (eD G S₁ 2 : ℤ) := by
    exact_mod_cast e1_le_e0_add_six_e2 (G := G) (S := S₁) h1.1 hf1
  have r2 : (eD G S₂ 1 : ℤ) ≤ (eD G S₂ 0 : ℤ) + 6 * (eD G S₂ 2 : ℤ) := by
    exact_mod_cast e1_le_e0_add_six_e2 (G := G) (S := S₂) h2.1 hf2
  have n11 : (0 : ℤ) ≤ (eD G S₁ 1 : ℤ) := Int.natCast_nonneg _
  have n12 : (0 : ℤ) ≤ (eD G S₁ 2 : ℤ) := Int.natCast_nonneg _
  have n21 : (0 : ℤ) ≤ (eD G S₂ 1 : ℤ) := Int.natCast_nonneg _
  have n22 : (0 : ℤ) ≤ (eD G S₂ 2 : ℤ) := Int.natCast_nonneg _
  constructor
  · rcases g1 with ⟨a3, _⟩ | s1 <;> rcases g2 with ⟨b3, _⟩ | s2
    · nlinarith
    · rw [s2.margin3] at m3
      rw [s2.2.1, s2.2.2.1, s2.2.2.2.1] at m3
      push_cast at m3
      nlinarith
    · rw [s1.margin3] at m3
      rw [s1.2.1, s1.2.2.1, s1.2.2.2.1] at m3
      push_cast at m3
      nlinarith
    · rw [s1.margin3, s2.margin3] at m3
      rw [s1.2.1, s1.2.2.1, s1.2.2.2.1, s2.2.1, s2.2.2.1, s2.2.2.2.1] at m3
      push_cast at m3
      linarith
  · rcases g1 with ⟨a3, a4⟩ | s1 <;> rcases g2 with ⟨b3, b4⟩ | s2
    · nlinarith
    · rw [s2.margin3, s2.margin4] at m4
      rw [s2.2.1, s2.2.2.1, s2.2.2.2.1] at m4
      push_cast at m4
      nlinarith
    · rw [s1.margin3, s1.margin4] at m4
      rw [s1.2.1, s1.2.2.1, s1.2.2.2.1] at m4
      push_cast at m4
      nlinarith
    · rw [s1.margin3, s1.margin4, s2.margin3, s2.margin4] at m4
      rw [s1.2.1, s1.2.2.1, s1.2.2.2.1, s2.2.1, s2.2.2.1, s2.2.2.2.1] at m4
      push_cast at m4
      linarith

end FivefoldForest
