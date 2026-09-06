import RequestProject.Codes

/-!
# Section 5: ultra-log-concavity of erasure-shadow profiles is false

The seven-bit code `Ω = {24, 47, 65, 67, 86, 93, 97, 99}` (integers encoding binary words)
has retained-coordinate profile `(p_0, …, p_7) = (1, 14, 80, 187, 220, 145, 52, 8)`.

Its profile *is* ordinary log-concave, but normalized (ultra) log-concavity fails at `k = 6`:
`52^2 * C(7,5) * C(7,7) = 56784 < 56840 = 145 * 8 * C(7,6)^2`.
-/

open scoped BigOperators

namespace MatchingBag

/-- The word of length 7 whose `i`-th bit is the `i`-th binary digit of `n`. -/
def wordOfNat (n : ℕ) : Fin 7 → Bool := fun i => n.testBit i.val

/-- The seven-bit code `Ω = {24, 47, 65, 67, 86, 93, 97, 99}` of Section 5. -/
def counterCode : Finset (Fin 7 → Bool) :=
  ({24, 47, 65, 67, 86, 93, 97, 99} : Finset ℕ).image wordOfNat

set_option maxRecDepth 100000

/-- The code really has eight words. -/
theorem counterCode_card : counterCode.card = 8 := by decide

/-- The retained-coordinate profile of the seven-bit code is `(1,14,80,187,220,145,52,8)`. -/
theorem counterCode_profile :
    codeP counterCode 0 = 1 ∧ codeP counterCode 1 = 14 ∧ codeP counterCode 2 = 80 ∧
      codeP counterCode 3 = 187 ∧ codeP counterCode 4 = 220 ∧ codeP counterCode 5 = 145 ∧
      codeP counterCode 6 = 52 ∧ codeP counterCode 7 = 8 := by
  refine ⟨by decide, by decide, by decide, by decide, by decide, by decide, by decide, by decide⟩

/-- Normalized (ultra) log-concavity fails at `k = 6`:
`p_6^2 / C(7,6)^2 < (p_5 / C(7,5)) * (p_7 / C(7,7))`. -/
theorem counterCode_not_ultraLogConcave :
    (codeP counterCode 6) ^ 2 * (Nat.choose 7 5 * Nat.choose 7 7)
      < codeP counterCode 5 * codeP counterCode 7 * (Nat.choose 7 6) ^ 2 := by
  decide

/-- Ordinary log-concavity nevertheless holds for this profile. -/
theorem counterCode_logConcave (k : ℕ) (hk : 1 ≤ k) (hk' : k ≤ 6) :
    codeP counterCode (k - 1) * codeP counterCode (k + 1) ≤ (codeP counterCode k) ^ 2 := by
  interval_cases k <;> decide

end MatchingBag
