import RequestProject.Codes

/-!
# Exhaustive checks for short binary codes

Section 5 of `matching_bag_poset_reduction_2026-08-20.md` reports that ordinary *and*
normalized (ultra) log-concavity hold for every nonempty binary code of length at most four.
This file verifies that claim by exhaustive enumeration.

The length-four check enumerates all `2^16 - 1` nonempty codes on `16` words and is
correspondingly expensive; it is discharged by `native_decide` (kernel reduction is not
feasible at that size), so it additionally depends on the `Lean.ofReduceBool` axiom.
-/

open scoped BigOperators
open Finset

namespace MatchingBag

/-- Ordinary log-concavity of the retained-coordinate profile of a code of length `m`. -/
def LogConcaveProfile {m : ℕ} (Ω : Finset (Fin m → Bool)) : Prop :=
  ∀ k ∈ Finset.Icc 1 (m - 1), codeP Ω (k - 1) * codeP Ω (k + 1) ≤ (codeP Ω k) ^ 2

/-- Normalized (ultra) log-concavity of the retained-coordinate profile of a code of
length `m`: `(p_k / C(m,k))^2 ≥ (p_{k-1} / C(m,k-1)) * (p_{k+1} / C(m,k+1))`. -/
def NormLogConcaveProfile {m : ℕ} (Ω : Finset (Fin m → Bool)) : Prop :=
  ∀ k ∈ Finset.Icc 1 (m - 1),
    codeP Ω (k - 1) * codeP Ω (k + 1) * (m.choose k) ^ 2
      ≤ (codeP Ω k) ^ 2 * (m.choose (k - 1) * m.choose (k + 1))

instance {m : ℕ} (Ω : Finset (Fin m → Bool)) : Decidable (LogConcaveProfile Ω) := by
  unfold LogConcaveProfile; infer_instance

instance {m : ℕ} (Ω : Finset (Fin m → Bool)) : Decidable (NormLogConcaveProfile Ω) := by
  unfold NormLogConcaveProfile; infer_instance

set_option maxRecDepth 100000

/-- Every nonempty binary code of length `1` has a log-concave and normalized log-concave
profile. -/
theorem logConcave_length_one (Ω : Finset (Fin 1 → Bool)) (hΩ : Ω.Nonempty) :
    LogConcaveProfile Ω ∧ NormLogConcaveProfile Ω := by
  revert hΩ; revert Ω; decide

/-- Every nonempty binary code of length `2` has a log-concave and normalized log-concave
profile. -/
theorem logConcave_length_two (Ω : Finset (Fin 2 → Bool)) (hΩ : Ω.Nonempty) :
    LogConcaveProfile Ω ∧ NormLogConcaveProfile Ω := by
  revert hΩ; revert Ω; decide

/-- Every nonempty binary code of length `3` has a log-concave and normalized log-concave
profile. -/
theorem logConcave_length_three (Ω : Finset (Fin 3 → Bool)) (hΩ : Ω.Nonempty) :
    LogConcaveProfile Ω ∧ NormLogConcaveProfile Ω := by
  revert hΩ; revert Ω; decide

/-- Every nonempty binary code of length `4` has a log-concave and normalized log-concave
profile. Checked by exhaustive enumeration of all `2^16 - 1` nonempty codes. -/
theorem logConcave_length_four (Ω : Finset (Fin 4 → Bool)) (hΩ : Ω.Nonempty) :
    LogConcaveProfile Ω ∧ NormLogConcaveProfile Ω := by
  revert hΩ; revert Ω; native_decide

end MatchingBag
