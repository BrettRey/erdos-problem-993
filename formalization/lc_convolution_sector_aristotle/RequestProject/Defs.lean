import Mathlib

/-!
# Log-concavity of sequences: definitions

Mathlib (as of the pinned checkout) has no notion of log-concavity for sequences,
so we introduce the four notions requested in the formalization packet.

We work with sequences `f : ℕ → ℝ` (the "sequence form"); polynomial corollaries
are derived from the sequence statements via `Polynomial.coeff`.
-/

open scoped BigOperators

namespace LC

/-- `Nonneg f` : every entry of the sequence is nonnegative. -/
def Nonneg (f : ℕ → ℝ) : Prop := ∀ k, 0 ≤ f k

/-- `LogConcave f` : `f (k+1) ^ 2 ≥ f k * f (k+2)` for every `k`. -/
def LogConcave (f : ℕ → ℝ) : Prop := ∀ k, f (k + 1) * f (k + 1) ≥ f k * f (k + 2)

/-- `NoInternalZeros f` : the support of `f` is "interval-like". -/
def NoInternalZeros (f : ℕ → ℝ) : Prop :=
  ∀ i j k, i ≤ j → j ≤ k → f i ≠ 0 → f k ≠ 0 → f j ≠ 0

/-- Discrete convolution of two sequences. -/
def Conv (f g : ℕ → ℝ) : ℕ → ℝ := fun n => ∑ i ∈ Finset.range (n + 1), f i * g (n - i)

/-- Reading a finite vector `![a, b, c]` as a sequence `ℕ → ℝ` padded with zeros.
This is only used to give a literal meaning to statements such as
`LogConcave ![w, 1]`, whose vectors are `Fin n`-indexed. -/
def ofVec {n : ℕ} (v : Fin n → ℝ) : ℕ → ℝ := fun k => if h : k < n then v ⟨k, h⟩ else 0

@[simp] lemma ofVec_apply {n : ℕ} (v : Fin n → ℝ) (k : ℕ) (h : k < n) :
    ofVec v k = v ⟨k, h⟩ := dif_pos h

@[simp] lemma ofVec_apply_of_le {n : ℕ} (v : Fin n → ℝ) (k : ℕ) (h : n ≤ k) :
    ofVec v k = 0 := dif_neg (by omega)

@[simp] lemma ofVec_two_zero (a b : ℝ) : ofVec ![a, b] 0 = a := by simp [ofVec]

@[simp] lemma ofVec_two_one (a b : ℝ) : ofVec ![a, b] 1 = b := by simp [ofVec]

@[simp] lemma ofVec_three_zero (a b c : ℝ) : ofVec ![a, b, c] 0 = a := by simp [ofVec]

@[simp] lemma ofVec_three_one (a b c : ℝ) : ofVec ![a, b, c] 1 = b := by simp [ofVec]

@[simp] lemma ofVec_three_two (a b c : ℝ) : ofVec ![a, b, c] 2 = c := by simp [ofVec]

end LC
