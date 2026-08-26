import RequestProject.Defs

/-!
# The `2 × 2` Cauchy–Binet identity and its positivity consequence

For four families `u v x y : ι → ℝ` and a finite index set `s`,

  `(∑ u p x p)(∑ v q y q) - (∑ u p y p)(∑ v q x q)
      = ∑_{p < q} (u p v q - u q v p) (x p y q - x q y p)`.

This is the `2 × 2` case of the Cauchy–Binet formula; it immediately gives that
a product of matrices with nonnegative `2 × 2` minors again has nonnegative
`2 × 2` minors.
-/

open scoped BigOperators

namespace LC

variable {ι : Type*} [LinearOrder ι]

/-- The `2 × 2` Cauchy–Binet identity. -/
theorem cauchyBinet_two (s : Finset ι) (u v x y : ι → ℝ) :
    (∑ p ∈ s, u p * x p) * (∑ q ∈ s, v q * y q)
      - (∑ p ∈ s, u p * y p) * (∑ q ∈ s, v q * x q)
    = ∑ p ∈ s, ∑ q ∈ s,
        (if p < q then (u p * v q - u q * v p) * (x p * y q - x q * y p) else 0) := by
  set T : ι → ι → ℝ := fun p q => u p * v q * (x p * y q - y p * x q) with hT
  have hdiag : ∀ p, T p p = 0 := by intro p; simp only [hT]; ring
  have hLHS :
      (∑ p ∈ s, u p * x p) * (∑ q ∈ s, v q * y q)
        - (∑ p ∈ s, u p * y p) * (∑ q ∈ s, v q * x q)
      = ∑ p ∈ s, ∑ q ∈ s, T p q := by
    rw [Finset.sum_mul_sum, Finset.sum_mul_sum, ← Finset.sum_sub_distrib]
    refine Finset.sum_congr rfl (fun p _ => ?_)
    rw [← Finset.sum_sub_distrib]
    exact Finset.sum_congr rfl (fun q _ => by simp only [hT]; ring)
  have hpt : ∀ p q : ι,
      (if p < q then (u p * v q - u q * v p) * (x p * y q - x q * y p) else 0)
        = (if p < q then T p q else 0) + (if p < q then T q p else 0) := by
    intro p q
    by_cases h : p < q
    · simp only [if_pos h, hT]; ring
    · simp [h]
  have hswap :
      ∑ p ∈ s, ∑ q ∈ s, (if p < q then T q p else 0)
        = ∑ p ∈ s, ∑ q ∈ s, (if q < p then T p q else 0) := by
    rw [Finset.sum_comm]
  rw [hLHS]
  calc ∑ p ∈ s, ∑ q ∈ s, T p q
      = ∑ p ∈ s, ∑ q ∈ s, ((if p < q then T p q else 0) + (if q < p then T p q else 0)) := by
        refine Finset.sum_congr rfl (fun p _ => Finset.sum_congr rfl (fun q _ => ?_))
        rcases lt_trichotomy p q with h | h | h
        · simp [h, asymm h]
        · subst h; simp [hdiag p]
        · simp [h, asymm h]
    _ = (∑ p ∈ s, ∑ q ∈ s, (if p < q then T p q else 0))
          + ∑ p ∈ s, ∑ q ∈ s, (if q < p then T p q else 0) := by
        rw [← Finset.sum_add_distrib]
        exact Finset.sum_congr rfl (fun p _ => by rw [← Finset.sum_add_distrib])
    _ = (∑ p ∈ s, ∑ q ∈ s, (if p < q then T p q else 0))
          + ∑ p ∈ s, ∑ q ∈ s, (if p < q then T q p else 0) := by rw [hswap]
    _ = ∑ p ∈ s, ∑ q ∈ s,
          (if p < q then (u p * v q - u q * v p) * (x p * y q - x q * y p) else 0) := by
        rw [← Finset.sum_add_distrib]
        refine (Finset.sum_congr rfl (fun p _ => ?_)).symm
        rw [← Finset.sum_add_distrib]
        exact Finset.sum_congr rfl (fun q _ => hpt p q)

/-- Positivity form of `2 × 2` Cauchy–Binet: if all `2 × 2` minors of the two
"factors" are nonnegative, so is the corresponding minor of the "product". -/
theorem cauchyBinet_two_nonneg (s : Finset ι) (u v x y : ι → ℝ)
    (h1 : ∀ p ∈ s, ∀ q ∈ s, p < q → 0 ≤ u p * v q - u q * v p)
    (h2 : ∀ p ∈ s, ∀ q ∈ s, p < q → 0 ≤ x p * y q - x q * y p) :
    (∑ p ∈ s, u p * y p) * (∑ q ∈ s, v q * x q)
      ≤ (∑ p ∈ s, u p * x p) * (∑ q ∈ s, v q * y q) := by
  have h := cauchyBinet_two s u v x y
  have hnn : 0 ≤ ∑ p ∈ s, ∑ q ∈ s,
      (if p < q then (u p * v q - u q * v p) * (x p * y q - x q * y p) else 0) := by
    refine Finset.sum_nonneg (fun p hp => Finset.sum_nonneg (fun q hq => ?_))
    by_cases hpq : p < q
    · simpa [hpq] using mul_nonneg (h1 p hp q hq hpq) (h2 p hp q hq hpq)
    · simp [hpq]
  linarith [h ▸ hnn]

end LC
