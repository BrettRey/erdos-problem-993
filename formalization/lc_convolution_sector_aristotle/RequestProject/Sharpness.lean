import RequestProject.Conv

/-!
# `NoInternalZeros` cannot be dropped from G0

The convolution theorem `LC.logConcave_conv` genuinely needs the
no-internal-zeros hypothesis.  Here is the smallest explicit witness:

* `f = (1, 0, 0, 1, 0, 0, …)`  (coefficients of `1 + x³`) is nonnegative and
  log-concave, but has an internal zero;
* `g = (1, 1, 0, 0, …)`        (coefficients of `1 + x`) is nonnegative,
  log-concave and internal-zero-free;
* `Conv f g = (1, 1, 0, 1, 1, 0, …)` fails log-concavity at `k = 1`, since
  `c₂ * c₂ = 0 < 1 = c₁ * c₃`.
-/

open scoped BigOperators

namespace LC

/-- The coefficient sequence of `1 + x³`. -/
def cexF : ℕ → ℝ := fun k => if k = 0 ∨ k = 3 then 1 else 0

/-- The coefficient sequence of `1 + x`. -/
def cexG : ℕ → ℝ := fun k => if k ≤ 1 then 1 else 0

lemma cexF_nonneg : Nonneg cexF := by
  intro k; unfold cexF; split <;> norm_num

lemma cexG_nonneg : Nonneg cexG := by
  intro k; unfold cexG; split <;> norm_num

lemma cexF_logConcave : LogConcave cexF := by
  intro k
  match k with
  | 0 => norm_num [cexF]
  | 1 => norm_num [cexF]
  | 2 => norm_num [cexF]
  | 3 => norm_num [cexF]
  | (n + 4) =>
    norm_num [cexF, show ¬(n + 4 = 3) by omega, show ¬(n + 4 + 1 = 3) by omega,
      show ¬(n + 4 + 2 = 3) by omega]

lemma cexG_logConcave : LogConcave cexG := by
  intro k
  match k with
  | 0 => norm_num [cexG]
  | (n + 1) =>
    norm_num [cexG, show ¬(n + 1 + 1 ≤ 1) by omega, show ¬(n + 1 + 2 ≤ 1) by omega]

lemma cexG_noInternalZeros : NoInternalZeros cexG := by
  intro i j k hij hjk _ hk
  have hk1 : k ≤ 1 := by
    by_contra h
    exact hk (by simp [cexG, show ¬ (k ≤ 1) by omega])
  simp [cexG, show j ≤ 1 by omega]

/-- `f = 1 + x³` has an internal zero. -/
lemma cexF_not_noInternalZeros : ¬ NoInternalZeros cexF := by
  intro h
  have := h 0 1 3 (by omega) (by omega) (by norm_num [cexF]) (by norm_num [cexF])
  exact this (by norm_num [cexF])

/-- **G0 is false without `NoInternalZeros`.**  With `f` and `g` as above, all
of the other hypotheses hold, yet `Conv f g` is not log-concave (it fails at
index `k = 1`). -/
theorem conv_not_logConcave_without_noInternalZeros :
    Nonneg cexF ∧ Nonneg cexG ∧ LogConcave cexF ∧ LogConcave cexG
      ∧ NoInternalZeros cexG ∧ ¬ NoInternalZeros cexF
      ∧ ¬ LogConcave (Conv cexF cexG) := by
  refine ⟨cexF_nonneg, cexG_nonneg, cexF_logConcave, cexG_logConcave, cexG_noInternalZeros,
    cexF_not_noInternalZeros, ?_⟩
  intro h
  have hc1 : Conv cexF cexG 1 = 1 := by
    simp only [Conv, Finset.sum_range_succ, Finset.sum_range_zero]
    norm_num [cexF, cexG]
  have hc2 : Conv cexF cexG 2 = 0 := by
    simp only [Conv, Finset.sum_range_succ, Finset.sum_range_zero]
    norm_num [cexF, cexG]
  have hc3 : Conv cexF cexG 3 = 1 := by
    simp only [Conv, Finset.sum_range_succ, Finset.sum_range_zero]
    norm_num [cexF, cexG]
  have h1 := h 1
  rw [show (1 : ℕ) + 1 = 2 from rfl, show (1 : ℕ) + 2 = 3 from rfl, hc1, hc2, hc3] at h1
  norm_num at h1

end LC
