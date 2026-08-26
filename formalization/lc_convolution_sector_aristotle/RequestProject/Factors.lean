import RequestProject.Conv

/-!
# G1 : the two real factor shapes

A monic real polynomial factors over `ℝ` into linear factors `x + w` and
quadratic factors `x² + 2ρ cos φ · x + ρ²`.  Here we determine exactly when the
coefficient sequences `![w, 1]` and `![ρ², 2ρ cos φ, 1]` are log-concave.

The vectors are `Fin n`-indexed, so they are read as sequences `ℕ → ℝ` through
`LC.ofVec` (padding with zeros).
-/

open scoped BigOperators

namespace LC

/-! ### Linear factors -/

/-- **G1, linear factor.** The coefficient sequence of `x + w` (with `w > 0`) is
nonnegative and log-concave. -/
theorem logConcave_linear (w : ℝ) (hw : 0 < w) :
    LogConcave (ofVec ![w, 1]) ∧ Nonneg (ofVec ![w, 1]) := by
  constructor
  · intro k
    match k with
    | 0 => simp
    | (n + 1) => simp
  · intro k
    match k with
    | 0 => simpa using hw.le
    | 1 => simp
    | (n + 2) => simp

/-- The coefficient sequence of `x + w` (with `w ≠ 0`) has no internal zeros. -/
theorem noInternalZeros_linear (w : ℝ) (hw : w ≠ 0) :
    NoInternalZeros (ofVec ![w, 1]) := by
  intro i j k hij hjk hi hk
  have hk2 : k ≤ 1 := by
    by_contra h
    exact hk (ofVec_apply_of_le _ _ (by omega))
  match j, (show j ≤ 1 by omega) with
  | 0, _ => simpa using hw
  | 1, _ => simp

/-! ### Quadratic factors -/

/-- **G1, quadratic factor.** The coefficient sequence of `x² + 2ρ cos φ · x + ρ²`
is log-concave *exactly* when `cos φ ^ 2 ≥ 1/4`, i.e. `|cos φ| ≥ 1/2`.
Since `cos (π/3) = 1/2`, `π/3` is precisely the threshold. -/
theorem logConcave_quadratic (rho phi : ℝ) (hrho : 0 < rho) :
    LogConcave (ofVec ![rho ^ 2, 2 * rho * Real.cos phi, 1]) ↔ Real.cos phi ^ 2 ≥ 1 / 4 := by
  have hr2 : (0 : ℝ) < rho ^ 2 := by positivity
  constructor
  · intro h
    have h0 := h 0
    simp only [ge_iff_le, zero_add, ofVec_three_zero, ofVec_three_one, ofVec_three_two] at h0
    nlinarith [h0, hr2]
  · intro hc k
    match k with
    | 0 =>
      simp only [ge_iff_le, zero_add, ofVec_three_zero, ofVec_three_one, ofVec_three_two]
      nlinarith [hc, hr2]
    | 1 => simp
    | (n + 2) => simp

/-- The coefficient sequence of `x² + 2ρ cos φ · x + ρ²` is nonnegative when
`ρ > 0` and `cos φ ≥ 0`. -/
theorem nonneg_quadratic (rho phi : ℝ) (hrho : 0 < rho) (hc : 0 ≤ Real.cos phi) :
    Nonneg (ofVec ![rho ^ 2, 2 * rho * Real.cos phi, 1]) := by
  intro k
  match k with
  | 0 => simpa using sq_nonneg rho
  | 1 => simpa using by positivity
  | 2 => simp
  | (n + 3) => simp

/-- The coefficient sequence of `x² + 2ρ cos φ · x + ρ²` has no internal zeros
when `ρ ≠ 0` and `cos φ ≠ 0`. -/
theorem noInternalZeros_quadratic (rho phi : ℝ) (hrho : rho ≠ 0) (hc : Real.cos phi ≠ 0) :
    NoInternalZeros (ofVec ![rho ^ 2, 2 * rho * Real.cos phi, 1]) := by
  intro i j k hij hjk hi hk
  have hk3 : k ≤ 2 := by
    by_contra h
    exact hk (ofVec_apply_of_le _ _ (by omega))
  match j, (show j ≤ 2 by omega) with
  | 0, _ => simpa using pow_ne_zero 2 hrho
  | 1, _ => simpa using ⟨hrho, hc⟩
  | 2, _ => simp

/-! ### The `π/3` sector bound for a single quadratic factor -/

lemma half_le_cos_of_abs_le_pi_div_three {phi : ℝ} (h : |phi| ≤ Real.pi / 3) :
    1 / 2 ≤ Real.cos phi := by
  have h1 : Real.cos (Real.pi / 3) ≤ Real.cos |phi| :=
    Real.cos_le_cos_of_nonneg_of_le_pi (abs_nonneg _) (by linarith [Real.pi_pos]) h
  rw [Real.cos_abs, Real.cos_pi_div_three] at h1
  linarith

/-- For a root at angle `φ` from the negative real axis with `|φ| ≤ π/3`, the
quadratic factor `x² + 2ρ cos φ · x + ρ²` has a nonnegative, log-concave
coefficient sequence with no internal zeros. -/
theorem logConcave_quadratic_of_abs_le_pi_div_three (rho phi : ℝ) (hrho : 0 < rho)
    (hphi : |phi| ≤ Real.pi / 3) :
    Nonneg (ofVec ![rho ^ 2, 2 * rho * Real.cos phi, 1])
      ∧ LogConcave (ofVec ![rho ^ 2, 2 * rho * Real.cos phi, 1])
      ∧ NoInternalZeros (ofVec ![rho ^ 2, 2 * rho * Real.cos phi, 1]) := by
  have hc : 1 / 2 ≤ Real.cos phi := half_le_cos_of_abs_le_pi_div_three hphi
  refine ⟨nonneg_quadratic rho phi hrho (by linarith), ?_,
    noInternalZeros_quadratic rho phi (ne_of_gt hrho) (by linarith)⟩
  rw [logConcave_quadratic rho phi hrho]
  nlinarith

end LC
