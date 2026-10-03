import Mathlib

/-!
# Endpoint-aware curvature propagation for Poisson-binomial reserve

This is the first formalization milestone for the universal finite
Poisson-binomial effective-drop theorem.  It deliberately starts after the
Hillion--Johnson cubic inequalities: `hstep` below is the normalized
consequence of that external input.

The target is a truthful, reusable Lean verification of the recurrence,
endpoint exclusion, first-crossing algebra, and raw-from-effective corollary.
-/

noncomputable section

namespace PBReserve

/-- The propagated upper bound for curvature at distance `r`. -/
def curvatureBound (d : ℝ) (r : ℕ) : ℝ :=
  d / (1 - (r : ℝ) * d)

theorem curvatureBound_zero (d : ℝ) :
    curvatureBound d 0 = d := by
  simp [curvatureBound]

/-- The exact algebraic induction step behind curvature propagation. -/
theorem curvatureBound_succ_identity
    (d : ℝ) (r : ℕ)
    (hden : 1 - (r : ℝ) * d ≠ 0)
    (hnext : 1 - ((r + 1 : ℕ) : ℝ) * d ≠ 0) :
    curvatureBound d r / (1 - curvatureBound d r) =
      curvatureBound d (r + 1) := by
  unfold curvatureBound
  push_cast at hnext ⊢
  have key :
      (1 : ℝ) - d / (1 - (r : ℝ) * d) =
        (1 - ((r : ℝ) + 1) * d) / (1 - (r : ℝ) * d) := by
    rw [eq_div_iff hden, sub_mul, one_mul, div_mul_cancel₀ _ hden]
    ring
  rw [key, div_div_div_cancel_right₀]
  exact hden

theorem curvatureBound_lt_one
    (d : ℝ) (r : ℕ)
    (hd : 0 ≤ d)
    (hr : (((r + 1 : ℕ) : ℝ) * d) < 1) :
    curvatureBound d r < 1 := by
  unfold curvatureBound
  push_cast at hr
  have hden : 0 < 1 - (r : ℝ) * d := by
    nlinarith
  rw [div_lt_one hden]
  nlinarith

/-- One normalized cubic step, with a positive denominator. -/
theorem curvature_step
    (current next : ℝ)
    (hcurrent_lt : current < 1)
    (hstep : next * (1 - current) ≤ current) :
    next ≤ current / (1 - current) := by
  rw [le_div_iff₀ (by linarith)]
  linarith

/--
If neighboring curvatures obey the normalized Hillion--Johnson recurrence,
then curvature at distance `r` is at most `d / (1-r*d)` for as long as
`(r+1)*d < 1`.
-/
theorem curvature_propagation
    (κ : ℕ → ℝ) (d : ℝ)
    (hd : 0 ≤ d)
    (hκ_nonneg : ∀ r, 0 ≤ κ r)
    (hbase : κ 0 = d)
    (hstep : ∀ r, κ (r + 1) * (1 - κ r) ≤ κ r)
    (r : ℕ)
    (hr : (((r + 1 : ℕ) : ℝ) * d) < 1) :
    κ r ≤ curvatureBound d r := by
  induction r with
  | zero => rw [curvatureBound_zero, hbase]
  | succ n ih =>
    push_cast at hr
    have hn1 : (((n + 1 : ℕ) : ℝ) * d) < 1 := by
      push_cast
      nlinarith
    have ihn := ih hn1
    have hBrlt : curvatureBound d n < 1 :=
      curvatureBound_lt_one d n hd hn1
    have hκn_lt : κ n < 1 := lt_of_le_of_lt ihn hBrlt
    have hstepn := curvature_step (κ n) (κ (n + 1)) hκn_lt (hstep n)
    have hmono :
        κ n / (1 - κ n) ≤
          curvatureBound d n / (1 - curvatureBound d n) := by
      rw [div_le_div_iff₀ (by linarith) (by linarith)]
      nlinarith [hκ_nonneg n]
    have hden : (1 : ℝ) - (n : ℝ) * d ≠ 0 := by
      nlinarith
    have hnext : (1 : ℝ) - ((n + 1 : ℕ) : ℝ) * d ≠ 0 := by
      push_cast
      nlinarith
    have hid := curvatureBound_succ_identity d n hden hnext
    calc
      κ (n + 1) ≤ κ n / (1 - κ n) := hstepn
      _ ≤ curvatureBound d n / (1 - curvatureBound d n) := hmono
      _ = curvatureBound d (n + 1) := hid

/--
An endpoint has curvature one, so it cannot occur inside the propagated
window.  This is the endpoint convention that an abstract truncated
recurrence would miss.
-/
theorem endpoint_exclusion
    (κ : ℕ → ℝ) (d : ℝ)
    (hd : 0 ≤ d)
    (hκ_nonneg : ∀ r, 0 ≤ κ r)
    (hbase : κ 0 = d)
    (hstep : ∀ r, κ (r + 1) * (1 - κ r) ≤ κ r)
    (r : ℕ)
    (hr : (((r + 1 : ℕ) : ℝ) * d) < 1)
    (hendpoint : κ r = 1) :
    False := by
  have h1 := curvature_propagation κ d hd hκ_nonneg hbase hstep r hr
  have h2 := curvatureBound_lt_one d r hd hr
  rw [hendpoint] at h1
  linarith

/-- The first-crossing ratio bound used to seed the right mass window. -/
theorem crossing_ratio_lower_bound
    (d qPrev qNext δPrev : ℝ)
    (hd : 0 ≤ d)
    (hd_quarter : d < 1 / 4)
    (hqPrev : 1 ≤ qPrev)
    (hδPrev_nonneg : 0 ≤ δPrev)
    (hδPrev : δPrev ≤ d / (1 - d))
    (hqNext : qNext = qPrev * (1 - δPrev)) :
    (1 - 2 * d) / (1 - d) ≤ qNext := by
  have hden : 0 < 1 - d := by
    linarith
  have h1 : d / (1 - d) < 1 := by
    rw [div_lt_one hden]
    linarith
  have hpos : 0 < 1 - δPrev := by
    linarith
  have hfac : (1 - 2 * d) / (1 - d) = 1 - d / (1 - d) := by
    field_simp
    ring
  rw [hqNext]
  calc
    (1 - 2 * d) / (1 - d) = 1 - d / (1 - d) := hfac
    _ ≤ 1 - δPrev := by linarith
    _ = 1 * (1 - δPrev) := by ring
    _ ≤ qPrev * (1 - δPrev) :=
      mul_le_mul_of_nonneg_right hqPrev (le_of_lt hpos)

/--
At a strict descent, the raw drop is at least the effective Turán drop.
The positivity and indexing hypotheses are explicit to prevent the common
off-by-one error.
-/
theorem raw_drop_ge_effective
    (fPrev fNow fNext : ℝ)
    (hfPrev : 0 < fPrev)
    (hfNow : 0 < fNow)
    (hfNext : 0 ≤ fNext)
    (hdescent : fNow < fPrev) :
    1 - fNext / fNow ≥
      1 - (fPrev * fNext) / fNow ^ 2 := by
  have h : fNext / fNow ≤ (fPrev * fNext) / fNow ^ 2 := by
    rw [div_le_div_iff₀ hfNow (by positivity)]
    nlinarith [
      mul_nonneg
        (mul_nonneg hfNow.le hfNext)
        (sub_nonneg.mpr hdescent.le)
    ]
  linarith

theorem raw_quarter_of_effective
    (V fPrev fNow fNext : ℝ)
    (hV : 0 ≤ V)
    (hfPrev : 0 < fPrev)
    (hfNow : 0 < fNow)
    (hfNext : 0 ≤ fNext)
    (hdescent : fNow < fPrev)
    (heffective :
      1 / 4 ≤ V * (1 - (fPrev * fNext) / fNow ^ 2)) :
    1 / 4 ≤ V * (1 - fNext / fNow) := by
  have h :=
    raw_drop_ge_effective fPrev fNow fNext hfPrev hfNow hfNext hdescent
  have hle :
      V * (1 - (fPrev * fNext) / fNow ^ 2) ≤
        V * (1 - fNext / fNow) :=
    mul_le_mul_of_nonneg_left h hV
  linarith

end PBReserve
