import RequestProject.Factors

/-!
# G2 : the π/3 root-sector criterion

A real polynomial with positive coefficients all of whose complex roots lie
within `π/3` of the negative real axis has a log-concave coefficient sequence.

The argument is the one sketched in the request: factor over `ℝ` into linear
factors `X + w` (`w > 0`) and quadratic factors `X² + 2ρ cos φ · X + ρ²`
(`|φ| ≤ π/3`), each of which is log-concave by G1, and conclude by G0 and
induction on the degree.
-/

open scoped BigOperators
open Polynomial

namespace LC

/-! ### Coefficients of the two factor shapes -/

lemma natDegree_linear_factor_le (w : ℝ) : (X + C w : ℝ[X]).natDegree ≤ 1 := by
  compute_degree

lemma natDegree_quad_factor (b c : ℝ) : (X ^ 2 + C b * X + C c : ℝ[X]).natDegree = 2 := by
  compute_degree!

lemma coeff_linear_factor (w : ℝ) : (fun k => (X + C w).coeff k) = ofVec ![w, 1] := by
  funext k
  match k with
  | 0 => simp
  | 1 => simp
  | (n + 2) =>
    rw [ofVec_apply_of_le _ _ (by omega)]
    exact Polynomial.coeff_eq_zero_of_natDegree_lt
      (lt_of_le_of_lt (natDegree_linear_factor_le w) (by omega))

lemma coeff_quad_factor (b c : ℝ) :
    (fun k => (X ^ 2 + C b * X + C c).coeff k) = ofVec ![c, b, 1] := by
  funext k
  match k with
  | 0 => simp
  | 1 => simp
  | 2 => simp
  | (n + 3) =>
    rw [ofVec_apply_of_le _ _ (by omega)]
    exact Polynomial.coeff_eq_zero_of_natDegree_lt
      (by rw [natDegree_quad_factor]; omega)

/-! ### Polynomials with positive coefficients -/

/-- All coefficients up to the degree are strictly positive. -/
def PosCoeffs (p : Polynomial ℝ) : Prop := ∀ k ≤ p.natDegree, 0 < p.coeff k

lemma posCoeffs_nonneg {p : Polynomial ℝ} (hp : PosCoeffs p) : Nonneg (fun k => p.coeff k) := by
  intro k
  show (0 : ℝ) ≤ p.coeff k
  by_cases hk : k ≤ p.natDegree
  · exact (hp k hk).le
  · rw [Polynomial.coeff_eq_zero_of_natDegree_lt (show p.natDegree < k by omega)]

lemma posCoeffs_noInternalZeros {p : Polynomial ℝ} (hp : PosCoeffs p) :
    NoInternalZeros (fun k => p.coeff k) := by
  intro i j k hij hjk _ hk
  have hkd : k ≤ p.natDegree := by
    by_contra h
    exact hk (Polynomial.coeff_eq_zero_of_natDegree_lt (by omega))
  exact ne_of_gt (hp j (by omega))

lemma posCoeffs_mul {p q : Polynomial ℝ} (hp0 : p ≠ 0) (hq0 : q ≠ 0)
    (hp : PosCoeffs p) (hq : PosCoeffs q) : PosCoeffs (p * q) := by
  intro k hk
  rw [Polynomial.natDegree_mul hp0 hq0] at hk
  rw [congrFun (coeff_mul_eq_conv p q) k]
  set i0 := min k p.natDegree with hi0
  have h1 : i0 ≤ p.natDegree := min_le_right _ _
  have h2 : k - i0 ≤ q.natDegree := by simp only [hi0]; omega
  have hterm : 0 < p.coeff i0 * q.coeff (k - i0) := mul_pos (hp i0 h1) (hq _ h2)
  have hle : p.coeff i0 * q.coeff (k - i0) ≤ Conv (fun n => p.coeff n) (fun n => q.coeff n) k :=
    Finset.single_le_sum (f := fun i => p.coeff i * q.coeff (k - i))
      (fun i _ => mul_nonneg (posCoeffs_nonneg hp i) (posCoeffs_nonneg hq _))
      (by simp only [Finset.mem_range, hi0]; omega)
  linarith

lemma posCoeffs_linear_factor {w : ℝ} (hw : 0 < w) : PosCoeffs (X + C w) := by
  intro k hk
  rw [Polynomial.natDegree_X_add_C] at hk
  rw [congrFun (coeff_linear_factor w) k]
  interval_cases k <;> simp [hw]

lemma posCoeffs_quad_factor {b c : ℝ} (hb : 0 < b) (hc : 0 < c) :
    PosCoeffs (X ^ 2 + C b * X + C c) := by
  intro k hk
  rw [natDegree_quad_factor] at hk
  rw [congrFun (coeff_quad_factor b c) k]
  interval_cases k <;> simp [hb, hc]

/-! ### The inductive step -/

/-- If `p` factors as `F * q` with both factors having positive, log-concave
coefficient sequences, then so does `p`. -/
lemma sector_step {p F q : Polynomial ℝ} (hpq : p = F * q)
    (hF0 : F ≠ 0) (hq0 : q ≠ 0)
    (hF : PosCoeffs F) (hFlc : LogConcave (fun k => F.coeff k))
    (hq : PosCoeffs q) (hqlc : LogConcave (fun k => q.coeff k)) :
    PosCoeffs p ∧ LogConcave (fun k => p.coeff k) := by
  subst hpq
  refine ⟨posCoeffs_mul hF0 hq0 hF hq, ?_⟩
  exact (logConcave_coeff_mul F q (posCoeffs_nonneg hF) (posCoeffs_nonneg hq) hFlc hqlc
    (posCoeffs_noInternalZeros hF) (posCoeffs_noInternalZeros hq)).2.1

/-! ### The main induction -/

theorem sector_aux : ∀ n : ℕ, ∀ p : Polynomial ℝ, p.natDegree = n → 0 < p.leadingCoeff →
    p.coeff 0 ≠ 0 →
    (∀ z : ℂ, (p.map (algebraMap ℝ ℂ)).IsRoot z → |Complex.arg (-z)| ≤ Real.pi / 3) →
    PosCoeffs p ∧ LogConcave (fun k => p.coeff k) := by
  intro n
  induction n using Nat.strong_induction_on with
  | _ n IH =>
    intro p hdeg hlead h0 hroots
    have hp0 : p ≠ 0 := fun h => h0 (by simp [h])
    rcases Nat.eq_zero_or_pos n with hn | hn
    · -- constant polynomial
      subst hn
      have hzero : ∀ m : ℕ, 1 ≤ m → p.coeff m = 0 := fun m hm =>
        Polynomial.coeff_eq_zero_of_natDegree_lt (by omega)
      refine ⟨fun k hk => ?_, fun k => ?_⟩
      · have : k = 0 := by omega
        subst this
        rw [show (0 : ℕ) = p.natDegree from hdeg.symm]
        exact hlead
      · show p.coeff (k + 1) * p.coeff (k + 1) ≥ p.coeff k * p.coeff (k + 2)
        rw [hzero (k + 1) (by omega), hzero (k + 2) (by omega)]
        simp
    · -- degree at least one: extract a root
      have hdegC : 0 < (p.map (algebraMap ℝ ℂ)).degree := by
        rw [Polynomial.degree_map_eq_of_injective (algebraMap ℝ ℂ).injective,
          Polynomial.degree_eq_natDegree hp0, hdeg]
        exact_mod_cast hn
      obtain ⟨z, hz⟩ := Complex.exists_root hdegC
      have hzne : z ≠ 0 := by
        rintro rfl
        apply h0
        have hc0 : (p.map (algebraMap ℝ ℂ)).coeff 0 = 0 := by
          rw [Polynomial.coeff_zero_eq_eval_zero]; exact hz
        rw [Polynomial.coeff_map] at hc0
        exact (map_eq_zero_iff _ (algebraMap ℝ ℂ).injective).mp hc0
      have harg := hroots z hz
      by_cases him : z.im = 0
      · -- real root
        have hzeq : z = ((z.re : ℝ) : ℂ) := by
          apply Complex.ext <;> simp [him]
        set r : ℝ := z.re with hr
        have hrroot : p.IsRoot r := by
          have h2 : ((p.eval r : ℝ) : ℂ) = 0 := by
            have := hz
            rw [hzeq, Polynomial.IsRoot, Polynomial.eval_map,
              show (((r : ℝ)) : ℂ) = algebraMap ℝ ℂ r from rfl,
              Polynomial.eval₂_at_apply] at this
            exact this
          exact_mod_cast h2
        have hrne : r ≠ 0 := by
          intro h
          exact hzne (by rw [hzeq, h]; simp)
        have hrneg : r < 0 := by
          rcases lt_trichotomy r 0 with h | h | h
          · exact h
          · exact absurd h hrne
          · exfalso
            have : -z = ((-r : ℝ) : ℂ) := by rw [hzeq]; simp
            rw [this, Complex.arg_ofReal_of_neg (by simpa using h)] at harg
            rw [abs_of_pos Real.pi_pos] at harg
            linarith [Real.pi_pos]
        set w : ℝ := -r with hw
        have hwpos : 0 < w := by simp only [hw]; linarith
        have hFeq : (X + C w : Polynomial ℝ) = X - C r := by
          simp only [hw, map_neg]; ring
        obtain ⟨q, hq⟩ := (Polynomial.dvd_iff_isRoot).mpr hrroot
        have hpq : p = (X + C w) * q := by rw [hFeq]; exact hq
        have hF0 : (X + C w : Polynomial ℝ) ≠ 0 := by
          intro h
          have := Polynomial.natDegree_X_add_C w
          rw [h] at this
          simp at this
        have hq0 : q ≠ 0 := by
          intro h; apply hp0; rw [hpq, h, mul_zero]
        have hdq : q.natDegree = n - 1 := by
          have := Polynomial.natDegree_mul hF0 hq0
          rw [← hpq, hdeg, Polynomial.natDegree_X_add_C] at this
          omega
        have hlt : q.natDegree < n := by omega
        have hqlead : 0 < q.leadingCoeff := by
          have := Polynomial.leadingCoeff_mul (X + C w : Polynomial ℝ) q
          rw [← hpq, (Polynomial.monic_X_add_C w).leadingCoeff, one_mul] at this
          rwa [this] at hlead
        have hq00 : q.coeff 0 ≠ 0 := by
          intro h
          apply h0
          rw [hpq, Polynomial.mul_coeff_zero, h, mul_zero]
        have hqroots : ∀ z' : ℂ, (q.map (algebraMap ℝ ℂ)).IsRoot z' →
            (p.map (algebraMap ℝ ℂ)).IsRoot z' := by
          intro z' hz'
          have hz'' : Polynomial.eval z' (q.map (algebraMap ℝ ℂ)) = 0 := hz'
          show Polynomial.eval z' (p.map (algebraMap ℝ ℂ)) = 0
          rw [hpq, Polynomial.map_mul, Polynomial.eval_mul, hz'', mul_zero]
        obtain ⟨hqpos, hqlc⟩ :=
          IH q.natDegree hlt q rfl hqlead hq00 (fun z' hz' => hroots z' (hqroots z' hz'))
        refine sector_step hpq hF0 hq0 (posCoeffs_linear_factor hwpos) ?_ hqpos hqlc
        rw [coeff_linear_factor w]
        exact (logConcave_linear w hwpos).1
      · -- non-real root: quadratic factor
        set rho : ℝ := ‖z‖ with hrho
        set phi : ℝ := Complex.arg (-z) with hphi
        have hrhopos : 0 < rho := by simpa [hrho] using hzne
        have hcos : Real.cos phi = -z.re / rho := by
          rw [hphi, Complex.cos_arg (neg_ne_zero.mpr hzne)]
          simp [hrho]
        have hbeq : 2 * rho * Real.cos phi = -(2 * z.re) := by
          rw [hcos]
          field_simp
        have hceq : rho ^ 2 = ‖z‖ ^ 2 := by rw [hrho]
        have haeval : (Polynomial.aeval z) p = 0 := by
          rw [Polynomial.aeval_def, ← Polynomial.eval_map]
          exact hz
        have hdvd := Polynomial.quadratic_dvd_of_aeval_eq_zero_im_ne_zero p haeval him
        obtain ⟨q, hq⟩ := hdvd
        set Q : Polynomial ℝ := X ^ 2 + C (2 * rho * Real.cos phi) * X + C (rho ^ 2) with hQ
        have hQeq : Q = X ^ 2 - C (2 * z.re) * X + C (‖z‖ ^ 2) := by
          rw [hQ, hbeq, hceq, map_neg]; ring
        have hpq : p = Q * q := by rw [hQeq]; exact hq
        have hcospos : 1 / 2 ≤ Real.cos phi :=
          half_le_cos_of_abs_le_pi_div_three (by rw [hphi]; exact harg)
        have hbpos : 0 < 2 * rho * Real.cos phi := by nlinarith
        have hcpos : 0 < rho ^ 2 := by positivity
        have hQpos : PosCoeffs Q := posCoeffs_quad_factor hbpos hcpos
        have hQ0 : Q ≠ 0 := by
          intro h
          have := natDegree_quad_factor (2 * rho * Real.cos phi) (rho ^ 2)
          rw [← hQ, h] at this
          simp at this
        have hq0 : q ≠ 0 := by
          intro h; apply hp0; rw [hpq, h, mul_zero]
        have hdq : q.natDegree = n - 2 ∧ 2 ≤ n := by
          have hmul := Polynomial.natDegree_mul hQ0 hq0
          rw [← hpq, hdeg, hQ, natDegree_quad_factor] at hmul
          omega
        have hlt : q.natDegree < n := by omega
        have hqlead : 0 < q.leadingCoeff := by
          have hl := Polynomial.leadingCoeff_mul Q q
          rw [← hpq] at hl
          have hQmonic : Q.leadingCoeff = 1 := by
            have : Q.coeff 2 = 1 := by
              rw [hQ, congrFun (coeff_quad_factor (2 * rho * Real.cos phi) (rho ^ 2)) 2]
              simp
            rw [Polynomial.leadingCoeff, hQ, natDegree_quad_factor, ← hQ, this]
          rw [hQmonic, one_mul] at hl
          rwa [hl] at hlead
        have hq00 : q.coeff 0 ≠ 0 := by
          intro h
          apply h0
          rw [hpq, Polynomial.mul_coeff_zero, h, mul_zero]
        have hqroots : ∀ z' : ℂ, (q.map (algebraMap ℝ ℂ)).IsRoot z' →
            (p.map (algebraMap ℝ ℂ)).IsRoot z' := by
          intro z' hz'
          have hz'' : Polynomial.eval z' (q.map (algebraMap ℝ ℂ)) = 0 := hz'
          show Polynomial.eval z' (p.map (algebraMap ℝ ℂ)) = 0
          rw [hpq, Polynomial.map_mul, Polynomial.eval_mul, hz'', mul_zero]
        obtain ⟨hqpos, hqlc⟩ :=
          IH q.natDegree hlt q rfl hqlead hq00 (fun z' hz' => hroots z' (hqroots z' hz'))
        refine sector_step hpq hQ0 hq0 hQpos ?_ hqpos hqlc
        rw [hQ, coeff_quad_factor]
        exact (logConcave_quadratic_of_abs_le_pi_div_three rho phi hrhopos
          (by rw [hphi]; exact harg)).2.1

/-- **G2.** A real polynomial with strictly positive coefficients whose complex
roots all lie within `π/3` of the negative real axis has a log-concave
coefficient sequence. -/
theorem logConcave_of_roots_in_sector (p : Polynomial ℝ)
    (hpos : ∀ k ≤ p.natDegree, 0 < p.coeff k)
    (hroots : ∀ z : ℂ, (p.map (algebraMap ℝ ℂ)).IsRoot z →
        |Complex.arg (-z)| ≤ Real.pi / 3) :
    LogConcave (fun k => p.coeff k) :=
  (sector_aux p.natDegree p rfl (hpos p.natDegree le_rfl) (ne_of_gt (hpos 0 (by omega)))
    hroots).2

end LC
