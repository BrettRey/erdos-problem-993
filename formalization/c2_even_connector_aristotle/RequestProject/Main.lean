/-
Source: `c2_even_connector_aristotle_input_20260820.md`
  ("Formalization request: the even-connector Laurent block").

Main file.  All statements below are the formal counterparts of the objects and
claims of that request; see `README.md` for the traceability table.

Over the Laurent polynomials `ℚ[T;T⁻¹]` (with `z := T 1`) we put

    X   = z + z⁻¹                      (`XL`)
    H n = z ^ n + z ^ (-n)             (`HL n`)

and, for `r, s ≥ 1` and rationals `c, d ≥ 1`,

    G(r,s,c,d) = c*d*X^(r+s+1) + d*X^s*H(r+1) + c*X^r*H(s+1) + H(r+s+1)   (`Gpoly`)

The main theorem `Gpoly_centrallyUnimodal` states that `G(r,s,c,d)` has centrally
unimodal coefficients, and `Gpoly_coeff_sub_nonneg` records the resulting
coefficient inequality

    coeff (G r s c d) m - coeff (G r s c d) (m+2) ≥ 0     for every `m ≥ 0`.
-/
import Mathlib
import RequestProject.Binomial
import RequestProject.Coeff

-- The option block below is the one supplied with the request.
open scoped BigOperators
open scoped Real
open scoped Nat
open scoped Classical
open scoped Pointwise

set_option maxHeartbeats 8000000
set_option maxRecDepth 4000
set_option synthInstance.maxHeartbeats 20000
set_option synthInstance.maxSize 128

set_option relaxedAutoImplicit false
set_option autoImplicit false

set_option pp.fullNames true
set_option pp.structureInstances true
set_option pp.coercions.types true
set_option pp.funBinderTypes true
set_option pp.letVarTypes true
set_option pp.piBinderTypes true

set_option grind.warning false

namespace EvenConnector

open LaurentPolynomial

/-! ### The even-connector block -/

/-- The even-connector block
`G(r,s,c,d) = c*d*X^(r+s+1) + d*X^s*H(r+1) + c*X^r*H(s+1) + H(r+s+1)`. -/
noncomputable def Gpoly (r s : ℕ) (c d : ℚ) : ℚ[T;T⁻¹] :=
  C c * (C d * XL ^ (r + s + 1))
    + C d * (XL ^ s * HL ((r : ℤ) + 1))
    + C c * (XL ^ r * HL ((s : ℤ) + 1))
    + HL ((r : ℤ) + (s : ℤ) + 1)

/-- The coefficient of `G(r,s,c,d)` at the exponent `2i - (r+s+1)`. -/
noncomputable def coefG (r s : ℕ) (c d : ℚ) (i : ℤ) : ℚ :=
  c * d * czz (r + s + 1) i
    + d * (czz s i + czz s (i - ((r : ℤ) + 1)))
    + c * (czz r i + czz r (i - ((s : ℤ) + 1)))
    + (if i = (r : ℤ) + (s : ℤ) + 1 then 1 else 0) + (if i = 0 then 1 else 0)

/-- The exact identity suggested in the request:
`c*d*A + d*B + c*C + D = (c-1)*(d-1)*A + (c-1)*(A+C) + (d-1)*(A+B) + (A+B+C+D)`,
where `A = X^(r+s+1)`, `B = X^s*H(r+1)`, `C = X^r*H(s+1)` and `D = H(r+s+1)`. -/
theorem Gpoly_decomposition (r s : ℕ) (c d : ℚ) :
    Gpoly r s c d =
      C ((c - 1) * (d - 1)) * XL ^ (r + s + 1)
      + C (c - 1) * (XL ^ (r + s + 1) + XL ^ r * HL ((s : ℤ) + 1))
      + C (d - 1) * (XL ^ (r + s + 1) + XL ^ s * HL ((r : ℤ) + 1))
      + (XL ^ (r + s + 1) + XL ^ s * HL ((r : ℤ) + 1) + XL ^ r * HL ((s : ℤ) + 1)
          + HL ((r : ℤ) + (s : ℤ) + 1)) := by
  unfold Gpoly
  simp only [map_sub, map_mul, map_one]
  ring

/-! ### Coefficients of the block -/

/-- The coefficients of `G(r,s,c,d)` at exponents of the same parity as `r+s+1`. -/
theorem lcoeff_Gpoly (r s : ℕ) (c d : ℚ) (i : ℤ) :
    lcoeff (Gpoly r s c d) (2 * i - ((r : ℤ) + (s : ℤ) + 1)) = coefG r s c d i := by
  have hN : ((r + s + 1 : ℕ) : ℤ) = (r : ℤ) + (s : ℤ) + 1 := by push_cast; ring
  set m : ℤ := 2 * i - ((r : ℤ) + (s : ℤ) + 1) with hm
  have e1 : lcoeff (XL ^ (r + s + 1)) m = czz (r + s + 1) i := by
    have := lcoeff_XL_pow (r + s + 1) i
    rw [hN] at this
    exact this
  have e2 : lcoeff (XL ^ s) (m - ((r : ℤ) + 1)) = czz s (i - ((r : ℤ) + 1)) := by
    have := lcoeff_XL_pow s (i - ((r : ℤ) + 1))
    have h' : 2 * (i - ((r : ℤ) + 1)) - (s : ℤ) = m - ((r : ℤ) + 1) := by rw [hm]; ring
    rw [h'] at this
    exact this
  have e3 : lcoeff (XL ^ s) (m + ((r : ℤ) + 1)) = czz s i := by
    have := lcoeff_XL_pow s i
    have h' : 2 * i - (s : ℤ) = m + ((r : ℤ) + 1) := by rw [hm]; ring
    rw [h'] at this
    exact this
  have e4 : lcoeff (XL ^ r) (m - ((s : ℤ) + 1)) = czz r (i - ((s : ℤ) + 1)) := by
    have := lcoeff_XL_pow r (i - ((s : ℤ) + 1))
    have h' : 2 * (i - ((s : ℤ) + 1)) - (r : ℤ) = m - ((s : ℤ) + 1) := by rw [hm]; ring
    rw [h'] at this
    exact this
  have e5 : lcoeff (XL ^ r) (m + ((s : ℤ) + 1)) = czz r i := by
    have := lcoeff_XL_pow r i
    have h' : 2 * i - (r : ℤ) = m + ((s : ℤ) + 1) := by rw [hm]; ring
    rw [h'] at this
    exact this
  unfold Gpoly coefG
  simp only [lcoeff_add, lcoeff_C_mul, lcoeff_mul_HL]
  rw [e1, e2, e3, e4, e5]
  unfold HL
  rw [lcoeff_add, lcoeff_T, lcoeff_T]
  have hi1 : (m = (r : ℤ) + (s : ℤ) + 1) ↔ (i = (r : ℤ) + (s : ℤ) + 1) := by
    rw [hm]; constructor <;> intro h <;> omega
  have hi2 : (m = -((r : ℤ) + (s : ℤ) + 1)) ↔ (i = 0) := by
    rw [hm]; constructor <;> intro h <;> omega
  rw [if_congr hi1 rfl rfl, if_congr hi2 rfl rfl]
  ring

/-- The coefficients of `G(r,s,c,d)` vanish at exponents of the wrong parity. -/
theorem lcoeff_Gpoly_odd (r s : ℕ) (c d : ℚ) (m : ℤ)
    (h : ¬ (2 ∣ (m + ((r : ℤ) + (s : ℤ) + 1)))) : lcoeff (Gpoly r s c d) m = 0 := by
  have hN : ((r + s + 1 : ℕ) : ℤ) = (r : ℤ) + (s : ℤ) + 1 := by push_cast; ring
  have e1 : lcoeff (XL ^ (r + s + 1)) m = 0 := by
    refine lcoeff_XL_pow_odd _ _ ?_
    rw [hN]
    exact h
  have e2 : lcoeff (XL ^ s) (m - ((r : ℤ) + 1)) = 0 := by
    refine lcoeff_XL_pow_odd _ _ ?_
    intro hc
    exact h (by omega)
  have e3 : lcoeff (XL ^ s) (m + ((r : ℤ) + 1)) = 0 := by
    refine lcoeff_XL_pow_odd _ _ ?_
    intro hc
    exact h (by omega)
  have e4 : lcoeff (XL ^ r) (m - ((s : ℤ) + 1)) = 0 := by
    refine lcoeff_XL_pow_odd _ _ ?_
    intro hc
    exact h (by omega)
  have e5 : lcoeff (XL ^ r) (m + ((s : ℤ) + 1)) = 0 := by
    refine lcoeff_XL_pow_odd _ _ ?_
    intro hc
    exact h (by omega)
  unfold Gpoly
  simp only [lcoeff_add, lcoeff_C_mul, lcoeff_mul_HL]
  rw [e1, e2, e3, e4, e5]
  unfold HL
  rw [lcoeff_add, lcoeff_T, lcoeff_T]
  rw [if_neg (by intro hc; exact h (by omega)), if_neg (by intro hc; exact h (by omega))]
  ring

/-! ### Properties of the coefficient function -/

theorem coefG_nonneg (r s : ℕ) (c d : ℚ) (hc : 0 ≤ c) (hd : 0 ≤ d) (i : ℤ) :
    0 ≤ coefG r s c d i := by
  have h1 := czz_nonneg (r + s + 1) i
  have h2 := czz_nonneg s i
  have h3 := czz_nonneg s (i - ((r : ℤ) + 1))
  have h4 := czz_nonneg r i
  have h5 := czz_nonneg r (i - ((s : ℤ) + 1))
  unfold coefG
  have i1 : (0 : ℚ) ≤ if i = (r : ℤ) + (s : ℤ) + 1 then 1 else 0 := by split <;> norm_num
  have i2 : (0 : ℚ) ≤ if i = 0 then (1 : ℚ) else 0 := by split <;> norm_num
  have : 0 ≤ c * d := mul_nonneg hc hd
  nlinarith [mul_nonneg hd (add_nonneg h2 h3), mul_nonneg hc (add_nonneg h4 h5),
    mul_nonneg this h1]

theorem coefG_eq_zero_of_gt (r s : ℕ) (c d : ℚ) (i : ℤ) (h : (r : ℤ) + (s : ℤ) + 1 < i) :
    coefG r s c d i = 0 := by
  have hN : ((r + s + 1 : ℕ) : ℤ) = (r : ℤ) + (s : ℤ) + 1 := by push_cast; ring
  unfold coefG
  rw [czz_of_gt (r + s + 1) i (by omega), czz_of_gt s i (by omega),
    czz_of_gt s _ (by omega), czz_of_gt r i (by omega), czz_of_gt r _ (by omega),
    if_neg (by omega), if_neg (by omega)]
  ring

/-- The coefficient function is symmetric about the centre `(r+s+1)/2`. -/
theorem coefG_symm (r s : ℕ) (c d : ℚ) (i : ℤ) :
    coefG r s c d ((r : ℤ) + (s : ℤ) + 1 - i) = coefG r s c d i := by
  have hN : ((r + s + 1 : ℕ) : ℤ) = (r : ℤ) + (s : ℤ) + 1 := by push_cast; ring
  have a1 : czz (r + s + 1) ((r : ℤ) + (s : ℤ) + 1 - i) = czz (r + s + 1) i := by
    rw [← hN]; exact czz_symm _ _
  have a2 : czz s ((r : ℤ) + (s : ℤ) + 1 - i) = czz s (i - ((r : ℤ) + 1)) := by
    have := czz_symm s (i - ((r : ℤ) + 1))
    rw [show (s : ℤ) - (i - ((r : ℤ) + 1)) = (r : ℤ) + (s : ℤ) + 1 - i by ring] at this
    exact this
  have a3 : czz s ((r : ℤ) + (s : ℤ) + 1 - i - ((r : ℤ) + 1)) = czz s i := by
    have := czz_symm s i
    rw [show (s : ℤ) - i = (r : ℤ) + (s : ℤ) + 1 - i - ((r : ℤ) + 1) by ring] at this
    exact this
  have a4 : czz r ((r : ℤ) + (s : ℤ) + 1 - i) = czz r (i - ((s : ℤ) + 1)) := by
    have := czz_symm r (i - ((s : ℤ) + 1))
    rw [show (r : ℤ) - (i - ((s : ℤ) + 1)) = (r : ℤ) + (s : ℤ) + 1 - i by ring] at this
    exact this
  have a5 : czz r ((r : ℤ) + (s : ℤ) + 1 - i - ((s : ℤ) + 1)) = czz r i := by
    have := czz_symm r i
    rw [show (r : ℤ) - i = (r : ℤ) + (s : ℤ) + 1 - i - ((s : ℤ) + 1) by ring] at this
    exact this
  unfold coefG
  rw [a1, a2, a3, a4, a5]
  have b1 : ((r : ℤ) + (s : ℤ) + 1 - i = (r : ℤ) + (s : ℤ) + 1) ↔ (i = 0) := by omega
  have b2 : ((r : ℤ) + (s : ℤ) + 1 - i = 0) ↔ (i = (r : ℤ) + (s : ℤ) + 1) := by omega
  rw [if_congr b1 rfl rfl, if_congr b2 rfl rfl]
  ring

/-- The central unimodality step for the coefficient function. -/
theorem coefG_step (r s : ℕ) (c d : ℚ) (hr : 1 ≤ r) (hs : 1 ≤ s) (hc : 1 ≤ c) (hd : 1 ≤ d)
    (i : ℤ) (hi : (r : ℤ) + (s : ℤ) + 1 ≤ 2 * i) :
    coefG r s c d (i + 1) ≤ coefG r s c d i := by
  have hc0 : (0 : ℚ) ≤ c := by linarith
  have hd0 : (0 : ℚ) ≤ d := by linarith
  rcases lt_or_ge ((r : ℤ) + (s : ℤ)) i with hbig | hsmall
  · rw [coefG_eq_zero_of_gt r s c d (i + 1) (by omega)]
    exact coefG_nonneg r s c d hc0 hd0 i
  · -- the interesting range : `i = r + s - k` with `2k + 2 ≤ r + s + 1`
    have hr1 : (1 : ℤ) ≤ (r : ℤ) := by exact_mod_cast hr
    have hs1 : (1 : ℤ) ≤ (s : ℤ) := by exact_mod_cast hs
    obtain ⟨k, hk⟩ : ∃ k : ℕ, i = ((r : ℤ) + (s : ℤ)) - (k : ℤ) :=
      ⟨((r : ℤ) + (s : ℤ) - i).toNat, by omega⟩
    have hkZ : (k : ℤ) ≤ (r : ℤ) + (s : ℤ) := by omega
    have hk2 : 2 * k + 2 ≤ r + s + 1 := by
      have : 2 * (k : ℤ) + 2 ≤ (r : ℤ) + (s : ℤ) + 1 := by omega
      exact_mod_cast this
    -- rewrite the six binomial values in terms of `k`
    have r1 : czz (r + s + 1) i = czz (r + s + 1) ((k : ℤ) + 1) := by
      have h := czz_symm (r + s + 1) ((k : ℤ) + 1)
      rw [show ((r + s + 1 : ℕ) : ℤ) - ((k : ℤ) + 1) = i by push_cast; omega] at h
      exact h
    have r2 : czz (r + s + 1) (i + 1) = czz (r + s + 1) (k : ℤ) := by
      have h := czz_symm (r + s + 1) (k : ℤ)
      rw [show ((r + s + 1 : ℕ) : ℤ) - (k : ℤ) = i + 1 by push_cast; omega] at h
      exact h
    have r3 : czz s (i - ((r : ℤ) + 1)) = czz s ((k : ℤ) + 1) := by
      have h := czz_symm s ((k : ℤ) + 1)
      rw [show (s : ℤ) - ((k : ℤ) + 1) = i - ((r : ℤ) + 1) by omega] at h
      exact h
    have r4 : czz s (i + 1 - ((r : ℤ) + 1)) = czz s (k : ℤ) := by
      have h := czz_symm s (k : ℤ)
      rw [show (s : ℤ) - (k : ℤ) = i + 1 - ((r : ℤ) + 1) by omega] at h
      exact h
    have r5 : czz r (i - ((s : ℤ) + 1)) = czz r ((k : ℤ) + 1) := by
      have h := czz_symm r ((k : ℤ) + 1)
      rw [show (r : ℤ) - ((k : ℤ) + 1) = i - ((s : ℤ) + 1) by omega] at h
      exact h
    have r6 : czz r (i + 1 - ((s : ℤ) + 1)) = czz r (k : ℤ) := by
      have h := czz_symm r (k : ℤ)
      rw [show (r : ℤ) - (k : ℤ) = i + 1 - ((s : ℤ) + 1) by omega] at h
      exact h
    -- the two "far" monotonicity facts
    have mA : czz s (i + 1) ≤ czz s i := czz_antitone s i (by omega)
    have mB : czz r (i + 1) ≤ czz r i := czz_antitone r i (by omega)
    -- the three master inequalities
    have P0 : czz (r + s + 1) (k : ℤ) ≤ czz (r + s + 1) ((k : ℤ) + 1) :=
      czz_monotone (r + s + 1) k (by omega)
    have PU : czz (r + s + 1) (k : ℤ) + czz r (k : ℤ)
        ≤ czz (r + s + 1) ((k : ℤ) + 1) + czz r ((k : ℤ) + 1) := by
      have h := key_two r s k hr hs hk2
      have e : ((k : ℤ) + 1) = ((k + 1 : ℕ) : ℤ) := by push_cast; ring
      rw [e, czz_natCast, czz_natCast, czz_natCast, czz_natCast]
      exact_mod_cast h
    have PV : czz (r + s + 1) (k : ℤ) + czz s (k : ℤ)
        ≤ czz (r + s + 1) ((k : ℤ) + 1) + czz s ((k : ℤ) + 1) := by
      have h := key_two s r k hs hr (by omega)
      rw [show s + r + 1 = r + s + 1 by omega] at h
      have e : ((k : ℤ) + 1) = ((k + 1 : ℕ) : ℤ) := by push_cast; ring
      rw [e, czz_natCast, czz_natCast, czz_natCast, czz_natCast]
      exact_mod_cast h
    have PUV : czz (r + s + 1) (k : ℤ) + czz r (k : ℤ) + czz s (k : ℤ)
          + (if k = 0 then 1 else 0)
        ≤ czz (r + s + 1) ((k : ℤ) + 1) + czz r ((k : ℤ) + 1) + czz s ((k : ℤ) + 1) := by
      have h := key_three r s k hr hs hk2
      have e : ((k : ℤ) + 1) = ((k + 1 : ℕ) : ℤ) := by push_cast; ring
      rw [e, czz_natCast, czz_natCast, czz_natCast, czz_natCast, czz_natCast, czz_natCast]
      exact_mod_cast h
    -- indicator values
    have j1 : ¬ (i = (r : ℤ) + (s : ℤ) + 1) := by omega
    have j2 : ¬ (i = 0) := by omega
    have j3 : ¬ (i + 1 = 0) := by omega
    have j4 : (i + 1 = (r : ℤ) + (s : ℤ) + 1) ↔ (k = 0) := by
      constructor
      · intro h; omega
      · intro h; simp [h] at hk; omega
    unfold coefG
    rw [r1, r2, r3, r4, r5, r6, if_neg j1, if_neg j2, if_neg j3, if_congr j4 rfl rfl]
    have t1 : 0 ≤ (c - 1) * (d - 1) *
        (czz (r + s + 1) ((k : ℤ) + 1) - czz (r + s + 1) (k : ℤ)) := by
      apply mul_nonneg (mul_nonneg (by linarith) (by linarith))
      linarith
    have t2 : 0 ≤ (c - 1) *
        ((czz (r + s + 1) ((k : ℤ) + 1) + czz r ((k : ℤ) + 1))
          - (czz (r + s + 1) (k : ℤ) + czz r (k : ℤ))) := by
      apply mul_nonneg (by linarith)
      linarith
    have t3 : 0 ≤ (d - 1) *
        ((czz (r + s + 1) ((k : ℤ) + 1) + czz s ((k : ℤ) + 1))
          - (czz (r + s + 1) (k : ℤ) + czz s (k : ℤ))) := by
      apply mul_nonneg (by linarith)
      linarith
    have t4 : 0 ≤ d * (czz s i - czz s (i + 1)) := mul_nonneg hd0 (by linarith)
    have t5 : 0 ≤ c * (czz r i - czz r (i + 1)) := mul_nonneg hc0 (by linarith)
    nlinarith [t1, t2, t3, t4, t5, PUV]

/-! ### Central unimodality -/

/-- A Laurent polynomial has *centrally unimodal coefficients* when it is symmetric
under `m ↦ -m`, has nonnegative coefficients, and its coefficients weakly decrease
as `m ≥ 0` increases within each parity class. -/
structure CentrallyUnimodal (f : ℚ[T;T⁻¹]) : Prop where
  symm : ∀ m : ℤ, lcoeff f (-m) = lcoeff f m
  nonneg : ∀ m : ℤ, 0 ≤ lcoeff f m
  decr : ∀ m : ℤ, 0 ≤ m → lcoeff f (m + 2) ≤ lcoeff f m

/-- **Main theorem.**  For all natural numbers `r, s ≥ 1` and rationals `c, d ≥ 1`,
the even-connector block `G(r,s,c,d)` has centrally unimodal coefficients. -/
theorem Gpoly_centrallyUnimodal (r s : ℕ) (c d : ℚ)
    (hr : 1 ≤ r) (hs : 1 ≤ s) (hc : 1 ≤ c) (hd : 1 ≤ d) :
    CentrallyUnimodal (Gpoly r s c d) := by
  have hc0 : (0 : ℚ) ≤ c := by linarith
  have hd0 : (0 : ℚ) ≤ d := by linarith
  refine ⟨?_, ?_, ?_⟩
  · intro m
    by_cases hpar : 2 ∣ (m + ((r : ℤ) + (s : ℤ) + 1))
    · obtain ⟨i, hi⟩ : ∃ i : ℤ, m = 2 * i - ((r : ℤ) + (s : ℤ) + 1) :=
        ⟨(m + ((r : ℤ) + (s : ℤ) + 1)) / 2, by omega⟩
      subst hi
      have e : -(2 * i - ((r : ℤ) + (s : ℤ) + 1))
          = 2 * ((r : ℤ) + (s : ℤ) + 1 - i) - ((r : ℤ) + (s : ℤ) + 1) := by ring
      rw [e, lcoeff_Gpoly, lcoeff_Gpoly, coefG_symm]
    · rw [lcoeff_Gpoly_odd r s c d m hpar, lcoeff_Gpoly_odd r s c d (-m) (by omega)]
  · intro m
    by_cases hpar : 2 ∣ (m + ((r : ℤ) + (s : ℤ) + 1))
    · obtain ⟨i, hi⟩ : ∃ i : ℤ, m = 2 * i - ((r : ℤ) + (s : ℤ) + 1) :=
        ⟨(m + ((r : ℤ) + (s : ℤ) + 1)) / 2, by omega⟩
      subst hi
      rw [lcoeff_Gpoly]
      exact coefG_nonneg r s c d hc0 hd0 i
    · rw [lcoeff_Gpoly_odd r s c d m hpar]
  · intro m hm
    by_cases hpar : 2 ∣ (m + ((r : ℤ) + (s : ℤ) + 1))
    · obtain ⟨i, hi⟩ : ∃ i : ℤ, m = 2 * i - ((r : ℤ) + (s : ℤ) + 1) :=
        ⟨(m + ((r : ℤ) + (s : ℤ) + 1)) / 2, by omega⟩
      subst hi
      have e : 2 * i - ((r : ℤ) + (s : ℤ) + 1) + 2
          = 2 * (i + 1) - ((r : ℤ) + (s : ℤ) + 1) := by ring
      rw [e, lcoeff_Gpoly, lcoeff_Gpoly]
      exact coefG_step r s c d hr hs hc hd i (by omega)
    · rw [lcoeff_Gpoly_odd r s c d m hpar,
        lcoeff_Gpoly_odd r s c d (m + 2) (by omega)]

/-- **Consequence.**  Every central difference of coefficients of `G(r,s,c,d)` is
nonnegative for `m ≥ 0`. -/
theorem Gpoly_coeff_sub_nonneg (r s : ℕ) (c d : ℚ)
    (hr : 1 ≤ r) (hs : 1 ≤ s) (hc : 1 ≤ c) (hd : 1 ≤ d) (m : ℤ) (hm : 0 ≤ m) :
    0 ≤ lcoeff (Gpoly r s c d) m - lcoeff (Gpoly r s c d) (m + 2) :=
  sub_nonneg.mpr ((Gpoly_centrallyUnimodal r s c d hr hs hc hd).decr m hm)

end EvenConnector
