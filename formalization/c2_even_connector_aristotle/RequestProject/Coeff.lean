/-
Source: derived for the even-connector Laurent block request
(`c2_even_connector_aristotle_input_20260820.md`).

Coefficient calculus for Laurent polynomials over `ℚ`:

* `lcoeff f m` is the coefficient of `z ^ m` in `f`;
* `czz n j` is the binomial coefficient `C(n, j)` with an integer lower index
  (zero for `j < 0`);
* `XL = z + z⁻¹` and `HL n = z ^ n + z ^ (-n)`;
* `lcoeff_XL_pow` computes the coefficients of `XL ^ n`.
-/
import Mathlib
import RequestProject.Binomial

namespace EvenConnector

open LaurentPolynomial

/-! ### Coefficients of Laurent polynomials -/

/-- The coefficient of `z ^ m` in the Laurent polynomial `f`. -/
noncomputable def lcoeff (f : ℚ[T;T⁻¹]) (m : ℤ) : ℚ := (f : ℤ →₀ ℚ) m

@[simp] theorem lcoeff_add (f g : ℚ[T;T⁻¹]) (m : ℤ) :
    lcoeff (f + g) m = lcoeff f m + lcoeff g m := by
  simp [lcoeff]

@[simp] theorem lcoeff_T (k m : ℤ) : lcoeff (T k) m = if m = k then 1 else 0 := by
  simp [lcoeff, LaurentPolynomial.T_apply, eq_comm]

@[simp] theorem lcoeff_one (m : ℤ) : lcoeff (1 : ℚ[T;T⁻¹]) m = if m = 0 then 1 else 0 := by
  rw [show (1 : ℚ[T;T⁻¹]) = T 0 from (LaurentPolynomial.T_zero).symm, lcoeff_T]

theorem lcoeff_mul_T (f : ℚ[T;T⁻¹]) (k m : ℤ) : lcoeff (f * T k) m = lcoeff f (m - k) := by
  show ((f * AddMonoidAlgebra.single k (1 : ℚ) : ℚ[T;T⁻¹]) : ℤ →₀ ℚ) m = _
  rw [AddMonoidAlgebra.mul_single_apply]
  simp [lcoeff, sub_eq_add_neg]

theorem lcoeff_C_mul (a : ℚ) (f : ℚ[T;T⁻¹]) (m : ℤ) :
    lcoeff (C a * f) m = a * lcoeff f m := by
  show ((AddMonoidAlgebra.single 0 a * f : ℚ[T;T⁻¹]) : ℤ →₀ ℚ) m = _
  rw [AddMonoidAlgebra.single_mul_apply]
  simp [lcoeff]

/-! ### Binomial coefficients with an integer lower index -/

/-- `C(n, j)` for an integer lower index `j`, as a rational number. -/
noncomputable def czz (n : ℕ) (j : ℤ) : ℚ := if 0 ≤ j then (n.choose j.toNat : ℚ) else 0

theorem czz_nonneg (n : ℕ) (j : ℤ) : 0 ≤ czz n j := by
  unfold czz; split <;> positivity

theorem czz_of_neg (n : ℕ) (j : ℤ) (h : j < 0) : czz n j = 0 := by
  unfold czz; rw [if_neg (by omega)]

theorem czz_natCast (n i : ℕ) : czz n (i : ℤ) = (n.choose i : ℚ) := by
  unfold czz
  rw [if_pos (by positivity)]
  simp

@[simp] theorem czz_zero (n : ℕ) : czz n 0 = 1 := by
  unfold czz; norm_num

theorem czz_of_gt (n : ℕ) (j : ℤ) (h : (n : ℤ) < j) : czz n j = 0 := by
  unfold czz
  rw [if_pos (by omega)]
  have : n < j.toNat := by omega
  simp [Nat.choose_eq_zero_of_lt this]

/-- Symmetry of the binomial row. -/
theorem czz_symm (n : ℕ) (j : ℤ) : czz n ((n : ℤ) - j) = czz n j := by
  rcases lt_or_ge j 0 with h | h
  · rw [czz_of_neg n j h, czz_of_gt n _ (by omega)]
  · obtain ⟨i, rfl⟩ : ∃ i : ℕ, j = (i : ℤ) := ⟨j.toNat, by omega⟩
    rcases le_or_gt i n with hi | hi
    · have e : (n : ℤ) - (i : ℤ) = ((n - i : ℕ) : ℤ) := by omega
      rw [e, czz_natCast, czz_natCast, Nat.choose_symm hi]
    · rw [czz_of_gt n (i : ℤ) (by exact_mod_cast hi), czz_of_neg n _ (by omega)]

/-- Pascal's rule. -/
theorem czz_succ (n : ℕ) (j : ℤ) : czz (n + 1) j = czz n (j - 1) + czz n j := by
  rcases lt_or_ge j 0 with h | h
  · rw [czz_of_neg _ _ h, czz_of_neg _ _ h, czz_of_neg _ _ (by omega)]
    ring
  · rcases eq_or_lt_of_le h with h0 | h0
    · rw [← h0, czz_of_neg n (0 - 1) (by norm_num), czz_zero, czz_zero]
      ring
    · obtain ⟨i, rfl⟩ : ∃ i : ℕ, j = ((i + 1 : ℕ) : ℤ) := ⟨j.toNat - 1, by omega⟩
      have h1 : ((i + 1 : ℕ) : ℤ) - 1 = (i : ℤ) := by push_cast; ring
      rw [h1, czz_natCast, czz_natCast, czz_natCast, Nat.choose_succ_succ]
      push_cast
      ring

/-- Monotonicity past the middle of the row. -/
theorem czz_antitone (n : ℕ) (j : ℤ) (h : (n : ℤ) ≤ 2 * j + 1) : czz n (j + 1) ≤ czz n j := by
  have hj : 0 ≤ j := by
    have : (0 : ℤ) ≤ (n : ℤ) := by positivity
    omega
  obtain ⟨i, rfl⟩ : ∃ i : ℕ, j = (i : ℤ) := ⟨j.toNat, by omega⟩
  have h1 : ((i : ℤ)) + 1 = ((i + 1 : ℕ) : ℤ) := by push_cast; ring
  rw [h1, czz_natCast, czz_natCast]
  have hn : n ≤ 2 * i + 1 := by exact_mod_cast h
  exact_mod_cast choose_succ_le_of_ge n i hn

/-- Monotonicity before the middle of the row. -/
theorem czz_monotone (n : ℕ) (i : ℕ) (h : 2 * i + 1 ≤ n) :
    czz n (i : ℤ) ≤ czz n ((i : ℤ) + 1) := by
  have h1 : ((i : ℤ)) + 1 = ((i + 1 : ℕ) : ℤ) := by push_cast; ring
  rw [h1, czz_natCast, czz_natCast]
  exact_mod_cast choose_le_succ_of_le n i h

/-! ### The Laurent polynomials `X` and `H` -/

/-- `X = z + z⁻¹`. -/
noncomputable def XL : ℚ[T;T⁻¹] := T 1 + T (-1)

/-- `H n = z ^ n + z ^ (-n)`. -/
noncomputable def HL (n : ℤ) : ℚ[T;T⁻¹] := T n + T (-n)

theorem lcoeff_mul_HL (f : ℚ[T;T⁻¹]) (b m : ℤ) :
    lcoeff (f * HL b) m = lcoeff f (m - b) + lcoeff f (m + b) := by
  unfold HL
  rw [mul_add, lcoeff_add, lcoeff_mul_T, lcoeff_mul_T]
  congr 2
  ring

/-- The coefficients of `X ^ n` are the binomial coefficients, at exponents `2t - n`. -/
theorem lcoeff_XL_pow (n : ℕ) (t : ℤ) : lcoeff (XL ^ n) (2 * t - (n : ℤ)) = czz n t := by
  induction n generalizing t with
  | zero =>
    simp only [pow_zero, Nat.cast_zero, sub_zero, lcoeff_one]
    rcases lt_or_ge t 0 with ht | ht
    · rw [if_neg (by omega), czz_of_neg _ _ ht]
    · obtain ⟨i, rfl⟩ : ∃ i : ℕ, t = (i : ℤ) := ⟨t.toNat, by omega⟩
      rw [czz_natCast]
      match i with
      | 0 => norm_num
      | (i + 1) =>
        rw [if_neg (by push_cast; omega)]
        simp
  | succ n ih =>
    have hpow : (XL : ℚ[T;T⁻¹]) ^ (n + 1) = XL ^ n * T 1 + XL ^ n * T (-1) := by
      rw [pow_succ]
      unfold XL
      ring
    rw [hpow, lcoeff_add, lcoeff_mul_T, lcoeff_mul_T]
    have e1 : 2 * t - ((n : ℤ) + 1) - 1 = 2 * (t - 1) - (n : ℤ) := by ring
    have e2 : 2 * t - ((n : ℤ) + 1) - (-1) = 2 * t - (n : ℤ) := by ring
    push_cast
    rw [e1, e2, ih, ih, czz_succ]

/-- Off-parity coefficients of `X ^ n` vanish. -/
theorem lcoeff_XL_pow_odd (n : ℕ) (m : ℤ) (h : ¬ (2 ∣ (m + (n : ℤ)))) :
    lcoeff (XL ^ n) m = 0 := by
  induction n generalizing m with
  | zero =>
    simp only [pow_zero, lcoeff_one]
    have hm : ¬ (m = 0) := by
      rintro rfl
      exact h (by norm_num)
    rw [if_neg hm]
  | succ n ih =>
    have hpow : (XL : ℚ[T;T⁻¹]) ^ (n + 1) = XL ^ n * T 1 + XL ^ n * T (-1) := by
      rw [pow_succ]
      unfold XL
      ring
    rw [hpow, lcoeff_add, lcoeff_mul_T, lcoeff_mul_T]
    push_cast at h
    rw [ih (m - 1) (by omega), ih (m - (-1)) (by omega)]
    ring

end EvenConnector
