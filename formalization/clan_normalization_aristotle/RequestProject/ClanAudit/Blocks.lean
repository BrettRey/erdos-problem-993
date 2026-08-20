import RequestProject.ClanAudit.Laurent

/-!
# The proposed imbalance blocks `A(r,c;z)` and `F(z)`

This file studies the Laurent polynomials proposed in the adjacent two-hub
target of the request:

```text
A(r,c;z) = c*(z + z^(-1))^r + z^r + z^(-r),   c > 0
F(z)     = A(r,c;z) * A(s,d;z) - z^(r+s) - z^(-(r+s)).
```

We prove the four-map identity, we *refute* the claim that the coefficients of
`F` are centrally unimodal for arbitrary `c, d > 0` (with an explicit small
witness), and we prove that they *are* centrally unimodal as soon as
`1 ≤ c` and `1 ≤ d`.
-/

namespace ClanAudit

open LaurentPolynomial Finset

/-- `z + z⁻¹`, the imbalance weight of a single unbalanced bipartite component. -/
noncomputable def zz : LaurentPolynomial ℚ := box 1

theorem zz_eq : zz = T 1 + T (-1) := by
  rw [zz, box, Finset.sum_range_succ, Finset.sum_range_succ, Finset.sum_range_zero]
  norm_num

/-- `z^t + z^(-t)`, the orientation weight of a bipartite component with imbalance `t`. -/
noncomputable def Pw (t : ℕ) : LaurentPolynomial ℚ := T (t : ℤ) + T (-(t : ℤ))

theorem coeffL_Pw (t : ℕ) (m : ℤ) :
    coeffL (Pw t) m = (if m = (t : ℤ) then 1 else 0) + (if m = -(t : ℤ) then 1 else 0) := by
  rw [Pw, coeffL_add, coeffL_T, coeffL_T]

theorem coeffL_mul_Pw (f : LaurentPolynomial ℚ) (t : ℕ) (m : ℤ) :
    coeffL (f * Pw t) m = coeffL f (m - t) + coeffL f (m + t) := by
  rw [Pw, mul_add, coeffL_add, coeffL_mul_T, coeffL_mul_T, sub_neg_eq_add]

/-- On the nonnegative side, `Pw t` has a single nonzero coefficient (for `t ≥ 1`). -/
theorem coeffL_Pw_of_nonneg {t : ℕ} (ht : 1 ≤ t) {x : ℤ} (hx : 0 ≤ x) :
    coeffL (Pw t) x = if x = (t : ℤ) then 1 else 0 := by
  rw [coeffL_Pw, if_neg (show ¬(x = -(t : ℤ)) by omega), add_zero]

/-! ### Coefficients of `(z + z⁻¹)^n` -/

/-- The (computable) coefficient sequence of `(z + z⁻¹)^n`. -/
def cbin : ℕ → ℤ → ℚ
  | 0, m => if m = 0 then 1 else 0
  | (n + 1), m => cbin n (m - 1) + cbin n (m + 1)

theorem coeffL_zz_pow_succ (n : ℕ) (m : ℤ) :
    coeffL (zz ^ (n + 1)) m = coeffL (zz ^ n) (m - 1) + coeffL (zz ^ n) (m + 1) := by
  rw [pow_succ, zz, coeffL_mul_box, Finset.sum_range_succ, Finset.sum_range_succ,
    Finset.sum_range_zero, zero_add]
  congr 2 <;> push_cast <;> ring

theorem coeffL_zz_pow (n : ℕ) (m : ℤ) : coeffL (zz ^ n) m = cbin n m := by
  induction n generalizing m with
  | zero => simp [cbin]
  | succ n ih => rw [coeffL_zz_pow_succ, ih, ih, cbin]

theorem isCU_zz_pow (n : ℕ) : IsCU (zz ^ n) := by
  induction n with
  | zero =>
      constructor
      · intro m; simp only [pow_zero, coeffL_one]; congr 1; simp only [eq_iff_iff]; omega
      · intro m; simp only [pow_zero, coeffL_one]; split <;> norm_num
      · intro m hm
        simp only [pow_zero, coeffL_one]
        rw [if_neg (show ¬(m + 2 = 0) by omega)]
        split <;> norm_num

  | succ n ih => rw [pow_succ, zz]; exact ih.mul_box 1

theorem coeffL_zz_pow_eq_zero (n : ℕ) (m : ℤ) (h : n < m.natAbs) : coeffL (zz ^ n) m = 0 := by
  induction n generalizing m with
  | zero => rw [pow_zero, coeffL_one, if_neg (by omega)]
  | succ n ih =>
      rw [coeffL_zz_pow_succ, ih (m - 1) (by omega), ih (m + 1) (by omega)]
      ring

theorem coeffL_zz_pow_top (n : ℕ) : coeffL (zz ^ n) (n : ℤ) = 1 := by
  induction n with
  | zero => rw [pow_zero, coeffL_one]; norm_num
  | succ n ih =>
      have h1 : ((n + 1 : ℕ) : ℤ) - 1 = (n : ℤ) := by push_cast; ring
      have h2 : coeffL (zz ^ n) (((n + 1 : ℕ) : ℤ) + 1) = 0 :=
        coeffL_zz_pow_eq_zero n _ (by omega)
      rw [coeffL_zz_pow_succ, h1, h2, ih]
      ring

theorem coeffL_zz_pow_sub_two (n : ℕ) (h : 1 ≤ n) : coeffL (zz ^ n) ((n : ℤ) - 2) = n := by
  induction n with
  | zero => omega
  | succ n ih =>
      rcases Nat.lt_or_ge n 1 with h1 | h1
      · have hn : n = 0 := by omega
        subst hn
        rw [coeffL_zz_pow]
        norm_num [cbin]
      · have e1 : ((n + 1 : ℕ) : ℤ) - 2 - 1 = (n : ℤ) - 2 := by push_cast; ring
        have e2 : ((n + 1 : ℕ) : ℤ) - 2 + 1 = (n : ℤ) := by push_cast; ring
        rw [coeffL_zz_pow_succ, e1, e2, ih h1, coeffL_zz_pow_top]
        push_cast; ring

/-- `(z+z⁻¹)^s + z^s + z^(-s)` is centrally unimodal for `s ≥ 1`. -/
theorem isCU_zz_pow_add_Pw (s : ℕ) (hs : 1 ≤ s) : IsCU (zz ^ s + Pw s) := by
  have hzz := isCU_zz_pow s
  constructor
  · intro m
    rw [coeffL_add, coeffL_add, hzz.symm, coeffL_Pw, coeffL_Pw]
    have h1 : (-m = (s:ℤ)) ↔ (m = -(s:ℤ)) := by omega
    have h2 : (-m = -(s:ℤ)) ↔ (m = (s:ℤ)) := by omega
    simp only [h1, h2]
    ring
  · intro m
    rw [coeffL_add, coeffL_Pw]
    have h0 := hzz.nonneg m
    have h1 : (0:ℚ) ≤ if m = (s:ℤ) then 1 else 0 := by split <;> norm_num
    have h2 : (0:ℚ) ≤ if m = -(s:ℤ) then 1 else 0 := by split <;> norm_num
    linarith
  · intro m hm
    have hd := hzz.decr m hm
    rw [coeffL_add, coeffL_add, coeffL_Pw_of_nonneg hs hm,
      coeffL_Pw_of_nonneg hs (show (0:ℤ) ≤ m + 2 by omega)]
    by_cases h1 : m = (s : ℤ) - 2
    · rw [if_neg (show ¬(m = (s:ℤ)) by omega), if_pos (show m + 2 = (s:ℤ) by omega)]
      have hs2 : 2 ≤ s := by omega
      have hv1 : coeffL (zz ^ s) m = (s : ℚ) := by rw [h1]; exact coeffL_zz_pow_sub_two s hs
      have hv2 : coeffL (zz ^ s) (m + 2) = 1 := by
        rw [show m + 2 = ((s : ℕ) : ℤ) by omega]; exact coeffL_zz_pow_top s
      have hq : (2:ℚ) ≤ (s:ℚ) := by exact_mod_cast hs2
      rw [hv1, hv2]
      linarith
    · by_cases h2 : m = (s : ℤ)
      · rw [if_pos h2, if_neg (show ¬(m + 2 = (s:ℤ)) by omega)]
        linarith
      · rw [if_neg h2, if_neg (show ¬(m + 2 = (s:ℤ)) by omega)]
        linarith

/-! ### The blocks -/

/-- The proposed one-hub normalized imbalance block `A(r,c;z) = c (z+z⁻¹)^r + z^r + z^(-r)`. -/
noncomputable def Ablock (r : ℕ) (c : ℚ) : LaurentPolynomial ℚ := c • zz ^ r + Pw r

/-- The proposed two-hub four-map block
`F(z) = A(r,c;z) * A(s,d;z) - z^(r+s) - z^(-(r+s))`. -/
noncomputable def Fblock (r s : ℕ) (c d : ℚ) : LaurentPolynomial ℚ :=
  Ablock r c * Ablock s d - Pw (r + s)

theorem Pw_mul_Pw {r s : ℕ} (h : s ≤ r) : Pw r * Pw s = Pw (r + s) + Pw (r - s) := by
  have hcast : (((r - s : ℕ)) : ℤ) = (r : ℤ) - s := by omega
  have hcast2 : (((r + s : ℕ)) : ℤ) = (r : ℤ) + s := by push_cast; ring
  simp only [Pw, hcast, hcast2]
  rw [add_mul, mul_add, mul_add, ← T_add, ← T_add, ← T_add, ← T_add,
    show (-(r:ℤ) + s) = -((r:ℤ) - s) by ring, show ((r:ℤ) + -s) = (r:ℤ) - s by ring,
    show (-(r:ℤ) + -s) = -((r:ℤ) + s) by ring]
  ring

/-- **The four-map identity.**  For `s ≤ r`,
`F(z) = c d (z+z⁻¹)^(r+s) + c (z+z⁻¹)^r (z^s+z^(-s)) + d (z+z⁻¹)^s (z^r+z^(-r))
        + z^(r-s) + z^(-(r-s))`.
The subtraction of `z^(r+s) + z^(-(r+s))` in `F` exactly converts the product
`(z^r+z^(-r))(z^s+z^(-s))` into the single connected imbalance term `|r-s|`. -/
theorem Fblock_expand {r s : ℕ} (h : s ≤ r) (c d : ℚ) :
    Fblock r s c d
      = (c * d) • zz ^ (r + s) + c • (zz ^ r * Pw s) + d • (zz ^ s * Pw r) + Pw (r - s) := by
  rw [Fblock, Ablock, Ablock]
  simp only [smul_eq_C_mul, map_mul]
  rw [show (C c * zz ^ r + Pw r) * (C d * zz ^ s + Pw s)
      = C c * C d * (zz ^ r * zz ^ s) + C c * (zz ^ r * Pw s) + C d * (zz ^ s * Pw r)
        + Pw r * Pw s from by ring, Pw_mul_Pw h, ← pow_add]
  ring

theorem coeffL_Fblock {r s : ℕ} (h : s ≤ r) (c d : ℚ) (m : ℤ) :
    coeffL (Fblock r s c d) m
      = c * d * cbin (r + s) m
        + c * (cbin r (m - s) + cbin r (m + s))
        + d * (cbin s (m - r) + cbin s (m + r))
        + ((if m = ((r - s : ℕ) : ℤ) then 1 else 0)
            + (if m = -((r - s : ℕ) : ℤ) then 1 else 0)) := by
  rw [Fblock_expand h]
  simp only [coeffL_add, coeffL_smul, coeffL_mul_Pw, coeffL_Pw, coeffL_zz_pow]

/-! ### Refutation of central unimodality for general `c, d > 0` -/

/-- With `r = 3`, `s = 1`, `c = 1/4`, `d = 1` the coefficient of `z^0` in `F` is `3`
while the coefficient of `z^2` is `4`. -/
theorem Fblock_coeffs_of_counterexample :
    coeffL (Fblock 3 1 (1/4) 1) 0 = 3 ∧ coeffL (Fblock 3 1 (1/4) 1) 2 = 4 := by
  constructor
  · rw [coeffL_Fblock (by norm_num)]
    norm_num [cbin]
  · rw [coeffL_Fblock (by norm_num)]
    norm_num [cbin]

/-- **Refutation.**  The claim that the Laurent coefficients of the four-map block
`F` are centrally unimodal for *arbitrary* `c, d > 0` is false: for `r = 3`, `s = 1`,
`c = 1/4`, `d = 1` one has `coeff F 0 = 3 < 4 = coeff F 2`, so the normalized two-row
coefficient `coeff F 0 - coeff F 2 = -1` is negative. -/
theorem Fblock_not_centrally_unimodal :
    ∃ (r s : ℕ) (c d : ℚ) (m : ℤ), 1 ≤ r ∧ 1 ≤ s ∧ 0 < c ∧ 0 < d ∧ 0 ≤ m ∧
      coeffL (Fblock r s c d) m < coeffL (Fblock r s c d) (m + 2) := by
  refine ⟨3, 1, 1/4, 1, 0, by norm_num, by norm_num, by norm_num, by norm_num, le_refl 0, ?_⟩
  obtain ⟨h0, h2⟩ := Fblock_coeffs_of_counterexample
  rw [h0, show (0:ℤ) + 2 = 2 by ring, h2]
  norm_num

/-- The same failure already happens for a single hub: `A(3, 1/4; z)` has
`coeff 1 = 3/4 < 5/4 = coeff 3`. -/
theorem Ablock_not_centrally_unimodal :
    coeffL (Ablock 3 (1/4)) 1 = 3/4 ∧ coeffL (Ablock 3 (1/4)) 3 = 5/4 := by
  constructor <;>
  · rw [Ablock]
    simp only [coeffL_add, coeffL_smul, coeffL_zz_pow, coeffL_Pw]
    norm_num [cbin]

/-! ### Central unimodality when `1 ≤ c` and `1 ≤ d` -/

private theorem Q_expand (r s : ℕ) : (zz ^ r + Pw r) * (zz ^ s + Pw s)
    = zz ^ (r + s) + zz ^ r * Pw s + zz ^ s * Pw r + Pw r * Pw s := by
  rw [pow_add]; ring

private theorem coeffL_Q_top {r s : ℕ} (hr : 1 ≤ r) (hs : 1 ≤ s) :
    coeffL ((zz ^ r + Pw r) * (zz ^ s + Pw s)) ((r : ℤ) + s) = 4 := by
  rw [Q_expand]
  simp only [coeffL_add, coeffL_mul_Pw, coeffL_Pw]
  have e1 : coeffL (zz ^ (r + s)) ((r:ℤ) + s) = 1 := by
    rw [show (r:ℤ) + s = ((r + s : ℕ) : ℤ) by push_cast; ring, coeffL_zz_pow_top]
  have e2 : coeffL (zz ^ r) ((r:ℤ) + s - s) = 1 := by
    rw [show (r:ℤ) + s - s = (r:ℤ) by ring, coeffL_zz_pow_top]
  have e3 : coeffL (zz ^ r) ((r:ℤ) + s + s) = 0 := coeffL_zz_pow_eq_zero r _ (by omega)
  have e4 : coeffL (zz ^ s) ((r:ℤ) + s - r) = 1 := by
    rw [show (r:ℤ) + s - r = (s:ℤ) by ring, coeffL_zz_pow_top]
  have e5 : coeffL (zz ^ s) ((r:ℤ) + s + r) = 0 := coeffL_zz_pow_eq_zero s _ (by omega)
  rw [e1, e2, e3, e4, e5,
    if_pos (show (r:ℤ) + s - s = (r:ℤ) by ring),
    if_neg (show ¬((r:ℤ) + s - s = -(r:ℤ)) by omega),
    if_neg (show ¬((r:ℤ) + s + s = (r:ℤ)) by omega),
    if_neg (show ¬((r:ℤ) + s + s = -(r:ℤ)) by omega)]
  norm_num

private theorem coeffL_Q_above {r s : ℕ} (hr : 1 ≤ r) (hs : 1 ≤ s) (m : ℤ)
    (hm : (r : ℤ) + s < m) : coeffL ((zz ^ r + Pw r) * (zz ^ s + Pw s)) m = 0 := by
  rw [Q_expand]
  simp only [coeffL_add, coeffL_mul_Pw, coeffL_Pw]
  have e1 : coeffL (zz ^ (r + s)) m = 0 := coeffL_zz_pow_eq_zero _ _ (by omega)
  have e2 : coeffL (zz ^ r) (m - s) = 0 := coeffL_zz_pow_eq_zero r _ (by omega)
  have e3 : coeffL (zz ^ r) (m + s) = 0 := coeffL_zz_pow_eq_zero r _ (by omega)
  have e4 : coeffL (zz ^ s) (m - r) = 0 := coeffL_zz_pow_eq_zero s _ (by omega)
  have e5 : coeffL (zz ^ s) (m + r) = 0 := coeffL_zz_pow_eq_zero s _ (by omega)
  rw [e1, e2, e3, e4, e5,
    if_neg (show ¬(m - s = (r:ℤ)) by omega),
    if_neg (show ¬(m - s = -(r:ℤ)) by omega),
    if_neg (show ¬(m + s = (r:ℤ)) by omega),
    if_neg (show ¬(m + s = -(r:ℤ)) by omega)]
  norm_num

/-- The base block `F` with `c = d = 1` has nonincreasing coefficients in steps of two
on `m ≥ 0`. -/
private theorem base_decr {r s : ℕ} (hr : 1 ≤ r) (hs : 1 ≤ s) (m : ℤ) (hm : 0 ≤ m) :
    coeffL ((zz ^ r + Pw r) * (zz ^ s + Pw s) - Pw (r + s)) (m + 2)
      ≤ coeffL ((zz ^ r + Pw r) * (zz ^ s + Pw s) - Pw (r + s)) m := by
  set Q : LaurentPolynomial ℚ := (zz ^ r + Pw r) * (zz ^ s + Pw s) with hQ
  have hQCU : IsCU Q := (isCU_zz_pow_add_Pw r hr).mul (isCU_zz_pow_add_Pw s hs)
  have hdec := hQCU.decr m hm
  simp only [coeffL_sub]
  rw [coeffL_Pw_of_nonneg (show 1 ≤ r + s by omega) hm,
    coeffL_Pw_of_nonneg (show 1 ≤ r + s by omega) (show (0:ℤ) ≤ m + 2 by omega)]
  by_cases h1 : m = ((r + s : ℕ) : ℤ)
  · rw [if_pos h1, if_neg (show ¬(m + 2 = ((r + s : ℕ) : ℤ)) by omega)]
    have h4 : coeffL Q m = 4 := by
      rw [h1, show ((r + s : ℕ) : ℤ) = (r:ℤ) + s by push_cast; ring]
      exact coeffL_Q_top hr hs
    have h0 : coeffL Q (m + 2) = 0 := coeffL_Q_above hr hs _ (by omega)
    rw [h4, h0]
    norm_num
  · rw [if_neg h1]
    by_cases h2 : m + 2 = ((r + s : ℕ) : ℤ)
    · rw [if_pos h2]; linarith
    · rw [if_neg h2]; linarith

/-- **Central unimodality of the four-map block for `1 ≤ c, d`.**  Every normalized
two-row coefficient of the block, `coeff F m - coeff F (m+2)` for `m ≥ 0`, is then
nonnegative. -/
theorem Fblock_decr {r s : ℕ} (hr : 1 ≤ r) (hs : 1 ≤ s) {c d : ℚ} (hc : 1 ≤ c) (hd : 1 ≤ d)
    (m : ℤ) (hm : 0 ≤ m) : coeffL (Fblock r s c d) (m + 2) ≤ coeffL (Fblock r s c d) m := by
  have hsplit : Fblock r s c d
      = ((zz ^ r + Pw r) * (zz ^ s + Pw s) - Pw (r + s))
        + ((c - 1) * (d - 1)) • zz ^ (r + s)
        + (c - 1) • (zz ^ r * (zz ^ s + Pw s))
        + (d - 1) • (zz ^ s * (zz ^ r + Pw r)) := by
    rw [Fblock, Ablock, Ablock, pow_add]
    simp only [smul_eq_C_mul, map_mul, map_sub, map_one]
    ring
  have h1 : IsCU (((c - 1) * (d - 1)) • zz ^ (r + s)) :=
    (isCU_zz_pow (r + s)).smul (by nlinarith)
  have h2 : IsCU ((c - 1) • (zz ^ r * (zz ^ s + Pw s))) :=
    ((isCU_zz_pow r).mul (isCU_zz_pow_add_Pw s hs)).smul (by linarith)
  have h3 : IsCU ((d - 1) • (zz ^ s * (zz ^ r + Pw r))) :=
    ((isCU_zz_pow s).mul (isCU_zz_pow_add_Pw r hr)).smul (by linarith)
  rw [hsplit]
  simp only [coeffL_add]
  have hb := base_decr hr hs m hm
  have a1 := h1.decr m hm
  have a2 := h2.decr m hm
  have a3 := h3.decr m hm
  linarith

end ClanAudit
