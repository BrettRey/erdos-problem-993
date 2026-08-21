import Mathlib
import RequestProject.Spider

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

/-!
# LAW-V when one hub is trivial: the final assembly

This file collects the final statements.

* `Cof_nil_eq`, `Cof_single_eq`, `Cof_pair_eq` : a spider with at most two arms *is* a path,
  so its independence polynomial is a path polynomial.
* `Cof_logConcave_of_length_le_two` : layer 4 for at most two arms (unconditional).
* `V_nonneg_of_length_le_two` : the LAW-V inequality, unconditional, for at most two arms —
  in particular the empty family and the one-arm family.
* `V_nonneg_conditional` : the LAW-V inequality for an arbitrary family of arms, taking the
  spider log-concavity layer as an explicit hypothesis (never as an axiom).
* `two_hub_cross_term_counterexample` : the optional adversarial check.

Source label: agent (final assembly).
-/

namespace LawV

open Polynomial

noncomputable section

/-! ### Small spiders are paths -/

/-- The spider with no arms is the path on one vertex. -/
lemma Cof_nil_eq : Cof [] = W 2 := by
  rw [Cof, Qof_nil, Rof_nil, W_succ_succ]
  simp

/-- A spider with one arm of length `a` is the path on `a+1` vertices. -/
lemma Cof_single_eq (a : ℕ) : Cof [a] = W (a + 2) := by
  rw [Cof, Qof_cons, Rof_cons, Qof_nil, Rof_nil, W_succ_succ]
  ring

/-- A spider with two arms of lengths `a` and `b` is the path on `a+b+1` vertices. -/
lemma Cof_pair_eq (a b : ℕ) : Cof [a, b] = W (a + b + 2) := by
  rw [Cof, Qof_cons, Qof_cons, Rof_cons, Rof_cons, Qof_nil, Rof_nil, W_add a b]
  ring

/-- **Layer 4 for at most two arms**: the spider is then a path, hence log-concave. -/
theorem Cof_logConcave_of_length_le_two (arms : List ℕ) (h : arms.length ≤ 2) :
    LogConcave (Cof arms) := by
  match arms, h with
  | [], _ => rw [Cof_nil_eq]; exact W_logConcave 2
  | [a], _ => rw [Cof_single_eq]; exact W_logConcave (a + 2)
  | [a, b], _ => rw [Cof_pair_eq]; exact W_logConcave (a + b + 2)

/-- **The LAW-V one-hub inequality, unconditionally, for families of at most two arms.**
This includes the empty family and the one-arm family. -/
theorem V_nonneg_of_length_le_two (arms : List ℕ) (h : arms.length ≤ 2) (k : ℤ) :
    0 ≤ Vco arms k :=
  V_nonneg_of_logConcave arms (Cof_logConcave_of_length_le_two arms h) k

/-- The empty family. -/
theorem V_nonneg_nil (k : ℤ) : 0 ≤ Vco [] k :=
  V_nonneg_of_length_le_two [] (by simp) k

/-- The one-arm family. -/
theorem V_nonneg_single (a : ℕ) (k : ℤ) : 0 ≤ Vco [a] k :=
  V_nonneg_of_length_le_two [a] (by simp) k

/-! ### Stars: layer 4, unconditionally, for an infinite family with arbitrarily many arms

For the star `K_{1,m}` (`m` arms of length `1`) we have `Q = (1+x)^m`, `R = 1`, hence
`C = (1+x)^m + x`.  Log-concavity of `C` then reduces to the binomial inequality
`choose m 3 * (m+1) ≤ (choose m 2)^2` at the single index `k = 2`; all other indices follow
from log-concavity of `Q`.  Note that no splitting of the likelihood-ratio order can give
this case: `Q ≤lr x²*R = x²` is already false for `m = 3`. -/

/-- Coefficient sequence of `x` itself, at integer indices. -/
lemma cf_X_eq (i : ℤ) : cf (X : ℤ[X]) i = if i = 1 then 1 else 0 := by
  have h : (X : ℤ[X]) = X * 1 := by ring
  rw [h, cf_X_mul]
  unfold cf
  split
  · rw [coeff_one]
    by_cases h0 : i = 1
    · subst h0; simp
    · rw [if_neg (by omega), if_neg h0]
  · rw [if_neg (by omega)]

/-- For the star, `Q = (1+x)^m`. -/
lemma Qof_star (m : ℕ) : Qof (List.replicate m 1) = (1 + X) ^ m := by
  induction m with
  | zero => simp
  | succ n ih =>
      rw [List.replicate_succ, Qof_cons, ih, pow_succ]
      have hW : W (1 + 1) = 1 + X := by rw [W_succ_succ]; simp
      rw [hW]; ring

/-- For the star, `R = 1`. -/
lemma Rof_star (m : ℕ) : Rof (List.replicate m 1) = 1 := by
  induction m with
  | zero => simp
  | succ n ih => rw [List.replicate_succ, Rof_cons, ih, W_one, one_mul]

/-- The independence polynomial of the star `K_{1,m}`. -/
lemma Cof_star (m : ℕ) : Cof (List.replicate m 1) = (1 + X) ^ m + X := by
  rw [Cof, Qof_star, Rof_star, mul_one]

/-- The one binomial inequality behind log-concavity of the star polynomial. -/
lemma star_choose_key (m : ℕ) : m.choose 3 * (m + 1) ≤ (m.choose 2) ^ 2 := by
  rcases lt_or_ge m 3 with h | h
  · rw [Nat.choose_eq_zero_of_lt h]; simp
  · have e3 : m.choose 3 * 3 = m.choose 2 * (m - 2) := Nat.choose_succ_right_eq m 2
    have e2 : m.choose 2 * 2 = m * (m - 1) := by
      have := Nat.choose_succ_right_eq m 1
      simpa [Nat.choose_one_right] using this
    have key : 2 * ((m - 2) * (m + 1)) ≤ 3 * (m * (m - 1)) := by
      obtain ⟨n, rfl⟩ : ∃ n, m = n + 3 := ⟨m - 3, by omega⟩
      have h1 : n + 3 - 2 = n + 1 := by omega
      have h2 : n + 3 - 1 = n + 2 := by omega
      rw [h1, h2]; nlinarith
    nlinarith [Nat.zero_le (m.choose 2)]

/-- **Layer 4 for stars**: the independence polynomial of `K_{1,m}` is log-concave, for
every number `m` of arms. -/
theorem Cof_logConcave_star (m : ℕ) : LogConcave (Cof (List.replicate m 1)) := by
  set Q : ℤ[X] := (1 + X) ^ m with hQdef
  have hQ : LogConcave Q := by
    have := Qof_logConcave (List.replicate m 1); rwa [Qof_star] at this
  have hnn : ∀ i : ℤ, 0 ≤ cf Q i := by
    intro i
    have := Qof_cf_nonneg (List.replicate m 1) i
    rwa [Qof_star] at this
  have hc : ∀ i : ℤ, cf (Cof (List.replicate m 1)) i = cf Q i + (if i = 1 then 1 else 0) := by
    intro i; rw [Cof_star, cf_add, cf_X_eq]
  have hcoeff : ∀ j : ℕ, cf Q (j : ℤ) = (m.choose j : ℤ) := by
    intro j; rw [hQdef, cf_ofNat, coeff_one_add_X_pow]
  intro k
  rcases lt_or_ge k 1 with hk | hk
  · have h1 : cf (Cof (List.replicate m 1)) (k - 1) = 0 := by
      rw [hc, cf_of_neg _ (by omega), if_neg (by omega)]
      simp
    rw [h1, mul_zero]
    exact mul_self_nonneg _
  rcases eq_or_lt_of_le hk with hk1 | hk1
  · have hk1' : k = 1 := hk1.symm
    subst hk1'
    rw [hc, hc, hc]
    norm_num
    have h := hQ 1
    norm_num at h
    nlinarith [h, hnn 1]
  · have hk2 : (2 : ℤ) ≤ k := hk1
    rcases eq_or_lt_of_le hk2 with hk3 | hk3
    · have hk3' : k = 2 := hk3.symm
      subst hk3'
      rw [hc, hc, hc]
      norm_num
      have e3 : cf Q 3 = (m.choose 3 : ℤ) := by simpa using hcoeff 3
      have e2 : cf Q 2 = (m.choose 2 : ℤ) := by simpa using hcoeff 2
      have e1 : cf Q 1 = (m : ℤ) := by simpa using hcoeff 1
      rw [e3, e2, e1]
      have hkey : (m.choose 3 : ℤ) * ((m : ℤ) + 1) ≤ ((m.choose 2 : ℤ)) ^ 2 := by
        exact_mod_cast star_choose_key m
      nlinarith [hkey]
    · have h1 : k + 1 ≠ 1 := by omega
      have h2 : k ≠ 1 := by omega
      have h3 : k - 1 ≠ 1 := by omega
      rw [hc, hc, hc, if_neg h1, if_neg h2, if_neg h3]
      simpa using hQ k

/-- **The LAW-V one-hub inequality, unconditionally, for every star `K_{1,m}`.** -/
theorem V_nonneg_star (m : ℕ) (k : ℤ) : 0 ≤ Vco (List.replicate m 1) k :=
  V_nonneg_of_logConcave _ (Cof_logConcave_star m) k

/-! ### The general theorem, conditional on the spider log-concavity layer

The following theorem is the requested LAW-V statement for an arbitrary family of arms.
Layers 1, 2, 3 and 5 are proved unconditionally in this development; the spider
log-concavity layer (layer 4) enters as the explicit hypothesis `hLC`.  It is a hypothesis,
never an axiom: discharging

```
theorem Cof_logConcave (arms : List ℕ) : LogConcave (Cof arms)
```

turns `V_nonneg_conditional` into an unconditional theorem. -/
theorem V_nonneg_conditional (hLC : ∀ arms : List ℕ, LogConcave (Cof arms))
    (arms : List ℕ) (k : ℤ) : 0 ≤ Vco arms k :=
  V_nonneg_of_logConcave arms (hLC arms) k

/-- **Sharpening of the remaining gap.**  Layer 4 follows from the single likelihood-ratio
inequality `Q ≤lr x*C`, i.e. from `q (k+1) * c (k-1) ≤ q k * c k` for all `k`: writing
`C = Q + x*R`, the left summand is handled by that hypothesis and the right one by the
unconditional `Rof_LR_Cof`.  So the whole development is unconditional except for

```
theorem Qof_LR_X_Cof (arms : List ℕ) : LR (Qof arms) (X * Cof arms)
```

which is equivalent to the coefficient inequality
`q (k+1) * r (k-2) ≤ q k * r (k-1) + (q k ^ 2 - q (k+1) * q (k-1))`,
a statement about products of path polynomials only. -/
theorem Cof_logConcave_of_K (arms : List ℕ) (hK : LR (Qof arms) (X * Cof arms)) :
    LogConcave (Cof arms) := by
  rw [logConcave_iff_LR]
  show LR (Qof arms + X * Rof arms) (X * Cof arms)
  exact LR.add_left hK (LR_X_mul (Rof_LR_Cof arms))

/-- The LAW-V inequality for an arbitrary family of arms, conditional on the sharper
hypothesis `Q ≤lr x*C` instead of log-concavity of `C`. -/
theorem V_nonneg_of_K (arms : List ℕ) (hK : LR (Qof arms) (X * Cof arms)) (k : ℤ) :
    0 ≤ Vco arms k :=
  V_nonneg_of_logConcave arms (Cof_logConcave_of_K arms hK) k

/-! ### Adversarial check: the naive two-hub cross term

For two spiders the naive cross term `u k * v (k-1) − u (k+1) * v (k-2)`, with `u`, `v` the
two spider polynomials, is *not* nonnegative: for arms `(1,1,5)` and `(1,1)` it equals `-3`
at `k = 3`.  This does not affect the one-trivial-hub theorem, where the second factor is
the arm product `Q`. -/
theorem two_hub_cross_term_counterexample :
    cf (Cof [1, 1, 5]) 3 * cf (Cof [1, 1]) 2 - cf (Cof [1, 1, 5]) 4 * cf (Cof [1, 1]) 1
      = -3 := by
  have h1 : Cof [1, 1] = 1 + 3 * X + X ^ 2 := by simp [Cof, Qof, Rof, W]; ring
  have h2 : Cof [1, 1, 5] = 1 + 8 * X + 21 * X ^ 2 + 21 * X ^ 3 + 8 * X ^ 4 + X ^ 5 := by
    simp [Cof, Qof, Rof, W]; ring
  have e1 : (1 : ℤ) = ((1 : ℕ) : ℤ) := by norm_num
  have e2 : (2 : ℤ) = ((2 : ℕ) : ℤ) := by norm_num
  have e3 : (3 : ℤ) = ((3 : ℕ) : ℤ) := by norm_num
  have e4 : (4 : ℤ) = ((4 : ℕ) : ℤ) := by norm_num
  rw [e1, e2, e3, e4, cf_ofNat, cf_ofNat, cf_ofNat, cf_ofNat, h1, h2]
  simp [coeff_one, coeff_X, coeff_X_pow]

/-! ### Axiom audit

`#print axioms` for the final theorems and the main closure lemmas.  Only the three
standard Lean axioms `propext`, `Classical.choice`, `Quot.sound` may appear. -/

section AxiomAudit

#print axioms path_choose_LR                -- layer 1
#print axioms LR.mul_left                   -- layer 2
#print axioms Rof_LR_Qof                    -- layer 3
#print axioms cross_term_nonneg             -- layer 5
#print axioms V_nonneg_of_logConcave        -- assembly
#print axioms V_nonneg_of_length_le_two     -- unconditional LAW-V, at most two arms
#print axioms V_nonneg_star                 -- unconditional LAW-V, stars K_{1,m}
#print axioms V_nonneg_conditional          -- LAW-V, conditional on layer 4
#print axioms V_nonneg_of_K                 -- LAW-V, conditional on the sharpened layer 4
#print axioms two_hub_cross_term_counterexample

end AxiomAudit

end

end LawV
