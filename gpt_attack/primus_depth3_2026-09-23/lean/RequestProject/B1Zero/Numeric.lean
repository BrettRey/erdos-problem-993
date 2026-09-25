import Mathlib.Tactic

/-!
# The finite parameter calculation for at most three forbidden vertices

`density_numeric` is a purely arithmetic statement: for the twelve admissible parameter
triples `(α, |Forced|, |Forbidden|)` in the window with at most three forbidden vertices, the
coefficientwise blocked bound is at most `(α-7)/(3(α-3))` times the extendable count, given
only the matching-bag deletion inequalities.  The proof is a finite case split followed by
integer linear arithmetic (a Farkas certificate for each case).
-/

namespace DepthThree.B1Zero

set_option maxHeartbeats 4000000 in
/-- The finite parameter calculation.  Here `p j` is the number of independent `j`-subsets of
the flexible vertices, `a` the independence number, `c` the number of forced vertices and `r`
the number of forbidden vertices. -/
theorem density_numeric (a c r : ℕ) (p : ℕ → ℕ)
    (ha1 : 17 ≤ a) (ha2 : a ≤ 19) (hca : c ≤ a) (hrc0 : r ≤ c)
    (hn1 : 33 + (c - r) ≤ 2 * a) (hn2 : 2 * a ≤ 38 + (c - r))
    (hrc : 1 ≤ r → 2 * r + 1 ≤ c) (hr3 : r ≤ 3)
    (hz : ∀ j, a - c < j → p j = 0)
    (hcone : ∀ l, (l + 1) * p (l + 1) ≤ 2 * (a - c - l) * p l) :
    3 * (a - 3) * (∑ t ∈ Finset.Icc 1 3, Nat.choose r t *
        ∑ i ∈ Finset.range (a - 4 - t + 1), Nat.choose (c - (2 * t + 1)) i * p (a - 4 - t - i))
      ≤ (a - 7) * ∑ i ∈ Finset.range (a - 4 + 1), Nat.choose c i * p (a - 4 - i) := by
  have hsum3 : ∀ f : ℕ → ℕ, ∑ t ∈ (Finset.Icc 1 3 : Finset ℕ), f t = f 1 + f 2 + f 3 := by
    intro f
    rw [show (Finset.Icc 1 3 : Finset ℕ) = {1, 2, 3} by decide]
    simp [Finset.sum_insert, add_assoc]
  rw [hsum3]
  have hc5 : 6 * p 6 ≤ 2 * (a - c - 5) * p 5 := hcone 5
  have hc6 : 7 * p 7 ≤ 2 * (a - c - 6) * p 6 := hcone 6
  have hc7 : 8 * p 8 ≤ 2 * (a - c - 7) * p 7 := hcone 7
  have hc8 : 9 * p 9 ≤ 2 * (a - c - 8) * p 8 := hcone 8
  have hc9 : 10 * p 10 ≤ 2 * (a - c - 9) * p 9 := hcone 9
  have hc10 : 11 * p 11 ≤ 2 * (a - c - 10) * p 10 := hcone 10
  have hc11 : 12 * p 12 ≤ 2 * (a - c - 11) * p 11 := hcone 11
  have hc12 : 13 * p 13 ≤ 2 * (a - c - 12) * p 12 := hcone 12
  have hc13 : 14 * p 14 ≤ 2 * (a - c - 13) * p 13 := hcone 13
  have hc14 : 15 * p 15 ≤ 2 * (a - c - 14) * p 14 := hcone 14
  have hc15 : 16 * p 16 ≤ 2 * (a - c - 15) * p 15 := hcone 15
  have hc16 : 17 * p 17 ≤ 2 * (a - c - 16) * p 16 := hcone 16
  have hc17 : 18 * p 18 ≤ 2 * (a - c - 17) * p 17 := hcone 17
  clear hcone
  interval_cases a <;> interval_cases c <;> interval_cases r <;>
    first
      | omega
      | (norm_num [Finset.sum_range_succ, hz, Nat.choose] <;> omega)

end DepthThree.B1Zero
