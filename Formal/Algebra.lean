/-
  Algebraic lemmas for Erdős Problem #993.

  1. Star+Star w_2 positivity: the key inequality Q(a,b) > 0
     via completing the square.
  2. Log-concavity of binomial coefficients: C(n,k)² ≥ C(n,k-1)·C(n,k+1).
-/
import Mathlib

open Nat

noncomputable section

namespace Formal.Algebra

/-! ### Star+Star w_2 positivity -/

/-- Q(a,b) = a² - a·b + a + b² + b + 9. -/
def Q (a b : ℤ) : ℤ := a ^ 2 - a * b + a + b ^ 2 + b + 9

/-- Completing the square: 4·Q(a,b) = (2a - b + 1)² + 3·(b + 1)² + 32. -/
theorem Q_complete_square (a b : ℤ) :
    4 * Q a b = (2 * a - b + 1) ^ 2 + 3 * (b + 1) ^ 2 + 32 := by
  unfold Q; ring;

/-- Q(a,b) > 0 for all integers a, b. -/
theorem Q_pos (a b : ℤ) : Q a b > 0 := by
  unfold Q; exact by nlinarith [ sq_nonneg ( 2 * a - b + 1 ), sq_nonneg ( b + 1 ) ] ;

/-- Star+Star w_2 positivity: Q(a,b) > 0 for all a, b ≥ 1.
    (Follows from Q_pos which proves it for all integers.) -/
theorem starstar_w2_pos (a b : ℕ) (ha : 1 ≤ a) (hb : 1 ≤ b) :
    0 < Q (↑a) (↑b) := by
  exact Q_pos _ _

/-! ### Log-concavity of binomial coefficients -/

/-- The recurrence (k+1) · C(n, k+1) = (n - k) · C(n, k) for k < n. -/
theorem choose_succ_mul (n k : ℕ) (hk : k < n) :
    (k + 1) * Nat.choose n (k + 1) = (n - k) * Nat.choose n k := by
  rw [ Nat.mul_comm, Nat.choose_succ_right_eq ];
  ring

/-- The recurrence k · C(n, k) = (n - k + 1) · C(n, k - 1) for 1 ≤ k ≤ n. -/
theorem mul_choose_eq (n k : ℕ) (hk : 1 ≤ k) (hkn : k ≤ n) :
    k * Nat.choose n k = (n - k + 1) * Nat.choose n (k - 1) := by
  rcases k with ( _ | k ) <;> simp_all +decide [ Nat.add_one_mul_choose_eq ];
  nlinarith [ Nat.add_one_mul_choose_eq n k, Nat.choose_succ_succ n k, Nat.sub_add_cancel ( by linarith : k + 1 ≤ n ) ]

/-- Cross-multiplication identity:
    C(n,k-1) · C(n,k+1) · (k+1) · (n-k+1) = k · (n-k) · C(n,k)²
    for 1 ≤ k and k+1 ≤ n. -/
theorem choose_cross_mul (n k : ℕ) (hk : 1 ≤ k) (hkn : k + 1 ≤ n) :
    Nat.choose n (k - 1) * Nat.choose n (k + 1) * (k + 1) * (n - k + 1) =
    k * (n - k) * (Nat.choose n k * Nat.choose n k) := by
  have h_choose_succ_mul : (k + 1) * Nat.choose n (k + 1) = (n - k) * Nat.choose n k :=
    choose_succ_mul n k (by linarith)
  have h_mul_choose_eq : k * Nat.choose n k = (n - k + 1) * Nat.choose n (k - 1) := by
    convert mul_choose_eq n k hk ( by linarith ) using 1
  simp_all +decide [ mul_comm, mul_assoc, mul_left_comm ];
  grind +ring

/-- (k+1)(n-k+1) ≥ k(n-k) for natural numbers. -/
theorem succ_factor_le (n k : ℕ) (hk : 1 ≤ k) (hkn : k + 1 ≤ n) :
    k * (n - k) ≤ (k + 1) * (n - k + 1) := by
  nlinarith [ Nat.sub_add_cancel ( by linarith : k ≤ n ) ]

/-- Log-concavity of binomial coefficients:
    C(n,k)² ≥ C(n,k-1) · C(n,k+1). -/
theorem choose_log_concave (n k : ℕ) (hk : 1 ≤ k) (hkn : k + 1 ≤ n) :
    Nat.choose n (k - 1) * Nat.choose n (k + 1) ≤ Nat.choose n k * Nat.choose n k := by
  have h_mul : (Nat.choose n (k - 1)) * (Nat.choose n (k + 1)) * (k + 1) * (n - k + 1) ≤ (Nat.choose n k) * (Nat.choose n k) * (k + 1) * (n - k + 1) := by
    rw [ choose_cross_mul _ _ hk hkn ];
    grind;
  nlinarith [ show 0 < ( k + 1 ) * ( n - k + 1 ) by positivity ]

end Formal.Algebra
