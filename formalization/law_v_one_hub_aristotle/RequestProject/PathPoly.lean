/-
# Layer 1: the path independence polynomials

`W n` is the independence polynomial of the path on `n - 1` vertices, i.e. `W n = P (n-1)`
in the notation of the request (`P 0 = P (-1) = 1`).  Working with `W : ℕ → ℤ[X]` makes the
shifted families `P (aᵢ)` and `P (aᵢ - 1)` unambiguous: an arm of length `a` contributes
`W (a+1)` to `Q` and `W a` to `R`.

Main results:

* `W_coeff` : the path coefficient formula `[x^j] W n = choose (n - j) j`, equivalently
  `[x^j] P n = choose (n+1-j) j`.
* `W_niceD` : `W n` has positive coefficients exactly on `0, …, ⌊n/2⌋`.
* `W_LR_succ` : the adjacent likelihood-ratio inequality `W n ≤lr W (n+1)`.
* `W_LR_X` : the shifted comparison `W (n+1) ≤lr X * W n`.
* `W_logConcave` : each `W n` is log-concave.
* `path_choose_LR` : layer 1 in the binomial-coefficient form requested.

Source label: agent (layer 1).
-/
import RequestProject.LR

namespace LawV

open Polynomial

noncomputable section

/-- `W n = P (n-1)`: the independence polynomial of the path on `n-1` vertices. -/
def W : ℕ → ℤ[X]
  | 0 => 1
  | 1 => 1
  | (n+2) => W (n+1) + X * W n

@[simp] lemma W_zero : W 0 = 1 := rfl
@[simp] lemma W_one : W 1 = 1 := rfl
lemma W_succ_succ (n : ℕ) : W (n+2) = W (n+1) + X * W n := rfl

/-- **Layer 1, coefficient formula.**  `[x^j] W n = choose (n-j) j`. -/
lemma W_coeff : ∀ n j : ℕ, (W n).coeff j = ((n - j).choose j : ℤ)
  | 0, j => by
      rcases Nat.eq_zero_or_pos j with rfl | hj
      · simp
      · have e : 0 - j = 0 := by omega
        rw [W_zero, coeff_one, if_neg (by omega), e, Nat.choose_eq_zero_of_lt hj]
        norm_num
  | 1, j => by
      rcases Nat.eq_zero_or_pos j with rfl | hj
      · simp
      · have e : 1 - j = 0 := by omega
        rw [W_one, coeff_one, if_neg (by omega), e, Nat.choose_eq_zero_of_lt hj]
        norm_num
  | (n+2), j => by
      rw [W_succ_succ, coeff_add]
      rcases Nat.eq_zero_or_pos j with rfl | hj
      · simp [W_coeff (n+1) 0]
      · obtain ⟨m, rfl⟩ : ∃ m, j = m + 1 := ⟨j - 1, by omega⟩
        rw [coeff_X_mul, W_coeff (n+1) (m+1), W_coeff n m]
        rcases le_or_gt m n with hmn | hmn
        · have e1 : n + 1 - (m + 1) = n - m := by omega
          have e2 : n + 2 - (m + 1) = (n - m) + 1 := by omega
          rw [e1, e2, Nat.choose_succ_succ (n - m) m]
          push_cast
          ring
        · have e1 : n + 1 - (m + 1) = 0 := by omega
          have e2 : n - m = 0 := by omega
          have hm : 0 < m := by omega
          have z1 : Nat.choose 0 (m + 1) = 0 := Nat.choose_eq_zero_of_lt (by omega)
          have z2 : Nat.choose 0 m = 0 := Nat.choose_eq_zero_of_lt hm
          rw [e1, e2, z1, z2]
          have e3 : n + 2 - (m + 1) = 0 := by omega
          rw [e3, z1]
          norm_num

/-- The splitting identity for path polynomials: cutting a path at an interior vertex.
(Combinatorially: the path on `a+b+1` vertices split at its `(a+1)`-st vertex.) -/
lemma W_add : ∀ a b : ℕ, W (a + b + 2) = W (a+1) * W (b+1) + X * (W a * W b)
  | 0, b => by
      simp only [W_zero, W_one, Nat.zero_add, one_mul]
      exact W_succ_succ b
  | 1, b => by
      have h1 : 1 + b + 2 = (b + 1) + 2 := by omega
      have hW2 : W 2 = 1 + X := by rw [W_succ_succ]; simp
      rw [h1, W_succ_succ (b+1), W_succ_succ b, hW2]
      simp only [W_one, one_mul]
      ring
  | (a+2), b => by
      have ih1 : W (a + b + 3) = W (a+2) * W (b+1) + X * (W (a+1) * W b) := by
        have h := W_add (a+1) b
        rwa [show a + 1 + b + 2 = a + b + 3 by omega, show a + 1 + 1 = a + 2 by omega] at h
      have ih0 : W (a + b + 2) = W (a+1) * W (b+1) + X * (W a * W b) := W_add a b
      have hrec : W (a + b + 4) = W (a + b + 3) + X * W (a + b + 2) := by
        have h := W_succ_succ (a + b + 2)
        rwa [show a + b + 2 + 2 = a + b + 4 by omega, show a + b + 2 + 1 = a + b + 3 by omega] at h
      have hW3 : W (a + 3) = W (a + 2) + X * W (a + 1) := by
        have h := W_succ_succ (a + 1)
        rwa [show a + 1 + 2 = a + 3 by omega, show a + 1 + 1 = a + 2 by omega] at h
      have hW2' : W (a + 2) = W (a + 1) + X * W a := by
        have h := W_succ_succ a
        rwa [show a + 1 = a + 1 from rfl] at h
      rw [show a + 2 + b + 2 = a + b + 4 by omega, show a + 2 + 1 = a + 3 by omega,
        hrec, ih1, ih0, hW3, hW2']
      ring

/-- `W n` is positive exactly on the coefficients `0, …, ⌊n/2⌋`. -/
lemma W_niceD (n : ℕ) : NiceD (n / 2) (W n) := by
  constructor
  · intro j hj
    have hle : j ≤ n - j := by omega
    rw [W_coeff]
    exact_mod_cast Nat.choose_pos hle
  · intro j hj
    have hlt : n - j < j := by omega
    rw [W_coeff]
    exact_mod_cast Nat.choose_eq_zero_of_lt hlt

lemma W_nice (n : ℕ) : Nice (W n) := (W_niceD n).nice

lemma W_cf_nonneg (n : ℕ) (i : ℤ) : 0 ≤ cf (W n) i := (W_nice n).cf_nonneg i

/-- A polynomial with a single (positive) coefficient is log-concave. -/
lemma NiceD.logConcave_zero {p : ℤ[X]} (hp : NiceD 0 p) : LogConcave p := by
  intro k
  have hnn := hp.nice.cf_nonneg
  rcases lt_or_ge k 0 with hk | hk
  · rw [cf_of_neg p (by omega : k - 1 < 0), mul_zero]
    exact mul_nonneg (hnn _) (hnn _)
  · rw [hp.cf_eq_zero (by omega : ((0:ℕ) : ℤ) < k + 1), zero_mul]
    exact mul_nonneg (hnn _) (hnn _)

/-- The two adjacent comparisons, proved by a simultaneous induction. -/
lemma W_LR_pair : ∀ n : ℕ, LR (W n) (W (n+1)) ∧ LR (W (n+1)) (X * W n)
  | 0 => by
      refine ⟨by rw [W_zero, W_one]; exact LR.refl _, ?_⟩
      rw [W_zero, W_one]
      exact logConcave_iff_LR.1 NiceD.one.logConcave_zero
  | (n+1) => by
      obtain ⟨hA, hC⟩ := W_LR_pair n
      have hB : LR (W (n+1)) (X * W (n+1)) :=
        LR.trans_X (W_nice (n+1)) (W_niceD n) (W_cf_nonneg (n+1)) hC hA
      refine ⟨?_, ?_⟩
      · show LR (W (n+1)) (W (n+2))
        rw [W_succ_succ]
        exact LR.add_right (LR.refl _) hC
      · show LR (W (n+2)) (X * W (n+1))
        rw [W_succ_succ]
        exact LR.add_left hB (LR_X_mul hA)

/-- **Layer 1, adjacent likelihood-ratio inequality.** -/
lemma W_LR_succ (n : ℕ) : LR (W n) (W (n+1)) := (W_LR_pair n).1

/-- The shifted comparison `W (n+1) ≤lr X * W n`. -/
lemma W_LR_X (n : ℕ) : LR (W (n+1)) (X * W n) := (W_LR_pair n).2

/-- Each path polynomial is log-concave. -/
lemma W_logConcave : ∀ n : ℕ, LogConcave (W n)
  | 0 => by rw [W_zero]; exact NiceD.one.logConcave_zero
  | (n+1) =>
      logConcave_iff_LR.2
        (LR.trans_X (W_nice (n+1)) (W_niceD n) (W_cf_nonneg (n+1)) (W_LR_X n) (W_LR_succ n))

/-- All `2×2` minors of two consecutive path polynomials are nonnegative. -/
lemma W_LRgen_succ (n : ℕ) : LRgen (W n) (W (n+1)) :=
  LR.gen (W_LR_succ n) (W_niceD n) (W_niceD (n+1))

/-- `W n ≤lr X^2 * W (n-1)`; used for the shifted cross terms. -/
lemma W_LR_X_X (n : ℕ) : LR (W (n+1)) (X * (X * W n)) :=
  LR.trans_X (W_nice (n+1)) (W_niceD n) (fun i => by
    rw [cf_X_mul]; exact W_cf_nonneg n _) (W_LR_X n)
    (logConcave_iff_LR.1 (W_logConcave n))

/-- **Layer 1 in binomial form.**  With `p (n, j) = choose (n + 1 - j) j` the coefficient of
`x^j` in the independence polynomial of the path on `n` vertices, this says
`p (n, j) * p (n-1, j-1) ≥ p (n, j-1) * p (n-1, j)`, with all support boundaries included. -/
theorem path_choose_LR (n j : ℕ) :
    (n - (j+1)).choose (j+1) * (n + 1 - j).choose j
      ≤ (n - j).choose j * (n - j).choose (j+1) := by
  have h := W_LR_succ n (j : ℤ)
  have e1 : ((j : ℤ) + 1) = ((j + 1 : ℕ) : ℤ) := by push_cast; ring
  rw [e1] at h
  simp only [cf_ofNat, W_coeff] at h
  have e2 : n + 1 - (j + 1) = n - j := by omega
  rw [e2] at h
  exact_mod_cast h

end

end LawV
