/-
Source: derived for the even-connector Laurent block request
(`c2_even_connector_aristotle_input_20260820.md`).

This file contains the purely combinatorial core of the proof: a family of
inequalities between binomial coefficients that expresses the unimodality of the
coefficient polynomial of the even-connector block.

The central statement is `key_three`:

  for `r, s ≥ 1`, `N = r + s + 1` and `2k + 2 ≤ N`,

      C(N,k) + C(r,k) + C(s,k) + [k = 0]  ≤  C(N,k+1) + C(r,k+1) + C(s,k+1).

The proof goes through the Catalan numbers:

* `choose_sub_succ_le_catalan` : `C(n,k) - C(n,k+1) ≤ catalan k` for every `n`
  (the maximum over `n` is attained at `n = 2k`);
* `catalan_succ_le_choose_diff` : `catalan (k+1) ≤ C(n,k+1) - C(n,k)`
  whenever `2k + 2 ≤ n`;
* `two_catalan_le_catalan_succ` : `2 * catalan k ≤ catalan (k+1)` for `k ≥ 1`.
-/
import Mathlib

namespace EvenConnector

/-! ### Monotonicity of a binomial row -/

/-- Binomial coefficients increase before the middle of the row. -/
theorem choose_le_succ_of_le (n k : ℕ) (h : 2 * k + 1 ≤ n) : n.choose k ≤ n.choose (k + 1) := by
  have hkey := Nat.choose_succ_right_eq n k
  have h1 : k + 1 ≤ n - k := by omega
  have h2 : n.choose k * (k + 1) ≤ n.choose (k + 1) * (k + 1) := by
    rw [hkey]; exact Nat.mul_le_mul_left _ h1
  exact Nat.le_of_mul_le_mul_right h2 (by omega)

/-- Binomial coefficients decrease after the middle of the row. -/
theorem choose_succ_le_of_ge (n k : ℕ) (h : n ≤ 2 * k + 1) : n.choose (k + 1) ≤ n.choose k := by
  have hkey := Nat.choose_succ_right_eq n k
  have h1 : n - k ≤ k + 1 := by omega
  have h2 : n.choose (k + 1) * (k + 1) ≤ n.choose k * (k + 1) := by
    rw [hkey]; exact Nat.mul_le_mul_left _ h1
  exact Nat.le_of_mul_le_mul_right h2 (by omega)

/-! ### Catalan numbers as the central descent -/

/-- `C(2k,k) - C(2k,k+1) = catalan k`. -/
theorem central_diff_eq_catalan (k : ℕ) :
    ((2 * k).choose k : ℤ) - (2 * k).choose (k + 1) = catalan k := by
  have h1 : (2 * k).choose k = (k + 1) * catalan k := by
    rw [succ_mul_catalan_eq_centralBinom]; rfl
  have h2 := Nat.choose_succ_right_eq (2 * k) k
  have h3 : 2 * k - k = k := by omega
  rw [h3, h1] at h2
  have h4 : (2 * k).choose (k + 1) = k * catalan k := by
    have h5 : (2 * k).choose (k + 1) * (k + 1) = k * catalan k * (k + 1) := by rw [h2]; ring
    exact Nat.eq_of_mul_eq_mul_right (by omega) h5
  rw [h1, h4]
  push_cast
  ring

/-- The descent `n ↦ C(n,k) - C(n,k+1)`. -/
private def gg (k n : ℕ) : ℤ := (n.choose k : ℤ) - n.choose (k + 1)

private theorem gg_step (m n : ℕ) :
    gg (m + 1) (n + 1) - gg (m + 1) n = (n.choose m : ℤ) - n.choose (m + 1) := by
  simp only [gg, Nat.choose_succ_succ n m, Nat.choose_succ_succ n (m + 1)]
  push_cast
  ring

private theorem gg_up (m n : ℕ) (h : n + 1 ≤ 2 * (m + 1)) : gg (m + 1) n ≤ gg (m + 1) (n + 1) := by
  have hstep := gg_step m n
  have h2 : n.choose (m + 1) ≤ n.choose m := choose_succ_le_of_ge n m (by omega)
  have h3 : ((n.choose (m + 1) : ℤ)) ≤ (n.choose m : ℤ) := by exact_mod_cast h2
  omega

private theorem gg_down (m n : ℕ) (h : 2 * (m + 1) ≤ n) :
    gg (m + 1) (n + 1) ≤ gg (m + 1) n := by
  have hstep := gg_step m n
  have h2 : n.choose m ≤ n.choose (m + 1) := choose_le_succ_of_le n m (by omega)
  have h3 : ((n.choose m : ℤ)) ≤ (n.choose (m + 1) : ℤ) := by exact_mod_cast h2
  omega

/-- The descent of a binomial row is bounded by a Catalan number, uniformly in `n`. -/
theorem choose_sub_succ_le_catalan (n k : ℕ) :
    (n.choose k : ℤ) - n.choose (k + 1) ≤ catalan k := by
  match k with
  | 0 => simp [Nat.choose_one_right, catalan_zero]
  | (m + 1) =>
    have hbase : gg (m + 1) (2 * (m + 1)) = catalan (m + 1) := central_diff_eq_catalan (m + 1)
    have hmax : ∀ d : ℕ, gg (m + 1) (2 * (m + 1) + d) ≤ gg (m + 1) (2 * (m + 1)) := by
      intro d
      induction d with
      | zero => simp
      | succ e ih =>
        have hd := gg_down m (2 * (m + 1) + e) (by omega)
        have h' : 2 * (m + 1) + (e + 1) = 2 * (m + 1) + e + 1 := by ring
        rw [h']
        exact le_trans hd ih
    have hmin : ∀ d : ℕ, gg (m + 1) (2 * (m + 1) - d) ≤ gg (m + 1) (2 * (m + 1)) := by
      intro d
      induction d with
      | zero => simp
      | succ e ih =>
        rcases Nat.lt_or_ge e (2 * (m + 1)) with he | he
        · have h1 : 2 * (m + 1) - e = 2 * (m + 1) - (e + 1) + 1 := by omega
          have hu := gg_up m (2 * (m + 1) - (e + 1)) (by omega)
          rw [← h1] at hu
          exact le_trans hu ih
        · have h1 : 2 * (m + 1) - (e + 1) = 2 * (m + 1) - e := by omega
          rw [h1]
          exact ih
    rcases Nat.lt_or_ge n (2 * (m + 1)) with hn | hn
    · have hle := hmin (2 * (m + 1) - n)
      have h2 : 2 * (m + 1) - (2 * (m + 1) - n) = n := by omega
      rw [h2, hbase] at hle
      exact hle
    · have hle := hmax (n - 2 * (m + 1))
      have h2 : 2 * (m + 1) + (n - 2 * (m + 1)) = n := by omega
      rw [h2, hbase] at hle
      exact hle

/-! ### The ascent `C(n,k+1) - C(n,k)` grows with `n` -/

/-- For `2k + 1 ≤ n`, the ascent at `k` is nondecreasing in `n`. -/
theorem ascent_mono (n k : ℕ) (h : 2 * k + 1 ≤ n) :
    (n.choose (k + 1) : ℤ) - n.choose k ≤ ((n + 1).choose (k + 1) : ℤ) - (n + 1).choose k := by
  match k with
  | 0 =>
    simp [Nat.choose_one_right]
  | (m + 1) =>
    have e1 : (n + 1).choose (m + 2) = n.choose (m + 1) + n.choose (m + 2) :=
      Nat.choose_succ_succ n (m + 1)
    have e2 : (n + 1).choose (m + 1) = n.choose m + n.choose (m + 1) :=
      Nat.choose_succ_succ n m
    have h2 : n.choose m ≤ n.choose (m + 1) := choose_le_succ_of_le n m (by omega)
    have h3 : ((n.choose m : ℤ)) ≤ (n.choose (m + 1) : ℤ) := by exact_mod_cast h2
    have e1' : ((n + 1).choose (m + 2) : ℤ) = n.choose (m + 1) + n.choose (m + 2) := by
      exact_mod_cast e1
    have e2' : ((n + 1).choose (m + 1) : ℤ) = n.choose m + n.choose (m + 1) := by
      exact_mod_cast e2
    show ((n.choose (m + 2) : ℤ)) - n.choose (m + 1) ≤ _
    rw [e1', e2']
    linarith

/-- If `2k + 2 ≤ n` then the ascent `C(n,k+1) - C(n,k)` is at least `catalan (k+1)`. -/
theorem catalan_succ_le_choose_diff (n k : ℕ) (h : 2 * k + 2 ≤ n) :
    (catalan (k + 1) : ℤ) ≤ (n.choose (k + 1) : ℤ) - n.choose k := by
  have base : (catalan (k + 1) : ℤ) =
      ((2 * k + 2).choose (k + 1) : ℤ) - (2 * k + 2).choose k := by
    have hs : (2 * k + 2).choose k = (2 * k + 2).choose (k + 2) := by
      have : (2 * k + 2).choose (2 * k + 2 - k) = (2 * k + 2).choose k := Nat.choose_symm (by omega)
      have h2 : 2 * k + 2 - k = k + 2 := by omega
      rw [h2] at this
      exact this.symm
    rw [hs]
    have := central_diff_eq_catalan (k + 1)
    have h2 : 2 * (k + 1) = 2 * k + 2 := by ring
    rw [h2] at this
    linarith [this]
  have main : ∀ d : ℕ, (catalan (k + 1) : ℤ) ≤
      ((2 * k + 2 + d).choose (k + 1) : ℤ) - (2 * k + 2 + d).choose k := by
    intro d
    induction d with
    | zero => simpa using base.le
    | succ e ih =>
      have hm := ascent_mono (2 * k + 2 + e) k (by omega)
      have h' : 2 * k + 2 + (e + 1) = 2 * k + 2 + e + 1 := by ring
      rw [h']
      linarith
  have := main (n - (2 * k + 2))
  have h2 : 2 * k + 2 + (n - (2 * k + 2)) = n := by omega
  rwa [h2] at this

/-! ### Doubling of Catalan numbers -/

/-- `2 * catalan k ≤ catalan (k+1)` for `k ≥ 1`. -/
theorem two_catalan_le_catalan_succ (k : ℕ) (hk : 1 ≤ k) : 2 * catalan k ≤ catalan (k + 1) := by
  have e1 : (k + 1) * catalan k = k.centralBinom := succ_mul_catalan_eq_centralBinom k
  have e2 : (k + 1 + 1) * catalan (k + 1) = (k + 1).centralBinom :=
    succ_mul_catalan_eq_centralBinom (k + 1)
  have e3 : (k + 1) * (k + 1).centralBinom = 2 * (2 * k + 1) * k.centralBinom :=
    Nat.succ_mul_centralBinom_succ k
  -- `(k+2) * catalan (k+1) = 2 * (2k+1) * catalan k`
  have e4 : (k + 1) * ((k + 2) * catalan (k + 1)) = (k + 1) * (2 * (2 * k + 1) * catalan k) := by
    calc (k + 1) * ((k + 2) * catalan (k + 1)) = (k + 1) * (k + 1).centralBinom := by
          rw [show k + 2 = k + 1 + 1 from rfl, e2]
      _ = 2 * (2 * k + 1) * k.centralBinom := e3
      _ = 2 * (2 * k + 1) * ((k + 1) * catalan k) := by rw [e1]
      _ = (k + 1) * (2 * (2 * k + 1) * catalan k) := by ring
  have e5 : (k + 2) * catalan (k + 1) = 2 * (2 * k + 1) * catalan k :=
    Nat.eq_of_mul_eq_mul_left (by omega) e4
  have h6 : (k + 2) * (2 * catalan k) ≤ (k + 2) * catalan (k + 1) := by
    rw [e5]
    have : (k + 2) * 2 ≤ 2 * (2 * k + 1) := by omega
    calc (k + 2) * (2 * catalan k) = ((k + 2) * 2) * catalan k := by ring
      _ ≤ (2 * (2 * k + 1)) * catalan k := Nat.mul_le_mul_right _ this
  exact Nat.le_of_mul_le_mul_left h6 (by omega)

theorem catalan_le_catalan_succ (k : ℕ) : catalan k ≤ catalan (k + 1) := by
  rcases Nat.eq_zero_or_pos k with hk | hk
  · subst hk; simp [catalan_zero, catalan_one]
  · have := two_catalan_le_catalan_succ k hk
    omega

/-! ### The two master inequalities -/

/-- Two-term version: `C(N,k) + C(r,k) ≤ C(N,k+1) + C(r,k+1)` for `N = r+s+1`, `2k+2 ≤ N`. -/
theorem key_two (r s k : ℕ) (hr : 1 ≤ r) (hs : 1 ≤ s) (hk : 2 * k + 2 ≤ r + s + 1) :
    ((r + s + 1).choose k : ℤ) + r.choose k ≤ ((r + s + 1).choose (k + 1) : ℤ) + r.choose (k + 1) := by
  rcases Nat.eq_zero_or_pos k with hk0 | hk0
  · subst hk0
    have c1 : (r + s + 1).choose (0 + 1) = r + s + 1 := by simp
    have c2 : r.choose (0 + 1) = r := by simp
    have hr' : (1 : ℤ) ≤ (r : ℤ) := by exact_mod_cast hr
    have hs' : (0 : ℤ) ≤ (s : ℤ) := by positivity
    simp only [Nat.choose_zero_right, c1, c2]
    push_cast
    linarith
  · obtain ⟨m, rfl⟩ : ∃ m, k = m + 1 := ⟨k - 1, by omega⟩
    have h1 := catalan_succ_le_choose_diff (r + s + 1) (m + 1) hk
    have h2 := choose_sub_succ_le_catalan r (m + 1)
    have h3 := catalan_le_catalan_succ (m + 1)
    have h3' : (catalan (m + 1) : ℤ) ≤ catalan (m + 2) := by exact_mod_cast h3
    linarith

/-- Three-term version. -/
theorem key_three (r s k : ℕ) (hr : 1 ≤ r) (hs : 1 ≤ s) (hk : 2 * k + 2 ≤ r + s + 1) :
    ((r + s + 1).choose k : ℤ) + r.choose k + s.choose k + (if k = 0 then 1 else 0) ≤
      ((r + s + 1).choose (k + 1) : ℤ) + r.choose (k + 1) + s.choose (k + 1) := by
  rcases Nat.eq_zero_or_pos k with hk0 | hk0
  · subst hk0
    have c1 : (r + s + 1).choose (0 + 1) = r + s + 1 := by simp
    have c2 : r.choose (0 + 1) = r := by simp
    have c3 : s.choose (0 + 1) = s := by simp
    have hr' : (1 : ℤ) ≤ (r : ℤ) := by exact_mod_cast hr
    have hs' : (1 : ℤ) ≤ (s : ℤ) := by exact_mod_cast hs
    simp only [Nat.choose_zero_right, c1, c2, c3]
    push_cast
    linarith
  · obtain ⟨m, rfl⟩ : ∃ m, k = m + 1 := ⟨k - 1, by omega⟩
    have h1 := catalan_succ_le_choose_diff (r + s + 1) (m + 1) hk
    have h2 := choose_sub_succ_le_catalan r (m + 1)
    have h3 := choose_sub_succ_le_catalan s (m + 1)
    have h4 := two_catalan_le_catalan_succ (m + 1) (by omega)
    have h4' : 2 * (catalan (m + 1) : ℤ) ≤ catalan (m + 2) := by exact_mod_cast h4
    have hif : (if m + 1 = 0 then (1 : ℤ) else 0) = 0 := by simp
    rw [hif]
    linarith

end EvenConnector
