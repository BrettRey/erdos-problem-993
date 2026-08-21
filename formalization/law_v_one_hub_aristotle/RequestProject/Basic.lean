/-
Copyright (c) 2026. Released under Apache 2.0 license.

# Basic machinery: integer-indexed coefficients, likelihood-ratio order, log-concavity

This file sets up the elementary framework used throughout the development:

* `cf p i` : the coefficient of `p : ℤ[X]` at the *integer* index `i` (zero for `i < 0`).
  Using integer indices makes all shift bookkeeping (`X * p`) painless.
* `LR p q` : the likelihood-ratio order `p ≤lr q`, i.e. `p (k+1) * q k ≤ p k * q (k+1)`
  for every integer `k`.
* `LogConcave p` : `p (k+1) * p (k-1) ≤ p k * p k` for every integer `k`.
* `Nice p` : `p` has strictly positive coefficients up to some index `d` and vanishes
  after `d`.  This is the "nonnegative, no internal zeros" class in
  which all our polynomials live; it is closed under products.

Source label: agent (framework).
-/
import Mathlib

namespace LawV

open Polynomial

noncomputable section

/-- Coefficient of `p` at an integer index; zero at negative indices. -/
def cf (p : ℤ[X]) (i : ℤ) : ℤ := if 0 ≤ i then p.coeff i.toNat else 0

@[simp] lemma cf_of_neg (p : ℤ[X]) {i : ℤ} (h : i < 0) : cf p i = 0 := by
  simp [cf, not_le.2 h]

@[simp] lemma cf_ofNat (p : ℤ[X]) (n : ℕ) : cf p (n : ℤ) = p.coeff n := by
  simp [cf]

lemma cf_eq_coeff_toNat (p : ℤ[X]) {i : ℤ} (h : 0 ≤ i) : cf p i = p.coeff i.toNat := by
  simp [cf, h]

@[simp] lemma cf_zero (i : ℤ) : cf (0 : ℤ[X]) i = 0 := by
  unfold cf; split <;> simp

@[simp] lemma cf_add (p q : ℤ[X]) (i : ℤ) : cf (p + q) i = cf p i + cf q i := by
  unfold cf; split <;> simp

lemma cf_X_mul (p : ℤ[X]) (i : ℤ) : cf (X * p) i = cf p (i - 1) := by
  unfold cf
  by_cases h : 0 ≤ i
  · by_cases h1 : 0 ≤ i - 1
    · have h0 : i.toNat = (i-1).toNat + 1 := by omega
      rw [if_pos h, if_pos h1, h0, coeff_X_mul]
    · have hi : i = 0 := by omega
      subst hi
      norm_num
  · have h2 : ¬ (0 ≤ i - 1) := by omega
    simp only [h, h2, if_false]

/-- Likelihood-ratio (LR) order: `p ≤lr q`. -/
def LR (p q : ℤ[X]) : Prop := ∀ k : ℤ, cf p (k+1) * cf q k ≤ cf p k * cf q (k+1)

/-- Log-concavity of the coefficient sequence. -/
def LogConcave (p : ℤ[X]) : Prop := ∀ k : ℤ, cf p (k+1) * cf p (k-1) ≤ cf p k * cf p k

/-- `p` has positive coefficients up to index `d` and vanishes after `d`. -/
def NiceD (d : ℕ) (p : ℤ[X]) : Prop :=
  (∀ j : ℕ, j ≤ d → 0 < p.coeff j) ∧ (∀ j : ℕ, d < j → p.coeff j = 0)

/-- `p` has positive coefficients on an initial segment and vanishes afterwards. -/
def Nice (p : ℤ[X]) : Prop := ∃ d : ℕ, NiceD d p

lemma NiceD.nice {d : ℕ} {p : ℤ[X]} (h : NiceD d p) : Nice p := ⟨d, h⟩

lemma Nice.cf_nonneg {p : ℤ[X]} (hp : Nice p) (i : ℤ) : 0 ≤ cf p i := by
  obtain ⟨d, hpos, hzero⟩ := hp
  unfold cf
  split
  · rename_i h
    rcases le_or_gt i.toNat d with h1 | h1
    · exact (hpos _ h1).le
    · simp [hzero _ h1]
  · exact le_rfl

lemma Nice.coeff_nonneg {p : ℤ[X]} (hp : Nice p) (j : ℕ) : 0 ≤ p.coeff j := by
  simpa using hp.cf_nonneg (j : ℤ)

/-- Products of nice polynomials are nice, with degrees adding. -/
lemma NiceD.mul {d e : ℕ} {p q : ℤ[X]} (hp : NiceD d p) (hq : NiceD e q) :
    NiceD (d + e) (p * q) := by
  have hp0 := hp.nice.coeff_nonneg
  have hq0 := hq.nice.coeff_nonneg
  obtain ⟨hdpos, hdzero⟩ := hp
  obtain ⟨hepos, hezero⟩ := hq
  refine ⟨?_, ?_⟩
  · intro j hj
    rw [coeff_mul]
    have hmem : (min j d, j - min j d) ∈ Finset.antidiagonal j := by
      simp only [Finset.mem_antidiagonal]; omega
    have hpos : 0 < p.coeff (min j d) * q.coeff (j - min j d) := by
      have h1 : 0 < p.coeff (min j d) := hdpos _ (min_le_right _ _)
      have h2 : 0 < q.coeff (j - min j d) := hepos _ (by omega)
      positivity
    refine lt_of_lt_of_le hpos (Finset.single_le_sum
      (f := fun x : ℕ × ℕ => p.coeff x.1 * q.coeff x.2) (fun x _ => ?_) hmem)
    exact mul_nonneg (hp0 _) (hq0 _)
  · intro j hj
    rw [coeff_mul]
    refine Finset.sum_eq_zero ?_
    rintro ⟨a, b⟩ hab
    simp only [Finset.mem_antidiagonal] at hab
    rcases lt_or_ge d a with h | h
    · simp [hdzero _ h]
    · have hb : e < b := by omega
      simp [hezero _ hb]

lemma Nice.mul {p q : ℤ[X]} (hp : Nice p) (hq : Nice q) : Nice (p * q) := by
  obtain ⟨d, hd⟩ := hp; obtain ⟨e, he⟩ := hq; exact ⟨d + e, hd.mul he⟩

lemma NiceD.one : NiceD 0 (1 : ℤ[X]) := by
  refine ⟨?_, ?_⟩
  · intro j hj
    interval_cases j
    norm_num
  · intro j hj
    rw [coeff_one, if_neg (by omega : ¬ (j = 0))]

lemma Nice.one : Nice (1 : ℤ[X]) := NiceD.one.nice

lemma NiceD.cf_pos {d : ℕ} {p : ℤ[X]} (hp : NiceD d p) {i : ℤ} (h0 : 0 ≤ i) (hi : i ≤ d) :
    0 < cf p i := by
  rw [cf_eq_coeff_toNat p h0]
  exact hp.1 _ (by omega)

lemma NiceD.cf_eq_zero {d : ℕ} {p : ℤ[X]} (hp : NiceD d p) {i : ℤ} (hi : (d : ℤ) < i) :
    cf p i = 0 := by
  have h0 : 0 ≤ i := by omega
  rw [cf_eq_coeff_toNat p h0]
  exact hp.2 _ (by omega)

/-- `p + X * q` is again positive on an initial segment: the shifted block of `q` fills in
exactly where `p` stops. -/
lemma NiceD.add_X_mul {d e : ℕ} {p q : ℤ[X]} (hp : NiceD d p) (hq : NiceD e q) :
    NiceD (max d (e+1)) (p + X * q) := by
  have hqn := hq.nice.cf_nonneg
  have hpn := hp.nice.coeff_nonneg
  have key : ∀ j : ℕ, (p + X * q).coeff j = p.coeff j + cf q ((j : ℤ) - 1) := by
    intro j
    rw [coeff_add]
    congr 1
    simpa using cf_X_mul q (j : ℤ)
  constructor
  · intro j hj
    rw [key]
    rcases le_or_gt j d with h | h
    · have h1 := hp.1 j h
      have h2 := hqn ((j : ℤ) - 1)
      omega
    · have hj1 : 1 ≤ j := by omega
      have hje : (j : ℤ) - 1 ≤ (e : ℤ) := by
        have : j ≤ e + 1 := by omega
        omega
      have h1 : 0 < cf q ((j : ℤ) - 1) := by
        rw [cf_eq_coeff_toNat q (by omega)]
        exact hq.1 _ (by omega)
      have h2 := hpn j
      omega
  · intro j hj
    rw [key, hp.2 j (by omega), hq.cf_eq_zero (by omega)]
    ring

end

end LawV
