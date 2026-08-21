/-
# The likelihood-ratio order: closure properties

Main results of this file (all "layer 2" material of the request):

* `LR.add_right`, `LR.add_left` : the LR order is preserved by sums in each argument.
* `LR_X_mul_iff` : shifting both sides by `X` does not change the LR order.
* `logConcave_iff_LR` : `p` is log-concave iff `p ≤lr X * p`.
* `LR.gen` : for polynomials with no internal zeros, the adjacent LR inequalities imply
  all the `2 × 2` minors `p j * q i ≤ p i * q j` (`i ≤ j`).
* `LogConcave.gen` : likewise all "Toeplitz" minors of a log-concave sequence are nonnegative.
* `LR.trans'` : transitivity.
* `LR.mul_left` : **common convolution preserves the LR order** — a Cauchy–Binet argument.
* `LogConcave.mul` : products of log-concave polynomials (no internal zeros) are log-concave.

Source label: agent (layer 2).
-/
import RequestProject.Basic

namespace LawV

open Polynomial

noncomputable section

variable {p q r h u v : ℤ[X]}

/-! ### Elementary closure properties -/

lemma LR.add_right (h1 : LR p q) (h2 : LR p r) : LR p (q + r) := by
  intro k
  have := h1 k
  have := h2 k
  simp only [cf_add]
  nlinarith [h1 k, h2 k]

lemma LR.add_left (h1 : LR p r) (h2 : LR q r) : LR (p + q) r := by
  intro k
  simp only [cf_add]
  nlinarith [h1 k, h2 k]

lemma LR_X_mul (hpq : LR p q) : LR (X * p) (X * q) := by
  intro k
  simp only [cf_X_mul]
  have := hpq (k - 1)
  have e1 : k + 1 - 1 = k := by ring
  rw [e1]
  have e2 : k - 1 + 1 = k := by ring
  rw [e2] at this
  exact this

lemma LR_of_X_mul (hpq : LR (X * p) (X * q)) : LR p q := by
  intro k
  have := hpq (k + 1)
  simp only [cf_X_mul] at this
  have e1 : k + 1 + 1 - 1 = k + 1 := by ring
  have e2 : k + 1 - 1 = k := by ring
  rw [e1, e2] at this
  exact this

lemma LR_X_mul_iff : LR (X * p) (X * q) ↔ LR p q := ⟨LR_of_X_mul, LR_X_mul⟩

lemma logConcave_iff_LR : LogConcave p ↔ LR p (X * p) := by
  constructor
  · intro hp k
    simpa only [cf_X_mul, add_sub_cancel_right] using hp k
  · intro hp k
    simpa only [cf_X_mul, add_sub_cancel_right] using hp k

/-! ### Degrees -/

/-- If `p ≤lr q` then the degree of `p` is at most the degree of `q`. -/
lemma NiceD.deg_le_of_LR {d e : ℕ} (hp : NiceD d p) (hq : NiceD e q) (hpq : LR p q) : d ≤ e := by
  by_contra hlt
  push_neg at hlt
  have h1 := hpq (e : ℤ)
  have h2 : cf q ((e : ℤ) + 1) = 0 := hq.cf_eq_zero (by omega)
  have h3 : 0 < cf p ((e : ℤ) + 1) := hp.cf_pos (by positivity) (by omega)
  have h4 : 0 < cf q (e : ℤ) := hq.cf_pos (by positivity) (by omega)
  rw [h2, mul_zero] at h1
  nlinarith

/-! ### All minors -/

/-- All `2 × 2` minors: from the adjacent LR inequalities to the general ones. -/
lemma LR.gen {d e : ℕ} (hpq : LR p q) (hp : NiceD d p) (hq : NiceD e q) :
    ∀ i j : ℤ, i ≤ j → cf p j * cf q i ≤ cf p i * cf q j := by
  have hde : d ≤ e := hp.deg_le_of_LR hq hpq
  have hpn : ∀ i, 0 ≤ cf p i := hp.nice.cf_nonneg
  have hqn : ∀ i, 0 ≤ cf q i := hq.nice.cf_nonneg
  have key : ∀ n : ℕ, ∀ i : ℤ, cf p (i + n) * cf q i ≤ cf p i * cf q (i + n) := by
    intro n
    induction n with
    | zero => intro i; simp
    | succ n ih =>
      intro i
      push_cast
      have hstep := hpq (i + n)
      have hIH := ih i
      by_cases hz : 0 < cf q (i + n)
      · have e1 : (i + (n+1) : ℤ) = (i + n) + 1 := by ring
        rw [e1]
        -- p_{i+n+1} q_i * q_{i+n} ≤ p_i q_{i+n+1} * q_{i+n}
        have h1 : cf p (i + n + 1) * cf q i * cf q (i + n)
            ≤ cf p i * cf q (i + n + 1) * cf q (i + n) := by
          calc cf p (i + n + 1) * cf q i * cf q (i + n)
              = (cf p (i + n + 1) * cf q (i + n)) * cf q i := by ring
            _ ≤ (cf p (i + n) * cf q (i + n + 1)) * cf q i := by
                exact mul_le_mul_of_nonneg_right hstep (hqn i)
            _ = (cf p (i + n) * cf q i) * cf q (i + n + 1) := by ring
            _ ≤ (cf p i * cf q (i + n)) * cf q (i + n + 1) := by
                exact mul_le_mul_of_nonneg_right hIH (hqn _)
            _ = cf p i * cf q (i + n + 1) * cf q (i + n) := by ring
        exact le_of_mul_le_mul_right (by linarith [h1]) hz
      · -- `q` vanishes at `i+n`, hence beyond
        have hq0 : cf q (i + n) = 0 := le_antisymm (not_lt.1 hz) (hqn _)
        rcases lt_or_ge (i + n) 0 with hneg | hpos
        · have hi : i < 0 := by omega
          have : cf p i = 0 := cf_of_neg p hi
          have h2 : cf q i = 0 := cf_of_neg q hi
          simp [h2, this]
        · -- `i + n > e` since `q` is positive up to `e`
          have hgt : (e : ℤ) < i + n := by
            by_contra hle
            push_neg at hle
            exact absurd (hq.cf_pos hpos hle) hz
          have hp1 : cf p (i + (n+1) : ℤ) = 0 := by
            refine hp.cf_eq_zero ?_
            have : (d : ℤ) ≤ (e : ℤ) := by exact_mod_cast hde
            omega
          rw [hp1, zero_mul]
          exact mul_nonneg (hpn i) (hqn _)
  intro i j hij
  have hn : ∃ n : ℕ, j = i + n := ⟨(j - i).toNat, by omega⟩
  obtain ⟨n, rfl⟩ := hn
  exact key n i

/-- All Toeplitz minors of a log-concave sequence without internal zeros are nonnegative. -/
lemma LogConcave.gen {d : ℕ} (hh : LogConcave h) (hn : NiceD d h) :
    ∀ s t : ℤ, s ≤ t → cf h s * cf h (t + 1) ≤ cf h (s + 1) * cf h t := by
  have hnn : ∀ i, 0 ≤ cf h i := hn.nice.cf_nonneg
  have key : ∀ n : ℕ, ∀ s : ℤ, cf h s * cf h (s + n + 1) ≤ cf h (s + 1) * cf h (s + n) := by
    intro n
    induction n with
    | zero => intro s; simp [mul_comm]
    | succ n ih =>
      intro s
      push_cast
      have hIH := ih s
      have hlc := hh (s + n + 1)
      have e1 : (s + (n+1) + 1 : ℤ) = (s + n + 1) + 1 := by ring
      have e2 : (s + (n+1) : ℤ) = (s + n) + 1 := by ring
      rw [e1, e2]
      by_cases hz : 0 < cf h (s + n + 1)
      · have hlc' : cf h (s + n + 1 + 1) * cf h (s + n) ≤ cf h (s + n + 1) * cf h (s + n + 1) := by
          have e3 : (s + n + 1 - 1 : ℤ) = s + n := by ring
          rw [e3] at hlc
          exact hlc
        have hA : cf h s * cf h (s + n + 1) * cf h (s + n + 1 + 1)
            ≤ cf h (s + 1) * cf h (s + n) * cf h (s + n + 1 + 1) :=
          mul_le_mul_of_nonneg_right hIH (hnn _)
        have hB : cf h (s + 1) * (cf h (s + n + 1 + 1) * cf h (s + n))
            ≤ cf h (s + 1) * (cf h (s + n + 1) * cf h (s + n + 1)) :=
          mul_le_mul_of_nonneg_left hlc' (hnn _)
        have h1 : cf h s * cf h (s + n + 1 + 1) * cf h (s + n + 1)
            ≤ cf h (s + 1) * cf h (s + n + 1) * cf h (s + n + 1) := by nlinarith [hA, hB]
        exact le_of_mul_le_mul_right (by linarith [h1]) hz
      · have h0 : cf h (s + n + 1) = 0 := le_antisymm (not_lt.1 hz) (hnn _)
        rcases lt_or_ge s 0 with hneg | hs0
        · rw [cf_of_neg h hneg, zero_mul]
          exact mul_nonneg (hnn _) (hnn _)
        · have hpos : 0 ≤ s + (n : ℤ) + 1 := by positivity
          have hgt : (d : ℤ) < s + n + 1 := by
            by_contra hle
            push_neg at hle
            exact absurd (hn.cf_pos hpos hle) hz
          have h1 : cf h (s + n + 1 + 1) = 0 := hn.cf_eq_zero (by omega)
          rw [h1, mul_zero]
          exact mul_nonneg (hnn _) (hnn _)
  intro s t hst
  obtain ⟨n, rfl⟩ : ∃ n : ℕ, t = s + n := ⟨(t - s).toNat, by omega⟩
  exact key n s

/-! ### The general-minor form of the LR order -/

/-- All `2 × 2` minors of the two coefficient sequences are nonnegative.  For sequences
without internal zeros this is equivalent to `LR` (see `LR.gen` and `LRgen.LR`). -/
def LRgen (u v : ℤ[X]) : Prop := ∀ i j : ℤ, i ≤ j → cf u j * cf v i ≤ cf u i * cf v j

lemma LRgen.LR (h : LRgen u v) : LR u v := fun k => h k (k+1) (by omega)

lemma LRgen.X_mul (h : LRgen u v) : LRgen (X * u) (X * v) := by
  intro i j hij
  simpa only [cf_X_mul] using h (i-1) (j-1) (by omega)

lemma NiceD.LRgen_of_LR {d e : ℕ} (hpq : LR p q) (hp : NiceD d p) (hq : NiceD e q) :
    LRgen p q := LR.gen hpq hp hq

lemma LogConcave.LRgen {d : ℕ} (hh : LogConcave h) (hn : NiceD d h) : LRgen h (X * h) := by
  intro i j hij
  simp only [cf_X_mul]
  rw [mul_comm]
  simpa using hh.gen hn (i-1) (j-1) (by omega)

/-! ### Zero patterns forced by an LR comparison -/

/-- If `p ≤lr q` with `q` positive exactly on `[0,d]`, then `p` vanishes above `d`. -/
lemma LR.cf_eq_zero_of_gt {d : ℕ} (hpq : LR p q) (hp : Nice p) (hq : NiceD d q) {i : ℤ}
    (hi : (d : ℤ) < i) : cf p i = 0 := by
  obtain ⟨dp, hdp⟩ := hp
  have : dp ≤ d := hdp.deg_le_of_LR hq hpq
  exact hdp.cf_eq_zero (by omega)

/-- If `p ≤lr X * q` with `q` positive exactly on `[0,d]`, then `p` vanishes above `d+1`. -/
lemma LR.cf_eq_zero_of_gt_X {d : ℕ} (hpq : LR p (X * q)) (hp : Nice p) (hq : NiceD d q) {i : ℤ}
    (hi : (d : ℤ) + 1 < i) : cf p i = 0 := by
  obtain ⟨dp, hdp⟩ := hp
  have hnn := hdp.nice.cf_nonneg
  have hkey : cf p ((d : ℤ) + 2) = 0 := by
    have h1 := hpq ((d : ℤ) + 1)
    simp only [cf_X_mul] at h1
    have e1 : (d : ℤ) + 1 + 1 - 1 = (d : ℤ) + 1 := by ring
    have e2 : (d : ℤ) + 1 - 1 = (d : ℤ) := by ring
    rw [e1, e2] at h1
    have h2 : cf q ((d : ℤ) + 1) = 0 := hq.cf_eq_zero (by omega)
    have h3 : 0 < cf q (d : ℤ) := hq.cf_pos (by positivity) le_rfl
    rw [h2, mul_zero] at h1
    have h4 : cf p ((d : ℤ) + 1 + 1) * cf q (d : ℤ) ≤ 0 := h1
    have h5 : 0 ≤ cf p ((d : ℤ) + 1 + 1) := hnn _
    have : cf p ((d : ℤ) + 1 + 1) = 0 := by nlinarith
    simpa only [show (d : ℤ) + 1 + 1 = (d : ℤ) + 2 by ring] using this
  have hdd : dp ≤ d + 1 := by
    by_contra hlt
    push_neg at hlt
    have : 0 < cf p ((d : ℤ) + 2) :=
      hdp.cf_pos (by positivity) (by exact_mod_cast (by omega : (d : ℤ) + 2 ≤ (dp : ℤ)))
    omega
  exact hdp.cf_eq_zero (by omega)

/-! ### Transitivity -/

lemma LR.refl (p : ℤ[X]) : LR p p := fun _ => le_of_eq (mul_comm _ _)

/-- Transitivity of the LR order through a middle term with no internal zeros. -/
lemma LR.trans' {d : ℕ} (hp : Nice p) (hq : NiceD d q) (hr : ∀ i, 0 ≤ cf r i)
    (h1 : LR p q) (h2 : LR q r) : LR p r := by
  intro k
  have hpn := hp.cf_nonneg
  have hqn := hq.nice.cf_nonneg
  rcases lt_or_ge k 0 with hk | hk
  · rw [cf_of_neg r hk, mul_zero]
    exact mul_nonneg (hpn _) (hr _)
  · by_cases hz : 0 < cf q k ∧ 0 < cf q (k+1)
    · obtain ⟨hz1, hz2⟩ := hz
      have ha := h1 k
      have hb := h2 k
      have hstep : (cf p (k+1) * cf r k) * (cf q k * cf q (k+1))
          ≤ (cf p k * cf r (k+1)) * (cf q k * cf q (k+1)) := by
        have h3 : (cf p (k+1) * cf q k) * (cf q (k+1) * cf r k)
            ≤ (cf p k * cf q (k+1)) * (cf q k * cf r (k+1)) := by
          have h4 : 0 ≤ cf q (k+1) * cf r k := mul_nonneg (hqn _) (hr _)
          have h5 : 0 ≤ cf p k * cf q (k+1) := mul_nonneg (hpn _) (hqn _)
          calc (cf p (k+1) * cf q k) * (cf q (k+1) * cf r k)
              ≤ (cf p k * cf q (k+1)) * (cf q (k+1) * cf r k) :=
                mul_le_mul_of_nonneg_right ha h4
            _ ≤ (cf p k * cf q (k+1)) * (cf q k * cf r (k+1)) :=
                mul_le_mul_of_nonneg_left hb h5
        nlinarith [h3]
      have hpos : 0 < cf q k * cf q (k+1) := mul_pos hz1 hz2
      exact le_of_mul_le_mul_right hstep hpos
    · have hzero : cf p (k+1) = 0 := by
        rcases not_and_or.1 hz with h | h
        · have hq0 : cf q k = 0 := le_antisymm (not_lt.1 h) (hqn _)
          have : (d : ℤ) < k := by
            by_contra hle
            push_neg at hle
            exact absurd (hq.cf_pos hk hle) (by omega)
          exact h1.cf_eq_zero_of_gt hp hq (by omega)
        · have hq0 : cf q (k+1) = 0 := le_antisymm (not_lt.1 h) (hqn _)
          have : (d : ℤ) < k + 1 := by
            by_contra hle
            push_neg at hle
            exact absurd (hq.cf_pos (by omega) hle) (by omega)
          exact h1.cf_eq_zero_of_gt hp hq (by omega)
      rw [hzero, zero_mul]
      exact mul_nonneg (hpn _) (hr _)

/-- Shifted transitivity: `p ≤lr X*u` and `u ≤lr v` give `p ≤lr X*v`. -/
lemma LR.trans_X {d : ℕ} (hp : Nice p) (hu : NiceD d u) (hv : ∀ i, 0 ≤ cf v i)
    (h1 : LR p (X * u)) (h2 : LR u v) : LR p (X * v) := by
  intro k
  simp only [cf_X_mul, add_sub_cancel_right]
  have hpn := hp.cf_nonneg
  have hun := hu.nice.cf_nonneg
  have ha := h1 k
  simp only [cf_X_mul, add_sub_cancel_right] at ha
  have hb := h2 (k - 1)
  have e2 : k - 1 + 1 = k := by ring
  rw [e2] at hb
  rcases lt_or_ge k 1 with hk | hk
  · rw [cf_of_neg v (by omega : k - 1 < 0), mul_zero]
    exact mul_nonneg (hpn _) (hv _)
  · by_cases hz : 0 < cf u (k-1) ∧ 0 < cf u k
    · obtain ⟨hz1, hz2⟩ := hz
      have hstep : (cf p (k+1) * cf v (k-1)) * (cf u (k-1) * cf u k)
          ≤ (cf p k * cf v k) * (cf u (k-1) * cf u k) := by
        have h3 : (cf p (k+1) * cf u (k-1)) * (cf u k * cf v (k-1))
            ≤ (cf p k * cf u k) * (cf u (k-1) * cf v k) := by
          have h4 : 0 ≤ cf u k * cf v (k-1) := mul_nonneg (hun _) (hv _)
          have h5 : 0 ≤ cf p k * cf u k := mul_nonneg (hpn _) (hun _)
          calc (cf p (k+1) * cf u (k-1)) * (cf u k * cf v (k-1))
              ≤ (cf p k * cf u k) * (cf u k * cf v (k-1)) :=
                mul_le_mul_of_nonneg_right ha h4
            _ ≤ (cf p k * cf u k) * (cf u (k-1) * cf v k) :=
                mul_le_mul_of_nonneg_left hb h5
        nlinarith [h3]
      exact le_of_mul_le_mul_right hstep (mul_pos hz1 hz2)
    · have hzero : cf p (k+1) = 0 := by
        rcases not_and_or.1 hz with h | h
        · have hu0 : cf u (k-1) = 0 := le_antisymm (not_lt.1 h) (hun _)
          have hgt : (d : ℤ) < k - 1 := by
            by_contra hle
            push_neg at hle
            exact absurd (hu.cf_pos (by omega) hle) (by omega)
          exact h1.cf_eq_zero_of_gt_X hp hu (by omega)
        · have hu0 : cf u k = 0 := le_antisymm (not_lt.1 h) (hun _)
          have hgt : (d : ℤ) < k := by
            by_contra hle
            push_neg at hle
            exact absurd (hu.cf_pos (by omega) hle) (by omega)
          exact h1.cf_eq_zero_of_gt_X hp hu (by omega)
      rw [hzero, zero_mul]
      exact mul_nonneg (hpn _) (hv _)

/-! ### Common convolution preserves the LR order (Cauchy–Binet) -/

/-- The convolution formula for `cf` of a product, with the summation range fixed by a
support bound on the second factor. -/
lemma cf_mul_eq_sum (h u : ℤ[X]) {N : ℕ} (hN : ∀ j : ℕ, N < j → u.coeff j = 0) (k : ℤ) :
    cf (h * u) k = ∑ i ∈ Finset.range (N+1), cf u (i : ℤ) * cf h (k - i) := by
  rcases lt_or_ge k 0 with hk | hk
  · rw [cf_of_neg _ hk]
    refine (Finset.sum_eq_zero ?_).symm
    intro i _
    rw [cf_of_neg h (by omega : k - (i:ℤ) < 0), mul_zero]
  · set n := k.toNat with hn
    have hkn : k = (n : ℤ) := by omega
    have hbase : cf (h * u) k = ∑ i ∈ Finset.range (n+1), cf u (i : ℤ) * cf h (k - i) := by
      rw [hkn, cf_ofNat, mul_comm, coeff_mul, Finset.Nat.sum_antidiagonal_eq_sum_range_succ_mk]
      refine Finset.sum_congr rfl ?_
      intro i hi
      simp only [Finset.mem_range] at hi
      rw [cf_ofNat]
      congr 1
      rw [cf_eq_coeff_toNat h (by omega)]
      congr 1
      omega
    set M := max N n with hM
    have h1 : ∑ i ∈ Finset.range (N+1), cf u (i : ℤ) * cf h (k - i)
        = ∑ i ∈ Finset.range (M+1), cf u (i : ℤ) * cf h (k - i) := by
      refine Finset.sum_subset (by intro x hx; simp only [Finset.mem_range] at *; omega) ?_
      intro x hx hx'
      simp only [Finset.mem_range] at hx hx'
      rw [cf_ofNat, hN x (by omega), zero_mul]
    have h2 : ∑ i ∈ Finset.range (n+1), cf u (i : ℤ) * cf h (k - i)
        = ∑ i ∈ Finset.range (M+1), cf u (i : ℤ) * cf h (k - i) := by
      refine Finset.sum_subset (by intro x hx; simp only [Finset.mem_range] at *; omega) ?_
      intro x hx hx'
      simp only [Finset.mem_range] at hx hx'
      rw [cf_of_neg h (by omega : k - (x:ℤ) < 0), mul_zero]
    rw [hbase, h2, h1]

/-- **Common convolution preserves the likelihood-ratio order.**  If `h` is log-concave with
no internal zeros, and all `2 × 2` minors of `(u,v)` are nonnegative, then `h*u ≤lr h*v`.
This is a Cauchy–Binet / TP2 argument. -/
lemma LR.mul_left {dh N : ℕ} (hh : NiceD dh h) (hlc : LogConcave h)
    (hu : ∀ j : ℕ, N < j → u.coeff j = 0) (hv : ∀ j : ℕ, N < j → v.coeff j = 0)
    (huv : LRgen u v) : LR (h * u) (h * v) := by
  intro k
  -- the four convolution expansions
  have eA := cf_mul_eq_sum h u hu (k+1)
  have eB := cf_mul_eq_sum h v hv k
  have eC := cf_mul_eq_sum h u hu k
  have eD := cf_mul_eq_sum h v hv (k+1)
  set T : ℕ → ℕ → ℤ := fun i j =>
    (cf u i * cf v j) * (cf h (k - i) * cf h (k + 1 - j) - cf h (k + 1 - i) * cf h (k - j))
    with hT
  have hsum : cf (h*u) k * cf (h*v) (k+1) - cf (h*u) (k+1) * cf (h*v) k
      = ∑ i ∈ Finset.range (N+1), ∑ j ∈ Finset.range (N+1), T i j := by
    rw [eA, eB, eC, eD, Finset.sum_mul_sum, Finset.sum_mul_sum, ← Finset.sum_sub_distrib]
    refine Finset.sum_congr rfl ?_
    intro i _
    rw [← Finset.sum_sub_distrib]
    refine Finset.sum_congr rfl ?_
    intro j _
    simp only [hT]
    ring
  have hswap : ∑ i ∈ Finset.range (N+1), ∑ j ∈ Finset.range (N+1), T i j
      = ∑ i ∈ Finset.range (N+1), ∑ j ∈ Finset.range (N+1), T j i := Finset.sum_comm
  have hnonneg : ∀ i ∈ Finset.range (N+1), ∀ j ∈ Finset.range (N+1), 0 ≤ T i j + T j i := by
    intro i _ j _
    have hD : T i j + T j i =
        (cf u i * cf v j - cf u j * cf v i) *
          (cf h (k - i) * cf h (k + 1 - j) - cf h (k + 1 - i) * cf h (k - j)) := by
      simp only [hT]; ring
    rw [hD]
    rcases le_total (i : ℤ) (j : ℤ) with hij | hij
    · have h1 : 0 ≤ cf u i * cf v j - cf u j * cf v i := by
        have := huv i j hij; omega
      have h2 : 0 ≤ cf h (k - i) * cf h (k + 1 - j) - cf h (k + 1 - i) * cf h (k - j) := by
        have := hlc.gen hh (k - j) (k - i) (by omega)
        have e1 : k - j + 1 = k + 1 - j := by ring
        have e2 : k - i + 1 = k + 1 - i := by ring
        rw [e1, e2] at this
        nlinarith [this]
      positivity
    · have h1 : cf u i * cf v j - cf u j * cf v i ≤ 0 := by
        have := huv j i hij; omega
      have h2 : cf h (k - i) * cf h (k + 1 - j) - cf h (k + 1 - i) * cf h (k - j) ≤ 0 := by
        have := hlc.gen hh (k - i) (k - j) (by omega)
        have e1 : k - j + 1 = k + 1 - j := by ring
        have e2 : k - i + 1 = k + 1 - i := by ring
        rw [e1, e2] at this
        nlinarith [this]
      nlinarith [mul_nonneg (neg_nonneg.2 h1) (neg_nonneg.2 h2)]
  have hpos : 0 ≤ ∑ i ∈ Finset.range (N+1), ∑ j ∈ Finset.range (N+1), (T i j + T j i) := by
    refine Finset.sum_nonneg ?_
    intro i hi
    refine Finset.sum_nonneg ?_
    intro j hj
    exact hnonneg i hi j hj
  have hsplit : ∑ i ∈ Finset.range (N+1), ∑ j ∈ Finset.range (N+1), (T i j + T j i)
      = (∑ i ∈ Finset.range (N+1), ∑ j ∈ Finset.range (N+1), T i j)
        + ∑ i ∈ Finset.range (N+1), ∑ j ∈ Finset.range (N+1), T j i := by
    rw [← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl ?_
    intro i _
    rw [← Finset.sum_add_distrib]
  rw [hsplit, ← hswap] at hpos
  linarith [hsum, hpos]

/-- Products of log-concave polynomials without internal zeros are log-concave. -/
lemma LogConcave.mul {dp dq : ℕ} (hp : NiceD dp p) (hlcp : LogConcave p)
    (hq : NiceD dq q) (hlcq : LogConcave q) : LogConcave (p * q) := by
  rw [logConcave_iff_LR]
  have hXq : X * (p * q) = p * (X * q) := by ring
  rw [hXq]
  refine LR.mul_left (N := dq + 1) hp hlcp ?_ ?_ (hlcq.LRgen hq)
  · intro j hj; exact hq.2 j (by omega)
  · intro j hj
    have h1 : cf (X * q) (j : ℤ) = cf q ((j : ℤ) - 1) := cf_X_mul q _
    have h2 : cf q ((j : ℤ) - 1) = 0 := hq.cf_eq_zero (by omega)
    simpa using h1.trans h2

end

end LawV
