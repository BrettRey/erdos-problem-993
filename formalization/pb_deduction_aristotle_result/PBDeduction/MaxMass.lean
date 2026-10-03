import PBDeduction.Defs

/-!
# The maximal-mass bound (G2)

`Var ≥ (M⁻² - 1)/12` for every mass function on `{0, …, n}` (Bobkov,
Marsiglietti and Melbourne, Corollary 3.2).

We give a purely discrete version of the paper's proof. With `μ` the mean and
`r = 1/(2M)`, put `Φ(y) = r² y - y³/3` (an antiderivative of `(r² - y²)₊` on
`[-r, r]`) and `Φc(y) = Φ(clip_{[-r,r]} y)`. For every real `y`,
`T(y) := Φc(y + 1/2) - Φc(y - 1/2)` satisfies `T(y) ≥ 0` and
`T(y) ≥ r² - y² - 1/12` (the latter is the midpoint form of
`∫_{y-1/2}^{y+1/2} (r² - x²) dx = r² - y² - 1/12`). Since `0 ≤ f k ≤ M`,
`f k ((k-μ)² + 1/12 - r²) ≥ -M T(k - μ)`; summing and telescoping
`∑ T(k-μ) ≤ Φ(r) - Φ(-r) = 4r³/3` gives `V + 1/12 - r² ≥ -4Mr³/3`, which for
`r = 1/(2M)` is exactly `V ≥ (M⁻² - 1)/12`.
-/

namespace PBDeduction

open Finset

noncomputable section

/-- Clipping to `[-r, r]`. -/
def clipR (r y : ℝ) : ℝ := max (-r) (min r y)

/-- `Φ(y) = r² y - y³/3`. -/
def PhiR (r y : ℝ) : ℝ := r ^ 2 * y - y ^ 3 / 3

/-- `Φc(y) = Φ(clip y)`. -/
def PhiC (r y : ℝ) : ℝ := PhiR r (clipR r y)

lemma PhiR_mono {r a b : ℝ} (ha : -r ≤ a) (hab : a ≤ b) (hb : b ≤ r) :
    PhiR r a ≤ PhiR r b := by
  unfold PhiR
  have h1 : a ^ 2 + a * b + b ^ 2 ≤ 3 * r ^ 2 := by nlinarith
  have : r ^ 2 * b - b ^ 3 / 3 - (r ^ 2 * a - a ^ 3 / 3) =
      (b - a) * (3 * r ^ 2 - (a ^ 2 + a * b + b ^ 2)) / 3 := by ring
  have : 0 ≤ (b - a) * (3 * r ^ 2 - (a ^ 2 + a * b + b ^ 2)) :=
    mul_nonneg (by linarith) (by linarith)
  linarith

lemma clipR_mem {r : ℝ} (hr : 0 ≤ r) (y : ℝ) : -r ≤ clipR r y ∧ clipR r y ≤ r := by
  unfold clipR
  constructor
  · exact le_max_left _ _
  · exact max_le (by linarith) (min_le_left _ _)

lemma clipR_mono {r a b : ℝ} (hab : a ≤ b) : clipR r a ≤ clipR r b := by
  unfold clipR
  exact max_le_max le_rfl (min_le_min le_rfl hab)

lemma PhiC_mono {r a b : ℝ} (hr : 0 ≤ r) (hab : a ≤ b) : PhiC r a ≤ PhiC r b := by
  unfold PhiC
  exact PhiR_mono (clipR_mem hr a).1 (clipR_mono hab) (clipR_mem hr b).2

lemma PhiC_le {r : ℝ} (hr : 0 ≤ r) (y : ℝ) : PhiC r y ≤ 2 * r ^ 3 / 3 := by
  have := PhiR_mono (clipR_mem hr y).1 (clipR_mem hr y).2 le_rfl
  unfold PhiC
  unfold PhiR at this ⊢
  linarith

lemma PhiC_ge {r : ℝ} (hr : 0 ≤ r) (y : ℝ) : -(2 * r ^ 3 / 3) ≤ PhiC r y := by
  have := PhiR_mono (le_refl (-r)) (clipR_mem hr y).1 (clipR_mem hr y).2
  unfold PhiC
  unfold PhiR at this ⊢
  linarith

/-- The midpoint comparison `T(y) ≥ r² - y² - 1/12`. -/
lemma T_ge {r : ℝ} (hr : 0 ≤ r) (y : ℝ) :
    r ^ 2 - y ^ 2 - 1 / 12 ≤ PhiC r (y + 1 / 2) - PhiC r (y - 1 / 2) := by
  have hT0 : 0 ≤ PhiC r (y + 1 / 2) - PhiC r (y - 1 / 2) := by
    have := PhiC_mono hr (show y - 1 / 2 ≤ y + 1 / 2 by linarith); linarith
  by_cases hy : r ^ 2 ≤ y ^ 2
  · linarith
  push_neg at hy
  have hyr : |y| < r := by
    rw [← abs_of_nonneg hr]; exact sq_lt_sq.mp hy
  rw [abs_lt] at hyr
  obtain ⟨hy1, hy2⟩ := hyr
  have hb : PhiR r (y + 1 / 2) ≤ PhiC r (y + 1 / 2) := by
    unfold PhiC clipR
    rw [max_eq_right (by apply le_min <;> linarith)]
    rcases le_total (y + 1 / 2) r with h | h
    · rw [min_eq_right h]
    · rw [min_eq_left h]
      unfold PhiR
      have : r ^ 2 * r - r ^ 3 / 3 - (r ^ 2 * (y + 1 / 2) - (y + 1 / 2) ^ 3 / 3) =
          (r - (y + 1 / 2)) ^ 2 * (2 * r + (y + 1 / 2)) / 3 := by ring
      have : 0 ≤ (r - (y + 1 / 2)) ^ 2 * (2 * r + (y + 1 / 2)) :=
        mul_nonneg (sq_nonneg _) (by linarith)
      linarith
  have ha : PhiC r (y - 1 / 2) ≤ PhiR r (y - 1 / 2) := by
    unfold PhiC clipR
    rw [min_eq_right (by linarith)]
    rcases le_total (-r) (y - 1 / 2) with h | h
    · rw [max_eq_right h]
    · rw [max_eq_left h]
      unfold PhiR
      have : r ^ 2 * (y - 1 / 2) - (y - 1 / 2) ^ 3 / 3 - (r ^ 2 * (-r) - (-r) ^ 3 / 3) =
          ((y - 1 / 2) + r) ^ 2 * (2 * r - (y - 1 / 2)) / 3 := by ring
      have : 0 ≤ ((y - 1 / 2) + r) ^ 2 * (2 * r - (y - 1 / 2)) :=
        mul_nonneg (sq_nonneg _) (by linarith)
      linarith
  have : PhiR r (y + 1 / 2) - PhiR r (y - 1 / 2) = r ^ 2 - y ^ 2 - 1 / 12 := by
    unfold PhiR; ring
  linarith

/-- Telescoping over `Icc 0 m` in `ℤ`. -/
lemma sum_Icc_telescope (F : ℤ → ℝ) (m : ℕ) :
    ∑ k ∈ Icc (0 : ℤ) m, (F (k + 1) - F k) = F (m + 1) - F 0 := by
  induction m with
  | zero => simp
  | succ m ih =>
    have : Icc (0 : ℤ) ((m + 1 : ℕ) : ℤ) = insert ((m : ℤ) + 1) (Icc (0 : ℤ) m) := by
      ext x; simp only [mem_Icc, mem_insert]; push_cast; omega
    rw [this, sum_insert (by simp), ih]
    push_cast; ring

/-- The pairwise variance equals `∑ f k (k - μ)²` when `∑ f = 1`. -/
lemma pairVar_eq_central (f : ℤ → ℝ) (n : ℕ)
    (hsum : ∑ k ∈ Finset.Icc (0 : ℤ) n, f k = 1) :
    pairVar f n = ∑ k ∈ Icc (0 : ℤ) n,
      f k * ((k : ℝ) - ∑ j ∈ Icc (0 : ℤ) n, f j * (j : ℝ)) ^ 2 := by
  set μ := ∑ j ∈ Icc (0 : ℤ) n, f j * (j : ℝ) with hμ
  have h0 : ∑ k ∈ Icc (0 : ℤ) n, f k * ((k : ℝ) - μ) = 0 := by
    simp only [mul_sub, sum_sub_distrib, ← sum_mul, hsum]; rw [← hμ]; ring
  unfold pairVar
  have this : ∀ i j : ℤ, f i * f j * ((i : ℝ) - (j : ℝ)) ^ 2 =
      f i * ((i : ℝ) - μ) ^ 2 * f j + f j * ((j : ℝ) - μ) ^ 2 * f i
        - 2 * (f i * ((i : ℝ) - μ)) * (f j * ((j : ℝ) - μ)) := by
    intro i j; ring
  have key : ∑ i ∈ Icc (0 : ℤ) n, ∑ j ∈ Icc (0 : ℤ) n, f i * f j * ((i : ℝ) - (j : ℝ)) ^ 2
      = 2 * ∑ k ∈ Icc (0 : ℤ) n, f k * ((k : ℝ) - μ) ^ 2 := by
    calc _ = ∑ i ∈ Icc (0 : ℤ) n, ∑ j ∈ Icc (0 : ℤ) n,
          (f i * ((i : ℝ) - μ) ^ 2 * f j + f j * ((j : ℝ) - μ) ^ 2 * f i
            - 2 * (f i * ((i : ℝ) - μ)) * (f j * ((j : ℝ) - μ))) := by
          simp only [this]
      _ = ∑ i ∈ Icc (0 : ℤ) n,
          (f i * ((i : ℝ) - μ) ^ 2 * ∑ j ∈ Icc (0 : ℤ) n, f j
            + (∑ j ∈ Icc (0 : ℤ) n, f j * ((j : ℝ) - μ) ^ 2) * f i
            - 2 * (f i * ((i : ℝ) - μ)) * ∑ j ∈ Icc (0 : ℤ) n, f j * ((j : ℝ) - μ)) := by
          apply sum_congr rfl; intro i _
          rw [sum_sub_distrib, sum_add_distrib, ← mul_sum, ← sum_mul, ← mul_sum]
      _ = ∑ i ∈ Icc (0 : ℤ) n,
          (f i * ((i : ℝ) - μ) ^ 2 + (∑ j ∈ Icc (0 : ℤ) n, f j * ((j : ℝ) - μ) ^ 2) * f i) := by
          rw [hsum, h0]; simp
      _ = _ := by rw [sum_add_distrib, ← mul_sum, hsum]; ring
  rw [key]; ring

/-- **G2** (the maximal-mass bound). -/
theorem max_mass_bound_proof (f : ℤ → ℝ) (n : ℕ) (hnn : ∀ k, 0 ≤ f k)
    (hsum : ∑ k ∈ Finset.Icc (0 : ℤ) n, f k = 1) :
    (1 / maxMass f n ^ 2 - 1) / 12 ≤ pairVar f n := by
  set M := maxMass f n with hMdef
  set S := Icc (0 : ℤ) n
  have hfM : ∀ k ∈ S, f k ≤ M := fun k hk => Finset.le_sup' f hk
  have hM : 0 < M := by
    by_contra hneg
    push_neg at hneg
    have : ∑ k ∈ S, f k ≤ 0 := sum_nonpos (fun k hk => (hfM k hk).trans hneg)
    linarith
  set μ := ∑ j ∈ S, f j * (j : ℝ)
  set r := 1 / (2 * M) with hr
  have hr0 : 0 ≤ r := by positivity
  have hpv := pairVar_eq_central f n hsum
  -- termwise bound
  have hterm : ∀ k ∈ S, f k * (((k : ℝ) - μ) ^ 2 + 1 / 12 - r ^ 2) ≥
      -(M * (PhiC r ((k : ℝ) - μ + 1 / 2) - PhiC r ((k : ℝ) - μ - 1 / 2))) := by
    intro k hk
    have hT := T_ge hr0 ((k : ℝ) - μ)
    have hT0 : 0 ≤ PhiC r ((k : ℝ) - μ + 1 / 2) - PhiC r ((k : ℝ) - μ - 1 / 2) := by
      have := PhiC_mono hr0 (show (k : ℝ) - μ - 1 / 2 ≤ (k : ℝ) - μ + 1 / 2 by linarith)
      linarith
    have h1 := hfM k hk
    have h2 := hnn k
    rcases le_total 0 (((k : ℝ) - μ) ^ 2 + 1 / 12 - r ^ 2) with hg | hg
    · have : 0 ≤ f k * (((k : ℝ) - μ) ^ 2 + 1 / 12 - r ^ 2) := mul_nonneg h2 hg
      nlinarith
    · have : M * (((k : ℝ) - μ) ^ 2 + 1 / 12 - r ^ 2) ≤
          f k * (((k : ℝ) - μ) ^ 2 + 1 / 12 - r ^ 2) :=
        mul_le_mul_of_nonpos_right h1 hg
      nlinarith
  have hsumT : ∑ k ∈ S, (PhiC r ((k : ℝ) - μ + 1 / 2) - PhiC r ((k : ℝ) - μ - 1 / 2)) ≤
      4 * r ^ 3 / 3 := by
    have := sum_Icc_telescope (fun k : ℤ => PhiC r ((k : ℝ) - μ - 1 / 2)) n
    simp only at this
    have e : ∀ k : ℤ, PhiC r (((k + 1 : ℤ) : ℝ) - μ - 1 / 2) = PhiC r ((k : ℝ) - μ + 1 / 2) := by
      intro k; push_cast; ring_nf
    simp only [e] at this
    rw [this]
    have := PhiC_le hr0 (((n : ℤ) : ℝ) - μ + 1 / 2)
    have := PhiC_ge hr0 (((0 : ℤ) : ℝ) - μ - 1 / 2)
    linarith
  have hmain : ∑ k ∈ S, f k * (((k : ℝ) - μ) ^ 2 + 1 / 12 - r ^ 2) ≥ -(M * (4 * r ^ 3 / 3)) := by
    calc ∑ k ∈ S, f k * (((k : ℝ) - μ) ^ 2 + 1 / 12 - r ^ 2)
        ≥ ∑ k ∈ S, -(M * (PhiC r ((k : ℝ) - μ + 1 / 2) - PhiC r ((k : ℝ) - μ - 1 / 2))) :=
          sum_le_sum hterm
      _ = -(M * ∑ k ∈ S, (PhiC r ((k : ℝ) - μ + 1 / 2) - PhiC r ((k : ℝ) - μ - 1 / 2))) := by
          rw [sum_neg_distrib, mul_sum]
      _ ≥ -(M * (4 * r ^ 3 / 3)) := by
          have := mul_le_mul_of_nonneg_left hsumT hM.le; linarith
  have hexp : ∑ k ∈ S, f k * (((k : ℝ) - μ) ^ 2 + 1 / 12 - r ^ 2) =
      pairVar f n + (1 / 12 - r ^ 2) := by
    rw [hpv]
    simp only [mul_add, mul_sub, sum_add_distrib, sum_sub_distrib, ← sum_mul, hsum]
    ring
  rw [hexp] at hmain
  have hr2 : r ^ 2 = 1 / (4 * M ^ 2) := by rw [hr]; field_simp; ring
  have hr3 : M * (4 * r ^ 3 / 3) = 1 / (6 * M ^ 2) := by rw [hr]; field_simp; ring
  rw [hr2, hr3] at hmain
  have : (1 / M ^ 2 - 1) / 12 = 1 / (4 * M ^ 2) - 1 / (6 * M ^ 2) - 1 / 12 := by
    field_simp; ring
  rw [this]
  linarith

end

end PBDeduction
