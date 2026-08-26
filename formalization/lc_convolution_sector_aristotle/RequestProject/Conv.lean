import RequestProject.ZBasic
import RequestProject.CauchyBinet

/-!
# G0 : log-concavity with no internal zeros is closed under convolution

The main result of this file is `LC.logConcave_conv`, together with the
polynomial corollary `LC.logConcave_coeff_mul`.

The proof is the total-positivity route: the convolution
`Conv f g` is the `(·, 0)`-column of the product of the Toeplitz matrices of
`f` and `g`, and the relevant `2 × 2` minor of the product is expanded by the
`2 × 2` Cauchy–Binet identity `LC.cauchyBinet_two` into a sum of products of
`2 × 2` minors of the factors, each of which is nonnegative by
`LC.inner_mul_ge`.
-/

open scoped BigOperators

namespace LC

/-! ### Rewriting the convolution as a sum over an integer interval -/

lemma sum_Icc_shift (H : ℤ → ℝ) (A B c : ℤ) :
    ∑ p ∈ Finset.Icc A B, H (p - c) = ∑ p ∈ Finset.Icc (A - c) (B - c), H p := by
  rw [show Finset.Icc A B = Finset.map (addRightEmbedding c) (Finset.Icc (A - c) (B - c)) by
    rw [Finset.map_add_right_Icc]; ring_nf]
  rw [Finset.sum_map]
  exact Finset.sum_congr rfl (fun p _ => by simp [addRightEmbedding])

/-- The convolution, written as a sum over any integer interval containing
`[0, m]`, using the zero-extensions of `f` and `g` to `ℤ`. -/
lemma conv_eq_sum_Icc (f g : ℕ → ℝ) (m : ℕ) {A B : ℤ} (hA : A ≤ 0) (hB : (m : ℤ) ≤ B) :
    Conv f g m = ∑ p ∈ Finset.Icc A B, ZExt f p * ZExt g ((m : ℤ) - p) := by
  have h1 : ∑ p ∈ Finset.Icc A B, ZExt f p * ZExt g ((m : ℤ) - p)
      = ∑ p ∈ Finset.Icc (0 : ℤ) (m : ℤ), ZExt f p * ZExt g ((m : ℤ) - p) := by
    refine (Finset.sum_subset ?_ ?_).symm
    · intro p hp
      simp only [Finset.mem_Icc] at hp ⊢; omega
    · intro p hp hp'
      simp only [Finset.mem_Icc] at hp hp'
      rcases lt_or_ge p 0 with h | h
      · rw [ZExt_of_neg f h, zero_mul]
      · rw [ZExt_of_neg g (by omega), mul_zero]
  rw [h1]
  rw [show Finset.Icc (0:ℤ) (m:ℤ)
      = Finset.map ⟨(Nat.cast : ℕ → ℤ), Nat.cast_injective⟩ (Finset.range (m+1)) by
    ext q
    simp only [Finset.mem_Icc, Finset.mem_map, Finset.mem_range, Function.Embedding.coeFn_mk]
    constructor
    · rintro ⟨h1, h2⟩; exact ⟨q.toNat, by omega, by omega⟩
    · rintro ⟨i, hi, rfl⟩; omega]
  rw [Finset.sum_map]
  refine Finset.sum_congr rfl (fun i hi => ?_)
  simp only [Finset.mem_range] at hi
  simp only [Function.Embedding.coeFn_mk, ZExt_natCast]
  rw [show ((m:ℤ) - (i:ℤ)) = ((m - i : ℕ) : ℤ) by omega, ZExt_natCast]

/-! ### Log-concavity of the convolution -/

/-- Log-concavity of a convolution of two nonnegative log-concave sequences with
no internal zeros. -/
theorem logConcave_of_conv {f g : ℕ → ℝ}
    (hf0 : Nonneg f) (hg0 : Nonneg g) (hf : LogConcave f) (hg : LogConcave g)
    (hfz : NoInternalZeros f) (hgz : NoInternalZeros g) :
    LogConcave (Conv f g) := by
  intro k
  set F : ℤ → ℝ := ZExt f with hF
  set G : ℤ → ℝ := ZExt g with hG
  -- basic properties of the zero extensions
  have hFpos : ∀ t, 0 ≤ F t := ZExt_nonneg hf0
  have hFlc : ∀ t : ℤ, F (t + 1) * F (t + 1) ≥ F t * F (t + 2) := ZExt_logConcave hf0 hf
  have hFnz : ∀ i j t : ℤ, i ≤ j → j ≤ t → F i ≠ 0 → F t ≠ 0 → F j ≠ 0 :=
    ZExt_noInternalZeros hfz
  have hGpos : ∀ t, 0 ≤ G t := ZExt_nonneg hg0
  have hGlc : ∀ t : ℤ, G (t + 1) * G (t + 1) ≥ G t * G (t + 2) := ZExt_logConcave hg0 hg
  have hGnz : ∀ i j t : ℤ, i ≤ j → j ≤ t → G i ≠ 0 → G t ≠ 0 → G j ≠ 0 :=
    ZExt_noInternalZeros hgz
  -- the four families entering Cauchy–Binet
  set u : ℤ → ℝ := fun p => G ((k : ℤ) + 1 - p) with hu
  set v : ℤ → ℝ := fun p => G ((k : ℤ) + 2 - p) with hv
  set x : ℤ → ℝ := fun p => F p with hx
  set y : ℤ → ℝ := fun p => F (p - 1) with hy
  set s : Finset ℤ := Finset.Icc (0 : ℤ) ((k : ℤ) + 2) with hs
  -- the minors of the two factors are nonnegative
  have h1 : ∀ p ∈ s, ∀ q ∈ s, p < q → 0 ≤ u p * v q - u q * v p := by
    intro p _ q _ hpq
    have := inner_mul_ge (F := G) hGpos hGlc hGnz
      (a := (k : ℤ) + 1 - q) (b := (k : ℤ) + 2 - q) (c := (k : ℤ) + 1 - p)
      (d := (k : ℤ) + 2 - p) (by omega) (by omega) (by omega) (by omega)
    simp only [hu, hv]
    nlinarith [this]
  have h2 : ∀ p ∈ s, ∀ q ∈ s, p < q → 0 ≤ x p * y q - x q * y p := by
    intro p _ q _ hpq
    have := inner_mul_ge (F := F) hFpos hFlc hFnz
      (a := p - 1) (b := p) (c := q - 1) (d := q) (by omega) (by omega) (by omega) (by omega)
    simp only [hx, hy]
    nlinarith [this]
  have main := cauchyBinet_two_nonneg s u v x y h1 h2
  -- identify the four sums with values of the convolution
  have e1 : ∑ p ∈ s, u p * x p = Conv f g (k + 1) := by
    rw [conv_eq_sum_Icc f g (k + 1) (A := 0) (B := (k : ℤ) + 2) (by norm_num) (by push_cast; omega)]
    refine Finset.sum_congr rfl (fun p _ => ?_)
    simp only [hu, hx]
    rw [show ((k + 1 : ℕ) : ℤ) - p = (k : ℤ) + 1 - p by push_cast; ring]
    ring
  have e4 : ∑ q ∈ s, v q * x q = Conv f g (k + 2) := by
    rw [conv_eq_sum_Icc f g (k + 2) (A := 0) (B := (k : ℤ) + 2) (by norm_num) (by push_cast; omega)]
    refine Finset.sum_congr rfl (fun p _ => ?_)
    simp only [hv, hx]
    rw [show ((k + 2 : ℕ) : ℤ) - p = (k : ℤ) + 2 - p by push_cast; ring]
    ring
  have e2 : ∑ q ∈ s, v q * y q = Conv f g (k + 1) := by
    have hcong : ∀ q ∈ s, v q * y q = (fun p => F p * G ((k : ℤ) + 1 - p)) (q - 1) := by
      intro q _
      simp only [hv, hy]
      rw [show (k : ℤ) + 1 - (q - 1) = (k : ℤ) + 2 - q by ring]
      ring
    rw [Finset.sum_congr rfl hcong, hs,
      sum_Icc_shift (fun p => F p * G ((k : ℤ) + 1 - p)) 0 ((k : ℤ) + 2) 1]
    rw [conv_eq_sum_Icc f g (k + 1) (A := (0 : ℤ) - 1) (B := (k : ℤ) + 2 - 1)
      (by norm_num) (by push_cast; omega)]
    refine Finset.sum_congr rfl (fun p _ => ?_)
    rw [show ((k + 1 : ℕ) : ℤ) - p = (k : ℤ) + 1 - p by push_cast; ring]
  have e3 : ∑ p ∈ s, u p * y p = Conv f g k := by
    have hcong : ∀ p ∈ s, u p * y p = (fun t => F t * G ((k : ℤ) - t)) (p - 1) := by
      intro p _
      simp only [hu, hy]
      rw [show (k : ℤ) - (p - 1) = (k : ℤ) + 1 - p by ring]
      ring
    rw [Finset.sum_congr rfl hcong, hs,
      sum_Icc_shift (fun t => F t * G ((k : ℤ) - t)) 0 ((k : ℤ) + 2) 1]
    rw [conv_eq_sum_Icc f g k (A := (0 : ℤ) - 1) (B := (k : ℤ) + 2 - 1)
      (by norm_num) (by omega)]
  rw [e1, e2, e3, e4] at main
  exact main

/-! ### No internal zeros for the convolution -/

lemma conv_ne_zero_iff {f g : ℕ → ℝ} (hf0 : Nonneg f) (hg0 : Nonneg g) (m : ℕ) :
    Conv f g m ≠ 0 ↔ ∃ a b : ℕ, a + b = m ∧ f a ≠ 0 ∧ g b ≠ 0 := by
  constructor
  · intro h
    obtain ⟨i, hi, hne⟩ := Finset.exists_ne_zero_of_sum_ne_zero h
    simp only [Finset.mem_range] at hi
    exact ⟨i, m - i, by omega, fun h' => hne (by rw [h', zero_mul]),
      fun h' => hne (by rw [h', mul_zero])⟩
  · rintro ⟨a, b, rfl, ha, hb⟩
    have hterm : 0 < f a * g (a + b - a) := by
      have : a + b - a = b := by omega
      rw [this]
      exact mul_pos (lt_of_le_of_ne (hf0 a) (Ne.symm ha)) (lt_of_le_of_ne (hg0 b) (Ne.symm hb))
    have hle : f a * g (a + b - a) ≤ Conv f g (a + b) :=
      Finset.single_le_sum (f := fun i => f i * g (a + b - i))
        (fun i _ => mul_nonneg (hf0 i) (hg0 _)) (by simp only [Finset.mem_range]; omega)
    exact ne_of_gt (lt_of_lt_of_le hterm hle)

/-- The convolution of two nonnegative sequences with no internal zeros again has
no internal zeros. -/
theorem noInternalZeros_of_conv {f g : ℕ → ℝ}
    (hf0 : Nonneg f) (hg0 : Nonneg g)
    (hfz : NoInternalZeros f) (hgz : NoInternalZeros g) :
    NoInternalZeros (Conv f g) := by
  intro i j k hij hjk hi hk
  obtain ⟨a₁, b₁, hab₁, ha₁, hb₁⟩ := (conv_ne_zero_iff hf0 hg0 i).1 hi
  obtain ⟨a₂, b₂, hab₂, ha₂, hb₂⟩ := (conv_ne_zero_iff hf0 hg0 k).1 hk
  refine (conv_ne_zero_iff hf0 hg0 j).2 ?_
  set A₀ := min a₁ a₂
  set A₁ := max a₁ a₂
  set B₀ := min b₁ b₂
  set a : ℕ := min A₁ (j - B₀) with ha
  refine ⟨a, j - a, by simp only [ha]; omega, ?_, ?_⟩
  · rcases le_total a₁ a₂ with h | h
    · exact hfz a₁ a a₂ (by simp only [ha, A₁, B₀]; omega) (by simp only [ha, A₁, B₀]; omega)
        ha₁ ha₂
    · exact hfz a₂ a a₁ (by simp only [ha, A₁, B₀]; omega) (by simp only [ha, A₁, B₀]; omega)
        ha₂ ha₁
  · rcases le_total b₁ b₂ with h | h
    · exact hgz b₁ (j - a) b₂ (by simp only [ha, A₁, B₀]; omega)
        (by simp only [ha, A₁, B₀]; omega) hb₁ hb₂
    · exact hgz b₂ (j - a) b₁ (by simp only [ha, A₁, B₀]; omega)
        (by simp only [ha, A₁, B₀]; omega) hb₂ hb₁

lemma nonneg_of_conv {f g : ℕ → ℝ} (hf0 : Nonneg f) (hg0 : Nonneg g) :
    Nonneg (Conv f g) := fun _ =>
  Finset.sum_nonneg (fun i _ => mul_nonneg (hf0 i) (hg0 _))

/-- **G0.** Log-concavity together with the no-internal-zeros condition is
preserved by convolution.

The finiteness hypotheses `hfin`, `hgin` are those requested in the statement;
the proof does not need them (each value of `Conv f g` is a finite sum anyway). -/
theorem logConcave_conv
    (f g : ℕ → ℝ)
    (hf0 : Nonneg f) (hg0 : Nonneg g)
    (hf : LogConcave f) (hg : LogConcave g)
    (hfz : NoInternalZeros f) (hgz : NoInternalZeros g)
    (hfin : f.support.Finite) (hgin : g.support.Finite) :
    LogConcave (Conv f g) ∧ NoInternalZeros (Conv f g) :=
  ⟨logConcave_of_conv hf0 hg0 hf hg hfz hgz, noInternalZeros_of_conv hf0 hg0 hfz hgz⟩

/-! ### Polynomial form -/

lemma coeff_mul_eq_conv (p q : Polynomial ℝ) :
    (fun n => (p * q).coeff n) = Conv (fun n => p.coeff n) (fun n => q.coeff n) := by
  funext n
  rw [Polynomial.coeff_mul, Finset.Nat.sum_antidiagonal_eq_sum_range_succ_mk]
  rfl

/-- **G0, polynomial form.** If `p` and `q` have nonnegative, log-concave
coefficient sequences with no internal zeros, then so does `p * q`. -/
theorem logConcave_coeff_mul (p q : Polynomial ℝ)
    (hp0 : Nonneg (fun n => p.coeff n)) (hq0 : Nonneg (fun n => q.coeff n))
    (hp : LogConcave (fun n => p.coeff n)) (hq : LogConcave (fun n => q.coeff n))
    (hpz : NoInternalZeros (fun n => p.coeff n)) (hqz : NoInternalZeros (fun n => q.coeff n)) :
    Nonneg (fun n => (p * q).coeff n) ∧ LogConcave (fun n => (p * q).coeff n)
      ∧ NoInternalZeros (fun n => (p * q).coeff n) := by
  rw [coeff_mul_eq_conv]
  exact ⟨nonneg_of_conv hp0 hq0, logConcave_of_conv hp0 hq0 hp hq hpz hqz,
    noInternalZeros_of_conv hp0 hq0 hpz hqz⟩

end LC
