import PBDeduction.PBLogConcave

/-!
# The Hillion–Johnson cubic inequalities (G4)

`C1(k) = g_{k-1} g_k² - 2 g_{k-1}² g_{k+1} + g_k g_{k+1} g_{k-2} ≥ 0`
(Hillion–Johnson, Theorem A.2, their (78)), by induction on the number of
Bernoulli summands, expanding in the cubic Bernstein basis in the new
parameter `q`:

`C1^{conv q g}(k) = (1-q)³ C1(k) + (1-q)² q X + (1-q) q² Y + q³ C1(k-1)`

with (writing `a, b, c, d, e = g_{k-3}, …, g_{k+1}`)

* `X = a d e - 3 b c e + 2 b d²` (their `D1(k)`), and `c X = 2 b C1(k) + e C1(k-1)`;
* `Y = a c e + a d² - 2 b² e - b c d + c³` (their Proposition A.4), and
  `c² d Y = b D_{k-1} C1(k) + e D_{k-1} C1(k-1) + D_{k-1}² (c d - b e)
            + a c d C1(k) + 2 c² e C1(k-1)`, where `D_{k-1} = c² - b d`.

The boundary cases `c = 0` or `d = 0` are handled by the support structure.
The mirrored inequality (79) (Corollary A.3) follows by the duality
`k ↦ n - k`, `p ↦ 1 - p`.
-/

namespace PBDeduction

open Finset

noncomputable section

/-- Hillion–Johnson's `C1(k)`, their (78). -/
def hjC1 (g : ℤ → ℝ) (k : ℤ) : ℝ :=
  g (k - 1) * g k ^ 2 - 2 * g (k - 1) ^ 2 * g (k + 1) + g k * g (k + 1) * g (k - 2)

/-- The cubic Bernstein expansion of `C1` for `conv q g`. -/
lemma hjC1_conv (q : ℝ) (g : ℤ → ℝ) (k : ℤ) :
    hjC1 (conv q g) k =
      (1 - q) ^ 3 * (g (k - 2) * g k * g (k + 1) - 2 * g (k - 1) ^ 2 * g (k + 1) +
          g (k - 1) * g k ^ 2) +
      (1 - q) ^ 2 * q * (g (k - 3) * g k * g (k + 1) - 3 * g (k - 2) * g (k - 1) * g (k + 1) +
          2 * g (k - 2) * g k ^ 2) +
      (1 - q) * q ^ 2 * (g (k - 3) * g (k - 1) * g (k + 1) + g (k - 3) * g k ^ 2 -
          2 * g (k - 2) ^ 2 * g (k + 1) - g (k - 2) * g (k - 1) * g k + g (k - 1) ^ 3) +
      q ^ 3 * (g (k - 3) * g (k - 1) * g k - 2 * g (k - 2) ^ 2 * g k +
          g (k - 2) * g (k - 1) ^ 2) := by
  unfold hjC1 conv
  rw [show k - 1 - 1 = k - 2 by ring, show k + 1 - 1 = k by ring,
    show k - 2 - 1 = k - 3 by ring]
  ring

/-- Nonnegativity of the two middle Bernstein coefficients. -/
lemma hj_middle_nonneg {a b c d e : ℝ} (ha : 0 ≤ a) (hb : 0 ≤ b) (hc : 0 ≤ c) (hd : 0 ≤ d)
    (he : 0 ≤ e)
    (hC1k : 0 ≤ b * d * e - 2 * c ^ 2 * e + c * d ^ 2)
    (hC1km1 : 0 ≤ a * c * d - 2 * b ^ 2 * d + b * c ^ 2)
    (hD : 0 ≤ c ^ 2 - b * d) (hX : 0 ≤ c * d - b * e)
    (hbd : (0 < c ∧ 0 < d) ∨ (a = 0 ∧ b = 0 ∧ c = 0) ∨ (d = 0 ∧ e = 0)) :
    0 ≤ a * d * e - 3 * b * c * e + 2 * b * d ^ 2 ∧
    0 ≤ a * c * e + a * d ^ 2 - 2 * b ^ 2 * e - b * c * d + c ^ 3 := by
  rcases hbd with ⟨hc0, hd0⟩ | ⟨rfl, rfl, rfl⟩ | ⟨rfl, rfl⟩
  · constructor
    · have key : c * (a * d * e - 3 * b * c * e + 2 * b * d ^ 2) =
          2 * b * (b * d * e - 2 * c ^ 2 * e + c * d ^ 2) +
            e * (a * c * d - 2 * b ^ 2 * d + b * c ^ 2) := by ring
      have : 0 ≤ c * (a * d * e - 3 * b * c * e + 2 * b * d ^ 2) := by
        rw [key]; positivity
      exact (mul_nonneg_iff_of_pos_left hc0).mp this
    · have key : c ^ 2 * d * (a * c * e + a * d ^ 2 - 2 * b ^ 2 * e - b * c * d + c ^ 3) =
          b * (c ^ 2 - b * d) * (b * d * e - 2 * c ^ 2 * e + c * d ^ 2) +
          e * (c ^ 2 - b * d) * (a * c * d - 2 * b ^ 2 * d + b * c ^ 2) +
          (c ^ 2 - b * d) ^ 2 * (c * d - b * e) +
          a * c * d * (b * d * e - 2 * c ^ 2 * e + c * d ^ 2) +
          2 * c ^ 2 * e * (a * c * d - 2 * b ^ 2 * d + b * c ^ 2) := by ring
      have h0 : 0 ≤ c ^ 2 * d * (a * c * e + a * d ^ 2 - 2 * b ^ 2 * e - b * c * d + c ^ 3) := by
        rw [key]; positivity
      have hpos : 0 < c ^ 2 * d := by positivity
      exact (mul_nonneg_iff_of_pos_left hpos).mp h0
  · constructor <;> simp
  · constructor
    · simp
    · simp only [mul_zero, ne_eq, OfNat.ofNat_ne_zero, not_false_eq_true, zero_pow, add_zero,
        sub_self, zero_add]
      positivity

/-- **Hillion–Johnson, Theorem A.2** (their (78)): `C1(k) ≥ 0` for every `k`. -/
theorem pbPmf_hjC1 (n : ℕ) (p : Fin n → ℝ) (hp : ∀ i, 0 < p i ∧ p i < 1) :
    ∀ k, 0 ≤ hjC1 (pbPmf n p) k := by
  induction n with
  | zero =>
    intro k
    unfold hjC1
    simp only [pbPmf_zero]
    split_ifs <;> first | omega | norm_num
  | succ n ih =>
    intro k
    set g := pbPmf n (Fin.tail p)
    have hp' : ∀ i, 0 < Fin.tail p i ∧ Fin.tail p i < 1 := fun i => hp i.succ
    have ih := ih (Fin.tail p) hp'
    obtain ⟨hnn, hpos⟩ := pbPmf_nonneg_pos n (Fin.tail p) hp'
    have hzero := pbPmf_zero_out n (Fin.tail p)
    have hlc := (pbPmf_lc n (Fin.tail p) hp').1
    have hcross := cross_ineq g n hnn hpos hzero hlc
    obtain ⟨q0, q1⟩ := hp 0
    rw [pbPmf_succ, hjC1_conv]
    have hC1k : 0 ≤ g (k - 2) * g k * g (k + 1) - 2 * g (k - 1) ^ 2 * g (k + 1) +
        g (k - 1) * g k ^ 2 := by
      have := ih k; unfold hjC1 at this; linarith
    have hC1km1 : 0 ≤ g (k - 3) * g (k - 1) * g k - 2 * g (k - 2) ^ 2 * g k +
        g (k - 2) * g (k - 1) ^ 2 := by
      have := ih (k - 1); unfold hjC1 at this
      rw [show k - 1 - 1 = k - 2 by ring, show k - 1 + 1 = k by ring,
        show k - 1 - 2 = k - 3 by ring] at this
      linarith
    have hD : 0 ≤ g (k - 1) ^ 2 - g (k - 2) * g k := by
      have := hlc (k - 1)
      rw [show k - 1 - 1 = k - 2 by ring, show k - 1 + 1 = k by ring] at this
      linarith
    have hX : 0 ≤ g (k - 1) * g k - g (k - 2) * g (k + 1) := by linarith [hcross k]
    have hbd : (0 < g (k - 1) ∧ 0 < g k) ∨ (g (k - 3) = 0 ∧ g (k - 2) = 0 ∧ g (k - 1) = 0) ∨
        (g k = 0 ∧ g (k + 1) = 0) := by
      by_cases h1 : 1 ≤ k
      · by_cases h2 : k ≤ n
        · exact Or.inl ⟨hpos _ (by omega) (by omega), hpos _ (by omega) h2⟩
        · exact Or.inr (Or.inr ⟨hzero _ (Or.inr (by omega)), hzero _ (Or.inr (by omega))⟩)
      · exact Or.inr (Or.inl ⟨hzero _ (Or.inl (by omega)), hzero _ (Or.inl (by omega)),
          hzero _ (Or.inl (by omega))⟩)
    obtain ⟨hXc, hYc⟩ := hj_middle_nonneg (hnn (k - 3)) (hnn (k - 2)) (hnn (k - 1)) (hnn k)
      (hnn (k + 1)) (by linarith) hC1km1 hD hX hbd
    have hs : 0 ≤ 1 - p 0 := by linarith
    have t1 := mul_nonneg (pow_nonneg hs 3) hC1k
    have t2 := mul_nonneg (mul_nonneg (pow_nonneg hs 2) q0.le) hXc
    have t3 := mul_nonneg (mul_nonneg hs (pow_nonneg q0.le 2)) hYc
    have t4 := mul_nonneg (pow_nonneg q0.le 3) hC1km1
    linarith

/-! ## Duality -/

/-- Reflection: `P(W = n - k) = P(W' = k)` where `W'` has parameters `1 - p i`. -/
lemma pbPmf_reflect (n : ℕ) (p : Fin n → ℝ) (k : ℤ) :
    pbPmf n p ((n : ℤ) - k) = pbPmf n (fun i => 1 - p i) k := by
  induction n generalizing k with
  | zero =>
    rw [pbPmf_zero, pbPmf_zero]
    simp only [Nat.cast_zero, zero_sub, neg_eq_zero]
  | succ n ih =>
    rw [pbPmf_succ, pbPmf_succ]
    have ht : Fin.tail (fun i => 1 - p i) = fun i => 1 - Fin.tail p i := rfl
    rw [ht]
    unfold conv
    have h1 := ih (Fin.tail p) (k - 1)
    have h2 := ih (Fin.tail p) k
    rw [show ((n : ℤ) - (k - 1)) = ((n + 1 : ℕ) : ℤ) - k by push_cast; ring] at h1
    rw [show ((n : ℤ) - k) = ((n + 1 : ℕ) : ℤ) - k - 1 by push_cast; ring] at h2
    rw [h1, h2]
    simp only [Fin.tail]
    ring

/-- **Hillion–Johnson, Corollary A.3** (their (79)), in `C1` form. -/
theorem pbPmf_hjC1_dual (n : ℕ) (p : Fin n → ℝ) (hp : ∀ i, 0 < p i ∧ p i < 1) (k : ℤ) :
    0 ≤ pbPmf n p (k + 1) * pbPmf n p k ^ 2 - 2 * pbPmf n p (k + 1) ^ 2 * pbPmf n p (k - 1) +
      pbPmf n p k * pbPmf n p (k - 1) * pbPmf n p (k + 2) := by
  have hp' : ∀ i, 0 < 1 - p i ∧ 1 - p i < 1 := fun i => ⟨by linarith [(hp i).2],
    by linarith [(hp i).1]⟩
  have := pbPmf_hjC1 n (fun i => 1 - p i) hp' ((n : ℤ) - k)
  have hr : ∀ j : ℤ, pbPmf n p j = pbPmf n (fun i => 1 - p i) ((n : ℤ) - j) := by
    intro j; rw [← pbPmf_reflect]; congr 1; ring
  unfold hjC1 at this
  rw [hr (k + 1), hr k, hr (k - 1), hr (k + 2)]
  rw [show (n : ℤ) - k - 1 = n - (k + 1) by ring, show (n : ℤ) - k + 1 = n - (k - 1) by ring,
    show (n : ℤ) - k - 2 = n - (k + 2) by ring] at this
  linarith

theorem pbPmf_hj_left (n : ℕ) (p : Fin n → ℝ) (hp : ∀ i, 0 < p i ∧ p i < 1) (k : ℤ) :
    (pbPmf n p (k - 1) ^ 2 - pbPmf n p (k - 2) * pbPmf n p k) * pbPmf n p (k + 1) ≤
      pbPmf n p (k - 1) * (pbPmf n p k ^ 2 - pbPmf n p (k - 1) * pbPmf n p (k + 1)) := by
  have := pbPmf_hjC1 n p hp k
  unfold hjC1 at this
  linarith

theorem pbPmf_hj_right (n : ℕ) (p : Fin n → ℝ) (hp : ∀ i, 0 < p i ∧ p i < 1) (k : ℤ) :
    (pbPmf n p (k + 1) ^ 2 - pbPmf n p k * pbPmf n p (k + 2)) * pbPmf n p (k - 1) ≤
      pbPmf n p (k + 1) * (pbPmf n p k ^ 2 - pbPmf n p (k - 1) * pbPmf n p (k + 1)) := by
  have := pbPmf_hjC1_dual n p hp k
  linarith

end

end PBDeduction
