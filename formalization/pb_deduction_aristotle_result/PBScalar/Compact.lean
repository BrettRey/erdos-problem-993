import PBScalar.Quadratic

/-!
# The compact range `1/17 ≤ δ < 1/4` (G1)

Method (exact, kernel-checked, different from the paper's Bernstein expansions but
covering the same cells): for fixed `K`, the weights `R r`, `L r` are nonnegative and
nonincreasing in `δ` as long as `(K+1) δ ≤ 1`, hence `A δ K` is nonincreasing in `δ`
(`A_antitone`). The target `(3+δ)/(4δ²)` is also nonincreasing in `δ`. Thus on a
sub-interval `[d_lo, d_hi]` it suffices to check the single rational inequality
`(3+d_lo)/(4 d_lo²) ≤ A d_hi K`, which we evaluate exactly in `ℚ` with `decide +kernel`.

The cells are `δ ∈ [1/(m+2), 1/(m+1)]` with `K = m`, i.e. `H = 1/δ - 1 ∈ [m, m+1]`.
For `K = m` the hypotheses `(K+1)δ < 1 ≤ (K+2)δ` place `δ` in `[1/(m+2), 1/(m+1))`,
i.e. `H ∈ (m, m+1]`, so each cell is used with the correct value of `K`.

* `m = 3`: 20 equal sub-intervals of `[1/5, 1/4]`;
* `m = 4, …, 15`: a single check per cell.
-/

namespace PBScalar

open Finset

/-- Exact rational model of `R δ r` (via its one-step recursion). -/
def Rq (d : ℚ) : ℕ → ℚ
  | 0 => 1
  | r + 1 => Rq d r * ((1 - 2 * d) / (1 - d) * (1 - (r : ℚ) * d))

/-- Exact rational model of `L δ r` (via its one-step recursion). -/
def Lq (d : ℚ) : ℕ → ℚ
  | 0 => 1
  | r + 1 => Lq d r * ((1 - ((r : ℚ) + 2) * d) / (1 - d))

def s0q (d : ℚ) (K : ℕ) : ℚ := 1 + ∑ r ∈ Finset.range K, (Rq d (r + 1) + Lq d (r + 1))

def s1q (d : ℚ) (K : ℕ) : ℚ :=
  ∑ r ∈ Finset.range K, ((r : ℚ) + 1) * (Rq d (r + 1) - Lq d (r + 1))

def s2q (d : ℚ) (K : ℕ) : ℚ :=
  ∑ r ∈ Finset.range K, ((r : ℚ) + 1) ^ 2 * (Rq d (r + 1) + Lq d (r + 1))

/-- Exact rational model of `A δ K`. -/
def Aq (d : ℚ) (K : ℕ) : ℚ := s0q d K * s2q d K - s1q d K ^ 2

/-- Exact rational model of the target `(3 + δ) / (4 δ²)`. -/
def tgq (d : ℚ) : ℚ := (3 + d) / (4 * d ^ 2)

lemma Rq_cast (d : ℚ) (r : ℕ) : ((Rq d r : ℚ) : ℝ) = R (d : ℝ) r := by
  induction r with
  | zero => simp [Rq, R_zero]
  | succ r ih => rw [Rq, R_succ, ← ih, a]; push_cast; ring

lemma Lq_cast (d : ℚ) (r : ℕ) : ((Lq d r : ℚ) : ℝ) = L (d : ℝ) r := by
  induction r with
  | zero => simp [Lq, L_zero]
  | succ r ih => rw [Lq, L_succ, ← ih]; push_cast; ring

lemma Aq_cast (d : ℚ) (K : ℕ) : ((Aq d K : ℚ) : ℝ) = A (d : ℝ) K := by
  rw [A_formula, Aq, s0q, s1q, s2q, s0, s1, s2]
  push_cast
  simp only [Rq_cast, Lq_cast]

/-- The target `(3+δ)/(4δ²)` is nonincreasing in `δ > 0`. -/
lemma target_antitone (d δ : ℝ) (hd : 0 < d) (hdδ : d ≤ δ) :
    (3 + δ) / (4 * δ ^ 2) ≤ (3 + d) / (4 * d ^ 2) := by
  have hδ : 0 < δ := lt_of_lt_of_le hd hdδ
  rw [div_le_div_iff₀ (by positivity) (by positivity)]
  nlinarith [mul_pos hd hδ, mul_nonneg (mul_pos hd hδ).le (sub_nonneg.mpr hdδ),
    mul_nonneg (add_pos hd hδ).le (sub_nonneg.mpr hdδ)]

/-- One sub-interval: the rational check at the endpoints gives the inequality
on the whole sub-interval. -/
theorem piece (K : ℕ) (dlo dhi : ℚ) (h0 : 0 < dlo) (hK : ((K : ℚ) + 1) * dhi ≤ 1)
    (hc : tgq dlo ≤ Aq dhi K) (δ : ℝ) (h1 : (dlo : ℝ) ≤ δ) (h2 : δ ≤ (dhi : ℝ)) :
    (3 + δ) / (4 * δ ^ 2) ≤ A δ K := by
  have h0' : (0 : ℝ) < dlo := by exact_mod_cast h0
  have hK' : ((K : ℝ) + 1) * (dhi : ℝ) ≤ 1 := by exact_mod_cast hK
  have hc' : ((tgq dlo : ℚ) : ℝ) ≤ ((Aq dhi K : ℚ) : ℝ) := by exact_mod_cast hc
  rw [Aq_cast, tgq] at hc'
  push_cast at hc'
  calc (3 + δ) / (4 * δ ^ 2) ≤ (3 + (dlo : ℝ)) / (4 * (dlo : ℝ) ^ 2) :=
        target_antitone _ _ h0' h1
    _ ≤ A (dhi : ℝ) K := hc'
    _ ≤ A δ K := A_antitone δ dhi K (lt_of_lt_of_le h0' h1) h2 hK'

/-- Uniform subdivision of `[a, a + (N+1) h]` into `N+1` pieces. -/
theorem cover (K : ℕ) (a h : ℚ) (ha : 0 < a) (hh : 0 < h) (N : ℕ)
    (hK : ((K : ℚ) + 1) * (a + ((N : ℚ) + 1) * h) ≤ 1)
    (hc : ∀ k ≤ N, tgq (a + (k : ℚ) * h) ≤ Aq (a + ((k : ℚ) + 1) * h) K)
    (δ : ℝ) (h1 : (a : ℝ) ≤ δ) (h2 : δ ≤ ((a + ((N : ℚ) + 1) * h : ℚ) : ℝ)) :
    (3 + δ) / (4 * δ ^ 2) ≤ A δ K := by
  induction N with
  | zero =>
    exact piece K a (a + h) ha (by simpa using hK) (by simpa using hc 0 le_rfl) δ h1
      (by simpa using h2)
  | succ N ih =>
    by_cases hδ : δ ≤ ((a + ((N : ℚ) + 1) * h : ℚ) : ℝ)
    · apply ih _ (fun k hk => hc k (by omega)) hδ
      push_cast at hK
      have : ((N : ℚ) + 1) * h ≤ ((N : ℚ) + 1 + 1) * h := by nlinarith
      have hK0 : (0 : ℚ) ≤ (K : ℚ) + 1 := by positivity
      nlinarith
    · push_neg at hδ
      apply piece K (a + ((N + 1 : ℕ) : ℚ) * h) (a + (((N + 1 : ℕ) : ℚ) + 1) * h)
        (by positivity) hK (hc (N + 1) le_rfl) δ
      · push_cast at hδ ⊢; exact hδ.le
      · exact h2

/-- **G1, cell `H ∈ [3,4]`** (`K = 3`, unsymmetrised weights):
the inequality holds for `1/5 ≤ δ ≤ 1/4`. -/
theorem cell_K3 (δ : ℝ) (h1 : 1 / 5 ≤ δ) (h2 : δ ≤ 1 / 4) :
    (3 + δ) / (4 * δ ^ 2) ≤ A δ 3 := by
  have hc : ∀ k : ℕ, k ≤ 19 → tgq (1 / 5 + (k : ℚ) * (1 / 400)) ≤
      Aq (1 / 5 + ((k : ℚ) + 1) * (1 / 400)) 3 := by decide +kernel
  refine cover 3 (1 / 5) (1 / 400) (by norm_num) (by norm_num) 19 (by norm_num) hc δ ?_ ?_
  · push_cast; linarith
  · push_cast; linarith

/-- **G1, cells `H ∈ [m, m+1]` for `m = 4, …, 15`** (with `K = m`):
the inequality holds for `1/(m+2) ≤ δ ≤ 1/(m+1)`. -/
theorem cell_K (m : ℕ) (hm4 : 4 ≤ m) (hm15 : m ≤ 15) (δ : ℝ)
    (h1 : 1 / ((m : ℝ) + 2) ≤ δ) (h2 : δ ≤ 1 / ((m : ℝ) + 1)) :
    (3 + δ) / (4 * δ ^ 2) ≤ A δ m := by
  have hc : ∀ m : ℕ, m ∈ Finset.Icc 4 15 →
      tgq (1 / ((m : ℚ) + 2)) ≤ Aq (1 / ((m : ℚ) + 1)) m := by decide +kernel
  apply piece m (1 / ((m : ℚ) + 2)) (1 / ((m : ℚ) + 1)) (by positivity)
    (by rw [mul_one_div_cancel (by positivity)]) (hc m (Finset.mem_Icc.mpr ⟨hm4, hm15⟩)) δ
  · push_cast; exact h1
  · push_cast; exact h2

/-- **G1** in its final form (same statement as `scalar_inequality_compact`). -/
theorem compact_main (δ : ℝ) (K : ℕ)
    (hδ17 : 1 / 17 ≤ δ) (hδ4 : δ < 1 / 4)
    (hK : ((K : ℝ) + 1) * δ < 1) (hK' : 1 ≤ ((K : ℝ) + 2) * δ) :
    (3 + δ) / (4 * δ ^ 2) ≤ A δ K := by
  have hδ0 : 0 < δ := by linarith
  have hKlt : (K : ℝ) < 16 := by nlinarith
  have hKgt : (2 : ℝ) < K := by nlinarith
  have hKlt' : K < 16 := by exact_mod_cast hKlt
  have hKgt' : 2 < K := by exact_mod_cast hKgt
  have hlo : 1 / ((K : ℝ) + 2) ≤ δ := by
    rw [div_le_iff₀ (by positivity)]; linarith
  have hhi : δ ≤ 1 / ((K : ℝ) + 1) := by
    rw [le_div_iff₀ (by positivity)]; linarith
  by_cases h3 : K = 3
  · subst h3
    norm_num at hlo
    exact cell_K3 δ (by linarith) hδ4.le
  · exact cell_K K (by omega) (by omega) δ hlo hhi

end PBScalar
