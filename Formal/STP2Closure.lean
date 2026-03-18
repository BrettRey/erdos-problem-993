/-
  STP2 Closure Under Convolution — THE OPEN PROBLEM

  This is the core gap in proving Erdős Problem #993 (unimodality of
  independent set sequences of trees).

  The conjecture: if two pairs of nonneg integer sequences each satisfy
  the STP2 (strongly log-concave ratio) condition, then their convolution
  products also satisfy STP2, under tree-realizable constraints.

  This is FALSE for arbitrary sequences (counterexample exists).
  It is TRUE for all tree-realizable pairs through n=22 (millions of checks).

  **Key discovery (this file):** the current hypotheses are too weak.
  Without a gap-freeness condition such as `noGaps`, sequences like
  `[1, 1, 0, 0, 5, 0, …]` satisfy the other listed hypotheses but their
  convolution violates STP2 at k=3.
-/
import Mathlib

open Finset

namespace Formal.STP2Closure

/-! ### Core definitions -/

/-- Convolution of two sequences (polynomial multiplication). -/
def conv (f g : ℕ → ℤ) (k : ℕ) : ℤ :=
  ∑ i ∈ Finset.range (k + 1), f i * g (k - i)

/-- Ladder minor: Λ(f,g)(k) = f(k)·g(k) - f(k-1)·g(k+1).
    Nonneg ladder minors ↔ f(k)/g(k) is nonincreasing (STP2). -/
def ladder (f g : ℕ → ℤ) (k : ℕ) : ℤ :=
  f k * g k - f (k - 1) * g (k + 1)

/-- STP2(f,g): the ratio f(k)/g(k) is nonincreasing.
    Equivalently, all ladder minors Λ(f,g)(k) ≥ 0. -/
def isSTP2 (f g : ℕ → ℤ) : Prop :=
  ∀ k, 0 ≤ ladder f g k

/-- Log-concavity: f(k)² ≥ f(k-1)·f(k+1) for all k.
    Note: at k = 0, `k - 1 = 0` in ℕ, so the condition becomes
    f(0)·f(1) ≤ f(0)², i.e. f(0)·(f(0) - f(1)) ≥ 0. -/
def isLC (f : ℕ → ℤ) : Prop :=
  ∀ k, f (k - 1) * f (k + 1) ≤ f k * f k

/-- Nonneg coefficients. -/
def isNonneg (f : ℕ → ℤ) : Prop :=
  ∀ k, 0 ≤ f k

/-- Finite support: f(k) = 0 for k > d. -/
def hasSupport (f : ℕ → ℤ) (d : ℕ) : Prop :=
  ∀ k, d < k → f k = 0

/-- The 2×2 minor: Δ(A,B)(i,j) = A(i)B(j) - A(j)B(i). -/
def minor₂ (A B : ℕ → ℤ) (i j : ℕ) : ℤ :=
  A i * B j - A j * B i

/-- **No internal zeros (gap-free support).**
    If f(k) = 0 for some k ≥ 1, then f(k+1) = 0.
    Combined with nonneg and f(0) > 0, this forces the support
    to be a contiguous initial segment {0, 1, …, d}.

    This is introduced here as a candidate extra hypothesis motivated by
    the internal-zeros counterexample below. The present file does not
    prove that it follows from the tree DP hypotheses. -/
def noGaps (f : ℕ → ℤ) : Prop :=
  ∀ k, 0 < k → f k = 0 → f (k + 1) = 0

/-! ### Basic convolution lemmas -/

theorem stp2_iff_ladder (f g : ℕ → ℤ) :
    isSTP2 f g ↔ ∀ k, f (k - 1) * g (k + 1) ≤ f k * g k := by
  simp only [isSTP2, ladder]
  constructor <;> intro h k <;> linarith [h k]

theorem conv_comm (f g : ℕ → ℤ) (k : ℕ) :
    conv f g k = conv g f k := by
  unfold conv
  rw [← Finset.sum_flip]
  exact Finset.sum_congr rfl fun i hi => by
    rw [Nat.sub_sub_self (Finset.mem_range_succ_iff.mp hi)]
    ring

theorem conv_nonneg (f g : ℕ → ℤ) (hf : isNonneg f) (hg : isNonneg g) :
    isNonneg (conv f g) :=
  fun k => Finset.sum_nonneg fun i _ => mul_nonneg (hf i) (hg _)

theorem conv_zero (f g : ℕ → ℤ) : conv f g 0 = f 0 * g 0 := by
  unfold conv
  norm_num

theorem conv_zero_right (f : ℕ → ℤ) (k : ℕ) :
    conv f (fun _ => 0) k = 0 :=
  Finset.sum_eq_zero fun i _ => mul_eq_zero_of_right _ rfl

theorem conv_zero_left (g : ℕ → ℤ) (k : ℕ) :
    conv (fun _ => 0) g k = 0 := by
  simp [conv]

/-! ### The Dirac delta: convolution identity -/

/-- The Dirac delta at zero: δ₀(0) = 1, δ₀(k) = 0 for k > 0.
    This is the identity element for convolution. -/
def δ0 : ℕ → ℤ
  | 0 => 1
  | _ => 0

/-- Convolution with δ₀ on the right is the identity. -/
theorem conv_delta_right (f : ℕ → ℤ) (k : ℕ) :
    conv f δ0 k = f k := by
  unfold conv
  have h_delta : ∀ i ∈ Finset.range (k + 1),
      δ0 (k - i) = if i = k then 1 else 0 := by
    unfold δ0
    grind
  rw [Finset.sum_congr rfl fun i hi => by rw [h_delta i hi]]
  aesop

/-- Convolution with δ₀ on the left is the identity. -/
theorem conv_delta_left (f : ℕ → ℤ) (k : ℕ) :
    conv δ0 f k = f k := by
  unfold conv δ0
  rw [Finset.sum_eq_single 0] <;> aesop

theorem delta_nonneg : isNonneg δ0 :=
  fun k => by
    cases k <;> norm_num [δ0]

theorem delta_lc : isLC δ0 := by
  intro k
  rcases k with (_ | _ | k) <;> norm_num [δ0]

theorem delta_noGaps : noGaps δ0 := by
  intro k hk_pos hk_zero
  cases k <;> aesop

theorem delta_support : hasSupport δ0 0 :=
  fun k hk => by
    cases k <;> trivial

/-- STP2(f, δ₀) holds for any nonneg f. -/
theorem stp2_delta_right (f : ℕ → ℤ)
    (hf : isNonneg f) :
    isSTP2 f δ0 := by
  intro k
  by_cases hk : k = 0
  · simp +decide [hk, ladder, δ0]
    exact hf 0
  · unfold ladder δ0
    cases k <;> aesop

/-! ### Leaf factors -/

/-- A single leaf contributes I(x) = 1 + x. -/
def leafI : ℕ → ℤ
  | 0 => 1
  | 1 => 1
  | _ => 0

/-- A single leaf contributes E(x) = 1 = δ₀. -/
def leafE : ℕ → ℤ
  | 0 => 1
  | _ => 0

theorem leafE_eq_delta : leafE = δ0 :=
  funext fun n => by
    rcases n with (_ | _ | n) <;> rfl

theorem leaf_stp2 : isSTP2 leafI leafE := by
  intro k
  rcases k with (_ | _ | k) <;> norm_num [ladder, leafI, leafE]

theorem leaf_lc_I : isLC leafI := by
  intro k
  rcases k with (_ | _ | k) <;> simp +arith +decide [*, leafI]

theorem leaf_lc_E : isLC leafE := by
  intro k
  rcases k with (_ | _ | k) <;> simp +decide [*]
  rfl

theorem leaf_nonneg_I : isNonneg leafI :=
  fun k => by
    rcases k with (_ | k | k | k) <;> norm_num [leafI]

theorem leaf_nonneg_E : isNonneg leafE :=
  fun k => by
    rcases k with (_ | _ | k) <;> norm_cast

theorem leaf_noGaps_I : noGaps leafI := by
  intro k hk_pos hk_zero
  simp [leafI] at *
  rcases k with (_ | _ | k) <;> tauto

theorem leaf_noGaps_E : noGaps leafE := by
  intro k hk
  aesop

/-! ### Structural lemmas -/

/-- Convolution preserves domination. -/
theorem conv_dom (I₁ E₁ I₂ E₂ : ℕ → ℤ)
    (h_nn_I1 : isNonneg I₁) (h_nn_E1 : isNonneg E₁)
    (h_nn_I2 : isNonneg I₂) (h_nn_E2 : isNonneg E₂)
    (h_dom_1 : ∀ k, E₁ k ≤ I₁ k)
    (h_dom_2 : ∀ k, E₂ k ≤ I₂ k) :
    ∀ k, conv E₁ E₂ k ≤ conv I₁ I₂ k :=
  fun k => Finset.sum_le_sum fun i _ =>
    mul_le_mul (h_dom_1 i) (h_dom_2 _) (h_nn_E2 _) (h_nn_I1 _)

theorem conv_zero_eq (f g : ℕ → ℤ) : conv f g 0 = f 0 * g 0 :=
  conv_zero f g

/-- Finite support is preserved by convolution. -/
theorem conv_hasSupport (f g : ℕ → ℤ) (d₁ d₂ : ℕ)
    (hf : hasSupport f d₁) (hg : hasSupport g d₂) :
    hasSupport (conv f g) (d₁ + d₂) := by
  intro k hk
  refine Finset.sum_eq_zero fun i hi => ?_
  by_cases hi' : i ≤ d₁ <;> simp_all +decide [hasSupport]
  exact Or.inr (hg _ (by omega))

/-- noGaps propagates: once a gap-free sequence hits 0, it stays 0. -/
theorem noGaps_tail (f : ℕ → ℤ) (hng : noGaps f)
    (k : ℕ) (hk : 0 < k) (hfk : f k = 0) :
    ∀ j, k ≤ j → f j = 0 := by
  intro j hj
  induction' hj with j _ ih
  · exact hfk
  · cases j <;> aesop

/-- For gap-free sequences with f(0) > 0, support is an initial segment. -/
theorem support_initial_segment (f : ℕ → ℤ)
    (hf : isNonneg f) (hng : noGaps f)
    (h0 : 0 < f 0) (d : ℕ) (hd : hasSupport f d) :
    ∀ k, k ≤ d → (f k = 0 → ∀ j, k ≤ j → f j = 0) := by
  intro k _ hfk j hj
  induction' hj with j _ ih
  · exact hfk
  · cases j <;> aesop

/-! ### The shifted-ladder lemma (key technical result) -/

/-- **Shifted ladder nonnegativity.**
    Under STP2 + LC of the second sequence + nonneg + noGaps:
    I(k-1)·E(k) - I(k-2)·E(k+1) ≥ 0 for all k ≥ 2.

    Proof: Case E(k-1) = 0 uses noGaps to propagate zeros.
    Case E(k-1) > 0 chains the STP2 inequality at k-1 with
    the LC inequality at k, then cancels the positive E(k-1). -/
theorem shifted_ladder_nonneg
    (I E : ℕ → ℤ)
    (h_stp2 : isSTP2 I E) (h_lc_E : isLC E)
    (h_nn_I : isNonneg I) (h_nn_E : isNonneg E)
    (h_ng_E : noGaps E)
    (k : ℕ) (hk : 2 ≤ k) :
    0 ≤ I (k - 1) * E k - I (k - 2) * E (k + 1) := by
  by_cases hk1 : E (k - 1) = 0
  · rcases k with (_ | _ | k) <;> simp_all +decide [noGaps]
  · have h_stp2_km1 : I (k - 1) * E (k - 1) ≥ I (k - 2) * E k := by
      rcases k with (_ | _ | k) <;> simp_all +decide [isSTP2]
      have := h_stp2 (k + 1)
      simp [isSTP2, ladder] at this
      linarith
    have h_lc_k : E k * E k ≥ E (k - 1) * E (k + 1) := h_lc_E k
    cases lt_or_gt_of_ne hk1 <;>
      nlinarith [h_nn_I (k - 1), h_nn_I (k - 2),
                 h_nn_E k, h_nn_E (k + 1), h_nn_E (k - 1)]

/-- Shifted ladder at k = 1: I(0)·E(1) ≥ 0 (trivial from nonneg). -/
theorem shifted_ladder_k1
    (I E : ℕ → ℤ) (h_nn_I : isNonneg I) (h_nn_E : isNonneg E) :
    0 ≤ I 0 * E 1 :=
  mul_nonneg (h_nn_I 0) (h_nn_E 1)

/-! ### Special case: one factor is δ₀ -/

/-- STP2 closure when one factor is δ₀ (trivial: conv with δ₀ is identity). -/
theorem stp2_conv_closure_delta
    (I₁ E₁ : ℕ → ℤ) (h_stp2 : isSTP2 I₁ E₁) :
    isSTP2 (conv I₁ δ0) (conv E₁ δ0) := by
  have h (k : ℕ) : conv I₁ δ0 k = I₁ k ∧ conv E₁ δ0 k = E₁ k :=
    ⟨conv_delta_right _ _, conv_delta_right _ _⟩
  intro k
  simpa only [funext fun i => (h i).1, funext fun i => (h i).2]
    using h_stp2 k

/-! ### Special case: one factor is a leaf pair -/

/-- Convolution with leafI formula for k ≥ 1. -/
theorem conv_leafI_succ (f : ℕ → ℤ) (k : ℕ) :
    conv f leafI (k + 1) = f (k + 1) + f k := by
  unfold conv
  have h_leafI : ∀ j, leafI j = if j = 0 ∨ j = 1 then 1 else 0 := by
    intro j
    rcases j with (_ | _ | j) <;> rfl
  simp +decide [Finset.sum_range_succ, h_leafI]
  rw [Finset.sum_eq_zero] <;> simp +arith +decide
  intros
  omega

theorem conv_leafI_zero (f : ℕ → ℤ) :
    conv f leafI 0 = f 0 := by
  unfold conv
  simp +decide [Finset.sum_range_succ', leafI]

/-- Convolution with leafE is the identity (since leafE = δ₀). -/
theorem conv_leafE (f : ℕ → ℤ) (k : ℕ) :
    conv f leafE k = f k := by
  convert conv_delta_right f k

/-- **STP2 closure when one factor is a leaf (I₂ = leafI, E₂ = leafE).**

    This is the base induction step: attaching a single leaf to a tree
    preserves the STP2 property of the IS/edge-cover polynomial pair.

    The proof decomposes the ladder of the convolution product into:
    - `ladder(I₁, E₁)(k)` (≥ 0 by STP2 hypothesis)
    - A shifted ladder term (≥ 0 by `shifted_ladder_nonneg`)

    Requires `noGaps E₁` to handle the shifted ladder term. -/
theorem stp2_conv_closure_leaf
    (I₁ E₁ : ℕ → ℤ)
    (h_stp2 : isSTP2 I₁ E₁)
    (h_lc_E1 : isLC E₁)
    (h_nn_I1 : isNonneg I₁)
    (h_nn_E1 : isNonneg E₁)
    (h_ng_E1 : noGaps E₁) :
    isSTP2 (conv I₁ leafI) (conv E₁ leafE) := by
  intro k
  by_cases hk : k = 0 ∨ k = 1
  · obtain rfl | rfl := hk
    · convert h_stp2 0 using 1
      unfold ladder
      norm_num [Finset.sum_range_succ', conv_leafI_zero, conv_leafE]
    · simp [ladder, conv_leafI_zero, conv_leafI_succ, conv_leafE]
      have := h_stp2 1
      unfold ladder at this
      nlinarith [h_nn_I1 0, h_nn_I1 1, h_nn_E1 0, h_nn_E1 1,
                 h_lc_E1 1]
  · have h_decomp :
        ladder (conv I₁ leafI) (conv E₁ leafE) k =
        ladder I₁ E₁ k +
        (I₁ (k - 1) * E₁ k - I₁ (k - 2) * E₁ (k + 1)) := by
      rcases k with (_ | _ | k) <;>
        simp_all +decide [ladder, conv_leafI_succ]
      ring!
      rw [show conv E₁ leafE = E₁ from funext fun n => conv_leafE E₁ n]
    rw [h_decomp]
    apply add_nonneg (h_stp2 k)
    exact shifted_ladder_nonneg I₁ E₁ h_stp2 h_lc_E1 h_nn_I1 h_nn_E1
      h_ng_E1 k (by omega)

/-! ### Notes on false lemmas

#### `lc_conv` (removed — false)

The statement `isLC f → isLC g → isNonneg f → isNonneg g → isLC (conv f g)`
is **FALSE** under the current definition of `isLC` which uses ℕ subtraction.

Counterexample: f = g = [3, 2, 0, 0, …] satisfies `isLC` but
conv(f,g)(0) = 9, conv(f,g)(1) = 12, violating `isLC` at k=0.

#### `stp2_closure_deg0` (removed — false)

For constant sequences, `conv (fun _ => a) (fun _ => b) k = (k+1)*a*b`,
which is increasing, violating STP2. This is an artifact of ℕ subtraction.

#### `cb_expansion_form1` (removed — incorrectly stated)

The Cauchy-Binet expansion had flawed support hypotheses.

#### Internal zeros counterexample (NEW — discovered in this file)

The sequence f = [1, 1, 0, 0, 5, 0, …] satisfies isLC, isNonneg,
isSTP2(f,f), and domination (f ≤ f). But
conv(f, f) = [1, 2, 1, 0, 10, 10, 0, 0, 25, …] has
ladder(conv, conv)(3) = 0² − 1·10 = −10 < 0, violating STP2.

This shows that some additional hypothesis, such as `noGaps`,
is necessary for the main theorem. -/

/-! ### THE OPEN PROBLEM -/

/-- **STP2 closure under convolution (2-child case).**

    Given two child factors (I₁, E₁) and (I₂, E₂) from the tree DP,
    each satisfying STP2(I, E), show the product pair also satisfies STP2.

    This is the core gap in proving unimodality of IS sequences of trees.

    **Known to be TRUE:** exhaustive verification through n=22.

    **Known to be FALSE without constraints:** the internal-zeros
    counterexample above shows the current hypothesis list is too weak.

    **Proved special cases in this file:**
    - One factor is δ₀ (`stp2_conv_closure_delta`)
    - One factor is a leaf pair (`stp2_conv_closure_leaf`) -/
theorem stp2_conv_closure
    (I₁ E₁ I₂ E₂ : ℕ → ℤ)
    (h_stp2_1 : isSTP2 I₁ E₁)
    (h_stp2_2 : isSTP2 I₂ E₂)
    (h_lc_I1 : isLC I₁) (h_lc_E1 : isLC E₁)
    (h_lc_I2 : isLC I₂) (h_lc_E2 : isLC E₂)
    (h_nn_I1 : isNonneg I₁) (h_nn_E1 : isNonneg E₁)
    (h_nn_I2 : isNonneg I₂) (h_nn_E2 : isNonneg E₂)
    (h_dom_1 : ∀ k, E₁ k ≤ I₁ k)
    (h_dom_2 : ∀ k, E₂ k ≤ I₂ k)
    (h_I1_0 : I₁ 0 = 1) (h_I2_0 : I₂ 0 = 1)
    (h_ng_I1 : noGaps I₁) (h_ng_E1 : noGaps E₁)
    (h_ng_I2 : noGaps I₂) (h_ng_E2 : noGaps E₂)
    (d₁ d₂ : ℕ)
    (h_supp_I1 : hasSupport I₁ d₁) (h_supp_E1 : hasSupport E₁ d₁)
    (h_supp_I2 : hasSupport I₂ d₂) (h_supp_E2 : hasSupport E₂ d₂) :
    isSTP2 (conv I₁ I₂) (conv E₁ E₂) := by
  sorry

/-! ### Multi-child closure

The multi-child case reduces to the 2-child case by induction,
folding convolution from the identity element δ₀.

**Important caveat:** For the induction to go through, the intermediate
convolution products must satisfy all hypotheses of `stp2_conv_closure`
(STP2, LC, nonneg, domination, noGaps). Some of these are themselves
nontrivial closure questions. -/

/-- Multi-child STP2 closure by iterated convolution from δ₀. -/
theorem stp2_multi_child_closure
    (n : ℕ) (I E : Fin n → ℕ → ℤ)
    (h_stp2 : ∀ i, isSTP2 (I i) (E i))
    (h_lc_I : ∀ i, isLC (I i)) (h_lc_E : ∀ i, isLC (E i))
    (h_nn_I : ∀ i, isNonneg (I i)) (h_nn_E : ∀ i, isNonneg (E i))
    (h_dom : ∀ i k, E i k ≤ I i k)
    (h_ng_I : ∀ i, noGaps (I i)) (h_ng_E : ∀ i, noGaps (E i))
    (PI PE : ℕ → ℤ)
    (hPI : PI = List.foldr conv δ0 (List.ofFn I))
    (hPE : PE = List.foldr conv δ0 (List.ofFn E)) :
    isSTP2 PI PE := by
  sorry

end Formal.STP2Closure
