/-
  STP2 Closure Under Convolution — THE OPEN PROBLEM

  This is the core gap in proving Erdős Problem #993 (unimodality of
  independent set sequences of trees).

  The conjecture: if two pairs of nonneg integer sequences each satisfy
  the STP2 (strongly log-concave ratio) condition, then their convolution
  products also satisfy STP2, under tree-realizable constraints.

  This is FALSE for arbitrary sequences (counterexample exists).
  It is TRUE for all tree-realizable pairs through n=22 (millions of checks).
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
    Equivalently, all ladder minors Λ(f,g)(k) ≥ 0.
    (Name: "Strongly log-concave Total Positivity order 2".) -/
def isSTP2 (f g : ℕ → ℤ) : Prop :=
  ∀ k, 0 ≤ ladder f g k

/-- Log-concavity: f(k)² ≥ f(k-1)·f(k+1) for all k.

    **Note on ℕ subtraction:** at k = 0, `k - 1 = 0` in ℕ, so the
    condition becomes f(0) · f(1) ≤ f(0)², i.e. f(0) · (f(0) - f(1)) ≥ 0.
    For nonneg sequences this means f(0) ≥ f(1) or f(0) = 0. -/
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

/-! ### Auxiliary lemmas (proved) -/

/-- STP2 diagonal form: STP2(f,g) ↔ f(k-1)·g(k+1) ≤ f(k)·g(k). -/
theorem stp2_iff_ladder (f g : ℕ → ℤ) :
    isSTP2 f g ↔ ∀ k, f (k - 1) * g (k + 1) ≤ f k * g k := by
  simp only [isSTP2, ladder]
  constructor <;> intro h k <;> linarith [h k]

/-- Convolution is commutative. -/
theorem conv_comm (f g : ℕ → ℤ) (k : ℕ) :
    conv f g k = conv g f k := by
  unfold conv
  rw [← Finset.sum_flip]
  exact Finset.sum_congr rfl fun i hi => by
    rw [Nat.sub_sub_self (Finset.mem_range_succ_iff.mp hi)]
    ring

/-- Convolution with nonneg sequences preserves nonnegativity. -/
theorem conv_nonneg (f g : ℕ → ℤ) (hf : isNonneg f) (hg : isNonneg g) :
    isNonneg (conv f g) := by
  exact fun k => Finset.sum_nonneg fun i _ => mul_nonneg (hf i) (hg _)

/-- Convolution at zero: conv(f,g)(0) = f(0) * g(0). -/
theorem conv_zero (f g : ℕ → ℤ) : conv f g 0 = f 0 * g 0 := by
  unfold conv
  norm_num

/-- Convolution with the zero function gives zero. -/
theorem conv_zero_right (f : ℕ → ℤ) (k : ℕ) :
    conv f (fun _ => 0) k = 0 := by
  exact Finset.sum_eq_zero fun i _ => mul_eq_zero_of_right _ rfl

/-- Convolution with the zero function on the left gives zero. -/
theorem conv_zero_left (g : ℕ → ℤ) (k : ℕ) :
    conv (fun _ => 0) g k = 0 := by
  simp [conv]

/-! ### STP2 for leaf factors -/

/-- A single leaf contributes I(x) = 1 + x, E(x) = 1. -/
def leafI : ℕ → ℤ
  | 0 => 1
  | 1 => 1
  | _ => 0

/-- A single leaf contributes E(x) = 1. -/
def leafE : ℕ → ℤ
  | 0 => 1
  | _ => 0

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

theorem leaf_nonneg_I : isNonneg leafI := by
  intro k
  rcases k with (_ | k | k | k) <;> norm_num [leafI]

theorem leaf_nonneg_E : isNonneg leafE := by
  exact fun k => by
    rcases k with (_ | _ | k) <;> norm_cast

/-! ### Notes on false lemmas

#### `lc_conv` (removed — false)

The statement `isLC f → isLC g → isNonneg f → isNonneg g → isLC (conv f g)`
is **FALSE** under the current definition of `isLC` which uses ℕ subtraction.

Counterexample: f = g = [3, 2, 0, 0, …] satisfies `isLC` (since at k=0,
f(0)·f(1) = 6 ≤ 9 = f(0)², and at k≥1 the standard LC holds).
But conv(f,g)(0) = 9, conv(f,g)(1) = 12, so at k=0:
  conv(0)·conv(1) = 108 > 81 = conv(0)², violating `isLC`.

The classical Newton inequality applies only to the standard LC condition
for k ≥ 1. The k=0 boundary artifact from ℕ subtraction breaks it.

#### `stp2_closure_deg0` (removed — false)

For constant sequences (degree 0), `conv (fun _ => a) (fun _ => b) k = (k+1)*a*b`,
so the convolution is *increasing*. The ladder at k=0 becomes
  (1·a₁a₂)·(1·b₁b₂) − (1·a₁a₂)·(2·b₁b₂) = −a₁a₂b₁b₂ < 0
for positive constants, violating `isSTP2`.
This is again an artifact of ℕ subtraction making `ladder` at k=0 use f(0) instead of f(−1).

#### `cb_expansion_form1` (removed — incorrectly stated)

The Cauchy-Binet expansion had support hypotheses coupling the support
bound to the evaluation point k, which is not how tree-factor supports work.
A correct formulation would decouple the support from k, but getting the
exact index arithmetic right with ℕ subtraction is part of the open problem.
-/

/-! ### THE OPEN PROBLEM -/

/-- **STP2 closure under convolution (2-child case).**

    Given two child factors (I₁, E₁) and (I₂, E₂) from the tree DP,
    each satisfying STP2(I, E), show the product pair also satisfies STP2.

    This is the core gap in proving unimodality of IS sequences of trees.

    **Known to be TRUE:** exhaustive verification through n=22
    (millions of tree-factor pairs, 0 failures).

    **Known to be FALSE without constraints:** generic sequence pairs
    can violate this (Round 13 counterexample).

    The constraints capture "tree-realizability": these aren't arbitrary
    sequences but arise from the tree DP recurrence. -/
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
    (d₁ d₂ : ℕ) (h_supp_I1 : hasSupport I₁ d₁) (h_supp_E1 : hasSupport E₁ d₁)
    (h_supp_I2 : hasSupport I₂ d₂) (h_supp_E2 : hasSupport E₂ d₂)
    : isSTP2 (conv I₁ I₂) (conv E₁ E₂) := by
  sorry

/-- **Corollary: multi-child closure by induction.**
    If the 2-child case holds, the s-child case follows by induction
    (since convolution is associative and the constraints are preserved). -/
theorem stp2_multi_child_closure
    (n : ℕ) (I E : Fin n → ℕ → ℤ)
    (h_stp2 : ∀ i, isSTP2 (I i) (E i))
    (h_lc_I : ∀ i, isLC (I i)) (h_lc_E : ∀ i, isLC (E i))
    (h_nn_I : ∀ i, isNonneg (I i)) (h_nn_E : ∀ i, isNonneg (E i))
    (h_dom : ∀ i k, E i k ≤ I i k)
    (PI PE : ℕ → ℤ)
    (hPI : PI = List.foldr conv (fun _ => 0) (List.ofFn I))
    (hPE : PE = List.foldr conv (fun _ => 0) (List.ofFn E))
    : isSTP2 PI PE := by
  sorry

end Formal.STP2Closure
