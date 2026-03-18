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

open Finset BigOperators

namespace Formal.STP2Closure

/-! ### Core definitions -/

/-- Convolution of two sequences (polynomial multiplication). -/
def conv (f g : ℕ → ℤ) (k : ℕ) : ℤ :=
  ∑ i in Finset.range (k + 1), f i * g (k - i)

/-- Ladder minor: Λ(f,g)(k) = f(k)·g(k) - f(k-1)·g(k+1).
    Nonneg ladder minors ↔ f(k)/g(k) is nonincreasing (STP2). -/
def ladder (f g : ℕ → ℤ) (k : ℕ) : ℤ :=
  f k * g k - f (k - 1) * g (k + 1)

/-- STP2(f,g): the ratio f(k)/g(k) is nonincreasing.
    Equivalently, all ladder minors Λ(f,g)(k) ≥ 0.
    (Name: "Strongly log-concave Total Positivity order 2".) -/
def isSTP2 (f g : ℕ → ℤ) : Prop :=
  ∀ k, 0 ≤ ladder f g k

/-- Log-concavity: f(k)² ≥ f(k-1)·f(k+1) for all k. -/
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

/-! ### Auxiliary lemmas (provable, included as building blocks) -/

/-- STP2 diagonal form: under LC(f), STP2(f,g) ↔ f(k)g(k) ≥ f(k-1)g(k+1). -/
theorem stp2_iff_ladder (f g : ℕ → ℤ) :
    isSTP2 f g ↔ ∀ k, f (k - 1) * g (k + 1) ≤ f k * g k := by
  simp [isSTP2, ladder]; intro h k; linarith [h k]

/-- Convolution is commutative. -/
theorem conv_comm (f g : ℕ → ℤ) (k : ℕ) :
    conv f g k = conv g f k := by
  sorry

/-- Convolution with nonneg sequences preserves nonnegativity. -/
theorem conv_nonneg (f g : ℕ → ℤ) (hf : isNonneg f) (hg : isNonneg g) :
    isNonneg (conv f g) := by
  sorry

/-- LC is preserved under convolution (Newton's inequality for products).
    This is a classical result: product of PF2 sequences is PF2. -/
theorem lc_conv (f g : ℕ → ℤ) (hf : isLC f) (hg : isLC g)
    (hfn : isNonneg f) (hgn : isNonneg g) :
    isLC (conv f g) := by
  sorry

/-! ### The Cauchy-Binet expansion -/

/-- **CB identity (Form 1):**
    Λ(P*A, Q*B)(k) = Σᵢ Σⱼ Δ(A,B)(i,j) · P(k-i) · Q(k+1-j)

    where Δ(A,B)(i,j) = A(i)B(j) - A(j)B(i).

    Terms with j ≥ i are nonneg (from STP2(A,B)).
    Terms with j < i can be negative. -/
theorem cb_expansion_form1 (P Q A B : ℕ → ℤ) (k : ℕ)
    (hPs : hasSupport P (k + 1)) (hAs : hasSupport A (k + 1))
    (hQs : hasSupport Q (k + 2)) (hBs : hasSupport B (k + 2)) :
    ladder (conv P A) (conv Q B) k =
    ∑ i in Finset.range (k + 1), ∑ j in Finset.range (k + 2),
      minor₂ A B i j * P (k - i) * Q (k + 1 - j) := by
  sorry

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
    -- STP2 for each child
    (h_stp2_1 : isSTP2 I₁ E₁)
    (h_stp2_2 : isSTP2 I₂ E₂)
    -- Log-concavity
    (h_lc_I1 : isLC I₁) (h_lc_E1 : isLC E₁)
    (h_lc_I2 : isLC I₂) (h_lc_E2 : isLC E₂)
    -- Nonneg coefficients
    (h_nn_I1 : isNonneg I₁) (h_nn_E1 : isNonneg E₁)
    (h_nn_I2 : isNonneg I₂) (h_nn_E2 : isNonneg E₂)
    -- Coefficientwise domination (I ≥ E, since I = E + xJ with J ≥ 0)
    (h_dom_1 : ∀ k, E₁ k ≤ I₁ k)
    (h_dom_2 : ∀ k, E₂ k ≤ I₂ k)
    -- Initial values (from tree DP: E(0) = 1 or 0, I(0) = 1)
    (h_I1_0 : I₁ 0 = 1) (h_I2_0 : I₂ 0 = 1)
    -- Finite support
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
    -- Product sequences
    (PI PE : ℕ → ℤ)
    (hPI : PI = List.foldr conv (fun _ => 0) (List.ofFn I))
    (hPE : PE = List.foldr conv (fun _ => 0) (List.ofFn E))
    : isSTP2 PI PE := by
  sorry

end Formal.STP2Closure
