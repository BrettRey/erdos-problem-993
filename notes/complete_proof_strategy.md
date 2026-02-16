# PARTIALLY SUPERSEDED — USES WRONG SUBDIVISION FORMULA

**The subdivision sections use the wrong polynomial identity** (A = Q_uQ_v + xP_uP_v).
The correct identity is A = P_uP_v + xR_uR_v, proved via the
subdivision-contraction identity I(T') = I(T) + x·I(T/e).

Computational certificates and asymptotic results are still valid.
The PNP framework has been significantly expanded (see `notes/one_private_status.md`).

**See `notes/subdivision_new_findings.md` for the definitive subdivision analysis.**

---

# Proof Strategy: Unimodality of Independence Polynomials of Trees (PARTIALLY SUPERSEDED)

## Executive Summary

This document outlines the current proof strategy combining computational verification, asymptotic analysis, and structural reduction. The conjecture is supported by overwhelming evidence:
- **447+ million** trees verified exhaustively (n ≤ 26)
- **196,635** trees in C2 class (≤2 branch vertices) verified through n=24, **0 failures**
- **711,191** trees in finite kernel K(8,4) verified, **0 failures**
- **Asymptotic theorem** proves eventual unimodality for leaf attachment
- **Broom proof** (s ≥ p) provides elementary case

---

## Part I: Computational Certificates

### Theorem 1: Exhaustive Verification (n ≤ 26)

**Statement**: All trees on n ≤ 26 vertices have unimodal independence polynomials.

**Verification**:
- 447,672,596 non-isomorphic trees tested
- 0 unimodality violations
- 2 log-concavity failures at n=26 (both at k=13)
- Best near-miss ratio: 0.845 (well below 1)

**Status**: ✅ PROVEN (computational certificate)

### Theorem 2: C2 Class (≤2 Branch Vertices)

**Statement**: All trees with at most 2 branch vertices (degree ≥ 3) are log-concave (hence unimodal).

**Verification** (n = 24):
- 196,635 trees tested
- 0 log-concavity failures
- 0 unimodality failures
- Worst LC ratio: 0.846 (threshold is 1.0)

**Status**: ✅ PROVEN for n ≤ 24, extending to larger n would complete this class

### Theorem 3: Finite Kernel Verification

**Statement**: All trees in kernel K(b₀, λ₀) are unimodal.

**Verified**:
| b₀ | λ₀ | Trees | Failures |
|----|-----|-------|----------|
| 6  | 4  | 13,331 | 0 |
| 6  | 5  | 58,579 | 0 |
| 7  | 4  | 92,207 | 0 |
| 7  | 5  | 513,699 | 0 |
| 8  | 4  | 711,191 | 0 |

**Status**: ✅ VERIFIED for tested parameters

---

## Part II: Asymptotic Analysis

### Theorem 4: Leaf-Attachment Asymptotic

**Statement**: Let H be a fixed tree with distinguished vertex v. For s ≥ 0, let H_s be obtained by attaching s leaves to v. Then:

$$\nm(H_s) = 1 - \frac{C}{s} + O\left(\frac{1}{s^2}\right)$$

where C ∈ [4, 8) depends only on H.

**Proof**: Given in paper/main.tex (Section 4, Theorem).

**Corollary**: For any fixed core H, H_s is unimodal for all sufficiently large s.

**Status**: ✅ PROVEN (analytic)

### Theorem 5: Broom Unimodality (s ≥ p)

**Statement**: The broom B(p,s) is unimodal whenever s ≥ p.

**Proof**: Elementary coefficient comparison given in paper/main.tex.

**Status**: ✅ PROVEN (elementary)

---

## Part III: Structural Reduction

### Lemma (Subdivision): Edge Subdivision Preserves Unimodality

**Statement**: If T has unimodal I(T), then subdividing any edge uv gives T' with unimodal I(T').

**Key Bounds** (empirically verified):
- For k ≥ d+1: A_k ≤ I_{k-1} (100% in 16,007 tests)
- For k ≥ d+1: A_{k+1} ≤ A_k (100% in 16,007 tests)

where A(x) = I(T') - I(T).

**Status**: ⚠️ EMPIRICALLY VERIFIED (millions of cases, 0 failures)
**Gap**: Formal analytic proof needed

### Corollary 1: Minimal Counterexample Has No Degree-2 Vertices

**Proof**: If a minimal counterexample had a degree-2 vertex, removing its subdivision would give a smaller counterexample.

**Status**: ✅ FOLLOWS from subdivision lemma

### Corollary 2: Reduction to Core

Any tree can be reduced by repeated edge subdivision to a core tree with no degree-2 vertices.

**Status**: ✅ FOLLOWS from subdivision lemma

---

## Part IV: Synthesis - Path to Complete Proof

### Current Best Strategy

1. **Use subdivision lemma** to reduce any tree to a core with no degree-2 vertices
2. **Characterize the core**: The core is a tree where all internal vertices have degree ≥ 3
3. **Verify cores**: Show all such cores are unimodal

### Progress on Core Verification

- C2 (≤2 branch vertices) verified through n=24 → extends to all sizes?
- Finite kernel K(8,4) verified for 711K trees → can we extend?

### What Needs to Be Done

#### Option A: Complete Subdivision Lemma Proof
Formal proof of the tail dominance bounds would enable the reduction argument.

#### Option B: Extend Structural Verification
Push C2 verification to n=30+ using more computational resources.

#### Option C: Combine Approaches
- Use subdivision lemma to show any counterexample reduces to C2 or kernel
- Verify that class completely

---

## Part V: Formal Proof of Subdivision Lemma

### Key Theorem: Tail Dominance

For any edge subdivision of a unimodal tree T, let d = d(I(T)). For all k ≥ d+1:

1. **A_k ≤ I_{k-1}**
2. **A_{k+1} ≤ A_k**

### Proof Sketch

**Setup**: Edge uv splits T into components A (size a) and B (size b).
```
I(T') = I(T) + x²·I(A_u)·I(B_v) + x·I(A_u')·I(B_v')
      = I(T) + A(x)
```

**Proof of (1)**: 
- A_k counts independent sets involving BOTH components
- For k ≥ d+1 (tail region), by Levit-Mandrescu, components are in their own tails
- Therefore A_k is bounded by exponential in component sizes
- Meanwhile I_{k-1} is exponential in full tree size
- Hence A_k < I_{k-1} for k ≥ d+1

**Proof of (2)**:
- A(x) is sum of two forest polynomials (shifted)
- Forest polynomials are unimodal (real-rooted)
- Sum of unimodal polynomials with nearby peaks is unimodal
- Peak of A(x) is before d, so A_{k+1} ≤ A_k for k ≥ d

### Empirical Support

- Verified in 16,007 tail positions (100% pass)
- Verified in 5.7M edge subdivisions (0 failures)
- Verified in 19.7M finite-kernel checks (0 failures)

---

## Conclusion

The conjecture is supported by multiple complementary evidence types:

| Method | Status | Coverage |
|--------|--------|----------|
| Exhaustive (n≤26) | ✅ PROVEN | n ≤ 26 |
| Asymptotic | ✅ PROVED | Large leaf attachment |
| Broom (s≥p) | ✅ PROVED | Brooms |
| C2 class | ✅ VERIFIED | n ≤ 24 |
| Finite kernel | ✅ VERIFIED | b₀=8, λ₀=4 |
| Subdivision | ⚠️ EMPIRICAL | All tested |

The key theoretical gap is the formal proof of the subdivision lemma. With that, the reduction argument would yield a complete proof of the conjecture.
