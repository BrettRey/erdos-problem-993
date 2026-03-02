# SCC ★ Identity Verification (2026-02-28)

## Summary

Verified Instance 1's ★ identity (bivariate cancellation) computationally and discovered
the correct decomposition for SCC induction.

## Key findings

### 1. P2 is NOT the right invariant

Instance 1's formula targets P2 (ratio dominance E ≽ J), but P2 FAILS at intermediate
product stages:

| n | P2 fails (final) | P2 fails (intermediate) | T_k < 0 (intermediate) |
|---|---|---|---|
| 8 | 22 | 1 | 1 |
| 14 | 11,323 | 775 | 841 |
| 18 | 745,941 | 58,646 | 66,876 |
| 20 | 6,142,205 | 505,322 | 584,528 |

The bridge lemma (T_k ≥ 0) is **FALSE**. Y_k monotonicity is **FALSE**.
P2 cannot be proved inductively on the number of factors.

### 2. SCC IS the right invariant

Strong Condition C holds at ALL intermediate stages (0 failures through n=22, 112M+ checks).
The product structure for SCC involves `e_new = e_old * Q + x*(1+x)*b_old * R`.

### 3. SCC 5-term decomposition (VERIFIED)

With α = e*Q, β = b*Q, γ = b*R (where e=(1+x)I, b=E accumulated; Q=E_c, R=J_c factor):

    Δ_k^new = Term1 + Term2 + Term3 + Term4 + Term5

| Term | Formula | Sign | Events < 0 |
|------|---------|------|------------|
| 1 | Δ_k(e*Q, b*Q) | **≥ 0 ALWAYS** | 0 |
| 2 | α_{k+1}·γ_{k-1} - α_k·γ_k | can be < 0 | 45.5M |
| 3 | γ_k·β_k - γ_{k-1}·β_{k+1} | can be < 0 | 16.1M |
| 4 | γ_{k-1}·β_k - γ_{k-2}·β_{k+1} | can be < 0 | 4.2M |
| 5 | γ_{k-1}² - γ_{k-2}·γ_k | **≥ 0 ALWAYS** | 0 |

Identity verified: 112M checks, 0 failures.

### 4. Term 1 dominates correction (12:1 ratio)

When correction = Terms 2+3+4+5 < 0 (1.23M events / 112M total ≈ 1.1%):
- **Min ratio Term1/|correction| = 12.09** (at n=20, k=10, stage 2/2)
- Most correction-negative events have ratio > 14

This means the TP2 main term provides 12x more positivity than needed.

### 5. Why Term 1 ≥ 0 is provable

Δ_k(e*Q, b*Q) ≥ 0 follows from:
- Old pair (e, b) satisfies SCC (Δ_k(e, b) ≥ 0) — inductive hypothesis
- Q = E_c is PF2 (LC, nonneg) — Q is a tree DP polynomial
- TP2 closure: convolution with PF2 preserves 2×2 Toeplitz minors (standard result)

### 6. Why Term 5 ≥ 0 is provable

γ_{k-1}² - γ_{k-2}·γ_k is the LC gap of γ = b*R = E_acc * J_c.
- E_acc is LC (verified at all stages)
- J_c = dp1s[c] is LC (product of LC tree polynomials)
- Convolution of LC sequences is LC (standard result)

## Proof landscape

**Provable:**
- Term 1 ≥ 0 (TP2 closure)
- Term 5 ≥ 0 (LC closure)
- SCC = Term1 + correction ≥ 0 (computationally verified)

**Open:** correction ≤ Term1 (i.e., terms 2+3+4+5 ≤ Term1)

The 12:1 ratio is huge slack. A proof needs to bound:
|correction| ≤ f(γ, β, α) ≤ Term1

Since γ ≤ β (because R ≤ Q coefficientwise), the correction terms are
"small perturbations" of the TP2 quantity.

## Relation to T1/T2/T3 decomposition

The T1/T2/T3 decomposition (S_k = b_{k-1}·d_k + b_k·d_{k-1} + a_{k-1}·c_k)
is at the TREE level (final product). The 5-term decomposition here is at the
INCREMENTAL PRODUCT level (one factor at a time). They capture different aspects
of the same phenomenon: curvature dominates ratio-dominance deficits.

The T1+T3 ratio (≥ 4.23 at n=22) and the Term1/|correction| ratio (≥ 12.09 at n=20)
are related but measure different decompositions.
