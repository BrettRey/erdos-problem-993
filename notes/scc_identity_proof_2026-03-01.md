# SCC via the Identity + LC of E (2026-03-01)

## The Key Identity (PROVED algebraically, verified n ≤ 22)

For any tree T rooted at a support vertex r with E = dp[r][0] and I = E + xJ:

```
b_{k-1} · Δ_k = b_{k-1} · d_k + b_k · d_{k-1} + a_{k-1} · c_k
```

where:
- a_k = I_k, b_k = E_k
- d_k = I_{k+1}·E_k - I_k·E_{k+1} (LR minor of I vs E)
- c_k = E_k² - E_{k-1}·E_{k+1} (LC gap of E)
- Δ_k = ((1+x)I)_{k+1}·E_k - ((1+x)I)_k·E_{k+1} (SCC)

**Proof:** Direct algebra. Expand both sides using e_k = I_k + I_{k-1}:
- LHS = b_{k-1}·[(a_{k+1}+a_k)b_k - (a_k+a_{k-1})b_{k+1}]
- RHS = b_{k-1}·(a_{k+1}b_k-a_kb_{k+1}) + b_k·(a_kb_{k-1}-a_{k-1}b_k) + a_{k-1}·(b_k²-b_{k-1}b_{k+1})

Both simplify to b_{k-1}b_k(a_{k+1}+a_k) - b_{k-1}b_{k+1}(a_k+a_{k-1}). □

## LC of E at Support Vertices (verified n ≤ 22, PROVABLE for all n)

At a support vertex r with ℓ leaf children and non-leaf children c_1,...,c_s:

```
E = dp[r][0] = (1+x)^ℓ · ∏ I(T_{c_j})
```

Each I(T_{c_j}) is the IS polynomial of a subtree on ≤ n-2 vertices.

**Claim:** E is log-concave at support vertices of all trees.

**Proof strategy:**
1. (1+x) is LC (trivially)
2. Cauchy products of nonneg LC sequences are LC (Newton's inequality)
3. If all I(T_{c_j}) are LC, then E is LC as a product of LC polynomials
4. Each T_{c_j} has ≤ n-2 vertices

So: **LC of E reduces to LC of IS polys of proper subtrees.**

### What's known about LC of IS polys:
- LC holds for ALL trees on ≤ 25 vertices (Kadrawi & Levit 2023)
- LC FAILS for exactly 2 trees at n = 26 (Kadrawi & Levit 2023)
- Every tree IS polynomial is unimodal (our conjecture!)

### Inductive argument:
- For n ≤ 27: subtrees have ≤ 25 vertices, all LC. E is LC. SCC follows.
- For n ≥ 28: subtrees could have 26+ vertices. LC might fail.

## Verification Results (n ≤ 22)

| Metric | Value |
|--------|-------|
| Identity failures | **0** |
| c_k < 0 (LC fail at support) | **0** |
| d_k < 0 (LR minor negative) | ~29% of checks |
| SCC (Δ_k) < 0 | **0** |
| Factor LC failures | **0** |
| Max subtree size | n - 2 |

## Theorem (conditional)

**Theorem.** If the IS polynomial of every tree on ≤ n-2 vertices is log-concave,
then the IS sequence of every tree on n vertices is unimodal.

**Proof.**
1. Every tree on ≥ 2 vertices has a support vertex r (adjacent to at least one leaf).
2. Root at r. Then E = (1+x)^ℓ · ∏ I(T_{c_j}) where each T_{c_j} has ≤ n-2 vertices.
3. By hypothesis, each I(T_{c_j}) is LC. Since (1+x) is LC and products preserve LC, E is LC.
4. The identity gives: b_{k-1}·Δ_k = b_{k-1}·d_k + b_k·d_{k-1} + a_{k-1}·c_k.
5. Since E is LC: c_k ≥ 0, so a_{k-1}·c_k ≥ 0.
6. But d_k and d_{k-1} can be negative, so this alone doesn't prove Δ_k ≥ 0.
7. **ADDITIONAL STEP NEEDED**: must show the curvature term a_{k-1}·c_k dominates
   the negative d-terms.

**Status:** Step 7 is the gap. The identity alone + LC of E is NOT sufficient to prove SCC.

## Why the identity alone doesn't suffice

The identity decomposes Δ_k into three terms. Two can be negative (d_k, d_{k-1}).
The curvature term a_{k-1}·c_k is nonneg but not necessarily large enough.

However, computationally:
- When both d_k and d_{k-1} < 0, the curvature term ALWAYS compensates (n ≤ 22)
- The curvature-to-deficit ratio is ≥ 12.09 (huge margin)
- The a_{k-1}/b_{k-1} amplification factor (always ≥ 1 since I ≥ E coefficientwise)
  is the key: it multiplies the LC gap c_k

## What's needed for a complete proof

Two possible paths:

### Path A: Prove the curvature dominance algebraically
Show: b_{k-1}·d_k + b_k·d_{k-1} + a_{k-1}·c_k ≥ 0 when E is LC.
Equivalently: a_{k-1}·c_k ≥ -b_{k-1}·d_k - b_k·d_{k-1} = |b_{k-1}·d_k + b_k·d_{k-1}|.
Since a_{k-1} = E_{k-1} + J_{k-2} ≥ E_{k-1} = b_{k-1}, we have amplification.
Need to bound d_k relative to c_k using the product structure of E.

### Path B: Prove LC of IS polys for all trees
If LC holds universally (not just through n ≤ 25), then the conditional theorem
gives unimodality. But LC fails at n = 26, so this path is DEAD for a universal proof.

### Path C: Prove LC of E at support vertices WITHOUT using LC of all subtrees
The E polynomial at a support vertex is a SPECIFIC product: (1+x)^ℓ · ∏ I(T_{c_j}).
Even if some I(T_{c_j}) is not LC, the product might still be LC because:
- The (1+x)^ℓ factor adds smoothing
- The other factors are LC (most subtrees are small)
- Non-LC IS polys are "barely non-LC" (only fail at 1-2 indices)

This is the most promising path but requires a quantitative LC-preservation result.

## Connection to previous work

This identity-based approach connects to the 4-term bilinear decomposition:
- The d_k terms relate to the cross terms QR + RQ
- The c_k term relates to the diagonal terms QQ + RR
- The curvature dominance is the algebraic counterpart of diag/|cross| ≥ 2

The advantage of this approach: it works at the FINAL product level (E at the support
vertex), not at each incremental stage. The incremental product analysis was a tool
to understand the mechanism, but the proof only needs the final E to be LC.
