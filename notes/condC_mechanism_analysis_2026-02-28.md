> **OBSOLETE (2026-03-01):** SCC / Condition C fails at n=28. This proof direction is dead. See `notes/scc_false_n28_2026-03-01.md`.

# Condition C Product Closure Mechanism Analysis (2026-02-28)

## Setup

Three parallel investigations into the mechanism of Strong Condition C product closure.

## Investigation 1: Decomposition of product-level compensation

Script: `diagnose_condC_product_mechanism.py`
2,000 factor pairs from trees n ≤ 13; 2,001,000 pairwise products; 25.99M (pair,k) checks.

### Results

**0 FAILURES** across all products.

When d_k^{prod} < 0 (1,939,383 events):

| Compensator | Count | Fraction |
|-------------|-------|----------|
| Curvature (T3) alone | 1,542,631 | 79.5% |
| Both T2 + T3 | 396,752 | 20.5% |
| Memory (T2) alone | 0 | 0.0% |
| T1 nonneg | 0 | 0.0% |

**Key finding: curvature is the sole rescue mechanism.** T3 = a_{k-1} · c_k participates in 100% of rescues. Memory (T2 = b_k · d_{k-1}) never suffices alone.

Curvature fraction of positive load: min = 0.4718, median = 1.0000, p10 = 0.9691.

### Where d_k < 0 occurs

- Concentrated in the **tail** (k ≥ mode), peaking at k = mode + 4 to mode + 5
- The off-diagonal Cauchy-Binet terms drive negativity (diagonal contribution = 0 at tightest margins)
- Factor-level d_k < 0 only at the penultimate index

### Tightest margins

All top 20 tightest margins occur at k = deg(product) - 1. Minimum margin = 0.159 (subtree-11 factors, both identical). Margins shrink as subtree sizes grow.

At the tightest pair (k=11): T1 = -7552, T2 = -224256, T3 = +319600, total = +87792. Curvature provides 319600 against combined negative of 231808 from T1+T2.

### Implication for proof strategy

Since curvature (LC of E_prod) is the primary rescue, and LC is preserved under Cauchy products of nonneg sequences, the key algebraic question is: **does the curvature bonus a_{k-1}·c_k grow fast enough to compensate the negative T1 + T2 contributions from off-diagonal cross-terms?**

The a_{k-1}/b_{k-1} ≥ 1 amplification is essential (it multiplies the curvature bonus). Any proof should focus on:
1. Bounding the negative excursion of d_k in the tail (off-diagonal terms)
2. Showing the curvature bonus c_k from the LC product is large enough
3. Using a_{k-1} ≥ b_{k-1} to amplify c_k

## Investigation 2: Synthetic failure profile

Script: `profile_condC_failures.py`
200,000 synthetic (I, E) pairs; 200,000 pairwise products.

### Results

| Category | Tested | Failures | Rate |
|----------|--------|----------|------|
| Both J ≤ E | 50,000 | 2 | 0.004% |
| One J > E | 75,000 | 9 | 0.012% |
| Both J > E | 75,000 | 31 | 0.041% |
| **Total** | **200,000** | **42** | **0.021%** |

### Critical correction

**J ≤ E is necessary but NOT sufficient** for product closure in the generic (non-tree) case. The previous 0/125K claim held because those were tree-derived pairs, not random synthetic pairs satisfying the same abstract constraints.

### Failure characterization

- No sharp threshold on relative excess (max(J_k - E_k)/max(E_k))
- Failures scattered across all excess levels (0.0% to 0.04% across bins)
- Failing pairs have **lower** mean relative excess (1.27) than passing pairs (2.02) -- counterintuitive
- Failures concentrate at k = 7-9 (upper tail) and product degree 9-10
- Violation magnitudes range from -14 to -782955

### Closest call for J ≤ E products

Min margin = 0 at k=4 (factors I=[1,2,1], E=[1,1] and I=[1,4,11,8], E=[1,4,8]). Margin exactly 0.

### Implication

The tree DP structure (E_c = product of subtree IS polys, each being a unimodal LC polynomial of a specific form) provides constraints beyond LC + J ≤ E + unimodality. **A proof must exploit the tree product structure**, not just abstract polynomial properties.

## Investigation 3: Inductive consistency (partial)

The conceptual question: does the induction "close"? At a support vertex r:
- Factor pair: (I_c, E_c) = (IS_poly(subtree_c), dp0[c])
- These are NOT the same as the (A_c, B_c) pair at c's own support vertex

The mismatch: at c's own support vertex, one works with (A_c, B_c) = (∏ I_{gc}, ∏ E_{gc}) over non-leaf grandchildren. But at r's factor level, one uses (I_c, E_c) = ((1+x)^{ℓ_c} A_c + x B_c, (1+x)^{ℓ_c} A_c).

This means: the induction must prove Condition C for (I_c, E_c), not for (A_c, B_c). Since I_c = (1+x)^ℓ A_c + x B_c and E_c = (1+x)^ℓ A_c, the factor-level Condition C involves the FULL support-vertex structure (with the (1+x)^ℓ leaves folded in).

**This is actually fine for the induction:** at each level, the factor pair (I_c, E_c) has the same structure as the overall (I, E) at the root (it IS the IS poly / exclude poly at vertex c). The induction hypothesis is: "for all subtrees on < n vertices, the (IS, dp0) pair satisfies Condition C at every vertex." The inductive step uses this for the factors.

The key subtlety: Condition C at the factor level (for (I_c, E_c)) follows from: (a) Condition C holds for (A_c, B_c) by product closure, (b) the (1+x)^ℓ step only helps (adds smoothing).

**Status: logically consistent.** The induction closes if product closure and the (1+x)^ℓ step are proved.

## Updated proof structure

The proof of unimodality via Condition C requires proving exactly **two algebraic facts**:

1. **Product closure**: If (I_1, E_1) and (I_2, E_2) are factor pairs from tree DP (each satisfying Condition C, E_i is LC, I_i ≥ E_i coefficientwise, and the tree-specific structure), then (I_1·I_2, E_1·E_2) satisfies Condition C.

2. **The (1+x) step**: If (A, B) satisfies Condition C and E = (1+x)·A, J = B, then P2 holds (Δ_k ≥ 0 for all k). [For ℓ = 1, this is exactly the Delta identity. For ℓ ≥ 2, extra slack.]

Both are verified computationally. Neither is proved algebraically. The curvature-dominance finding (Investigation 1) suggests the proof should focus on bounding the LC bonus c_k relative to the negative off-diagonal d_k terms.
