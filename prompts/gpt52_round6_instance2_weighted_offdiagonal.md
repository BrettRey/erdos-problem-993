# Task: Weighted off-diagonal control via shifted-minor recursion

## Context from your Round 5 output

Your Cauchy-product decomposition of d_k and the Cauchy-Binet expansion of c_k^{prod} are clean. Your framework correctly identifies the off-diagonal contribution Off(k) as the sole obstruction. Your proposed "tree-only bound" — controlling long minors |L_{p,q}| by telescoping sums of curvature gaps via the shifted-minor recursion — was the right idea.

## New result that invalidates the naive bound

The simple bound |L_{p,q}| ≤ Σ_{t=p+1}^{q} c_t is **FALSE** for tree-derived factors starting at n=15:

| n | max |L_{p,q}| / Σ c_t | Violation tree family |
|---|---------------------|-----------------------|
| ≤14 | 0.500 | (none) |
| 15 | 0.656 | deg seq [4,3,3,3,2,2,...] |
| 17 | 0.848 | same family |
| 19 | **1.067** | same family, L_{8,10}=-192, Σc=180 |

The violation always occurs at span 2 (q-p=2), at the tail (near the independence number), in subtrees of size n-1. Zero "infinite ratio" cases (curv_sum is always positive when L < 0).

But: **product-level Strong Condition C holds at all 89.9M intermediate stages (n≤22), zero failures.** So the control comes from the product level, not the factor level.

## The shifted-minor recursion at the product level

Your recursion (Rec) is: b_{k-1} M_k^{(t)} = b_k M_{k-1}^{(t-1)} + a_{k-t} c_k.

Applied at the PRODUCT level (with product coefficients a_k^{prod}, b_k^{prod}, c_k^{prod}), this recursion telescopes the product-level Δ_k into curvature contributions. The key advantage: product-level curvature c_k^{prod} is given by the Cauchy-Binet expansion (your eq. (2)), which aggregates curvature from ALL factors.

### The question

**Can the shifted-minor recursion at the product level absorb the off-diagonal debt, even though the factor-level recursion cannot?**

Specifically, in the two-factor decomposition

d_k^{prod} = D_k^{diag} + D_k^{off}

you showed that the diagonal part D_k^{diag} can be regrouped via factor-level SCC. The off-diagonal part D_k^{off} involves long minors of the factors.

Now apply the PRODUCT-LEVEL recursion to d_k^{prod}:

b_{k-1}^{prod} · M_k^{(1),prod} = b_k^{prod} · d_{k-1}^{prod} + a_{k-1}^{prod} · c_k^{prod}

This gives:
b_{k-1}^{prod} · Δ_k^{prod} = b_{k-1}^{prod} · d_k^{prod} + b_k^{prod} · d_{k-1}^{prod} + a_{k-1}^{prod} · c_k^{prod}

The curvature bonus a_{k-1}^{prod} · c_k^{prod} involves:
- a_{k-1}^{prod} ≥ b_{k-1}^{prod} (amplification from I ≥ E)
- c_k^{prod} = Σ (2×2 Toeplitz minors from both factors) (Cauchy-Binet, all nonneg)

### What I want from you

1. **Quantify the product curvature bonus.** Can you show that a_{k-1}^{prod} · c_k^{prod} ≥ |b_{k-1}^{prod} · d_k^{prod} + b_k^{prod} · d_{k-1}^{prod}| whenever the LHS has negativity from off-diagonal terms? The Cauchy-Binet sum for c_k^{prod} has terms from both factors; the d_k^{prod} negativity comes from off-diagonal cross-minors. Can you match them term-by-term?

2. **Iterated recursion at the product level.** Instead of applying Rec once (t=0→1, giving SCC), apply it multiple times (t=0→1→2→...→s). Each step harvests one more c_k^{prod} contribution. For the long minors of the product, the recursion depth needed is exactly the span (q-p). Does iterating the product-level Rec telescope the off-diagonal debt into a nonneg sum of product curvature terms?

3. **Cauchy-Binet matching.** In your expansion:

   c_k^{prod} = Σ_{0≤u<v} (2×2 minor of T(b^{(1)}))(2×2 minor of T(b^{(2)}))

   And the off-diagonal contribution from factor 1 involves L^{(1)}_{p,q} weighted by factor-2 coefficients. Can you show that the specific Cauchy-Binet term involving u=p, v=q in the c_k^{prod} sum dominates the corresponding L_{p,q} contribution in D_k^{off}?

## Key computational data

For the **product-level** SCC (not factor-level):
- 89,888,083 intermediate stages checked (all trees n≤22)
- Zero failures
- Min nontrivial margin = 6 (at n=6, g6=E?qo)
- Margin grows ~⌊n/2⌋, stabilizes quickly
- Worst-case intermediate stages are STARS K_{1,n-1} (margin = 2 at stage 1)

The product curvature c_k^{prod} is always large enough to absorb the off-diagonal negativity at the product level, even when factor-level curvature is insufficient.

## Dead ends

- Factor-level |L_{p,q}| ≤ Σ c_t: FALSE at n≥15. Do not pursue.
- HWZZ partial synchronicity: FALSE at n=12+
- Generic product closure: FALSE (non-tree-realizable counterexample)

## Notation (matching your Round 5)

| Symbol | Definition |
|--------|-----------|
| d_k | a_{k+1}b_k - a_kb_{k+1} (adjacent minor of (I,E)) |
| M_k^{(t)} | a_{k+1-t}b_k - a_{k-t}b_{k+1} (shifted minor) |
| c_k | b_k² - b_{k-1}b_{k+1} (LC gap of E) |
| L_{p,q} | b_pa_q - a_pb_q (long minor, same as your D_k^{off} source) |
| Δ_k | (a_{k+1}+a_k)b_k - (a_k+a_{k-1})b_{k+1} = d_k + M_k^{(1)} |
| S_k | b_{k-1}Δ_k = b_{k-1}d_k + b_kd_{k-1} + a_{k-1}c_k (Strong C) |
