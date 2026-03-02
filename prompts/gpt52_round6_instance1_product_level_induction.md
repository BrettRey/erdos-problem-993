# Task: Product-level induction for Strong Condition C

## Context from your Round 5 output

Your leaf-factor proof is clean and complete. Your LC brick wall observation is correct and critical: we cannot assume LC of the "denominator" E at intermediate stages, because arbitrary (non-LC) tree IS polynomials can appear as dp[v][0] via one-vertex attachment.

Your bivariate identity (3.5) isolates the obstruction cleanly:

F_{eP, bQ}(x,y) = F_{e,b}(x,y) · P(x)Q(y) + e(y)b(x) · F_{P,Q}(x,y)

where e = (1+x)I, b = E, P = I_c, Q = E_c. The correction term e(y)b(x)·F_{P,Q} is the problem.

## New result since your Round 5 output

### The factor-level curvature budget FAILS

Your proposed "tree-only bound" — that |L_{p,q}| ≤ Σ_{t=p+1}^{q} c_t for tree-derived factors — is **FALSE** starting at n=15. By n=19, the ratio exceeds 1:

| n | max |L_{p,q}| / Σ c_t |
|---|---------------------|
| ≤14 | 0.500 |
| 15 | 0.656 |
| 17 | 0.848 |
| 19 | **1.067** |

Concrete n=19 violation:
- P = [1, 18, 136, 565, 1417, 2211, 2124, 1188, 337, 34, 1]
- Q = [1, 17, 121, 470, 1087, 1530, 1284, 604, 145, 18, 1]
- L_{8,10} = Q_8·P_{10} - P_8·Q_{10} = 145·1 - 337·1 = -192
- Σ c_t (t=9,10) = 179 + 1 = 180
- ratio = 192/180 = 1.067

So the factor-level curvature alone cannot control the off-diagonal long minors.

### But product-level SCC holds universally

We verified: across all 9.1M trees n≤22, at every support vertex, at every intermediate product stage (89.9M total stages), Strong Condition C holds. Zero failures. Minimum nontrivial margin = 6 (at n=6), growing roughly as ⌊n/2⌋.

**So the control must come from the product level, not the factor level.**

## The question

In your kernel form (3.3):

Δ⁺_k = Σ_{u,v} e_{k-u} b_{k-v} · K_{u,v}

the weights e_{k-u} b_{k-v} come from the **accumulated product** (already-verified good pair). These weights are nonneg. The kernel K_{u,v} = p_{u+1}q_v - p_u q_{v+1} can be negative off-diagonal.

The product-level control must work through one of these mechanisms:

**(A) Weight concentration:** The e_{k-u}b_{k-v} weights might concentrate on indices (u,v) where K_{u,v} ≥ 0 (the diagonal and near-diagonal), and be small where K_{u,v} < 0 (far off-diagonal).

**(B) Product-level induction:** Instead of bounding K_{u,v} from factor data alone, use the fact that the accumulated pair (I, E) already satisfies Strong C. The old SCC gives Δ_k ≥ 0, which is a constraint on the e, b coefficients. These constraints might force the weighted sum to be nonneg.

**(C) Self-reinforcing induction:** The product-level SCC Δ⁺_k ≥ 0 might follow from the PRODUCT-LEVEL SCC at index k-1 (Δ⁺_{k-1} ≥ 0) plus factor-level data, giving a forward induction in k.

### What I want from you

1. **Investigate mechanism (B) or (C).** Can you formulate a product-level induction that avoids needing any factor-level curvature bound? The key hypothesis would be: "If the accumulated pair (I, E) satisfies Strong C and E is LC, and we multiply by a tree-derived factor (P, Q) with P = Q + xR, R ≤ Q, and Q is LC, then the new pair (IP, EQ) satisfies Strong C."

2. **Use the bivariate identity (3.5) in the PRODUCT direction.** In your identity, F_{eP,bQ} = F_{e,b}·PQ + eb·F_{P,Q}, the first term has diagonal coefficients related to Δ_k (old SCC), weighted by Cauchy products of P and Q. Can you bound the diagonal coefficients of the second term (eb·F_{P,Q}) using the old SCC + LC of E, without needing to bound K_{u,v} individually?

3. **Consider the recursion in k.** The SCC quantity Δ⁺_k at index k relates to Δ⁺_{k-1} via the shifted-minor recursion. At the product level, this recursion involves product-level curvature c_k^{prod} (which is ≥ 0 by Cauchy-Binet). Can you set up a forward induction: Δ⁺_0 ≥ 0 (trivially), and Δ⁺_k ≥ 0 implies Δ⁺_{k+1} ≥ 0 via the recursion + product curvature?

## Dead ends (do not pursue)

- Factor-level curvature budget |L_{p,q}| ≤ Σ c_t: FALSE at n≥15
- HWZZ partial synchronicity: FALSE at n=12+
- Generic product closure with J≤E + LC + Cond C: FALSE (counterexample exists)
- Global LC of I(T): FALSE at n=26

## Notation (consistent with your Round 5)

| Symbol | Definition |
|--------|-----------|
| e_k | [x^k](1+x)I = a_k + a_{k-1} |
| b_k | [x^k]E |
| Δ_k | e_{k+1}b_k - e_kb_{k+1} (Strong C target) |
| K_{u,v} | p_{u+1}q_v - p_u q_{v+1} (factor cross-minor) |
| c_k | b_k² - b_{k-1}b_{k+1} (LC gap of E) |
