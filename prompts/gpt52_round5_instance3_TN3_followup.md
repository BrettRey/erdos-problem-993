# Task: Construct a TN₃ / block-Toeplitz closure proof for Strong Condition C

## Context (what you already know)

You are continuing from your Round 4 output on Condition C product closure. Your counterexample was correct (verified computationally). Your diagnosis that tree-realizability is the essential missing constraint was spot-on. You offered to construct the TN₃/block-Toeplitz framework if given the exact formal product closure statement. Here it is.

## New results since your last output

### 1. Your counterexample is verified and correctly diagnosed

Factor 1: E₁ = [1,5,6,5,3], J₁ = [1,4,4], I₁ = [1,6,10,9,3]
Factor 2: E₂ = [1,5,4,2,1], J₂ = [1,3,4], I₂ = [1,6,7,6,1]
Product: Δ₇ = -6. CONFIRMED.

Both factors satisfy Strong Condition C, J ≤ E, E is LC, J is LC. But neither is tree-realizable (E₁ = [1,5,6,5,3] cannot be a product of IS polynomials of subtrees).

### 2. HWZZ partial synchronicity is FALSE for tree-derived (I, E) pairs

We scanned all trees n ≤ 18 (all rootings): 319 failures of the HWZZ partial sync condition a_m·b_n + a_n·b_m ≥ a_{m+1}·b_{n-1} + a_{n-1}·b_{m+1}. All failures are diagonal (m = n). First failure at n = 12.

**Partial synchronicity is dead as a product-closure route.** Tree-derived (I, E) pairs are NOT partially synchronized in Hu's sense.

### 3. GMTW gives the wrong number of (1+x) factors

Condition C is equivalent to (1+x)·A ≽ B (ratio dominance of (1+x)A over B). GMTW Cor. 2.21 says ratio dominance is preserved under convolution with an LC sequence. But applied factor by factor, GMTW gives (1+x)^s · ∏I_c ≽ ∏E_c (one (1+x) per factor), while we need (1+x) · ∏I_c ≽ ∏E_c (single (1+x)). The s-1 extra factors are wasted. GMTW doesn't apply directly.

### 4. Mechanism analysis (2M products, n ≤ 13 factors)

When d_k^{prod} < 0 at the product level (1.94M events out of 26M checks):
- **T3 (curvature = a_{k-1}·c_k) participates in 100% of rescues.** Alone in 79.5%.
- **T2 (memory = b_k·d_{k-1}) alone: 0%.** Never sufficient by itself.
- d_k < 0 concentrated in the tail (k ≥ mode), driven by off-diagonal Cauchy-Binet cross-terms.
- Tightest margin: 0.159 at k = 11 (subtree-11 factors). Margins shrink with subtree size.
- At the tightest point: T1 = -7552, T2 = -224256, T3 = +319600. Curvature provides 319600 against combined negative of 231808.

### 5. J ≤ E is necessary but NOT sufficient (even with LC + Cond C)

- Generic LC pairs without J ≤ E: 166 failures / 125K products (as you diagnosed).
- Generic LC pairs WITH J ≤ E + LC + Cond C on factors: 2 failures / 50K products.
- Tree-derived pairs: 0 failures / 701K exhaustive + 100K random + 2M mechanism products.

**So even J ≤ E + LC + Cond C is not enough for generic polynomials.** The tree product-of-IS-polys structure provides additional constraints that abstract axioms miss.

## The exact formal product closure statement

### Definition: Tree-realizable factor triple

A triple (I, E, J) of polynomials with nonneg integer coefficients is **tree-realizable** if there exists a rooted tree T with root r such that:
- E = dp[r][0] = ∏_{c child of r} I(subtree(c); x)
- J = dp[r][1]/x = ∏_{c child of r} dp[c][0]
- I = E + x·J

In other words, E is a product of subtree IS polynomials, and J is a product of the corresponding exclude-root polynomials, and the **same children c** govern both products. This is the "shared-factor" constraint.

### Properties of tree-realizable triples (all computationally verified)

1. **I = E + x·J** (by DP definition)
2. **E[0] = 1** (constant term)
3. **J ≤ E coefficientwise** (verified 61.8M subtree checks, 0 failures)
4. **E is log-concave** (verified 11.9M factor checks, 0 failures)
5. **Strong Condition C holds**: b_{k-1}·d_k + b_k·d_{k-1} + a_{k-1}·c_k ≥ 0 for all k ≥ 1, where a_k = [x^k]I, b_k = [x^k]E, d_k = a_{k+1}·b_k - a_k·b_{k+1}, c_k = b_k² - b_{k-1}·b_{k+1} (verified 11.9M factor checks, 0 failures)

### The product operation

Given two tree-realizable triples (I₁, E₁, J₁) and (I₂, E₂, J₂), their product is:
- E' = E₁ · E₂
- J' = J₁ · E₂ + E₁ · J₂ + x · J₁ · J₂
- I' = E' + x · J' = I₁ · I₂

This corresponds to attaching two independent subtrees as children of a common root.

**Note:** J' ≤ E' does NOT hold in general (even for leaf factors: J₁=J₂=1, E₁=E₂=1 gives J'=[2,1] > E'=[1,2,1] at degree 0). But Strong Condition C is still preserved.

### The product closure claim (to be proved)

**Theorem (Product Closure).** If (I₁, E₁, J₁) and (I₂, E₂, J₂) are triples satisfying properties 1-5 above, then the product triple (I₁·I₂, E₁·E₂, J') also satisfies Strong Condition C (property 5).

**Moreover**, E' = E₁·E₂ is LC (product of LC nonneg sequences), and I' = I₁·I₂ ≥ E₁·E₂ = E' coefficientwise (product of nonneg sequences preserves dominance). So properties 1, 2, 4 are automatic. Property 3 (J' ≤ E') fails but is not needed at the product level; it serves as a hypothesis on the factors.

**The hard part is: property 5 propagates.** This is what needs proving.

## Computational status

| Check | Result |
|-------|--------|
| Product closure, tree-derived pairs | 0 fails / 2.8M products |
| Product closure, generic J ≤ E + LC + Cond C pairs | 2 fails / 50K products |
| Product closure, generic without J ≤ E | 166 fails / 125K products |

**Conclusion:** Product closure holds for tree-realizable triples but NOT for the abstract axiom set {J ≤ E, LC, Cond C}. The proof must use the product-of-IS-polys structure.

## What I need from you

### Primary: TN₃ / block-Toeplitz construction

You suggested encoding each tree-realizable triple as a planar-network or block-Toeplitz object whose order-3 minors encode the Strong Condition C minors. The tree recursion (product of children) would correspond to composition/multiplication of these objects, and TN₃ closure under multiplication would give product closure.

**Please construct this explicitly:**

1. **The encoding.** Given a tree-realizable triple (I, E, J) with coefficients a_k, b_k, j_k, define the explicit matrix M(I, E, J) whose minors capture Strong Condition C.

   Strong Condition C in integer form is: b_{k-1}·d_k + b_k·d_{k-1} + a_{k-1}·c_k ≥ 0, which equals b_{k-1} · [(a_{k+1}+a_k)·b_k - (a_k+a_{k-1})·b_{k+1}] = b_{k-1} · Δ_k.

   Equivalently (when b_{k-1} > 0): the ratio sequence (a_k+a_{k-1})/(b_k) is nonincreasing, i.e., the polynomial (1+x)·I ratio-dominates E. This is a statement about the interleaved Toeplitz matrices of (1+x)I and E.

2. **The product map.** Show that M(I₁·I₂, E₁·E₂, J') can be expressed in terms of M(I₁, E₁, J₁) and M(I₂, E₂, J₂) (Cauchy product, Hadamard product, matrix multiplication, or other composition).

3. **TN₃ closure.** Invoke the appropriate TP/TN closure theorem to conclude that the relevant minors stay nonneg.

### Key structural constraints to use

- **Shared-factor constraint**: E and J are products of the SAME list of polynomials, just different ones (I_c vs E_c for each child c). This links the ratio structure of E and J.
- **J ≤ E on factors**: Each individual factor has J_c ≤ E_c. This constrains the "gap" I_c - E_c = x·J_c ≤ x·E_c.
- **Factor-level Condition C**: Each factor already satisfies Strong Condition C (the base of the induction).
- **E is a product of LC polys**: LC is preserved under nonneg convolution, so E is always LC. The LC gaps c_k provide the curvature rescue.
- **a_{k-1} ≥ b_{k-1}**: The curvature bonus has amplification factor a_{k-1}/b_{k-1} ≥ 1. This is essential (the weak version with bare c_k fails at n=17).

### Secondary: Alternative approaches

If TN₃ doesn't work out cleanly, consider:

1. **Direct Cauchy-Binet decomposition.** Write d_k^{prod} = Σ_{i+j=k} [diagonal + off-diagonal], bound the off-diagonal negativity using factor-level Condition C + J ≤ E, and show the curvature bonus a_{k-1}^{prod} · c_k^{prod} absorbs it. The mechanism analysis shows curvature handles 100% of rescues.

2. **Induction on number of factors.** Instead of two-factor closure, prove: if (I₁, E₁, J₁) satisfies properties 1-5 and (I₂, E₂, J₂) = (1+x, 1, 1) (leaf factor), then the product satisfies Condition C. Then iterate. The leaf factor is special: E₂ = 1 (LC trivially), J₂ ≤ E₂ (equality), d_k^{(2)} = 0 for all k. This might be tractable.

3. **Weighted Cauchy product identity.** Find an algebraic identity expressing b_{k-1}^{prod} · Δ_k^{prod} as a nonneg combination of factor-level quantities. This is what the binomial-minor expansion (from Instance 2) partially achieves.

## Dead ends (do not pursue)

- **HWZZ partial synchronicity**: FALSE for tree-derived (I, E) pairs starting at n = 12.
- **Factor-level ratio dominance**: FALSE in ~31% of factors. Cannot assume d_k ≥ 0.
- **Global LC of I(T)**: FALSE at n = 26. Cannot use.
- **Weak Condition C** (with b_{k-1} instead of a_{k-1}): FALSE at n = 17.
- **J ≤ E as product-closed invariant**: FALSE (J' > E' even for leaf factors).
- **3×3 determinant from Instance 1**: Trivially block-triangular (b_{k-1} × 2×2 minor). Not genuine TP₃.

## Notation summary

| Symbol | Definition |
|--------|-----------|
| I(T; x) | Independence polynomial of tree T |
| E = dp[r][0] | Exclude-root polynomial at root r |
| J = dp[r][1]/x | Include-root polynomial (shifted) |
| a_k = [x^k]I | Coefficients of I |
| b_k = [x^k]E | Coefficients of E |
| d_k = a_{k+1}·b_k - a_k·b_{k+1} | LR minors of I vs E |
| c_k = b_k² - b_{k-1}·b_{k+1} | LC gaps of E |
| Δ_k = (a_{k+1}+a_k)·b_k - (a_k+a_{k-1})·b_{k+1} | P2 target |
| Strong Condition C | b_{k-1}·d_k + b_k·d_{k-1} + a_{k-1}·c_k ≥ 0 |
