# Condition C: the formulation and what needs proving

Here is the complete self-contained formulation. Ignore the earlier SV1 prompt (SV1 is false).

## The problem

Prove that the independence polynomial of every tree is unimodal.

## Reduction to P2 at support vertices

Every tree has a support vertex r (a vertex adjacent to at least one leaf). At r with ℓ leaf children:

    E(x) = (1+x)^ℓ · A(x),    J(x) = B(x),
    I(T; x) = E(x) + x·J(x)

where A = ∏_c I_c and B = ∏_c E_c, products over the non-leaf children c of r. Here I_c = IS polynomial of subtree(c), E_c = dp0 (root-excluded polynomial of subtree(c)).

P3 (tail domination, e_k ≥ j_{k-1} for k ≥ mode) is PROVED at all support vertices.

P2 (prefix ratio dominance) requires: for all k = 0, ..., m-1 where m = mode(I(T)):

    Δ_k := e_{k+1}·j_k - e_k·j_{k+1} ≥ 0

With E = (1+x)^ℓ · A and J = B, and writing a_k = [x^k]A, b_k = [x^k]B, for ℓ = 1 this becomes:

    Δ_k = (a_{k+1} + a_k)·b_k - (a_k + a_{k-1})·b_{k+1}

## The key identity (PROVED, trivial algebra)

Define:
- d_k := a_{k+1}·b_k - a_k·b_{k+1}  (LR minors of A vs B)
- c_k := b_k² - b_{k-1}·b_{k+1}      (LC gaps of B)

Then for k ≥ 1 with b_{k-1} > 0, clearing denominators:

    b_{k-1}·Δ_k = b_{k-1}·d_k + b_k·d_{k-1} + a_{k-1}·c_k

This is the **integer form**. The three terms are:
1. T1 = b_{k-1}·d_k: current LR minor (can be negative)
2. T2 = b_k·d_{k-1}: memory from previous minor (can be negative)
3. T3 = a_{k-1}·c_k: curvature bonus from B's log-concavity (always ≥ 0 when B is LC)

Note: a_{k-1} ≥ b_{k-1} (since A ≥ B coefficientwise), so T3 ≥ b_{k-1}·c_k. This amplification is ESSENTIAL.

## Strong Condition C

**Strong Condition C for (I, E):** For all k ≥ 1 with b_{k-1} > 0:

    b_{k-1}·d_k + b_k·d_{k-1} + a_{k-1}·c_k ≥ 0

where a_k = [x^k]I, b_k = [x^k]E, d_k = a_{k+1}·b_k - a_k·b_{k+1}, c_k = b_k² - b_{k-1}·b_{k+1}.

By the identity, Strong Condition C for (A, B) implies Δ_k ≥ 0, hence P2.

**WARNING:** There is a WEAKER version that replaces a_{k-1} with b_{k-1} in the curvature term. The weak version FAILS at n=17 (3 failures). You MUST use the strong version with a_{k-1}.

## Computational verification

| Check | Result |
|-------|--------|
| Strong Cond C at support vertices | 0 fails / 59.9M checks (all 9.1M trees n ≤ 22) |
| Strong Cond C at factor level (I_c, E_c) | 0 fails / 11.9M checks (all trees n ≤ 20) |
| Factor-level ratio dominance (d_k ≥ 0) | FAILS in ~31% of factors |
| Product closure (exhaustive, n ≤ 12) | 0 fails / 701K pairwise products |
| Product closure (random, n ≤ 15) | 0 fails / 100K products |
| Factor E_c is LC | 0 fails / 11.9M checks |

## Mechanism analysis (2M products tested)

When d_k^{prod} < 0 at the product level (1.94M events out of 26M checks):
- **T3 (curvature) participates in 100% of rescues.** Alone in 79.5%. Memory alone: 0%.
- d_k < 0 concentrated in the tail (k ≥ mode), driven by off-diagonal Cauchy-Binet cross-terms.
- Tightest margin = 0.159 (subtree-11 factors). Margins shrink with subtree size.

## What needs proving

### Step 1: Product closure of Strong Condition C

**Claim:** If (I_1, E_1) and (I_2, E_2) are factor pairs from the tree DP, each satisfying:
- I_i = E_i + x·J_i with J_i ≤ E_i coefficientwise
- E_i is log-concave
- Strong Condition C holds for (I_i, E_i)

Then (I_1·I_2, E_1·E_2) satisfies Strong Condition C.

**Key constraints on this step:**
- J ≤ E is NOT product-closed (J' = J_1E_2 + E_1J_2 + xJ_1J_2 exceeds E_1E_2 even for leaf factors). But it holds at each individual factor from tree DP (verified 61.8M checks).
- Strong Condition C is NOT product-closed for generic polynomial pairs (166 failures / 125K synthetic pairs without J ≤ E). With J ≤ E on factors: 0 failures for tree-derived pairs. BUT even with J ≤ E, 2 failures / 50K for random synthetic LC polynomials. So the tree product-of-IS-polys structure provides constraints beyond J ≤ E + LC.
- A proof must exploit the tree product structure, not just abstract polynomial axioms.

### Step 2: The (1+x)^ℓ step

For ℓ = 1: the identity gives Δ_k = (1/b_{k-1}) · [Strong Cond C expression], so Cond C directly implies P2.
For ℓ ≥ 2: extra smoothing provides additional slack. 95.6% of trees have max ℓ ≥ 2.

### The induction

By strong induction on n:
1. T has a support vertex r.
2. P3 holds at r (proved).
3. Each non-leaf child c gives a factor pair (I_c, E_c) with the same structure.
4. By induction, each factor satisfies Cond C + E_c is LC + J_c ≤ E_c.
5. By product closure, (A, B) = (∏I_c, ∏E_c) satisfies Cond C.
6. By the identity, Δ_k ≥ 0. P2 holds. Combined with P3, unimodality follows.

## Possible approaches to product closure

### Approach A: Total positivity / Toeplitz-block formulation

The three terms in Strong Condition C involve 2×2 minors:
- d_k is a 2×2 minor of the matrix [A; B]
- c_k is a 2×2 minor of the Toeplitz matrix of B

If Condition C can be written as nonnegativity of a 3×3 minor of a block matrix, then TP_3 (total positivity of order 3) is preserved under Cauchy convolution, giving product closure for free.

**Can you find the explicit 3×3 minor formulation?**

### Approach B: Direct Cauchy product analysis

The product-level LR minor:

    d_k^{prod} = Σ_{i+j=k} [a1_{i+1}·a2_{j+1}·b1_i·b2_j - a1_i·a2_j·b1_{i+1}·b2_{j+1}]

This is NOT Σ d_i^{(1)}·d_j^{(2)} (there are cross-terms). But the cross-terms involve the factor-level d_k, c_k, and the J ≤ E constraint. Can you bound the cross-term negativity using factor-level Condition C?

### Approach C: GMTW / Hu-Wang-Zhao-Zhao frameworks

- Gross-Mansour-Tucker-Wang: ratio dominance preserved under convolution with LC sequence (Cor. 2.21). But factor-level ratio dominance FAILS.
- Hu-Wang-Zhao-Zhao (arXiv:1507.08430): partial synchronicity preserved under convolution for LC sequences. Potentially relevant but not directly applicable.

## Hard constraints

- Do NOT assume global log-concavity of I(T). Fails at n = 26.
- Do NOT assume A ratio-dominates B (d_k ≥ 0). False in ~31% of factors.
- Products preserve: nonnegativity, I ≥ E coefficientwise, log-concavity of E.
- Products do NOT preserve: ratio dominance, J ≤ E, the weak version of Condition C.
