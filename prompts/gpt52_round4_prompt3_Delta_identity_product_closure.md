# Task: Prove Condition C (local deficit bound) and its product closure

## Context

I am working on Erdos Problem #993: prove that the independence polynomial of every tree is unimodal. P2 (prefix ratio dominance at support vertices) is the remaining open piece. It has been verified computationally on 59.9M support-vertex checks (all trees n <= 22), zero failures.

**Your task: prove the local deficit bound (Condition C) for the (A, B) pair at support vertices, and establish that it is preserved under polynomial products.**

## Setup

At a support vertex r with leaf child u:

    E(x) = (1 + x) * A(x),    J(x) = B(x),

where A = prod_{c != u} I_c(x), B = prod_{c != u} E_c(x).

P2 requires Delta_k >= 0 for k < m = mode(I(T)), where:

    Delta_k := (a_{k+1} + a_k) * b_k - (a_k + a_{k-1}) * b_{k+1}

## The key identity

Define:
- d_k := a_{k+1} * b_k - a_k * b_{k+1}  (LR minors of A vs B; measures ratio dominance)
- c_k := b_k^2 - b_{k-1} * b_{k+1}       (LC gaps of B; nonneg when B is LC)

Then for k >= 1 with b_{k-1} > 0:

    Delta_k = d_k + (b_k / b_{k-1}) * d_{k-1} + (a_{k-1} / b_{k-1}) * c_k    (dagger)

This decomposes the P2 target into three terms:
1. d_k: the "current" ratio-dominance minor (can be negative)
2. (b_k/b_{k-1}) * d_{k-1}: a "memory" term from the previous minor
3. (a_{k-1}/b_{k-1}) * c_k: a "curvature bonus" from B's log-concavity

Since A >= B coefficientwise, we have a_{k-1}/b_{k-1} >= 1 (when b_{k-1} > 0), so the curvature bonus is at least c_k.

## Condition C

**Condition C(A, B; M):** For all 1 <= k <= M:

    d_k + (b_k / b_{k-1}) * d_{k-1} + c_k >= 0

(using the weaker bound a_{k-1}/b_{k-1} >= 1).

By identity (dagger), Condition C implies Delta_k >= 0, hence P2.

**Computational status (COMPLETED):**
- Condition C verified: 0 failures across 59,916,124 support-vertex checks, ALL 9.1M trees n <= 22.
- Factor-level Condition C: 0 failures across 11,941,358 factor pairs (I_c, E_c), all trees n <= 20.
- Factor-level ratio dominance (d_k^(c) >= 0): FAILS in ~31% of factors. But Condition C is STRICTLY STRONGER and always holds.
- Product closure: 701,520 exhaustive pairwise products of 1,184 unique factors (n <= 12), 0 Condition C failures. 100,000 random products of 19,535 unique factors (n <= 15), 0 failures.
- Factor-level LC of E_c: 0 failures (all factors are LC).

## What I need from you

### Part 1: Prove the identity (dagger)

Give a clean algebraic proof that

    Delta_k = d_k + (b_k / b_{k-1}) * d_{k-1} + (a_{k-1} / b_{k-1}) * c_k

for k >= 1 with b_{k-1} > 0. (This should be straightforward algebra.)

### Part 2: Product closure of Condition C

This is the hard part. A and B arise as products:
- A = prod_c I_c where each I_c has nonneg coefficients
- B = prod_c E_c where each E_c has nonneg coefficients
- I_c >= E_c coefficientwise (for each factor c)
- Each E_c is log-concave (empirically always true; product of subtree polynomials)

The question is: **if each factor pair (I_c, E_c) satisfies some local condition, does the product pair (A, B) satisfy Condition C?**

Specifically, define for each factor:
- d_k^{(c)} := (I_c)_{k+1} * (E_c)_k - (I_c)_k * (E_c)_{k+1}
- c_k^{(c)} := (E_c)_k^2 - (E_c)_{k-1} * (E_c)_{k+1}

When we multiply two pairs (I_1, E_1) and (I_2, E_2) to get (I_1*I_2, E_1*E_2), the LR minors of the product involve the Cauchy product:

    d_k^{prod} = sum_{i+j=k} [(I_1)_{i+1}(I_2)_{j+1}(E_1)_i(E_2)_j - (I_1)_i(I_2)_j(E_1)_{i+1}(E_2)_{j+1}]
    (this is NOT simply d_i^{(1)} * d_j^{(2)})

The interaction terms make direct product closure nontrivial. But there is hope:

**GMTW framework (Gross-Mansour-Tucker-Wang):** They proved that ratio dominance IS preserved under convolution with a log-concave sequence (their Cor. 2.21). Since our factors E_c are LC, and ratio dominance of the factors (I_c ≽ E_c) may hold at the factor level...

Wait, but I_c ≽ E_c is FALSE in general (the 3-vertex star K_{1,2} has I = 1+3x+x^2, E = 1+2x+x^2, and d_1 = 3*1 - 1*2 = 1 > 0 but one should check all examples). Actually it might be true at the factor level. Let me know what you find.

**Hu-Wang-Zhao-Zhao:** Proved that "partial synchronicity" is preserved under convolution for LC sequences. Their partial synchronicity involves the two-index inequality a_m*b_n + a_n*b_m >= a_{m+1}*b_{n-1} + a_{n-1}*b_{m+1}. This is stronger than ratio dominance but weaker than log-concavity of the interleaved sequence.

### Part 3: Factor-level analysis

For each factor pair (I_c, E_c) where c is a non-leaf child of r:
- I_c = E_c + x * J_c (DP decomposition at c)
- E_c = product of (I_{grandchild}) for grandchildren of r through c
- J_c = product of (E_{grandchild})

So at the factor level, (I_c, E_c) has the SAME structure as (A, B) at the global level, just one level down in the tree. This is what makes induction plausible.

**ANSWERS (from computation):**
1. Does I_c ratio-dominate E_c? **NO.** Fails in ~31% of factors (e.g., P_3: I=1+3x+x^2, E=1+2x+x^2, d_1=-1).
2. N/A -- factor-level ratio dominance fails.
3. **The right factor-level condition IS Condition C itself.** It holds for every factor pair (0 failures, 11.9M checks) and appears to be preserved under products (0 failures, 701K exhaustive pairs + 100K samples).

**So the inductive structure is: Condition C at factors ⟹ (via product closure) Condition C at the product ⟹ (via identity) P2.**

### Part 4: Alternative - Toeplitz-block TP formulation

The identity (dagger) involves adjacent LR minors d_k, d_{k-1} and the LC gap c_k. These are all 2x2 minors of small matrices:

- d_k is a 2x2 minor of the matrix [A; B] (rows are coefficient sequences)
- c_k is a 2x2 minor of the Toeplitz matrix of B

Condition C says: a specific linear combination of these 2x2 minors is nonneg. This looks like it could be rewritten as nonnegativity of a 3x3 minor of an explicitly constructed block matrix.

If so, the condition "all such 3x3 minors are nonneg" is a total positivity condition (TP_3), and TP_r is preserved under convolution of Toeplitz matrices. This would give product closure for free.

**Can you find the explicit 3x3 minor formulation of Condition C?** Even if it requires additional positivity conditions on the entries, this would connect the problem to a well-developed theory.

## CRITICAL: Weak vs Strong Condition C

There are two versions. You MUST work with the STRONG version:

- **Weak (WRONG):** d_k + (b_k/b_{k-1}) * d_{k-1} + c_k >= 0  (bare c_k)
- **Strong (CORRECT):** b_{k-1}*d_k + b_k*d_{k-1} + a_{k-1}*c_k >= 0  (integer form with a_{k-1})

The weak version FAILS at n=17. The strong version uses a_{k-1} >= b_{k-1} (from A >= B) to amplify the curvature bonus. This amplification is ESSENTIAL.

## CRITICAL: J <= E is NOT product-closed

In the tree DP, I_c = E_c + x*J_c where J_c <= E_c coefficientwise (VERIFIED at every single subtree, 61.8M checks). But under the product operation:

J' = J_1*E_2 + E_1*J_2 + x*J_1*J_2

This EXCEEDS E_1*E_2 even for leaf factors (J=E=1 gives J'=[2,1] > E'=[1]).

So J <= E is a FACTOR-LEVEL HYPOTHESIS, not a product-closed invariant. The product closure claim for Condition C uses J_c <= E_c on the inputs but does NOT require J' <= E' on the output.

## Constraints

- Do NOT assume global log-concavity of I(T). Fails at n=26.
- Do NOT assume A ≽ B (ratio dominance). It's false in general.
- Products preserve: nonnegativity, coefficientwise dominance (I >= E), log-concavity.
- Products do NOT preserve: ratio dominance, interlacing, J <= E.
- The mode m depends on I(T) = (1+x)A + xB, not on A or B alone.
- Strong Condition C is NOT product-closed for generic (I,E) pairs. It IS product-closed when J <= E holds on the factors (0 fails, 701K tree pairs + 125K synthetic pairs).
