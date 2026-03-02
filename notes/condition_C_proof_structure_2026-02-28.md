> **OBSOLETE (2026-03-01):** SCC / Condition C fails at n=28. This proof direction is dead. See `notes/scc_false_n28_2026-03-01.md`.

# Condition C Proof Structure (2026-02-28)

## The Goal

Prove that the independence polynomial of every tree is unimodal.

## Reduction Chain

1. **P⋆ ⟹ unimodality**: If there exists a rooting where the ratio j_k/e_k is nonincreasing on the prefix (P2) and e_k ≥ j_{k-1} on the tail (P3), then the coefficient sequence is unimodal. (Standard convexity argument.)

2. **Support vertices suffice**: Every tree has at least one support vertex (vertex adjacent to a leaf). P3 is PROVED at all support vertices (leaf-swap injection). P2 holds at all support vertices (59.9M checks, 0 failures).

3. **P2 ⟹ Condition C**: At a support vertex r with ℓ leaf children:
   - E = (1+x)^ℓ · A, J = B
   - P2 requires Δ_k ≥ 0 for k < mode
   - The algebraic identity (PROVED): b_{k-1}·Δ_k = b_{k-1}·d_k + b_k·d_{k-1} + a_{k-1}·c_k
   - Condition C: this sum is ≥ 0 for all k ≥ 1

## The Identity

For k ≥ 1 with b_{k-1} > 0:

    b_{k-1} · Δ_k = b_{k-1} · d_k + b_k · d_{k-1} + a_{k-1} · c_k

where:
- Δ_k = (a_{k+1}+a_k)·b_k - (a_k+a_{k-1})·b_{k+1}  (the P2 target)
- d_k = a_{k+1}·b_k - a_k·b_{k+1}  (LR minor of A vs B)
- c_k = b_k² - b_{k-1}·b_{k+1}  (LC gap of B)

The identity is trivial algebra (sympy-verified). The content is in the sign of the RHS.

## What's Known

| Statement | Status |
|-----------|--------|
| Identity is correct | **PROVED** (algebra) |
| Δ_0 ≥ 0 | **PROVED** (Δ_0 = a_1-b_1+1 ≥ 1) |
| Condition C at root level | **VERIFIED** 59.9M checks, 0 fails (n≤22) |
| Condition C at factor level | **VERIFIED** 11.9M checks, 0 fails (n≤20) |
| Product closure of Condition C | **VERIFIED** 701K exhaustive + 100K random, 0 fails |
| Factor E_c is always LC | **VERIFIED** 11.9M checks, 0 fails |
| Factor-level ratio dominance | **FALSE** (~31% of factors) |

## What Needs Proving

### Condition C for the base case

At a leaf, I = 1+x, E = 1. d_0 = 1, d_k = 0 for k ≥ 1. Condition C is trivially satisfied.

At a subtree of depth 1 (star): A and B are both products of (1+x) factors. Everything is nonneg.

### Product closure (THE KEY STEP)

**Claim:** If (I_1, E_1) and (I_2, E_2) each satisfy:
- I_i = E_i + x·J_i with J_i ≤ E_i coefficientwise (tree DP structure)
- E_i is log-concave
- Condition C holds for (I_i, E_i)

Then (I_1·I_2, E_1·E_2) also satisfies Condition C.

**IMPORTANT:** The J ≤ E constraint is NOT itself product-closed (J' = J_1E_2 + E_1J_2 + xJ_1J_2 > E_1E_2 even for leaf factors). But Condition C still holds at the product level. The J ≤ E constraint is needed as a HYPOTHESIS on the factors, not as an INVARIANT to maintain. In the tree DP, J_c ≤ E_c always holds at each individual subtree (verified 61.8M checks, n≤18).

Strong Condition C fails for generic (I, E) pairs without J ≤ E (166 failures / 125K). With J ≤ E on factors: 0 failures (125K synthetic + 701K tree pairs).

This is the hard algebraic step. The Cauchy product of the d_k and c_k sequences involves cross terms. One approach (from the GPT 5.2 prompt): rewrite Condition C as nonnegativity of a 3×3 minor of a block Toeplitz matrix, then use the fact that TP is preserved under Cauchy products.

### The (1+x)^ℓ step

After Condition C is established for (A, B), P2 follows because E = (1+x)^ℓ · A, and the (1+x)^ℓ factor only helps (adds more smoothing). For ℓ = 1, the identity gives exactly Condition C. For ℓ ≥ 2, there's extra slack.

## Proof Outline (if product closure is established)

**Theorem.** The independence polynomial of every tree is unimodal.

*Proof.* By strong induction on n. Let T be a tree on n ≥ 2 vertices.

1. T has a support vertex r (standard fact).
2. At r, P3 holds (leaf-swap injection, proved).
3. For P2: let the non-leaf children be c_1, ..., c_s. Each gives a factor pair (I_{c_i}, E_{c_i}).
4. By induction (applied to the subtrees rooted at c_i), each factor pair satisfies Condition C, I_{c_i} ≥ E_{c_i}, E_{c_i} is LC, and J_{c_i} ≤ E_{c_i} (from tree DP).
5. By product closure (using J_c ≤ E_c on factors), (A, B) = (∏I_{c_i}, ∏E_{c_i}) satisfies Condition C.
6. The algebraic identity gives Δ_k ≥ 0 for k ≥ 1. The k=0 case is trivial.
7. P2 holds at r. Combined with P3, the IS polynomial is unimodal. □

**Note:** Step 4 requires Condition C to hold for the (I_c, E_c) pair where I_c = dp0[c] + x·dp1s[c] and E_c = dp0[c]. This is the (A', B') pair at vertex c, which has the same structure -- so the induction applies.

## The a_{k-1}/b_{k-1} ratio is essential

The WEAKER version of Condition C (replacing a_{k-1} with b_{k-1}) FAILS at n=17 (3 failures by n=18). So the coefficientwise dominance A ≥ B provides a nontrivial amplification of the curvature bonus that is strictly necessary. Any proof must use A ≥ B, not just the structure of B alone.

## Open Questions

1. Can product closure be proved via total positivity of a block Toeplitz matrix?
2. Is there a simpler direct proof that avoids products (e.g., using the tree structure directly)?
3. Does the GMTW framework (Cor. 2.21: convolution with LC preserves ratio dominance) help?
4. J ≤ E is NOT product-closed. How does it enter the product closure proof of Condition C? The factors satisfy J_c ≤ E_c, which constrains the d_k and c_k cross-terms. Possibly: J ≤ E implies d_k ≥ -(k-th coeff of x·J·E_tail) or similar bound on the negative excursion of d_k.
