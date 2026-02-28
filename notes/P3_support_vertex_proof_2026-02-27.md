# P3 Tail Domination at Support Vertices: Proof

**Status:** PROVED. Ready to write into paper when the time comes.

## Theorem (Support-Root Domination)

Let T be a tree on n ≥ 2 vertices, and let r be a vertex adjacent to at least one leaf ℓ. With the standard DP notation, E(x) = dp[r][0] and J(x) = dp[r][1]/x. Then

  e_k ≥ j_{k-1}  for all k ≥ 0.

In particular, P3 holds at this rooting for every choice of mode threshold m.

## Proof (Injection)

Define φ on size-k independent sets S with r ∈ S by

  φ(S) = (S \ {r}) ∪ {ℓ}.

Since ℓ is a leaf with N(ℓ) = {r}, and r ∉ φ(S), the set φ(S) is independent in T. The map is injective: given φ(S), recover S = (φ(S) \ {ℓ}) ∪ {r}. So the number of size-k IS containing r is at most the number avoiding r, i.e. j_{k-1} ≤ e_k for all k ≥ 1. At k = 0: e_0 = 1 ≥ 0 = j_{-1}. ∎

## Proof (Algebraic, via DP structure)

Let the children of r be ℓ = c_1, c_2, ..., c_d. For the leaf child: dp[ℓ][0] = 1, dp[ℓ][1] = x. Define

  A(x) = Π_{i=2}^d (dp[c_i][0] + dp[c_i][1])   (product of subtree IS polys, excluding leaf)
  B(x) = Π_{i=2}^d dp[c_i][0]                    (product of exclude-root polys)

Then:
  E(x) = (1 + x) · A(x)
  xJ(x) = x · B(x)

Since dp[c_i][0] + dp[c_i][1] ≥ dp[c_i][0] coefficientwise (dp[c_i][1] has nonneg coefficients), and products of coefficientwise-dominant nonneg polynomials preserve dominance: A(x) ≥ B(x) coefficientwise.

Therefore:
  E(x) - xJ(x) = (1+x)A(x) - xB(x) = A(x) + x(A(x) - B(x))

Both A(x) and A(x) - B(x) have nonneg coefficients, so E(x) - xJ(x) ≥ 0 coefficientwise. ∎

## Corollary

Every tree on ≥ 2 vertices has a leaf; the neighbor of any leaf is a support vertex. Rooting at that neighbor gives P3 for all k. Hence every nontrivial tree admits at least one rooting where P3 holds.

## Computational verification

- 9,114,283 trees n ≤ 22
- 59,916,124 support-vertex checks
- 0 P3 failures (as expected from the proof)
- P3 also holds at many non-support vertices, but fails at ~60% of all rootings

## Context

P3 is one of two conditions in the P⋆ invariant. Combined with P2 (prefix TP2, computationally verified but OPEN), P⋆ implies unimodality of I(T; x). The complete proof reduces to proving P2 at support vertices.

## Source

GPT 5.2 Pro Round 2, Instances 1 and 3 independently arrived at the same proof. Verified by direct computation.
