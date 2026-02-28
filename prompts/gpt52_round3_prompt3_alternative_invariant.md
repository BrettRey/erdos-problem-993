# Task: Find an inductive invariant for tree IS polynomial unimodality

## Context

I am working on Erdos Problem #993: prove the independence polynomial of every tree is unimodal. I have reduced the problem to a ratio-dominance condition (P2) at support vertices, which holds computationally but resists a direct proof. I want you to step back and consider whether a different inductive invariant might work better.

## What I have

### Proved:
- **P3 (tail domination):** For every tree T and every support vertex r (adjacent to a leaf u), the injection S -> (S\{r}) union {u} gives e_k >= j_{k-1} for all k. This means: among IS of size k, at least half avoid the root. Holds for ALL k, not just k >= mode.

### Verified but unproved:
- **P2 (prefix ratio dominance):** e_{k+1}*j_k >= e_k*j_{k+1} for k = 0,...,m-1. This says j_k/e_k is nonincreasing up to the mode. Holds at every support vertex of every tree n <= 22 (59,916,124 checks, 0 failures). Fails at non-support vertices.

### The algebraic structure at support vertices:
Root at r with leaf child u. Then:
  E = (1+x)*A, where A = Product of I(T_c) over non-leaf children c
  J = B, where B = Product of E_c over non-leaf children c

P2 becomes: (1+x)A is ratio-dominant over B up to mode.

Known: B is always fully log-concave. A >= B (ratio-dominant) fails at ~10% of support vertices, but (1+x)A >= B always holds.

### The obstacle:
The binomial upgrade lemma (A >= B + LC of B implies (1+x)A >= B) has hypotheses that are too strong. A >= B fails. No weaker sufficient condition is known that covers all cases.

## What I want you to explore

Maybe P2 is not the right inductive invariant. The DP for trees has the form:

  At a leaf v: dp[v][0] = 1, dp[v][1] = x.
  At internal v with children c_1,...,c_d:
    dp[v][0] = Product_i (dp[c_i][0] + dp[c_i][1])
    dp[v][1] = x * Product_i dp[c_i][0]

The IS polynomial is I(T) = dp[root][0] + dp[root][1].

I need an invariant on the pair (dp[v][0], dp[v][1]) that:
1. Holds at leaves (base case)
2. Is preserved by the product-and-sum steps of the DP
3. Implies unimodality of dp[v][0] + dp[v][1] at the top level

P* (P2 + P3) satisfies (1) and (3) but I can't prove (2).

## Specific directions to explore

### Direction 1: Interlacing or common-interlacing invariants

Two polynomials f, g with positive coefficients "interlace" if their zeros alternate on the negative real axis. If f and g interlace, then f + g is log-concave (hence unimodal) and f/g has monotone coefficients.

IS polynomials of trees are NOT always real-rooted, so global interlacing fails. But maybe a WEAKER notion works:
- "Prefix interlacing" (zeros of truncations interlace)?
- "Approximate interlacing" (zeros are close to alternating)?
- Common interlacing: there exists h that interlaces both f and g?

If dp[v][0] and dp[v][1] "weakly interlace" in some sense preserved by the DP, that might imply unimodality of their sum.

### Direction 2: Ratio-monotonicity invariants

Instead of P2 (j_k/e_k nonincreasing), maybe the right invariant is:

  e_{k+1}/e_k >= j_k/j_{k-1}    for k in some range

This says E's growth ratios dominate J's growth ratios. If both E and J are LC (so their ratios are individually nonincreasing), this is related to synchronicity.

Or maybe the invariant should be on the ratios within each polynomial:

  For E = dp[v][0]: the ratio e_{k+1}/e_k is nonincreasing (i.e., E is LC)
  For J = dp[v][1]/x: some condition on its LC or on E/J

Products of LC sequences preserve LC (GMTW), so if we could maintain LC of both E and J through the DP, plus some cross-condition, we might get unimodality.

### Direction 3: Entropy/information-theoretic invariants

The hard-core model at fugacity lambda on a tree gives a probability distribution on IS. The partition function Z(T; lambda) = I(T; lambda). Properties like:
- The variance of |S| under the hard-core measure
- The Fisher information of the size distribution
- Log-concavity of the size distribution (which IS unimodality)

Trees have special structure in the hard-core model: belief propagation is exact, and marginals can be computed by message-passing. Maybe an invariant on the BP messages (which correspond to the ratios dp[v][1]/dp[v][0]) is the right object.

### Direction 4: Total positivity (TP) invariants

A pair of sequences (a_k, b_k) is TP2 if a_{k+1}b_k >= a_kb_{k+1} for all k. This is ratio dominance. Higher-order TP conditions (TP3, etc.) involve determinantal inequalities on larger minors.

Maybe the right invariant is a TP condition on the 2×(alpha+1) matrix:

  M = [ e_0  e_1  e_2  ...  e_alpha ]
      [ j_0  j_1  j_2  ...  j_alpha ]

Being TP2 on the prefix {0,...,m} is exactly P2. But maybe a weaker condition (like "all 2×2 minors involving non-adjacent columns are nonneg") suffices for unimodality and is better preserved?

### Direction 5: Completely different decomposition

Instead of the root-based DP (exclude root / include root), consider:
- Edge contraction/deletion decomposition
- Vertex elimination ordering (not necessarily rooted)
- Subdivision-based decomposition (splitting at degree-2 vertices)

I have extensive work on a degree-2 bridge decomposition (see below) that led to a different set of conditions. Maybe there's a decomposition where the invariant IS nicely preserved.

## Background on my other proof attempts (for context)

### PNP framework:
- Everything reduces to "Conjecture A": for trees with all leaf-depths >= 2, the mode is at most floor(n/3) + 1.
- Steiner peeling proves the mean mu(T, lambda_m) < n/3 for such trees.
- The gap is mode <= ceil(mu), which is open.

### Degree-2 bridge decomposition:
When a vertex s has degree 2 (adjacent to l and u), with l a leaf:
  I(T) = (1+2x)P + (1+x)Q
where P = dp_B[u][0], Q = dp_B[u][1], and B is the tree with the pendant {l,s} removed.
This led to four proof "lanes," all converging on "prove mode(P) >= m-1," which is verified but unproved.

### What has been ruled out:
- Real-rootedness: fails for general trees
- Global LC: fails at n=26
- Mode superadditivity for products of LC polys: FALSE
- Various matching/LP/SDR approaches to Conjecture A: all fail
- Subdivision-based approaches: blocked

## Computational resources

If you propose a new invariant, I can test it on all 9.1 million trees n <= 22 (with all rootings or all support vertices). The test takes about 10-15 minutes. So propose something concrete and testable.

## What would constitute success

A. An invariant that provably implies unimodality AND is provably preserved by the tree DP. (Full proof.)

B. An invariant that provably implies unimodality, is conjectured to be preserved by the DP, and can be computationally verified. (Testable conjecture.)

C. A clear identification of exactly what algebraic/combinatorial property is needed for the DP closure step, even if you can't prove it holds. (Sharpened problem statement.)

## Standards

Be rigorous about what you prove vs. conjecture. If you propose an invariant, state it precisely enough to be tested computationally. If you cite a theorem (e.g., from GMTW or the TP literature), give the exact statement. Check your proposals against:
- P_5 rooted at vertex 1: E = [1,4,4,1], J = [1,2]. Mode 2.
- K_{1,4} rooted at center: E = [1,4,6,4,1], J = [1]. Mode 2.
- Path P_7 rooted at vertex 1: E = [1,6,12,10,3], J = [1,4,4,1]. Mode 3.
