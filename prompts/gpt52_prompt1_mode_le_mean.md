# Task: Prove that mode ≤ ⌈μ⌉ for independence polynomials of trees

## The Erdős conjecture (background)

Erdős Problem #993 (Alavi, Malde, Schwenk, and Erdős 1987) conjectures that the independence polynomial of every tree is unimodal. This has been verified for all 1.2 billion trees on n ≤ 27 vertices, but no proof exists.

## Definitions

**Independence polynomial.** For a tree T on n vertices, an independent set is a set of vertices with no two adjacent. Define

  I(T; x) = Σ_{k=0}^{α(T)} i_k x^k

where i_k = i_k(T) counts independent sets of size k, with i_0 = 1 and α(T) = independence number. Example: the path P_4 has I = 1 + 4x + 3x^2.

**Mode.** The mode of I(T) is the smallest index m such that i_m = max_k i_k.

**Mean.** The mean of the coefficient distribution is

  μ = I'(1)/I(1) = Σ_k k·i_k / Σ_k i_k.

**Log-concavity (LC).** A sequence (a_0, a_1, ...) with positive terms is log-concave if a_k^2 ≥ a_{k-1}·a_{k+1} for all k.

**Unimodality.** A sequence is unimodal if it increases then decreases (no interior valley). For positive sequences, LC implies unimodality.

## Tree DP

For a tree rooted at any vertex r, we can compute I(T) by dynamic programming:

  dp[v][0] = Π_{c child of v} (dp[c][0] + dp[c][1])     (v excluded from IS)
  dp[v][1] = x · Π_{c child of v} dp[c][0]               (v included in IS)

  I(T; x) = dp[r][0] + dp[r][1]

This runs in O(n^2) time per tree (polynomial multiplication at each node). Note that dp[c][0] + dp[c][1] = I(T_c), the IS polynomial of the subtree rooted at c.

## Known results

1. **Darroch's theorem (1964).** For a Poisson Binomial Distribution (PBD), where f(x) = Π_{i=1}^n (1 + p_i·x) with 0 < p_i ≤ 1, the mode is in {⌊μ⌋, ⌈μ⌉} where μ = Σ p_i. For tree IS polynomials, the DP gives I(T) as a sum of products, not a pure product, so Darroch does not directly apply.

2. **Strict LC does NOT imply mode ∈ {⌊μ⌋, ⌈μ⌉}.** The sequence (1, 2, 3, 4) is strictly LC (4 > 3, 9 > 8) with μ = 2.0 and mode = 3 > ⌈μ⌉. So any proof for trees must use structure beyond LC.

3. **LC status of trees.** Tree IS polynomials are LC for all trees with n ≤ 25 (hundreds of millions of trees). At n = 26, exactly 2 out of 279,793,450 trees fail LC (matching Kadrawi & Levit 2023). Galvin (2025) constructed infinite families with LC failures. But ALL trees through n = 27 are unimodal.

4. **d_leaf ≤ 1 trees are strictly LC.** All 931,596 trees through n = 23 where every vertex is a leaf or adjacent to a leaf have strictly LC independence sequences (0 failures).

5. **Real-rootedness fails.** Tree IS polynomials are not always real-rooted (fails for trees containing K_{1,3} as induced subgraph, which is almost all trees). Chudnovsky-Seymour (2007) proved real-rootedness only for claw-free graphs, which excludes most trees. So "real roots ⟹ LC ⟹ unimodal" is blocked.

6. **Stars and paths.** For the star K_{1,n-1}: I = (1+x)^{n-1} + x. This is a PBD plus a small perturbation. For paths, I satisfies a Fibonacci-type recurrence.

## The specific conjecture (your target)

**For every tree T: mode(I(T)) ≤ ⌈μ(I(T))⌉**, where μ = I'(1)/I(1).

**Computational evidence:** Verified for ALL 9,114,283 trees with n ≤ 22, with 0 violations. This includes trees that fail LC at larger n (the two LC-failing trees at n = 26 have μ ≈ 7.83 and mode = 8, so mode = ⌈μ⌉ holds).

## Why this matters

I have proved ("Steiner peeling" theorem) that for all trees where every vertex is a leaf or adjacent to a leaf (the "d_leaf ≤ 1" class): μ(T) < n/3. Combined with mode ≤ ⌈μ⌉, this immediately gives mode ≤ ⌊n/3⌋ + 1 for these trees, which is exactly "Conjecture A" -- the key open step in my reduction of the Erdős unimodality conjecture. So proving mode ≤ ⌈μ⌉ would close a major gap.

## Potential approaches to consider

- **Extensions of Darroch.** Tree IS polynomials are "almost" PBDs (products of (1+x) and (1+2x) factors at the DP level, combined by addition). Can Darroch's proof technique extend to sums of products?

- **The hard-core model.** At fugacity λ > 0, the probability distribution π_λ(S) ∝ λ^|S| over independent sets has mean μ(λ) = λI'(λ)/I(λ). The mode at λ = 1 is the coefficient mode. Can statistical mechanics (e.g., correlation decay, log-Sobolev inequalities, FKG-type results on trees) control the mode-mean relationship?

- **Recursive structure.** I(T) = dp[r][0] + dp[r][1] where dp[r][0] is a product of subtree IS polynomials. Can an inductive argument on the tree build mode ≤ ⌈μ⌉ step by step?

- **Stochastic dominance or coupling.** Can you couple the distributions at different fugacities to control where the mode falls relative to the mean?

## What NOT to try

- Do not try to prove LC for all trees (it fails at n = 26).
- Do not try real-rootedness (fails for trees with K_{1,3}).
- Do not try ultra-log-concavity (fails at n = 8).

## Standards

Be rigorous. Flag any step that relies on an unverified assumption. I would much rather have an honest "I can reduce it to X" or "here's where the argument breaks down" than a claimed complete proof with a gap. If you find a counterexample, that's equally valuable -- describe it precisely.
