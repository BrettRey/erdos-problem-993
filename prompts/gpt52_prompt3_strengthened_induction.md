# Task: Find the right strengthened inductive hypothesis to prove unimodality of tree independence polynomials

## The problem

Erdős Problem #993 (1987): For every tree T, the independence polynomial I(T; x) = Σ_{k=0}^{α} i_k x^k has unimodal coefficients (the sequence i_0, i_1, ..., i_α increases then decreases, with no interior valley).

Verified for all 1,198,738,056 trees on n ≤ 27 vertices. No proof exists.

## Definitions

An independent set in a graph is a set of vertices with no two adjacent. i_k counts independent sets of size k, with i_0 = 1. A sequence is unimodal if there's a peak: a_0 ≤ a_1 ≤ ... ≤ a_m ≥ a_{m+1} ≥ ... ≥ a_d.

Log-concavity (LC): a_k^2 ≥ a_{k-1}·a_{k+1} for all k. For positive sequences, LC ⟹ unimodal.

## The recursive structure of trees

### Vertex elimination
For any vertex v in tree T:

  I(T; x) = I(T - v; x) + x · I(T - N[v]; x)

where T - v deletes v and its edges, T - N[v] deletes v and all its neighbors. For a tree, T - v is a forest, and T - N[v] is a smaller forest.

### Leaf elimination (the most natural induction)
If v is a leaf with unique neighbor u:

  I(T; x) = I(T - v; x) + x · I(T - {v, u}; x)

Here T - v is a tree on n-1 vertices, and T - {v,u} is a forest (connected components of T after removing the pendant edge). Both are smaller and have known IS polynomials by induction.

### Tree DP (bottom-up product structure)
Root T at any vertex r:

  dp[v][0] = Π_{c child of v} (dp[c][0] + dp[c][1])     (v excluded)
  dp[v][1] = x · Π_{c child of v} dp[c][0]               (v included)
  I(T; x) = dp[r][0] + dp[r][1]

Note: dp[c][0] + dp[c][1] = I(T_c), the IS polynomial of the subtree at c. So dp[v][0] is a product of subtree IS polynomials, and dp[v][1] is x times a product of "exclusion polynomials."

### Forest factorization
For a forest F with components T_1, ..., T_k:

  I(F; x) = I(T_1; x) · I(T_2; x) · ... · I(T_k; x)

Products of LC polynomials with positive coefficients are LC (Cauchy product preserves LC). So IS polynomials of forests with LC components are LC.

## Why "unimodal" alone fails as an inductive hypothesis

The sum of a unimodal polynomial and x times another unimodal polynomial is NOT always unimodal. Simple example: f = 1 + 10x + x^2 (unimodal, mode 1) and g = 1 + x + 10x^2 (unimodal, mode 2). Then f + x·g = 1 + 11x + 2x^2 + 10x^3 (unimodal). But with different coefficients this can fail. So to make induction work via I(T) = I(T-v) + x·I(T-N[v]), we need the inductive hypothesis to be STRONGER than unimodality.

## What's known about the coefficients

1. **Mode ≈ n/3.** For most trees, the mode is near n/3. In the hard-core model at fugacity λ = 1, the mean is μ = I'(1)/I(1), and computationally mode ≤ ⌈μ⌉ for all trees n ≤ 22 (9.1M trees, 0 failures).

2. **LC holds for almost all trees.** Only 2 out of 279M trees at n = 26 fail LC. All trees n ≤ 25 are LC. So the inductive invariant need not be LC exactly, but something close might work.

3. **Products are well-behaved.** Products of tree IS polynomials are LC (since each factor is LC for small trees, and products preserve LC). The DP product step (computing dp[v][0]) preserves LC. The problematic step is the SUM dp[v][0] + dp[v][1].

4. **Ratio bounds.** For LC sequences, the ratios r_k = a_{k+1}/a_k are decreasing. Trees almost always satisfy this, with the ratios decreasing smoothly from r_0 ≈ n down to r_α ≈ 0.

5. **The mode is the smallest index where i_k = max.** For trees, there are typically no ties at the mode (the peak is strict).

## Approaches that have been tried and FAILED

- **LC as inductive hypothesis:** Fails because LC fails at n = 26.
- **Ultra-log-concavity (ULC):** Fails at n = 8.
- **Real-rootedness:** Fails for most trees (those containing K_{1,3}).
- **Mode superadditivity** (mode(f·g) ≥ mode(f) + mode(g) for LC polynomials): FALSE.
- **Per-coefficient dominance bounds** (e.g., i_m ≥ i_{m±1} plus specific corrections): Synthetic counterexamples exist showing these are insufficient.

## Your task

Find a property P of polynomials (with positive coefficients) such that:

1. **P implies unimodality.** (Or at least implies unimodality for the specific coefficient patterns of tree IS polynomials.)

2. **P is preserved (or nearly preserved) by the tree DP.** Specifically:
   - **Product preservation:** If f and g satisfy P, then f·g satisfies P.
   - **Sum preservation:** If A satisfies P and B satisfies P, then A + x·B satisfies P (or satisfies P under conditions that tree IS polynomials always meet).

3. **Base cases hold.** The single-vertex polynomial 1 + x satisfies P. The single-edge polynomial 1 + 2x + x^2 satisfies P. Small trees satisfy P.

4. **P is not too strong.** It must not fail at n = 26 where LC fails (since unimodality still holds there). If P implies LC, it can only do so for specific tree structures, not universally.

### Suggestions for the shape of P

- **Ratio-based:** A condition on the ratios r_k = a_{k+1}/a_k (e.g., the ratios satisfy some convexity or Lipschitz condition).
- **Mode-mean relationship:** mode ≤ f(μ) for some function f.
- **Weighted LC:** a_k^2 ≥ c_k · a_{k-1} · a_{k+1} where c_k depends on k and the degree of the polynomial.
- **Tail control:** Bounds on how fast coefficients decay away from the mode.
- **Generating function property:** A condition on I(T; x) as an analytic function (e.g., log-convexity of |I(T; re^{iθ})| in some variable, or a bound on the argument of I on a circle).
- **Something entirely different** that I haven't thought of.

### A concrete test case

The star S_5 = K_{1,4} has I = 1 + 5x + 6x^2 + x^3. Verify your property P on this. Then check the path P_7, which has I = 1 + 7x + 15x^2 + 10x^3 + x^4 (wait, let me not state this -- compute it yourself from the DP).

A useful stress test: the tree that is a path of length 3 with a leaf attached to one interior vertex (sometimes called a "T-shape" or fork). This has 5 vertices and I = 1 + 5x + 6x^2 + 2x^3.

## What I value

- **The key insight I need is what P should be.** Even a conjectural P that you can verify on small examples and that has the right structural properties for the DP is extremely valuable.
- **Be explicit about where the argument breaks down.** If you can show product preservation but not sum preservation, say so clearly.
- **Partial results are welcome.** A proof for a restricted class (e.g., caterpillars, trees of max degree ≤ 3, trees of diameter ≤ d) is still valuable.
- **Do NOT claim a complete proof unless every step is rigorous.** Flag assumptions clearly. A correct partial result is worth infinitely more than an incorrect complete proof.
