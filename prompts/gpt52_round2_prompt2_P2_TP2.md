# Task: Prove prefix TP2 (P2) for independence polynomials of trees

## Context

I have computationally verified a two-state invariant ("P⋆") on rooted trees that implies unimodality of the independence polynomial. It holds for all 9.1 million trees on n ≤ 22 vertices, at every possible rooting (~140 million checks), with zero failures. I need a proof.

This prompt focuses on one of the two conditions: **prefix TP2 (P2)**.

## Setup

For a tree T with n vertices, root it at any vertex r. The standard DP gives:

  dp[r][0] = Π_{c child of r} (dp[c][0] + dp[c][1])     (root excluded from IS)
  dp[r][1] = x · Π_{c child of r} dp[c][0]               (root included in IS)

Define:
  E(x) = dp[r][0] = Σ e_k x^k    (the "exclude-root" polynomial)
  J(x) = dp[r][1] / x = Π_{c child of r} dp[c][0] = Σ j_k x^k    (the "include-root" polynomial, shifted)

Then I(T; x) = E(x) + x · J(x), so the kth coefficient of I(T) is:

  i_k(T) = e_k + j_{k-1}    (with j_{-1} = 0)

## Combinatorial meaning

- e_k = number of independent sets of size k in T that do NOT contain root r.
- j_k = number of independent sets of size k in T - N[r] (the tree minus the closed neighborhood of r).

So j_k counts IS of the smaller forest T - N[r], while e_k counts IS of T that avoid r.

## The property P2 (prefix TP2)

Let m = mode(I(T)) = the smallest index where i_k(T) is maximized.

**P2:** For k = 0, 1, ..., m-1, we have:

  e_{k+1} · j_k ≥ e_k · j_{k+1}

Equivalently: the ratio sequence j_k / e_k is nonincreasing for k = 0, 1, ..., m.

In probability language: for IS of size k drawn from T, the conditional probability P(r ∈ S | |S| = k) = j_{k-1} / (e_k + j_{k-1}) is nonincreasing for k = 1, ..., m (up to the mode).

**Verified existentially:** For all 9,114,283 trees on n ≤ 22, there EXISTS at least one rooting where both P2 and P3 hold, 0 failures. However, P2 does NOT hold at all rootings — it fails at roughly 15-20% of rootings. P3 fails more often (~60% of rootings). For example, P_5 passes P2 at most rootings but fails P3 at endpoints and center.

## Why P2 matters

P2 is one of two conditions in my candidate inductive invariant for proving unimodality of tree IS polynomials. Together with P3 (tail domination), P2 ensures that I(T) = E + xJ is unimodal. Specifically, P2 ensures the mixture E + xJ has a single peak up to and including the mode, because the "J contribution" is front-loaded (higher relative weight at small k, lower at large k).

The combination P2 + P3 gives unimodality by this argument: P2 (nonincreasing j_k/e_k) means the shifted component xJ doesn't create a new peak before m, and P3 (e_k ≥ j_{k-1} for k ≥ m) means xJ doesn't create a new peak after m.

## Key structural facts

1. **Product structure.** If r has children c_1, ..., c_d with subtrees T_1, ..., T_d:
   - E = Π_{i=1}^d I(T_i; x)  where I(T_i) = E_i + x·J_i is the IS polynomial of subtree T_i
   - J = Π_{i=1}^d E_i(x)     where E_i = dp[c_i][0]

   So each factor of E is I(T_i) = E_i + x·J_i, while the corresponding factor of J is just E_i. The ratio of factors is I(T_i)/E_i = 1 + x·J_i/E_i.

2. **TP2 and products.** A key result of Gross, Mansour, Tucker, and Wang (SIDMA ~2015) is: if f(x) and g(x) are polynomials with positive coefficients such that (f_k) and (g_k) are "synchronised" (i.e., f and g satisfy a TP2-like condition), then their Hadamard product and ordinary product preserve this property.

   More precisely: if two pairs (A, B) and (C, D) each satisfy a_{k+1}b_k ≥ a_k b_{k+1}, then under convolution (product of generating functions), the resulting pair may or may not preserve this. The question is whether the TP2 condition on (E, J) = (Π I(T_i), Π E_i) follows from TP2 on each factor pair (I(T_i), E_i).

3. **Each factor pair.** For a single subtree T_i rooted at c_i:
   - I(T_i) = E_i + x·J_i
   - The ratio J_i,k / I(T_i)_k = j_{i,k} / (e_{i,k} + j_{i,k-1})

   If P2 holds inductively for T_i at root c_i, we have j_{i,k}/e_{i,k} nonincreasing up to mode(I(T_i)).

4. **For leaves.** If T_i is a single vertex (leaf), then E_i = 1 (constant), J_i = 1 (constant), I(T_i) = 1 + x. The ratio J_i/I(T_i) = 1/(1+x), which has coefficients alternating in sign — but as polynomials with positive coefficients, e_{i,k} = [k=0] and j_{i,k} = [k=0], so the ratio j_0/e_0 = 1 and all higher terms vanish. P2 is vacuously true.

5. **For stars K_{1,d} rooted at center.** E = (1+x)^d, J = 1. The ratio j_k/e_k = 1/C(d,k) which is strictly decreasing. P2 holds easily.

## The DP induction step

The core challenge is proving P2 is preserved when we combine subtrees at a root. Suppose r has children c_1, ..., c_d. We need:

  [Π I(T_i)]_{k+1} · [Π E_i]_k ≥ [Π I(T_i)]_k · [Π E_i]_{k+1}    for k < m

given that P2 holds for each (T_i, c_i).

Since I(T_i) = E_i + x·J_i and the corresponding J-factor is E_i, each "factor ratio" is I(T_i)/E_i = 1 + x·(J_i/E_i). The overall ratio is:

  E(x)/J(x) = Π I(T_i) / Π E_i = Π (1 + x · J_i/E_i)

P2 says the coefficients of J/E are nonincreasing up to mode m. Equivalently, the coefficients of E/J = Π(1 + x · J_i/E_i) form a sequence whose "ratio inverse" is nonincreasing — this is related to the coefficients of E/J being log-concave or having a related monotonicity property.

## Your task

Prove that P2 holds for all trees at all rootings. Approaches to consider:

**Approach A (TP2 product closure).** Show that if (E_i, J_i) satisfies the TP2 ratio condition for each subtree, then the products (Π I(T_i), Π E_i) satisfy it. The key reference is Gross-Mansour-Tucker-Wang's theory of synchronised sequences. In their framework: if sequences (a_n) and (b_n) are "synchronised" (a_{n+1}/a_n ≥ b_{n+1}/b_n), then their convolutions (â_n) and (b̂_n) obtained from multiplying the generating functions of two such pairs preserve synchronisation. Can you apply this iteratively?

**Approach B (Ratio analysis).** The ratio R(x) = J(x)/E(x) = Π E_i / Π(E_i + xJ_i) = Π 1/(1 + x·J_i/E_i). If each J_i/E_i has nonincreasing coefficients (which would follow from a stronger form of P2 for subtrees), does 1/(1 + xf) have nonincreasing coefficients when f does? And do products of such preserve the property?

**Approach C (Combinatorial injection).** P2 says: among IS of size k ≤ mode in the whole tree, the fraction containing root r is nonincreasing in k. Combinatorially: for k < m, the "k+1 to k transfer" of root-inclusion probability doesn't increase. Can you construct an injection that witnesses this?

**Approach D (Hard-core model).** In the hard-core lattice gas on T at fugacity λ, the occupation probability of r conditioned on |S| = k is j_{k-1}/i_k. P2 says this is nonincreasing up to mode. The unconditional occupation probability at fugacity λ is P(r ∈ S; λ) = λJ(λ)/I(T;λ), and the "conditional on size" version conditions on the Poisson-like size distribution. The monotonicity might follow from properties of the pressure function or its derivatives.

## Important warnings

1. **Global LC fails for trees.** Do NOT assume the IS polynomial is log-concave. It fails at n = 26 (Kadrawi & Levit 2023). Your proof must work without LC.

2. **P2 is only a PREFIX condition.** The ratio j_k/e_k may fail to be nonincreasing past the mode. This is critical — don't try to prove it for all k.

3. **The mode m depends on the WHOLE tree,** not just the subtrees. So the prefix range changes as you combine subtrees. Any inductive argument needs to handle this carefully.

## Standards

Be rigorous. If you find a proof, check it on these concrete examples:
- Path P_5 = [1,2,3,4,5], root at vertex 3.
  E = I(P_2) · I(P_2) where P_2 = {1,2} and {4,5}. So E = (1+2x)^2 = 1+4x+4x^2.
  J = E_1 · E_2 where E_i = dp[c_i][0] for each P_2 subtree. E_i = 1+x. So J = (1+x)^2 = 1+2x+x^2.
  I = E + xJ = 1+4x+4x^2 + x+2x^2+x^3 = 1+5x+6x^2+x^3. Mode m = 2.
  Ratios j_k/e_k: j_0/e_0 = 1, j_1/e_1 = 2/4 = 0.5. Nonincreasing up to m-1=1. ✓

- Star K_{1,4}, root at center. E = (1+x)^4, J = 1.
  Ratios: j_0/e_0 = 1/1 = 1, and j_k = 0 for k ≥ 1. Trivially nonincreasing. ✓

- The caterpillar on 7 vertices: path [1,2,3,4,5] with leaves 6,7 attached to vertex 3. Root at 3.

If you get stuck, describe precisely where the argument breaks down.
