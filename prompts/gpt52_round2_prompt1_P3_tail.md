# Task: Prove tail domination (P3) for independence polynomials of trees

## Context

I have computationally verified a two-state invariant ("P⋆") on rooted trees that implies unimodality of the independence polynomial. It holds for all 9.1 million trees on n ≤ 22 vertices, at every possible rooting (~140 million checks), with zero failures. I need a proof.

This prompt focuses on one of the two conditions: **tail domination (P3)**.

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
- j_{k-1} = number of independent sets of size k in T that DO contain root r.
  (Each such set S has r ∈ S, so S \ {r} is an IS of size k-1 in T - N[r], and |S \ {r}| = k-1.)

So i_k(T) = e_k + j_{k-1} splits each IS of size k by whether it contains the root.

## The property P3 (tail domination)

Let m = mode(I(T)) = the smallest index where i_k(T) is maximized.

**P3:** For all k ≥ m, we have e_k ≥ j_{k-1}.

Equivalently: for IS of size k ≥ mode, at least half do not contain root r. In probability language:

  P(r ∈ S | |S| = k) ≤ 1/2    for all k ≥ mode(I(T))

**Verified existentially:** For all 9,114,283 trees on n ≤ 22, there EXISTS at least one rooting where P3 holds, 0 failures. However, P3 does NOT hold at all rootings — it fails at roughly 60% of rootings. For example, P_5 rooted at an endpoint or at the center vertex fails P3 (the unique maximum IS contains the root, so j_{α-1} > 0 = e_α).

## Why P3 matters

P3 is one of two conditions in my candidate inductive invariant for proving unimodality of tree IS polynomials. If both P2 (a monotone likelihood ratio condition) and P3 hold, then I(T) = E + xJ is unimodal. Specifically, P3 prevents the shifted J from creating a secondary bump in the tail of the coefficient sequence.

## Structural observations

1. **Including r is costly.** If r is in the IS, then all deg(r) neighbors of r are excluded, "losing" deg(r) vertices. For large IS (size k ≥ mode ≈ n/3), this loss is increasingly constraining. Intuitively, P3 says this cost makes root-inclusion a minority strategy for large IS.

2. **At the independence number α(T).** e_α counts IS of size α not containing r. j_{α-1} counts IS of size α containing r. For many trees, only the maximum IS contain r (or don't), but generically both are positive.

3. **For stars K_{1,d} rooted at the center.** E = (1+x)^d (all IS of d leaves), J = [1] (only the empty set of T-N[r]). So j_k = 0 for k ≥ 1, and P3 is trivially true.

4. **For paths P_n rooted at an endpoint.** The root has degree 1, so removing N[r] = {r, neighbor} leaves a path of length n-2. E and J are both IS polynomials of shorter paths. P3 holds (verified) but the proof would need to use properties of path IS polynomials.

5. **For a general rooting.** E = Π_{c child of r} I(T_c; x) where T_c is the subtree at child c. J = Π_{c child of r} E_c(x) where E_c is the "exclude c" polynomial of subtree T_c. So E is a product of subtree IS polynomials, and J is a product of "exclude-root" polynomials from the same subtrees.

## Your task

Prove that P3 holds for all trees at all rootings. Approaches to consider:

**Approach A (Injection).** For each IS S of size k ≥ m containing r, find an injection φ: S ↦ S' where S' is an IS of size k not containing r. The map could: remove r from S, then add some vertex not in N[S \ {r}]. The challenge is making φ injective.

**Approach B (Generating function).** Use the product structure E = Π I(T_c), J = Π E_c. Since I(T_c) = E_c + x·J_c, each factor of E "dominates" the corresponding factor of J (because I(T_c) = E_c + x·J_c ≥ E_c coefficientwise, with the shift). Can you show this termwise domination survives the product?

**Approach C (Hard-core model).** In the hard-core model at fugacity λ on a tree, the probability of including any vertex v is P(v) = λ · Z_{T-N[v]}(λ) / Z_T(λ). The conditional probability given |S| = k is P(r ∈ S | |S| = k) = j_{k-1} / i_k. Show this is ≤ 1/2 for k ≥ mode using properties of the partition function on trees.

**Approach D (Induction on tree size).** If v is a leaf of T with neighbor u, then I(T) = I(T-v) + x·I(T-{v,u}). Does P3 for T-v and T-{v,u} imply P3 for T? What inductive hypothesis would be needed?

## Standards

Be rigorous. If you find a proof, check it on these concrete examples:
- Path P_5 = [1,2,3,4,5], root at vertex 3. I(P_5) = 1 + 5x + 6x^2 + x^3, mode m = 2.
- Star K_{1,4}, root at center. I = 1 + 5x + 6x^2 + 4x^3 + x^4, mode m = 2. (Wait, compute this yourself.)
- The fork: vertex 1 adjacent to 2, 2 adjacent to 3 and 4. I = 1 + 4x + 3x^2, mode m = 1. Root at vertex 2.

If you get stuck, describe precisely where the argument breaks down.
