# Task: Prove P2 (ratio dominance) at support vertices of trees

## Context

I am working on Erdos Problem #993: prove that the independence polynomial of every tree is unimodal. I have a two-state invariant P* that implies unimodality, and I have proved half of it (P3, tail domination) at support vertices via a simple injection. The remaining open piece is P2 (ratio dominance) at support vertices, which has been computationally verified on 59,916,124 support-vertex checks across all 9,114,283 trees on n <= 22 vertices, with zero failures.

**Your task: prove P2 at support vertices.**

## Setup

For a tree T on n vertices, let r be a support vertex (a vertex adjacent to at least one leaf). The standard tree DP gives:

  E(x) = dp[r][0]    (IS polynomial counting sets avoiding root r)
  J(x) = dp[r][1]/x  (IS polynomial of T - N[r])
  I(T; x) = E(x) + x*J(x)

Let m = mode(I(T)) = smallest index maximizing the coefficient sequence.

**P2 (prefix ratio dominance):** For k = 0, 1, ..., m-1:

  e_{k+1} * j_k >= e_k * j_{k+1}

Equivalently: the ratio j_k / e_k is nonincreasing on {0, ..., m}.

## The support-vertex factorization

Let u be a leaf child of r. The key structural fact is:

  E(x) = (1 + x) * A(x)
  J(x) = B(x)

where:
  A(x) = Product_{c != u} I(T_c; x)    (product of IS polynomials of non-leaf child subtrees)
  B(x) = Product_{c != u} E_c(x)       (product of exclude-root polynomials of those subtrees)

Here T_c is the subtree rooted at child c, and E_c = dp[c][0] for that subtree.

So P2 becomes: (1+x)A is ratio-dominant over B up to mode m.

## Known results

### Binomial upgrade lemma (PROVED):
If A is ratio-dominant over B on {0,...,m} AND B is log-concave on that range, then (1+x)A is ratio-dominant over B on that range.

Proof: e_{k+1}b_k - e_kb_{k+1} = (a_{k+1}b_k - a_kb_{k+1}) + (a_kb_k - a_{k-1}b_{k+1}). First bracket >= 0 by A >= B. Second bracket >= 0 by combining ratio dominance at k-1 with LC of B.

### But A >= B fails!
A is NOT always ratio-dominant over B. It fails at roughly 10% of support vertices. Examples:

- n=5, path P_5, root at vertex 1: A = [1,3,1], B = [1,2,1]. At k=1: a_2*b_1 = 1*2 = 2 < 3*1 = 3 = a_1*b_2.
- n=7: A = [1,5,6,4,1], B = [1,4,6,4,1]. At k=1: a_2*b_1 = 6*4 = 24 < 5*6 = 30 = a_1*b_2.

So the binomial upgrade lemma's hypothesis is too strong.

### B is always log-concave (verified):
B = J = Product of E_c polynomials is always fully log-concave (59,061 support-vertex checks at n <= 15, zero failures). This is expected: each E_c is itself an IS polynomial (of T_c - c), and products of LC polynomials with positive coefficients preserve LC.

### P2 always holds despite A >= B failing (verified):
(1+x)A is ratio-dominant over B at all 59,916,124 support vertices checked. The (1+x) factor provides essential margin beyond what the binomial upgrade lemma guarantees.

## What needs to be proved

Show that for every tree T and every support vertex r:

  e_{k+1} * j_k >= e_k * j_{k+1}    for k = 0, ..., m-1

where E = (1+x)A, J = B, and m = mode(E + xB).

## Approaches to consider

### Approach A: Strengthen the upgrade lemma

The binomial upgrade lemma requires A >= B (ratio-dominant). But (1+x)A >= B always holds even when A >= B fails. What weaker condition on (A, B) suffices?

Key observation: e_{k+1}b_k - e_kb_{k+1} = (a_{k+1} + a_k)b_k - (a_k + a_{k-1})b_{k+1} = a_{k+1}b_k - a_{k-1}b_{k+1} + a_k(b_k - b_{k+1}).

If b_k >= b_{k+1} (B is nonincreasing in this range, i.e., we're past B's mode), then the a_k(b_k - b_{k+1}) term is nonneg. So only the INCREASING part of B matters.

Can you decompose the range {0,...,m-1} into "B increasing" and "B decreasing" parts and handle each separately?

### Approach B: Direct coefficient analysis

Since E = (1+x)A and B are both products of IS polynomials of tree substructures, their coefficients have combinatorial meaning.

  e_k = a_k + a_{k-1} = #{IS of size k in T - r}
  b_k = j_k = #{IS of size k in T - N[r]}

P2 at index k says: (a_{k+1} + a_k) * b_k >= (a_k + a_{k-1}) * b_{k+1}.

Rearranging: a_{k+1}*b_k - a_{k-1}*b_{k+1} >= a_k*(b_{k+1} - b_k).

If the mode of B is at most m-1 (so b_{k+1} <= b_k for k >= mode(B)), the RHS is <= 0 and we need the LHS to be nonneg. Is there a combinatorial argument?

### Approach C: Hard-core model / probabilistic

In the hard-core model at fugacity lambda on tree T, the partition function is Z(T; lambda) = I(T; lambda). The ratio j_k/e_k = #{IS of size k in T-N[r]} / #{IS of size k in T-r}.

The denominator counts IS of T-r (a forest), and the numerator counts IS of T-N[r] (a smaller forest). Since T-N[r] is obtained from T-r by deleting deg(r)-1 additional vertices (the non-leaf neighbors' neighborhoods), the ratio should decrease as k grows (bigger IS are more constrained, so the deletion of vertices hurts more).

Can this be formalized using the hard-core model on forests?

### Approach D: Use the product structure

A = Product of S_c = Product of I(T_c) over non-leaf children.
B = Product of E_c over non-leaf children.

Each factor satisfies S_c = E_c + x*J_c. So A/B = Product(1 + x*J_c/E_c).

P2 for (E, B) = ((1+x)A, B) can be rewritten in terms of the factor ratios. Since each J_c/E_c involves the DP values of subtrees, and the product structure converts ratio comparisons into sums of log-ratios, maybe there's a telescoping or convexity argument.

## Concrete test cases

1. Path P_5, root at vertex 1 (support vertex):
   E = [1,4,4,1], J = B = [1,2], A = [1,3,1]. Mode m = 2.
   P2 at k=0: e_1*b_0 = 4 >= 1*2 = e_0*b_1 = 2. OK (slack 2).
   P2 at k=1: e_2*b_1 = 4*2 = 8 >= 4*0 = e_1*b_2 = 0. OK.

2. Path P_7, root at vertex 1 (support):
   E = [1,6,12,10,3], J = B = [1,4,4,1]. Mode m = 3.
   P2: check k=0,1,2.

3. Star K_{1,4}, root at center (support):
   E = (1+x)^4 = [1,4,6,4,1], J = B = [1]. Mode m = 2.
   P2: trivially OK (j_k = 0 for k >= 1).

4. Caterpillar: path [0,1,2,3,4] with extra leaves 5,6 at vertex 2. Root at vertex 2 (support, degree 4).

## Standards

Be rigorous. If you find a proof, verify it on the test cases. If you find a weaker sufficient condition than A >= B that always holds and implies (1+x)A >= B, state it precisely. If you get stuck, say exactly where.
