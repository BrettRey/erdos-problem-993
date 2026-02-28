# Task: Product closure for ratio dominance of IS polynomials

## Context

I am working on a proof that the independence polynomial I(T; x) of every tree T is unimodal. I have reduced the problem to a single algebraic statement at support vertices (vertices adjacent to a leaf). The statement involves ratio dominance of products of IS polynomials, and I need to understand when ratio dominance is preserved under products (convolution of coefficient sequences).

## The specific question

Let T be a tree, r a support vertex with leaf child u. Define:

  A(x) = Product_{c in C} I(T_c; x)     (product of IS polynomials of subtrees at non-leaf children)
  B(x) = Product_{c in C} E_c(x)         (product of exclude-root polynomials at same children)

where C = {children of r} \ {u}, T_c is the subtree rooted at child c, and E_c = dp[c][0].

I need to show: (1+x)A(x) is ratio-dominant over B(x) up to the mode m of I(T; x) = (1+x)A + xB.

That is: for k = 0, ..., m-1,

  [(1+x)A]_{k+1} * B_k >= [(1+x)A]_k * B_{k+1}

## What I know about the factors

For each non-leaf child c with subtree T_c:
- S_c(x) = I(T_c; x) = E_c(x) + x * J_c(x), where E_c, J_c are the DP polynomials at c
- S_c has nonneg coefficients and is unimodal (by induction on n, assuming the conjecture holds for smaller trees)
- E_c has nonneg coefficients
- S_c >= E_c coefficientwise (since x*J_c >= 0)
- S_c and E_c are both IS polynomials of related graphs: S_c = I(T_c), E_c = I(T_c - c)

The key structural fact: each factor of A dominates the corresponding factor of B coefficientwise (S_c >= E_c), but this does NOT imply A >= B in the ratio-dominance sense. Factorwise ratio dominance S_c >= E_c (in the TP2 sense) is FALSE: for the 3-vertex star K_{1,2} rooted at center, S = 1+3x+x^2 and E = (1+x)^2 = 1+2x+x^2, and s_2*e_1 = 1*2 = 2 < 3*1 = 3 = s_1*e_2.

## Relevant theory

### Gross-Mansour-Tucker-Wang (SIDMA ~2015)
Two nonneg sequences (a_n) and (b_n) are "synchronised" if a_{n+1}/a_n >= b_{n+1}/b_n for all valid n (equivalently, a_{n+1}b_n >= a_nb_{n+1}).

**Theorem (GMTW):** If (a_n), (b_n) are synchronised and log-concave, and (c_n), (d_n) are synchronised and log-concave, then the convolution pairs (a*c, b*d) are also synchronised.

This means: if each factor pair (S_c, E_c) were synchronised AND log-concave, then the products (A, B) = (Product S_c, Product E_c) would be synchronised. But the factorwise synchronisation S_c >= E_c fails (as above).

### Hu-Wang-Zhao-Zhao (arXiv:1507.08430, 2015)
They introduce "partial synchronicity" as an intermediate between synchronicity and weak synchronicity, and prove it is preserved under convolution for LC sequences. Could this apply?

### The obstacle
Even if each factor pair had the right structure, the final step needs (1+x)A >= B, not just A >= B. The (1+x) multiplication is a specific convolution with (1, 1). The binomial upgrade lemma says: if A >= B (ratio-dominant) AND B is LC, then (1+x)A >= B. But A >= B fails at ~10% of support vertices.

## Computational data

### A and B are individually well-behaved:
- B is always fully log-concave (verified 59,061 support vertices, n <= 15, 0 failures)
- A is always fully log-concave (same check needed, but expected since A = Product of IS polynomials of trees, and products of LC polys preserve LC)
- Both A and B are unimodal

### A >= B fails but (1+x)A >= B always holds:
Example: n=7, A = [1,5,6,4,1], B = [1,4,6,4,1].
- A >= B fails at k=1: a_2*b_1 = 6*4 = 24 < 5*6 = 30 = a_1*b_2
- But (1+x)A = [1,6,11,10,5,1] and (1+x)A >= B holds:
  k=0: 6*1 = 6 >= 1*4 = 4. OK.
  k=1: 11*4 = 44 >= 6*6 = 36. OK.
  k=2: 10*6 = 60 >= 11*4 = 44. OK.

The (1+x) factor contributes enough at each index to compensate for A's deficit.

## Your task

Find a condition C on pairs (A, B) that:
1. Is satisfied by (Product S_c, Product E_c) at support vertices (i.e., is preserved under convolution of factor pairs)
2. Together with LC of B, implies (1+x)A >= B (ratio-dominant)
3. Is weaker than A >= B (which fails)

Alternatively: prove (1+x)A >= B directly without going through A >= B, by using the product structure and properties of IS polynomials.

## Specific sub-questions

### Q1: What IS-polynomial-specific properties could help?

IS polynomials of trees are:
- Unimodal (by assumption / induction)
- Log-concave up to n=25 (fails at n=26, but only for 2 trees out of 280 million)
- Real-rooted for paths, stars, caterpillars (but NOT for all trees)
- Have the "ultra-log-concave" property for forests (products of path IS polys)

Is there a property between "LC" and "ratio-dominant" that holds for all factor pairs (S_c, E_c) and is preserved under convolution?

### Q2: Does partial synchronicity (Hu et al.) apply?

In Hu-Wang-Zhao-Zhao's framework, partial synchronicity is defined as: there exists an index N such that a_{n+1}b_n >= a_nb_{n+1} for n < N and a_{n+1}b_n <= a_nb_{n+1} for n >= N. That is, the ratio a_n/b_n is first nondecreasing then nonincreasing.

Do the factor pairs (S_c, E_c) satisfy partial synchronicity? If so, does convolution preserve it? And does partial synchronicity of (A, B) suffice to prove (1+x)A >= B up to mode?

### Q3: Can we exploit the specific structure S_c = E_c + x*J_c?

Since A/B = Product(S_c/E_c) = Product(1 + x*J_c/E_c), and each J_c/E_c involves the DP of a subtree, maybe we can bound the coefficient ratios of A/B directly.

The ratio A_k/B_k measures how much extra "density" the inclusion option (x*J_c) adds at IS-size k. Can we show this ratio is nonincreasing up to mode (after the (1+x) adjustment)?

## Standards

Be rigorous. Cite specific theorems from GMTW or Hu et al. if you use them. If you propose a new condition, verify it on the examples above. If a full proof is out of reach, identify the sharpest sufficient condition you can and say precisely what remains.
