# Task: Prove P⋆ is preserved by the tree DP (full closure argument)

## Context

I have computationally verified a two-state invariant ("P⋆") on rooted trees that implies unimodality of the independence polynomial. For all 9.1 million trees on n ≤ 22 vertices, there EXISTS at least one rooting where P⋆ holds (zero failures in the existential sense).

**Important correction:** P⋆ does NOT hold at all rootings. In fact, P3 (tail domination) fails at roughly 60% of rootings, and P2 (prefix TP2) fails at roughly 15-20%. For example, the path P_5 passes P⋆ only when rooted at its degree-2 interior vertices (roots 1 and 3 in 0-indexed path [0,1,2,3,4]), but fails at the endpoints and the center vertex.

This prompt asks for help understanding which rootings satisfy P⋆ and whether the invariant can be proved inductively, given that it's not universal over rootings.

## Setup

For a tree T with n vertices, root it at any vertex r. The standard DP gives:

  dp[r][0] = Π_{c child of r} (dp[c][0] + dp[c][1])     (root excluded from IS)
  dp[r][1] = x · Π_{c child of r} dp[c][0]               (root included in IS)

Define for each vertex v:
  E_v(x) = dp[v][0] = Σ e_k x^k
  J_v(x) = dp[v][1] / x = Π_{c child of v} dp[c][0] = Σ j_k x^k
  S_v(x) = E_v + x·J_v = dp[v][0] + dp[v][1] = I(T_v; x)    (IS polynomial of subtree T_v)

## The P⋆ invariant

For a rooted tree (T, r) with S = I(T; x) = E + xJ, let m = mode(S) (smallest maximizer).

P⋆(E, J; m) consists of:

1. **P2 (prefix TP2):** e_{k+1} · j_k ≥ e_k · j_{k+1}  for k = 0, 1, ..., m-1.
   Equivalently: the ratio j_k/e_k is nonincreasing on {0, ..., m}.

2. **P3 (tail domination):** e_k ≥ j_{k-1}  for all k ≥ m.
   Equivalently: among IS of size k ≥ mode, at most half contain the root.

## Computational data on which rootings satisfy P⋆

### P_5 = path on 5 vertices [0,1,2,3,4]

I(P_5) = 1 + 5x + 6x² + x³, mode m = 2.

| Root | deg | E | J | P2 | P3 | P⋆ |
|------|-----|---|---|----|----|-----|
| 0 | 1 | [1,4,3] | [1,3,1] | ✓ | ✗ (k=3) | ✗ |
| 1 | 2 | [1,4,4,1] | [1,2] | ✓ | ✓ | ✓ |
| 2 | 2 | [1,4,4] | [1,2,1] | ✓ | ✗ (k=3) | ✗ |
| 3 | 2 | [1,4,4,1] | [1,2] | ✓ | ✓ | ✓ |
| 4 | 1 | [1,4,3] | [1,3,1] | ✗ | ✗ (k=3) | ✗ |

Pattern: P3 fails when J extends to high degrees (when T - N[r] is large enough to have big IS). Roots 1 and 3 pass because removing their closed neighborhood leaves a short path, making J low-degree.

### Aggregate statistics (n = 5 to 10)

| n | total rootings | P⋆ failures | P3 failures | P2 failures |
|---|---------------|-------------|-------------|-------------|
| 5 | 15 | 9 (60%) | 9 | 4 |
| 6 | 36 | 21 (58%) | 21 | 8 |
| 7 | 77 | 47 (61%) | 47 | 10 |
| 8 | 184 | 110 (60%) | 110 | 37 |
| 9 | 423 | 250 (59%) | 250 | 84 |
| 10 | 1060 | 624 (59%) | 624 | 168 |

P3 is the binding constraint: every P⋆ failure is a P3 failure. P2 failures are a strict subset.

### The existential result

Despite these per-rooting failures, for every tree on n ≤ 22 vertices (9,114,283 trees), there exists at least one rooting where P⋆ holds. This is a weaker but still nontrivial statement.

## Why P⋆ implies unimodality

If P⋆ holds at ANY rooting of T, then I(T) is unimodal. Here's why:

Since I(T) = E + xJ and the mode m is the peak of I(T)'s coefficients:
- P2 ensures the "J contribution" doesn't create a pre-mode dip (the mixture stays increasing up to mode).
- P3 ensures E dominates the tail, so the sum stays decreasing after mode.

More precisely: i_k = e_k + j_{k-1}. For k ≥ m, P3 gives e_k ≥ j_{k-1}, so at least half the IS of size k avoid the root. Combined with the decreasing tail of E (which follows from E being a product of unimodal polynomials, plus structural properties), this prevents a secondary bump.

## The DP structure

Suppose root r has children c_1, ..., c_d, each heading subtree T_i. Then:

  E_r = Π_{i=1}^d S_i       where S_i = I(T_i; x) = E_i + x·J_i
  J_r = Π_{i=1}^d E_i
  S_r = E_r + x·J_r = Π S_i + x · Π E_i

## Three questions for you

### Question 1: Characterize good rootings

What structural property of the root vertex r makes P⋆ hold? Hypotheses:
- **Degree constraint:** Maybe r needs degree ≤ some threshold?
- **Subtree balance:** Maybe the subtrees at r need to be "balanced" in some sense?
- **Neighborhood size:** P3 fails when T - N[r] is too big relative to T. Maybe we need |N[r]| = deg(r) + 1 to be large enough that T - N[r] has independence number < m?

From the P_5 data: root 1 (degree 2) passes, root 0 (degree 1, leaf) fails, root 2 (degree 2, center) fails. So it's not purely about degree. Root 1 has N[1] = {0,1,2}, so T - N[1] = {3,4}, a path of length 2 with α = 1. Root 2 has N[2] = {1,2,3}, so T - N[2] = {0,4}, two isolated vertices with α = 2. The higher α of T - N[r] causes J to extend to higher degrees, violating P3.

Conjecture: P⋆ holds at root r if and only if α(T - N[r]) ≤ m - 1 (the independence number of the complement is strictly less than the mode).

### Question 2: Does the existential version suffice?

For proving unimodality, we only need P⋆ at ONE rooting. Can we show: for every tree T, there exists a vertex r such that P⋆(E_r, J_r; m) holds?

One approach: show that P⋆ holds when r is chosen to maximize some structural quantity (like centrality, or maximum degree, or "centroid"). If there's a canonical choice of r that always works, the proof might go through.

### Question 3: Can P⋆ be made universal?

Is there a modified version of P3 that holds at ALL rootings? For example:
- P3': e_k ≥ j_k (instead of j_{k-1}) for k ≥ m.
  This removes the index shift. Since e_k counts IS of size k in T avoiding r, and j_k counts IS of size k in T - N[r], and T - N[r] ⊂ T \ {r}, every IS of T - N[r] is also an IS of T avoiding r. So e_k ≥ j_k ALWAYS holds, trivially! But this is too weak — it doesn't constrain the mixture E + xJ enough.

- Maybe combine P3' with a stronger P2 condition?
- Or use a different base for the domination: replace j_{k-1} with something between j_k and j_{k-1}?

## Relevant literature

- **Gross, Mansour, Tucker, Wang (~2015, SIDMA):** Synchronised sequences and products. If (a_n)↔(b_n) are synchronised (a_{n+1}/a_n ≥ b_{n+1}/b_n) and (c_n)↔(d_n) are synchronised, then convolutions (a*c)↔(b*d) are synchronised.
- **Brenti (1989):** Unimodal, log-concave and Pólya frequency sequences.
- **Bendjeddou & Hardiman (2024, arXiv:2405.00511):** Pre-Lorentzian gluing for tree IS polynomials.

## Concrete examples for testing

- Path P_3 rooted at vertex 1 (center): E = (1+x)², J = 1, S = 1+3x+x². Mode 1. P2: ✓. P3: e_1=2≥j_0=1 ✓, e_2=1≥j_1=0 ✓. ✓

- Star K_{1,4} rooted at center: E = (1+x)⁴, J = 1. All j_k = 0 for k ≥ 1. P3 trivially holds.

- Star K_{1,4} rooted at a leaf: E = (1+x)³ + x, J = (1+x)³. Compute this yourself and check P2 and P3.

- Double star S(3,3): center edge with 3 leaves on each side. Root at one center, root at a leaf.

## Standards

Be rigorous. If you make a conjecture (like "P⋆ holds iff α(T - N[r]) < m"), test it against the examples above. If you propose a modified invariant, check it computationally on the P_5 data.

If you get stuck, describe precisely where the argument breaks down and what additional data would help.
