# Tight Mode Analysis: mode(I(T-w)) = d(I(T))

## Setup

For a tree T with independence polynomial I(T) = (i_0, i_1, ..., i_alpha),
define:
- mode(I(T)) = index of the (last) maximum coefficient
- d(I(T)) = first descent index (first i with i_{i} < i_{i-1})
- The conjecture: for every vertex w, mode(I(T-w)) <= d(I(T))

This analysis characterizes the **tight cases** where equality holds:
mode(I(T-w)) = d(I(T)).

Data: exhaustive enumeration for n = 5 through 16.

## Summary Statistics

| n  | Trees   | Vertex-cases | Tight | Fraction | Tight trees |
|----|---------|-------------|-------|----------|-------------|
| 5  | 3       | 15          | 0     | 0        | 0           |
| 6  | 6       | 36          | 0     | 0        | 0           |
| 7  | 11      | 77          | 2     | 2.60%    | 2           |
| 8  | 23      | 184         | 0     | 0        | 0           |
| 9  | 47      | 423         | 1     | 0.24%    | 1           |
| 10 | 106     | 1,060       | 10    | 0.94%    | 9           |
| 11 | 235     | 2,585       | 30    | 1.16%    | 22          |
| 12 | 551     | 6,612       | 4     | 0.06%    | 4           |
| 13 | 1,301   | 16,913      | 149   | 0.88%    | 92          |
| 14 | 3,159   | 44,226      | 521   | 1.18%    | 356         |
| 15 | 7,741   | 116,115     | 280   | 0.24%    | 141         |
| 16 | 19,320  | 309,120     | 2,011 | 0.65%    | 1,008       |

Total: 3,008 tight cases out of 497,366 vertex-deletions (0.60%).

Note the parity effect: even n tends to produce more tight cases (0.70% vs
0.34% for odd n), and the counts oscillate with a period-3 or period-4
substructure. This likely reflects the interplay between floor functions in
the mode position and the parity of alpha(T).

## Invariant Properties (100% of tight cases)

These properties hold in **every** tight case observed (n = 5..16):

1. **w is never a leaf.** deg(w) >= 2 in all 3,008 cases.

2. **alpha(T) = alpha(T-w).** The independence number is preserved: w is
   never in every maximum independent set. Equivalently, there exists a
   maximum independent set that excludes w.

3. **d(I(T)) = mode(I(T)) + 1.** In every tight case, the first descent
   comes exactly one step after the mode of I(T). That is, I(T) has a
   "sharp peak": i_{mode} > i_{mode+1} with mode+1 = d(I(T)).

4. **mode(I(T-w)) - mode(I(T-N[w])) >= 2.** Removing the open neighborhood
   (T-w) has its mode at d(I(T)), while removing the closed neighborhood
   (T-N[w]) has its mode at least 2 lower. The distribution is:

   | mode(g) - mode(h) | Count | Fraction |
   |--------------------|-------|----------|
   | 2                  | 2,393 | 79.6%    |
   | 3                  | 601   | 20.0%    |
   | 4                  | 14    | 0.5%     |

## Vertex Properties

### deg(w) distribution

| deg(w) | Count | Fraction |
|--------|-------|----------|
| 2      | 13    | 0.4%     |
| 3      | 782   | 26.0%    |
| 4      | 1,195 | 39.7%    |
| 5      | 694   | 23.1%    |
| 6      | 262   | 8.7%     |
| 7      | 54    | 1.8%     |
| 8      | 7     | 0.2%     |
| 9      | 1     | 0.0%     |

The modal degree of w is 4, and deg(w) >= 3 accounts for 99.6% of cases.
The 13 cases with deg(w) = 2 all occur at n = 16, suggesting they become
more common at larger n.

### Is w the maximum-degree vertex?

| n range | w = max-deg | Fraction |
|---------|-------------|----------|
| n <= 12 | 45/47       | 95.7%    |
| n = 13  | 101/149     | 67.8%    |
| n = 14  | 339/521     | 65.1%    |
| n = 15  | 112/280     | 40.0%    |
| n = 16  | 1068/2011   | 53.1%    |

w is usually but not always the highest-degree vertex. The fraction declines
with n, suggesting that at larger n, tight cases can occur at vertices with
only moderately high degree, not necessarily the global maximum.

### w in center of tree

About 19.6% of tight vertices are center vertices (minimum eccentricity).
This is higher than the base rate (most trees have 1 or 2 center vertices
out of n), but far from universal.

## Tree Structure

### Type distribution

| Type       | Count | Fraction |
|------------|-------|----------|
| other      | 1,949 | 64.8%    |
| caterpillar| 1,039 | 34.5%    |
| spider     | 20    | 0.7%     |

No paths or stars appear, which makes sense: paths have no vertex with
deg >= 3, and stars have trivial independence polynomials.

### Two subclasses of tight cases (n <= 12 analysis)

**Type A (palindromic I(T-w)):** 23/47 small cases have I(T-w) palindromic.
When I(T-w) is palindromic, the mode is at floor(alpha(T-w)/2), and the
equality mode(g) = d_f holds because d_f happens to equal this midpoint.

Palindromic I(T-w) arises when the components of T-w form a forest whose
independence polynomial is self-reciprocal. This happens, for example, when
the components are paths (I(P_n) is unimodal with known symmetry properties)
or disjoint unions with matching structure.

**Type B (non-palindromic I(T-w)):** 24/47 small cases are not palindromic.
Some have a flat top (g[d_f] = g[d_f-1]) while others have a strict
maximum at d_f. In the strict cases, the large component of T-w
contributes enough independent sets at index d_f to push the mode there.

At larger n, the palindromic fraction declines (12.8% at n=13, 4.8% at
n=16), while the "flat-top" fraction (g[d_f] = g[d_f-1], not necessarily
palindromic) also declines from ~50% to ~9%. So at larger n, most tight
cases have a genuine strict maximum at d_f.

## d(I(T)) / n ratio

In tight cases, d(I(T))/n ranges from 0.33 to 0.47, with mean 0.375.
This means the first descent consistently occurs around n/3 to n/2.7,
which is consistent with the known behavior of independence polynomial
modes for trees (the mode is typically around alpha/3 to alpha/2, and
alpha is typically around n/2 to 2n/3).

## Implications for Proof Strategy

### The key constraint: d_f = mode_f + 1

The invariant d(I(T)) = mode(I(T)) + 1 means tight cases only arise when
I(T) has a **sharp peak** (an immediate descent after the mode). Trees where
I(T) has a plateau at the peak (multiple indices achieving the maximum) or
a gentle descent cannot produce tight cases, because d_f > mode_f + 1 would
give more room.

This suggests the proof should focus on trees with sharp peaks. These tend
to be trees with moderate branching (not too path-like, not too star-like).

### Why w is never a leaf

When w is a leaf, T-w is a tree on n-1 vertices, and its independence
polynomial is very similar to I(T). The mode can shift by at most 1 (it
typically stays the same or drops by 1), which is not enough to reach d_f
= mode_f + 1 from below.

More precisely: if w is a leaf attached to v, then I(T-w) = I(T-v-w) +
I(T-v-w-N(v)\w) * x (roughly). The mode of this sum is controlled by the
mode of I(T-w), which for a tree of size n-1 is close to mode_f.

### The role of closed-neighborhood mode drop

The gap mode(I(T-w)) - mode(I(T-N[w])) >= 2 means that each neighbor of w
contributes at least one unit of "mode boost." This is consistent with the
DP recurrence: I(T-w) = product of (I(T_i) where T_i are components), and
including w in an independent set forces exclusion of neighbors, reducing
the product. The mode of I(T-w) incorporates contributions from all deg(w)
components, while I(T-N[w]) loses those contributions.

### Possible proof approach

The data suggests a proof by considering the DP structure at the vertex w:

I(T) = dp0[w] + dp1[w]
     = I(T-w evaluated with w excluded) + x * I(T-N[w])

The mode of I(T-w) can reach d_f = mode_f + 1 only when the polynomial
product over components of T-w has its peak shifted just enough relative
to I(T). The constraint is that each factor in the product (one per
component of T-w) must align their modes so the convolution peaks at
exactly d_f.

For the inequality mode(I(T-w)) <= d_f, it suffices to show that the
convolution of component polynomials cannot push the mode past d_f.
Since I(T) = I(T-w) + x * I(T-N[w]), the coefficients satisfy
i_k(T) = i_k(T-w) + i_{k-1}(T-N[w]). At k = d_f, we need
i_{d_f}(T-w) <= i_{d_f}(T) = i_{d_f-1}(T) - margin, which relates
the mode constraint to the descent margin of I(T).

## Files

- `results/tight_mode_cases.json`: full enumeration results
- `scripts/tight_mode_analysis.py`: main enumeration script
- `scripts/tight_mode_deeper.py`: detailed case-by-case analysis (n <= 12)
- `scripts/tight_mode_classify.py`: palindromic vs non-palindromic classification
- `scripts/tight_mode_largn.py`: property verification for n = 13..16
- `scripts/tight_mode_verify_df.py`: confirmation that d_f - mode_f = 1 for n = 13..16
