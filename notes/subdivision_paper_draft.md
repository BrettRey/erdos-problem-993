# Subdivision Lemma: Draft for Paper

## Key formula

When edge uv of tree T is subdivided (inserting vertex w between u and v) to form T':

I(T'; x) = I_u(x) · I_v(x) + x · R_u(x) · R_v(x)

where T_u, T_v are the subtrees on either side of edge uv, I_u = IS polynomial of T_u,
and R_u = IS polynomial of T_u with u excluded.

Equivalently: I(T') = I(T) + A(x) where A(x) = P_u P_v + x R_u R_v
(P_u = IS polynomial of T_u with u included).

The term P_u P_v counts IS of T' that include both u and v (impossible in T since uv was an edge);
the term x R_u R_v counts IS that include w (forcing u, v out).

## Conditional theorem

**Theorem (Subdivision Lemma, conditional).** Let T be a tree with unimodal IS polynomial
having first descent at position d. If

(C1) A(x) is unimodal, and
(C2) first_descent(A) ≥ d,

then I(T') is unimodal.

**Proof.** Write d_A = first_descent(A). I(T') = I + A with I unimodal (mode d) and A unimodal (mode d_A ≥ d).

*Case d_A = d:* Both I and A are nondecreasing for k < d and nonincreasing for k ≥ d.
Their sum inherits both properties. Unimodal with mode d.

*Case d_A = d + j for j ≥ 1:* For k ≤ d-1, both ascending; for k ≥ d+j, both descending.
The uncertain zone is {d, d+1, ..., d+j-1} — at most j positions.

For j = 1: one uncertain position (k = d). The sign pattern of Δ(I+A) is
+...+, ?, -...-. Regardless of ?, at most one descent-to-ascent transition. Unimodal.

For j ≥ 2: the uncertain zone has ≥ 2 positions. A valley at (?_i, ?_{i+1}) = (-, +)
would violate unimodality. This must be ruled out separately. □

## Empirical verification

| Check | Edges | n range | Failures |
|-------|-------|---------|----------|
| A(x) unimodal (C1) | 24,xxx,xxx | 3–20 | 0 |
| first_descent(A) ≥ d (C2) | 24,xxx,xxx | 3–20 | 0 |
| I(T') unimodal | 24,xxx,xxx | 3–20 | 0 |
| I(T') LC | 24,xxx,xxx | 3–20 | 0 |
| Combined tail (I' desc. from d+1) | 24,xxx,xxx | 3–20 | 0 |

mode(A) - d distribution: +0 (26.6%), +1 (73.3%), +2 (0.09%).

For j = 1 cases (73.3%): the conditional theorem applies directly.
For j = 0 cases (26.6%): trivial.
For j = 2 cases (0.09% = 8,405 edges through n=19):
  - Sign pattern always (+, -) — never a valley. The combined tail
    (I+A nonincreasing from d+1) holds with large margin.
  - Max rise/drop ratio: 0.164.

## HI-Hub Lemma

**Lemma.** Every homeomorphically irreducible tree on n ≥ 3 vertices has max d_leaf ≥ 2.

**Proof.** Let T be HI with n ≥ 3. Every internal vertex has degree ≥ 3.

Case 1: T has only one internal vertex v. Then all other vertices are leaves adjacent to v,
so d_leaf(v) = n-1 ≥ 2.

Case 2: T has k ≥ 2 internal vertices. They form a connected subtree T_I with k ≥ 2 vertices.
T_I has at least 2 leaves. Let w be a leaf of T_I. Then w has exactly 1 internal neighbor
(in T_I) and d_leaf(w) leaf children. Since T is HI, deg(w) = 1 + d_leaf(w) ≥ 3,
so d_leaf(w) ≥ 2. □

**Corollary.** There are no HI trees with d_leaf ≤ 1 on n ≥ 3 vertices.

Verified: among 23,055 HI trees through n = 22, exactly 0 have d_leaf ≤ 1 (except P_2 at n = 2).

## Synthesis with PNP framework

The subdivision lemma and PNP framework attack unimodality from complementary directions:

1. **Subdivision lemma (if proved):** Any minimal counterexample is HI.
2. **HI-Hub lemma (proved):** Every HI tree on n ≥ 3 has a vertex with d_leaf ≥ 2.
3. **Hub Exclusion + Transfer (proved):** For any tree with a d_leaf ≥ 2 vertex,
   the PNP reduction applies, reducing to Conjecture A on residual d_leaf ≤ 1 components.

Together: a minimal counterexample would need to be an HI tree where the PNP reduction's
residual components violate Conjecture A. Since every HI tree has a d_leaf ≥ 2 hub,
the reduction always starts. The residual has fewer vertices and d_leaf ≤ 1 everywhere.

**If Conjecture A holds** (verified for 528,196 d_leaf ≤ 1 trees through n = 23, 0 violations):
No HI tree can be a counterexample. Hence no tree is a counterexample.

## What remains for a full proof

1. **Prove the subdivision lemma:** Requires proving C1 (A unimodal) and C2 (mode(A) ≥ d)
   algebraically, plus closing the 0.09% gap (j = 2 cases). The combined tail condition at d+1
   handles the gap empirically but needs a proof.

2. **Prove Conjecture A:** mode ≤ floor(n/3)+1 for all d_leaf ≤ 1 trees. The hard-core model
   gives μ < n/3 (verified), and for LC sequences mode ≤ ceil(μ), but LC of d_leaf ≤ 1 trees
   is itself unproved (though verified through n = 22).

Either result alone would be significant. Both together close the conjecture.
