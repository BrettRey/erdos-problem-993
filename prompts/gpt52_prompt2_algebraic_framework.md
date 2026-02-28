# Task: Assess which algebraic/combinatorial frameworks can prove unimodality of tree independence polynomials

## Background

Erdős Problem #993 (1987): Is the sequence of independent set counts (i_0, i_1, ..., i_α) of every tree unimodal? Verified computationally for all ~1.2 billion trees on n ≤ 27 vertices, but no proof exists. This is one of the oldest open problems in algebraic graph theory.

## Definitions

**Independence polynomial.** For a graph G on n vertices:

  I(G; x) = Σ_{k=0}^{α(G)} i_k x^k

where i_k counts independent sets of size k (i_0 = 1). A sequence is unimodal if it has no interior valley.

**Log-concavity (LC).** A positive sequence is LC if a_k^2 ≥ a_{k-1}·a_{k+1} for all k. LC ⟹ unimodal.

**Multivariate independence polynomial.** For a graph G = (V, E):

  I(G; x_1, ..., x_n) = Σ_{S independent} Π_{v ∈ S} x_v

The univariate I(G; x) is the specialization x_v = x for all v.

## Tree-specific structure

For a tree T rooted at any vertex r, the DP gives:

  dp[v][0] = Π_{c child of v} (dp[c][0] + dp[c][1])     (v not in IS)
  dp[v][1] = x · Π_{c child of v} dp[c][0]               (v in IS)
  I(T; x) = dp[r][0] + dp[r][1]

Key structural facts:
- dp[c][0] + dp[c][1] = I(T_c), the IS polynomial of the subtree at child c.
- I(T) is built from **products** (over children) and **sums** (dp[v][0] + dp[v][1]).
- Each factor I(T_c) is itself a tree IS polynomial.
- For a leaf v: dp[v][0] = 1, dp[v][1] = x, so I(leaf) = 1 + x.

## What is known (status as of Feb 2027)

### Results that WORK:
1. **Computational verification:** All 1,198,738,056 trees n ≤ 27 are unimodal. Zero violations.
2. **LC for small trees:** Tree IS polynomials are LC for ALL trees n ≤ 25 (hundreds of millions). At n = 26, exactly 2 trees fail LC. Galvin (2025) constructed infinite families with LC failures. But unimodality holds regardless.
3. **d_leaf ≤ 1 trees are strictly LC:** Trees where every vertex is a leaf or adjacent to a leaf are strictly LC through n = 23 (931,596 trees, 0 failures).
4. **Mode-mean bound:** mode ≤ ⌈μ⌉ where μ = I'(1)/I(1), verified for 9.1 million trees n ≤ 22 (0 failures). This would close the main open conjecture if proved.
5. **Edge bound (proved):** For any tree edge uv, the hard-core occupation probabilities satisfy P(u) + P(v) < 2/3.
6. **Steiner peeling (proved):** μ(T) < n/3 for all d_leaf ≤ 1 trees (via belief propagation cavity fields).

### Results that FAIL:
1. **Real-rootedness:** FAILS for trees. Chudnovsky-Seymour (2007) proved real-rootedness only for claw-free graphs, but any tree with a vertex of degree ≥ 3 contains K_{1,3}. So "real roots ⟹ LC ⟹ unimodal" is dead.
2. **Ultra-log-concavity:** FAILS at n = 8.
3. **Stable polynomials:** Choe-Oxley-Sokal-Wagner (2004) proved the multivariate IS polynomial of claw-free graphs has the half-plane property. Trees contain claws, so this doesn't apply.
4. **LC for all trees:** FAILS at n = 26 (2 trees). So unimodality must be proved by something weaker than LC.
5. **Mode superadditivity for products:** mode(f·g) ≥ mode(f) + mode(g) is FALSE for products of LC polynomials.

## Your task

Survey the landscape of algebraic and combinatorial techniques and give an honest assessment of which could prove unimodality (or an equivalent structural property) of tree IS polynomials. For each technique below, I want:

**(a)** Does it have a realistic chance of working for trees?
**(b)** What's the specific obstacle or gap?
**(c)** If it could work, sketch how the proof might go.

### Techniques to assess:

1. **Lorentzian polynomials (Brändén-Huh 2020).** The basis generating polynomial of a matroid is Lorentzian, which implies LC. Can the multivariate IS polynomial of a tree (or some transform of it) be shown to be Lorentzian? Note: IS polynomials of trees are NOT matroid basis generating polynomials, so this requires a new connection. Can homogenization or some change of variables help?

2. **Interlacing families (Marcus-Spielman-Srivastava 2015).** The IS polynomial of a tree has a recursive structure (vertex elimination: I(T) = I(T-v) + x·I(T-N[v])). Can common interlacing along this recursion prove anything? Note that real-rootedness fails, so classical interlacing won't directly give LC, but partial interlacing or interlacing of specific families might still give unimodality.

3. **Sector-stable polynomials and related notions.** Are there stability-type conditions weaker than half-plane property that tree IS polynomials satisfy? The Lee-Yang theorem gives zero-free regions for certain partition functions. What's known about the zero-free region of tree IS polynomials, and can it imply coefficient properties?

4. **Hodge theory / Lefschetz type.** Adiprasito-Huh-Katz (2018) used Hodge theory on the Chow ring of a matroid. Is there a graded algebra naturally associated to trees where the IS counts are dimensions of graded pieces, and where hard Lefschetz could give LC or unimodality?

5. **Combinatorial Atlas / convex geometry.** Adiprasito's combinatorial atlases generalize the Hodge-theoretic approach. Do they apply to independence complexes of trees?

6. **Recent work (2023-2026).** Are there new results on IS polynomials of trees or closely related families that I may have missed? Specifically:
   - Bencs, Csikvári, or Regts on zeros of IS polynomials
   - Extensions of Lorentzian framework beyond matroids
   - New results on the hard-core model partition function on trees
   - Any progress on Erdős Problem #993 itself

7. **Any other framework** you think has a realistic shot that I haven't listed.

## What I value

- **Honest negatives.** "This doesn't work because X" is far more valuable than optimistic handwaving. I've spent months on this problem and need to allocate effort wisely.
- **Specific obstacles.** Not "it might be hard to show stability" but "the specific issue is that the polynomial (1+x)^2 + x violates [property] because [reason]."
- **Actionable next steps.** If a framework looks promising, what's the concrete first thing to check or prove?
