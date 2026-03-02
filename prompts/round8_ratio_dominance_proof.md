# Round 8: Proving E ≽ J (Ratio Dominance) at Support Vertices

## The Problem

We're proving that the independence polynomial of every tree is unimodal (Erdős Problem #993). We've reduced the full problem to a single algebraic statement:

**At every support vertex r of every tree T, E ≽ J (ratio dominance) holds.**

Here E = dp[r][0], J = dp[r][1]/x, and E ≽ J means E_{k+1}·J_k ≥ E_k·J_{k+1} for ALL k ≥ 0.

## Why E ≽ J suffices for unimodality

**Proved chain:**

1. **P3 (tail domination)**: e_k ≥ j_{k-1} for all k. PROVED via leaf-swap injection at support vertices.

2. **SCC from E ≽ J + LC(E)**: Define Δ_k = ((1+x)I)_{k+1}·E_k - ((1+x)I)_k·E_{k+1}. Then:
   - Δ_k = c_k + LR_k (2-term decomposition, PROVED algebraically)
   - c_k = E_k² - E_{k-1}·E_{k+1} ≥ 0 (LC gap of E, PROVED: products of LC polynomials are LC)
   - LR_k = E_k·(J_k+J_{k-1}) - E_{k+1}·(J_{k-1}+J_{k-2}) = Δ_k(E, (1+x)J) (LR minor)
   - If E ≽ J, then by Karlin's TP2 theorem (since (1+x) is PF2), E ≽ (1+x)J... **NO, this FAILS** (~12% of checks).
   - However, SCC = c_k + LR_k ≥ 0 is verified independently (907M+ checks, 0 failures).

3. **SCC ⟹ unimodality**: Via the identity b_{k-1}·Δ_k = b_{k-1}·d_k + b_k·d_{k-1} + a_{k-1}·c_k.

So the proof strategy is: prove BOTH E ≽ J AND SCC ≥ 0 by induction on the product structure.

## Product Structure at Support Vertices

Root T at support vertex r with ℓ leaf neighbors and non-leaf children c₁,...,c_s.

```
E = (1+x)^ℓ · ∏ I(T_{c_j})    (exclude-root: each subtree contributes its IS poly)
J = ∏ E(T_{c_j})                (include-root: each subtree contributes its exclude-root poly)
```

**Incremental product:**
```
Stage 0:  E^{(0)} = (1+x)^ℓ,  J^{(0)} = [1]
Stage t:  E^{(t)} = E^{(t-1)} · I_t,  J^{(t)} = J^{(t-1)} · E_t
```
where I_t = E_t + x·J_t is the IS poly of the t-th non-leaf subtree.

## What's PROVED

### 1. LC(E) at every stage
E^{(t)} is a product of PF2 (nonneg + LC) polynomials. Products of PF2 polynomials are PF2 (Karlin). So E is LC at every stage. **PROVED.**

### 2. J ≤ E coefficientwise at every stage
J^{(t)} = J^{(t-1)} · E_t ≤ E^{(t-1)} · E_t (by induction: J^{(t-1)} ≤ E^{(t-1)}) and E^{(t-1)} · E_t ≤ E^{(t-1)} · I_t = E^{(t)} (since E_t ≤ I_t). So J^{(t)} ≤ E^{(t)}. **PROVED.**

### 3. Karlin main part of E ≽ J
At stage t, define A = E^{(t-1)} · E_t. Then:
- A ≽ J^{(t)} by Karlin's TP2 theorem: E^{(t-1)} ≽ J^{(t-1)} (induction) and E_t is PF2.
- E^{(t)} = A + x · E^{(t-1)} · J_t (since I_t = E_t + x·J_t)
- So Δ_k(E^{(t)}, J^{(t)}) = Δ_k(A, J^{(t)}) + Δ_k(x·E^{(t-1)}·J_t, J^{(t)})
- First term ≥ 0: **PROVED by Karlin.**
- Second term can be negative (~8% of stages).
- Sum is ALWAYS ≥ 0 (11.9M stages n ≤ 20, 0 failures).

### 4. P3 (tail domination)
e_k ≥ j_{k-1} for all k at support vertices. **PROVED via leaf-swap injection.**

## Computational Verification

| Property | Checks | Failures | Status |
|----------|--------|----------|--------|
| E ≽ J at ALL k | 907M+ (n ≤ 22) | **0** | VERIFIED |
| E ≽ (1+x)J | 907M+ | ~12% fail | Known to fail |
| SCC Δ_k ≥ 0 | 907M+ | **0** | VERIFIED |
| LC(E) | 59.9M support verts | **0** | PROVED |
| J ≤ E | 59.9M | **0** | PROVED |
| E^{(t)} ≽ J^{(t)} every stage | 11.9M stages (n ≤ 20) | **0** | VERIFIED |
| Karlin A ≽ J^{(t)} every stage | 11.9M | **0** | PROVED |
| Correction term negative | 11.9M | ~8% neg | Compensated |

## The Algebraic Gap

At each stage, we need:

**Δ_k(A, J^{(t)}) + Δ_k(x·B, J^{(t)}) ≥ 0**

where A = E^{(t-1)}·E_t (PROVED ≥ 0) and B = E^{(t-1)}·J_t.

Since E^{(t)} = A + x·B, this is just E^{(t)} ≽ J^{(t)}, the inductive claim.

The decomposition shows: the "Karlin part" (from E_t) is always nonneg. The "correction" (from x·J_t) can be negative but is compensated. The question is: WHY is it compensated?

## Structural Constraints (available for the proof)

1. **E_t is LC (PF2)**: the IS exclude-root polynomial of the subtree.
2. **J_t ≤ E_t**: coefficientwise, at the factor level. PROVED.
3. **I_t = E_t + x·J_t is LC**: the IS polynomial of the subtree is always LC? NO, LC can fail for IS polynomials. But I_t is always the IS polynomial of a tree, which has additional structure.
4. **J_t/E_t is nonincreasing**: i.e., E_t ≽ J_t at the factor level? This FAILS (~14% of non-leaf factors). So we CANNOT use factor-level ratio dominance.
5. **SCC at factor level**: Δ_k((1+x)I_t, E_t) ≥ 0. This is the SCC of the subtree. YES, always holds (verified 930M+ checks). This means (1+x)I_t ≽ E_t.
6. **Product-level E ≽ J always holds**: even when factor-level fails. The product structure rescues it.

## s=1 Reduction (NEW, handles 63% of support vertices)

When the support vertex r has exactly one non-leaf child c (and ℓ ≥ 1 leaves), E≽J **follows from SCC of the subtree** via a clean chain:

```
E = (1+x)^ℓ · I_c,  J = E_c

1. SCC at c:  (1+x)·I_c ≽ E_c                    [VERIFIED, 0 fails]
2. Karlin:    (1+x)^{ℓ-1} is PF2, so
              (1+x)^ℓ · I_c ≽ (1+x)^{ℓ-1} · E_c  [by Karlin on step 1]
3. Trivial:   (1+x)^{ℓ-1} ≽ [1]                   [obvious]
4. Karlin:    E_c is PF2, so
              (1+x)^{ℓ-1} · E_c ≽ E_c              [by Karlin on step 3]
5. Transitivity: E ≽ J                             [chain of 2 and 4]
```

**Statistics (n ≤ 18, 1,103,584 support vertices):**

| s | Count | % | Status |
|---|-------|---|--------|
| 0 | 16 | 0.0% | Trivial (star centers) |
| 1 | 692,736 | 62.8% | **PROVED** (SCC reduction) |
| 2 | 325,087 | 29.5% | Key battleground |
| 3 | 73,510 | 6.7% | Follows from s=2 by iteration |
| 4+ | 11,835 | 1.1% | Rare |

**Implication:** The s=1 case (63%) is PROVED conditional on SCC holding universally (verified, 0 fails). The real challenge is **s=2** (30%), which is the simplest non-trivial incremental step. If s=2 is proved, s≥3 follows by repeating the same argument.

**The pendant-star is s=1!** Its E≽J (the tightest case) is provable from SCC of the star K_{1,m}, which itself reduces to binomial coefficient identities.

## Key Observation: The x·B Correction

The correction term Δ_k(x·B, J^{(t)}) where B = E^{(t-1)}·J_t can be written:

Δ_k(x·B, J^{(t)}) = B_k · J^{(t)}_k - B_{k-1} · J^{(t)}_{k+1}

This uses J_t (the "smaller" part of I_t) multiplied by E^{(t-1)} (the accumulated product). When J_t is small relative to E_t, the correction is small relative to the Karlin main part. The constraint J_t ≤ E_t ensures this.

## Extremal Analysis: Pendant-Star

The tightest tree for the E ≽ J correction ratio is the **pendant-star**: one leaf u₁ attached to support vertex r, which is attached to star center c, which has m = n-3 pendant leaves.

At the single incremental stage (root at r, one non-leaf child c = star K_{1,m}):
- E_c = (1+x)^m, J_c = [1], I_c = (1+x)^m + x
- A = (1+x)^{m+1}, B = (1+x), J^{(1)} = (1+x)^m

At k=1:
```
main = C(m+1,2) · m - (m+1) · C(m,2) = m(m+1)/2
corr = 1·m - 1·C(m,2) = m(3-m)/2           (negative for m ≥ 4)
total = main + corr = 2m
ratio = (m+1)/(m-3) = (n-2)/(n-6)
```

**EXACT match with computational data** at every n from 7 to 20:
- n=10: 8/4 = 2.0 ✓, n=15: 13/9 = 1.444 ✓, n=20: 18/14 = 9/7 = 1.286 ✓

**Key properties:**
- **Ratio → 1** as n → ∞ (infimum = 1, not achieved). Ratio-based approaches DEAD.
- **Absolute margin = 2m = 2(n-3) → ∞**. Linearly growing gap.
- **Same extremal family as SCC**: pendant-star is tightest for both E ≽ J and SCC.
- The ratio (n-2)/(n-6) ≈ 1 + 4/(n-6) decays polynomially, NOT exponentially.

Any proof must exploit the absolute margin growth (Ω(n)), not the vanishing ratio.

## Diagonal Convolution Structure

The incremental step acts on (E, J) by:
```
E_new = I_t * E_old    (I_t = E_t + x·J_t)
J_new = E_t * J_old
```

This is a DIAGONAL convolution matrix M_t = [[I_t, 0], [0, E_t]]. Both diagonal entries are PF2 (E_t is always PF2; I_t is PF2 when the subtree's IS poly is LC, which holds for most but not all trees).

**Key observation:** If I_t ≽ E_t (subtree SCC) AND E_old is PF2, then by Karlin:
- E_old * I_t ≽ E_old * E_t (Karlin with PF2 kernel E_old)
- E_t * E_old ≽ E_t * J_old (Karlin with PF2 kernel E_t, using E_old ≽ J_old)

By transitivity of ≽: E_new = I_t * E_old ≽ E_old * E_t ≽ E_t * J_old = J_new.

**BUT**: I_t ≽ E_t FAILS at ~30% of factors! So the simple transitivity argument fails.

**Identity:** Δ_k(I_t, E_t) = (J_t)_k·(E_t)_k - (J_t)_{k-1}·(E_t)_{k+1}. This is the "J ascending relative to E" condition (ratio J_t/E_t nondecreasing), OPPOSITE to E_t ≽ J_t (ratio nonincreasing). Both can be true simultaneously (the subtree-level ratio increases while the root-level ratio decreases), but the subtree-level one fails for ~30% of individual factors.

**Despite this**, E_new ≽ J_new holds universally. The product structure of multiple factors must rescue the property. This is why the proof is hard: no single-factor argument works.

## What I Want from You

### Priority 1: Prove E ≽ J inductively

Show that if E^{(t-1)} ≽ J^{(t-1)} and LC(E^{(t-1)}) and J^{(t-1)} ≤ E^{(t-1)}, then after multiplying by factor (I_t, E_t):

E^{(t)} = E^{(t-1)} · I_t  ≽  J^{(t-1)} · E_t = J^{(t)}

Given: E_t is PF2 (LC), J_t ≤ E_t, and (1+x)I_t ≽ E_t (SCC of subtree).

The Karlin part (from E_t) is proved ≥ 0. Need: the x·J_t correction is bounded by the Karlin part.

### Priority 2: Alternative characterizations

- Is E ≽ J equivalent to some total positivity condition on a matrix built from E and J?
- Is there a direct combinatorial proof that the ratio J_k/E_k is nonincreasing?
- Combinatorially: J_k/E_k = #{IS of size k+1 containing r} / #{IS of size k not containing r}. Why does this ratio decrease in k?

### Priority 3: Unified proof of E ≽ J and SCC

Both properties have the same structure at each stage:
- Karlin main part (always ≥ 0, PROVED)
- Correction from x·J_t (can be negative, needs compensation)

Can they be proved simultaneously using a combined invariant?

## Notation Summary

- T: tree on n vertices
- r: support vertex (adjacent to at least one leaf)
- E = dp[r][0]: exclude-root polynomial (IS poly of T excluding r from all sets)
- J = dp[r][1]/x: include-root polynomial divided by x (IS poly of T including r, shifted)
- I = E + xJ: the full IS polynomial of the rooted tree
- E_k, J_k: coefficient of x^k
- E ≽ J: E_{k+1}·J_k ≥ E_k·J_{k+1} for all k (ratio dominance / likelihood ratio order)
- PF2: polynomial with nonneg, log-concave coefficients (Pólya frequency of order 2)
- Karlin's theorem: if A is PF2 and f ≽ g, then A·f ≽ A·g (convolution preserves LR order)

## Data for Verification

Tree independence polynomial DP:
```
dp[leaf][0] = [1], dp[leaf][1] = [0,1]
dp[v][0] = ∏_children (dp[c][0] + dp[c][1])
dp[v][1] = x · ∏_children dp[c][0]
```
At support vertex r: E = dp[r][0], J_k = dp[r][1]_{k+1} (shift by 1).

Example: Path P_5 = 0-1-2-3-4, root at vertex 1 (support vertex, adjacent to leaf 0):
- E = [1, 4, 4, 1], J = [1, 2]
- I = E + xJ = [1, 5, 6, 1]
- Check E ≽ J: E_1·J_0 - E_0·J_1 = 4-2=2 ≥ 0 ✓; E_2·J_1 - E_1·J_2 = 8-0=8 ≥ 0 ✓

Example: Pendant-star (n=7): leaf-r-center-{4 leaves}, root at r:
- E = (1+x)·[(1+x)^4 + x] = [1, 6, 11, 10, 5, 1], J = (1+x)^4 = [1, 4, 6, 4, 1]
- At k=1: main = 4·5/2 = 10, corr = 4·(-1)/2 = -2, total = 8 = 2m, ratio = 5.0
