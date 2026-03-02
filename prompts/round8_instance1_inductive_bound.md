# Round 8, Instance 1: Bounding the correction term in E ≽ J induction

## The Problem

We're proving that the independence polynomial of every tree is unimodal (Erdős Problem #993). We've reduced the full problem to a single algebraic statement at support vertices.

**At every support vertex r of every tree T, E ≽ J (ratio dominance) holds.**

E = dp[r][0], J = dp[r][1]/x, and E ≽ J means E_{k+1}·J_k ≥ E_k·J_{k+1} for ALL k ≥ 0.

Verified: 907M+ checks across 9.1M trees n ≤ 22, 0 failures.

## Reduction chain (all PROVED)

1. P3 (tail domination): e_k ≥ j_{k-1}. PROVED (leaf-swap injection).
2. SCC ≥ 0: Δ_k = c_k + LR_k where c_k = LC gap of E ≥ 0. VERIFIED (907M+ checks, 0 fails).
3. SCC ⟹ unimodality via 3-term identity.

## Product structure at support vertices

Root T at support vertex r with ℓ leaf neighbors and non-leaf children c₁,...,c_s.

```
E = (1+x)^ℓ · ∏ I(T_{c_j}),   J = ∏ E(T_{c_j})
```

Incremental:
```
Stage 0: E^{(0)} = (1+x)^ℓ,  J^{(0)} = [1]
Stage t: E^{(t)} = E^{(t-1)} · I_t,  J^{(t)} = J^{(t-1)} · E_t
```
where I_t = E_t + x·J_t.

## The inductive step

At stage t, define A = E^{(t-1)} · E_t (Karlin main part). Then:

Δ_k(E^{(t)}, J^{(t)}) = Δ_k(A, J^{(t)}) + Δ_k(x·B, J^{(t)})

where B = E^{(t-1)} · J_t.

- **First term ≥ 0: PROVED** by Karlin's TP2 theorem (E^{(t-1)} ≽ J^{(t-1)} by induction, E_t is PF2).
- **Second term can be negative** (~8% of stages).
- **Sum is ALWAYS ≥ 0** (11.9M stages n ≤ 20, 0 failures).

## s=1 reduction (PROVED, handles 63%)

When s = 1 (one non-leaf child c, ℓ ≥ 1 leaves):

```
E = (1+x)^ℓ · I_c,  J = E_c

1. SCC at c:   (1+x)·I_c ≽ E_c
2. Karlin:     (1+x)^ℓ · I_c ≽ (1+x)^{ℓ-1} · E_c
3. Trivial:    (1+x)^{ℓ-1} ≽ [1]
4. Karlin:     (1+x)^{ℓ-1} · E_c ≽ E_c
5. Transitivity: E ≽ J
```

| s | % of vertices | Status |
|---|---------------|--------|
| 0 | 0.0% | Trivial |
| 1 | 62.8% | **PROVED** (SCC reduction) |
| 2 | 29.5% | Key battleground |
| 3+ | 7.8% | Follows from s=2 |

The pendant-star (tightest case, ratio = (n-2)/(n-6) → 1) is s = 1 and thus PROVED.

## Your task: bound the correction for s ≥ 2

At stage t of the incremental product, the correction is:

Δ_k(x·B, J^{(t)}) = B_k · J^{(t)}_k - B_{k-1} · J^{(t)}_{k+1}

where B = E^{(t-1)} · J_t. Available constraints:

1. **J_t ≤ E_t** coefficientwise (PROVED)
2. **E_t is PF2** (LC with nonneg coefficients, PROVED)
3. **(1+x)I_t ≽ E_t** (SCC of subtree, VERIFIED 0 failures)
4. **E^{(t-1)} is PF2** (products of PF2 are PF2, PROVED)
5. **J^{(t-1)} ≤ E^{(t-1)}** coefficientwise (PROVED)
6. **E^{(t-1)} ≽ J^{(t-1)}** (inductive hypothesis)

Since B = E^{(t-1)} · J_t ≤ E^{(t-1)} · E_t = A coefficientwise (by constraint 1), we have B ≤ A.

**Key question:** Can you prove |Δ_k(x·B, J^{(t)})| ≤ Δ_k(A, J^{(t)}) when the former is negative?

### Pendant-star data (extremal for s=1)

At n=7 (m=4): E = [1, 6, 11, 10, 5, 1], J = [1, 4, 6, 4, 1].
A = (1+x)·(1+x)^4 = (1+x)^5 = [1,5,10,10,5,1]. B = (1+x)·[1] = [1,1].
- k=1: main = 5·1-1·4 = 1. Wait, J^{(1)} = (1+x)^4.
- Actually: main = A_2·J_1 - A_1·J_2 = 10·4 - 5·6 = 10, corr = B_1·J_1 - B_0·J_2 = 1·4 - 1·6 = -2.
- Total = 8 = 2m. Ratio = 5.0.

At large n (m=n-3): main = m(m+1)/2, corr = m(3-m)/2, total = 2m, ratio = (m+1)/(m-3).

### What makes this hard

- **Ratio → 1**: The correction grows as m²/2 while the main part grows as m²/2. They nearly cancel. But the absolute margin 2m grows linearly. Any bound must exploit absolute margin growth, not the vanishing ratio.

- **Factor-level E_t ≽ J_t FAILS** (~14% of factors). So we cannot argue factorwise.

- **I_t ≽ E_t FAILS** (~30% of factors). So transitivity E_new ≽ E_old·E_t ≽ J_new fails.

### Suggested approaches

1. **Absolute margin induction**: Show that the margin Δ_k(E^{(t)}, J^{(t)}) grows by at least some positive quantity at each stage, rather than trying to bound the ratio.

2. **Use J ≤ E + LC to control correction**: Since B ≤ A and J^{(t)} ≤ E^{(t)}, the correction involves "smaller" polynomials multiplied in the "wrong" (x-shifted) position. The LC gap of E^{(t)} provides quadratic surplus at each index. Can the LC surplus absorb the correction?

3. **Use SCC of subtree (constraint 3)**: (1+x)I_t ≽ E_t means I_t has a specific relationship to E_t. Since E^{(t)} = E^{(t-1)}·I_t and B relates to J_t, the SCC constraint bounds how much J_t can "distort" the product.

## Notation

| Symbol | Definition |
|--------|-----------|
| E ≽ J | E_{k+1}·J_k ≥ E_k·J_{k+1} for all k (ratio dominance) |
| PF2 | nonneg + log-concave coefficients |
| Karlin | if A is PF2 and f ≽ g, then A·f ≽ A·g |
| I_t = E_t + x·J_t | IS polynomial of t-th subtree |
| A = E^{(t-1)}·E_t | Karlin main part (always ≥ 0) |
| B = E^{(t-1)}·J_t | correction source (x·B is the x-shifted correction) |

## Verification data

Path P_5 = 0-1-2-3-4, root at vertex 1:
- E = [1, 4, 4, 1], J = [1, 2]
- Check: E_1·J_0 - E_0·J_1 = 4-2=2 ≥ 0; E_2·J_1 - E_1·J_2 = 8-0=8 ≥ 0. ✓
