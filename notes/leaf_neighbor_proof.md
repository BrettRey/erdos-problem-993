# Leaf-Neighbor Property: Analysis and Correction

## Date: 2026-02-14

## Summary

The Leaf-Neighbor Property (LNP) as originally stated **fails at n=18**.
The correct statement is the **Private Neighbor Property (PNP)**, which
holds through n=18 (99,147 cases, 0 failures).

## The Original Claim (LNP) — FALSE

**Claim:** For any tree T and any maximal IS S of size k < mode(I(T)),
some u in S has >= 2 leaf-neighbors (degree-1 vertices adjacent to u,
outside S).

**Status:** True for n <= 16 (13,042 cases). **Fails at n=18.**

### Counterexample

Tree: Spider with center vertex 17 (degree 5), four arms of length 2
(vertices 4-0, 5-1, 6-2, 7-3 through the center), and pendant leaf 8.

```
9--0--13--4--17--5--14--1--10
                 |
          11--2--15--6--17  (same 17)
                 |
          12--3--16--7--17  (same 17)
                 |
                 8
```

graph6: `Q??????_A?C?C?a?G_@C?CO?N_?`

- n=18, mode=6 (i_6 = 2316)
- Independence polynomial: [1, 18, 136, 566, 1430, 2281, 2316, 1462, 525, 82]
- Maximal IS S = {0, 1, 2, 3, 17}, size 5 < 6 = mode
- Every vertex in S has exactly 1 leaf-neighbor:
  - u=0: leaf-neighbor 9 only (neighbor 13 has degree 2)
  - u=1: leaf-neighbor 10 only
  - u=2: leaf-neighbor 11 only
  - u=3: leaf-neighbor 12 only
  - u=17: leaf-neighbor 8 only (neighbors 4,5,6,7 have degree 2)

### Why the LNP Held Through n=16

The tree has no MLSV (no vertex adjacent to >= 2 leaves). Through n=16,
all trees either (a) have no maximal IS below mode, or (b) have at
least one MLSV in every maximal IS below mode. The n=18 spider is the
first tree where (a) fails for a no-MLSV tree.

### Route C' Analysis (Corrected)

Route C' (S contains a support vertex with >= 2 leaf-children) covers
100% of cases through n=16, but this is because n=16 is too small
for the counterexample structure to appear.

## The Correct Statement (PNP) — Empirically Verified

**Conjecture (Private Neighbor Property):** For any tree T and any
maximal IS S of size k < mode(I(T)), some u in S has priv(u) >= 2.

A vertex v not in S is *private for u in S* if u is v's only S-neighbor.

### Verification

| n range | Cases | PNP failures | Pigeonhole covers | n < 3k cases |
|---------|-------|-------------|-------------------|-------------|
| 5-13    | 669   | 0           | 669 (100%)        | 0           |
| 14      | 1,481 | 0           | 1,471 (99.3%)     | 10          |
| 15      | 2,601 | 0           | 2,601 (100%)      | 0           |
| 16      | 8,291 | 0           | 8,281 (99.9%)     | 10          |
| 17      | 28,959| 0           | 28,584 (98.7%)    | 375         |
| 18      | 57,146| 0           | 57,136 (100%)     | 10          |
| **Total** | **99,147** | **0** | **98,742 (99.6%)** | **405** |

### For the n=18 Counterexample Tree

S = {0, 1, 2, 3, 17}:
- u=0: priv=2 (private neighbors 9, 13)
- u=1: priv=2 (private neighbors 10, 14)
- u=2: priv=2 (private neighbors 11, 15)
- u=3: priv=2 (private neighbors 12, 16)
- u=17: priv=5 (private neighbors 4, 5, 6, 7, 8)

max_priv = 5. The swap is trivially available: remove u=17, add any
two of {4,5,6,7,8} to get an IS of size 6 = mode.

## Proof Structure

### What Is Proved

1. **Private Neighbor Bound:** For any dominating IS S of size k in a
   tree on n vertices: P >= n - 2k + 1, where P = total private neighbors.
   (Proved via edge-counting.)

2. **Pigeonhole (n >= 3k):** If P >= k+1, some u has priv(u) >= 2.
   This holds when n >= 3k (since P >= n-2k+1 >= k+1).

3. **Tree Swap Lemma:** If u in S has 2 private neighbors v, w in a tree,
   then v, w are non-adjacent (acyclicity). So (S \ {u}) ∪ {v, w} is
   independent of size k+1.

### What Remains Open

The **n < 3k gap**: 405 cases (0.41%) where pigeonhole doesn't apply.
All satisfy PNP empirically but lack a proof. These are predominantly
double-star trees where one hub vertex concentrates many private neighbors.

To close the gap, need to prove: when n < 3k and k < mode(I(T)), some
u in S has deg(u) >= k+1 (which gives priv(u) >= deg(u) - (k-1) >= 2
by the shared-neighbor bound).

Alternatively: characterize trees where n < 3k and k < mode is possible,
and show PNP holds for this restricted class.

## Proof Routes Tested

### Route A: Pigeonhole via |L_out| > k
- Coverage: 98.8% (12,886/13,042 through n=16)
- Fails for 156 cases where |L_out| <= k
- Not sufficient alone

### Route B: Degree bound max deg >= k+1
- Coverage: 84.3% (10,993/13,042 through n=16)
- Fails for 2,049 cases
- Not sufficient alone

### Route A ∪ B: Combined
- Coverage: 99.96% (13,037/13,042)
- 5 uncovered cases (all covered by Route C')

### Route C': Support vertex with >= 2 leaf-children in S
- Coverage: 100% through n=16 (13,042/13,042)
- **FAILS at n=18** (counterexample above)
- Not the correct universal property

### PNP (Private Neighbor Property)
- Coverage: 100% through n=18 (99,147/99,147)
- The correct conjecture
- Proved for 99.6% of cases (pigeonhole when n >= 3k)
- Open for 0.4% (n < 3k cases)

## Implications for the SRI

The sign-reversing involution (containment-first strategy) works through
n=16 using the LNP for the swap step. At n=18, the LNP fails, so the
SRI as designed (canonical leaf-swap triple) would need modification.

The PNP guarantees a swap exists (via private neighbors), but the SRI
needs a *canonical* choice. Options:
1. Use max-degree vertex in S (as before) but select private neighbors
   instead of leaf-neighbors
2. The canonical triple becomes (u, v, w) where u is highest-degree
   in S, and v, w are the two smallest-label private neighbors of u
3. Need to reverify the SRI with this modified canonical triple

## Files

- `prove_leaf_neighbor.py`: Comprehensive route analysis (n <= 16)
- `verify_pnp_extended.py`: PNP verification through n=18
- `investigate_n18_counter.py`: Detailed analysis of LNP counterexample
- `analyze_support_vertex.py`: MLSV avoidance hypothesis (wrong direction)
- `analyze_lnp_mechanism.py`: Forced leaves vs mode analysis
- `analyze_lnp_deep.py`: Gap tree analysis
- `analyze_lnp_gap.py`: Detailed gap tree cascade analysis
- `analyze_lnp_alpha.py`: Min maximal IS avoiding MLSVs
- `analyze_no_mlsv.py`: No-MLSV tree analysis (found n=18 counterexample)
