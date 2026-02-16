# Core Average Verification Summary

## Task

Verify the hypothesis: For trees with d_leaf(v) ≤ 1 (every vertex has at most one leaf-child), the core vertices C (vertices with no leaf children and degree > 1) satisfy:

**Σ_{v∈C} P(v) < |C|/3**

where P(v) = occupation probability = P(v is in a random independent set).

## Method

1. Enumerate all trees n = 3 to 18 using nauty's `geng`
2. Filter to d_leaf ≤ 1 trees (at most one leaf-child per vertex)
3. For each tree:
   - Classify vertices into L (leaves), S (support: exactly 1 leaf-child), C (core: 0 leaf-children, not leaves)
   - Compute exact occupation probabilities P(v) via independence polynomial
   - Check: is core_sum = Σ_{v∈C} P(v) < |C|/3?
   - Check: does every heavy core vertex (P > 1/3) have a light core neighbor (P ≤ 1/3)?

## Results

### Summary Table

| n | d_leaf≤1 trees | Trees with core | Core violations | Max core avg | Heavy-light violations |
|---|----------------|-----------------|-----------------|--------------|------------------------|
| 3 | 0 | 0 | 0 | - | 0 |
| 4 | 1 | 0 | 0 | - | 0 |
| 5 | 1 | 1 | 0 | 0.307692 | 0 |
| 6 | 2 | 1 | 0 | 0.285714 | 0 |
| 7 | 3 | 3 | **1** | **0.333333** | 0 |
| 8 | 6 | 4 | **1** | **0.338983** | **1** |
| 9 | 10 | 10 | **2** | **0.360000** | **3** |
| 10 | 21 | 18 | **2** | **0.357143** | **3** |
| 11 | 39 | 39 | **6** | **0.380282** | **11** |
| 12 | 82 | 76 | **8** | **0.371795** | **19** |
| 13 | 167 | 167 | **18** | **0.400990** | **52** |
| 14 | 360 | 349 | **31** | **0.386503** | **106** |
| 15 | 766 | 766 | **65** | **0.416096** | **256** |
| 16 | 1692 | 1669 | **124** | **0.399577** | **571** |
| 17 | 3726 | 3726 | **255** | **0.431361** | **1362** |
| 18 | 8370 | 8323 | **470** | **0.410608** | **3139** |

**Core violation:** core_avg ≥ 1/3
**Heavy-light violation:** At least one heavy core vertex has no light core neighbor

### Conclusion

**Both approaches FAIL:**

1. **Core average bound (Σ_C P(v) < |C|/3)**: FAILS starting at n=7
   - First violation: n=7, core_avg = 1/3 (exactly, not strictly less)
   - By n=17: max core_avg reaches 0.431361 (significantly above 1/3)
   - Violation rate grows: 5.6% at n=18

2. **Heavy-light matching**: FAILS starting at n=8
   - Many trees have ALL core vertices heavy (P > 1/3)
   - Core vertices often form independent sets (no edges between them)
   - Violation rate: 37.7% at n=18

## Worst-Case Structure

The trees achieving maximum core average have a **remarkably simple structure**:

### Odd n (|C| = 1)
- Single degree-2 core vertex connecting two support vertices
- Each support vertex has one leaf-child and connects to high-degree support subtrees
- Example (n=17): deg = [5, 5, 2, 2, 2, 2, 2, 2, 2, 1, ...]
  - Two degree-5 support hubs connected by a degree-2 core vertex

### Even n (|C| = 2)
- Two degree-2 core vertices (NOT adjacent to each other)
- Each connects exactly two support vertices
- Core vertices have 0 core neighbors (isolated in core subgraph)
- Example (n=18): deg = [5, 4, 3, 2, 2, 2, 2, 2, 2, 2, 1, ...]

### Key Properties
1. All worst-case core vertices have degree exactly 2
2. Core vertices connect ONLY to support vertices (never to other core)
3. Core vertices form independent sets in the core-induced subgraph
4. Structure: high-degree support "hubs" (with many leaves) bridged by degree-2 core vertices

## Implications

The edge bound P(u) + P(v) ≤ 2/3 successfully constrains leaf+support vertices:
- **Σ_{L∪S} P(v) ≤ |L∪S|/3** ✓ (always holds)

But this is **insufficient** to force core vertices below average 1/3:
- **Σ_C P(v) < |C|/3** ✗ (FAILS; can reach ~0.43·|C|)

## Next Research Directions

1. **Compute asymptotic limit**: What is lim_{n→∞} max_v P(v) for degree-2 core vertices?

2. **Upper bound on core vertices**: Can we prove P(v) < 1/2 for all core vertices?
   - If yes: Σ_C < |C|/2, combined with Σ_{L∪S} ≤ |L∪S|/3, gives bounds

3. **Degree constraints**: In d_leaf ≤ 1 trees, can we prove core vertices must have small degree?
   - Observation: worst cases all have deg(core) = 2
   - Is deg(core) ≤ 2 always? Or just for worst cases?

4. **Alternative proof strategy**:
   - Instead of bounding μ directly, use mode-mean gap
   - Or use different decomposition (not L/S/C split)

5. **Structural analysis**:
   - Core vertices being isolated (independent set) is interesting
   - Can we exploit tree structure to bound individual P(v) better?

## Files

- **Script**: `verify_core_avg.py` (exhaustive verification)
- **Analysis**: `analyze_worst_core.py` (structural analysis of extremal trees)
- **Notes**: `notes/core_avg_analysis.md` (detailed findings)
- **This summary**: `VERIFICATION_CORE_AVG.md`

Run time: ~19 seconds for n ≤ 18 (8370 d_leaf ≤ 1 trees verified)

---

**Date**: 2026-02-15
**Verification**: Complete enumeration via `geng` + exact IS polynomial computation
**Status**: Hypothesis REFUTED; alternative approaches needed
