# Core Average Analysis for d_leaf ≤ 1 Trees

## Summary

The hypothesis that **Σ_{v∈C} P(v) < |C|/3** (core average strictly below 1/3) **FAILS**.

Starting at n=7, there are d_leaf ≤ 1 trees where the core vertices have average occupation probability ≥ 1/3.

## Key Results (n ≤ 18)

| n | Total d_leaf≤1 | With core | Core violations | Max core avg | Matching violations |
|---|----------------|-----------|-----------------|--------------|---------------------|
| 3-6 | 4 | 2 | 0 | 0.307692 | 0 |
| 7 | 3 | 3 | 1 | 0.333333 | 0 |
| 8 | 6 | 4 | 1 | 0.338983 | 1 |
| 9 | 10 | 10 | 2 | 0.360000 | 3 |
| 10 | 21 | 18 | 2 | 0.357143 | 3 |
| 11 | 39 | 39 | 6 | 0.380282 | 11 |
| 12 | 82 | 76 | 8 | 0.371795 | 19 |
| 13 | 167 | 167 | 18 | 0.400990 | 52 |
| 14 | 360 | 349 | 31 | 0.386503 | 106 |
| 15 | 766 | 766 | 65 | 0.416096 | 256 |
| 16 | 1692 | 1669 | 124 | 0.399577 | 571 |
| 17 | 3726 | 3726 | 255 | 0.431361 | 1362 |
| 18 | 8370 | 8323 | 470 | 0.410608 | 3139 |

**Core violations:** Trees where core_avg ≥ 1/3
**Matching violations:** Trees where at least one heavy core vertex (P > 1/3) has no light core neighbor (P ≤ 1/3)

## Observations

### 1. Core Average Bound Fails

The first violation occurs at n=7:
- Tree with deg sequence [3, 2, 2, 2, 1, 1, 1]
- Core size |C| = 1
- Core sum = 0.333333 (exactly 1/3)
- Core average = 1/3 (not strictly less than 1/3)

By n=18:
- 470/8323 trees with core (5.6%) violate core_avg < 1/3
- Maximum core average reaches 0.410608 (significantly above 1/3)

### 2. Maximum Core Average Grows with n

The worst-case core average appears to be growing toward some limit > 1/3:

| n | Max core avg | Champion tree deg sequence |
|---|--------------|----------------------------|
| 7 | 0.333333 | [3, 2, 2, 2, 1, 1, 1] |
| 11 | 0.380282 | [4, 3, 2, 2, 2, 2, 1, 1] |
| 13 | 0.400990 | [4, 4, 2, 2, 2, 2, 2, 1] |
| 15 | 0.416096 | [5, 4, 2, 2, 2, 2, 2, 2] |
| 17 | 0.431361 | [5, 5, 2, 2, 2, 2, 2, 2] |

The champion trees appear to be "dumbbells" with two high-degree core vertices connected by a chain of degree-2 core vertices, with leaves attached to the ends.

### 3. Heavy-Light Matching Also Fails

The alternative approach (match each heavy core vertex to a distinct light core neighbor) also fails:

- First matching violation at n=8
- By n=18: 3139/8323 trees (37.7%) have at least one heavy core vertex without a light core neighbor

This suggests that in many d_leaf ≤ 1 trees, **all core vertices are heavy** (all have P > 1/3), or the heavy vertices form an independent set in the core subgraph.

### 4. Pattern of Violations

Looking at sample violations:
- Many have core vertices all at degree 2
- Core vertices appear to form a path or chain structure
- High-degree vertices tend to be at the periphery (connecting core to leaves/support)

Example (n=8, max_core_avg=0.338983):
- deg = [3, 2, 2, 2, 2, 1, 1, 1]
- Core size = 2
- Both core vertices have P = 0.338983 (both heavy)
- They are neighbors (no light core vertices to match to)

## Implications

The edge bound approach P(u) + P(v) ≤ 2/3 successfully bounds the leaf+support contribution to Σ_{L∪S} ≤ |L∪S|/3, but this is **not sufficient** to force the core average below 1/3.

The core vertices can collectively exceed average 1/3 by a significant margin (~0.43 observed at n=17).

## Worst-Case Tree Structure

Detailed analysis reveals a **strikingly simple pattern**:

### Pattern 1: Odd n (worst case has |C| = 1)

The worst-case tree has exactly ONE core vertex:
- |L| = |S| = (n-1)/2 (half are leaves, half support)
- |C| = 1 (single core vertex of degree 2)
- The core vertex connects exactly 2 support vertices
- All core vertices have degree 2 and connect ONLY to support vertices (0 core neighbors)

Example (n=17, core_avg = 0.431361):
```
deg: [5, 5, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1]
     ↑  ↑   └─────────┘
   2 high-deg    single
   support      core v
     trees
```

Structure: Two high-degree support vertices (each with multiple leaves) connected by a single degree-2 core vertex.

### Pattern 2: Even n (worst case has |C| = 2)

The worst-case tree has exactly TWO core vertices:
- |L| = |S| = n/2 (perfectly balanced)
- |C| = 2 (two degree-2 core vertices)
- The two core vertices are **NOT adjacent** to each other (deg dist: {0: 2} means 0 core neighbors each)
- Each core vertex connects exactly 2 support vertices

Example (n=18, core_avg = 0.410608):
```
deg: [5, 4, 3, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1]
     ↑  ↑  ↑        └──┘
   3 high-deg      2 core vertices
   support         (non-adjacent)
     trees
```

### Key Observations

1. **All worst-case core vertices**:
   - Have degree exactly 2
   - Connect ONLY to support vertices (never to other core vertices)
   - Have identical occupation probabilities (by symmetry)

2. **The maximizing structure** appears to be:
   - Take k high-degree support "hubs" (each with many leaves)
   - Connect them via degree-2 core vertices
   - The core vertices act as "bridges" between leaf-heavy subtrees

3. **Core vertices are isolated** in the core subgraph:
   - For |C| = 1: trivially isolated
   - For |C| = 2: both have 0 core neighbors (connected only via support)
   - General pattern: core vertices form an independent set in the core-induced subgraph

4. **Asymptotic pattern**:
   - Odd n: max core avg grows from 1/3 (n=7) toward ~0.43 (n=17)
   - Even n: max core avg stays somewhat lower (~0.34–0.42)
   - The limit appears to approach some value around 0.43–0.45

## Next Steps

1. **Prove the upper bound analytically**:
   - For a degree-2 core vertex v connecting two support vertices s₁, s₂
   - Each support vertex connects to exactly one leaf
   - What is P(v) as a function of the subtree structures beyond s₁ and s₂?

2. **Compute asymptotic limit**:
   - As n → ∞, what is lim max P(v) for a degree-2 core vertex in a balanced tree?
   - Does it converge to some value < 1/2?

3. **Alternative approach**:
   - Since core vertices are ISOLATED in the core subgraph, can we bound them individually?
   - For degree-2 vertex v with two support neighbors: P(v) ≤ f(deg(s₁), deg(s₂))?
   - Can we prove core vertices must have small degree (≤ 2?) in d_leaf ≤ 1 trees?

4. **Sufficient condition**:
   - If we can prove max P(v) < 1/2 for all core vertices, then even if |C| large,
     we'd have Σ_C P(v) < |C|/2, which combined with Σ_{L∪S} ≤ |L∪S|/3 gives useful bounds

## Verification

All computations done via exhaustive enumeration with `geng` and exact IS polynomial computation.
Script: `verify_core_avg.py`
Run time: ~19 seconds for n ≤ 18 (8370 trees with d_leaf ≤ 1)
