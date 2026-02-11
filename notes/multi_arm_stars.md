# Multi-Arm Stars: A New Extremal Family for nm

**Finding:** Standard brooms are *not* the extremal family for the near-miss ratio.
"Multi-arm stars" (stars with multiple short paths emanating from the hub)
consistently achieve higher nm than any standard broom at the same vertex count.

## What is a multi-arm star?

A k-arm star `T(a₁, ..., aₖ; s)` has:
- A central hub vertex
- k paths of lengths a₁ ≥ a₂ ≥ ... ≥ aₖ ≥ 2 emanating from the hub
- s pendant leaves attached directly to the hub

Total vertices: n = 1 + (a₁ + ... + aₖ) + s.

A standard broom `broom(p, s)` is the special case k = 1 with a₁ = p - 1.

## Champion configurations

| n range | Best config | Arms |
|---------|------------|------|
| 50 | 2-arm(4, 2) | [4, 2] |
| 75 | 2-arm(6, 2) | [6, 2] |
| 100--150 | 3-arm(6, 3, 2) | [6, 3, 2] |
| ≥ 200 | 4-arm(5, 5, 4, 2) | [5, 5, 4, 2] |

## Comparison with brooms

| n | best broom nm | best multi-arm nm | config | gap |
|---|---|---|---|---|
| 50 | 0.9066 | **0.9117** | 2-arm(4,2) | +0.0051 |
| 75 | 0.9412 | **0.9437** | 2-arm(6,2) | +0.0025 |
| 100 | 0.9551 | **0.9575** | 3-arm(6,3,2) | +0.0025 |
| 200 | 0.9779 | **0.9792** | 4-arm(5,5,4,2) | +0.0013 |
| 500 | 0.9913 | **0.9918** | 4-arm(5,5,4,2) | +0.0006 |
| 1000 | 0.9957 | **0.9959** | 4-arm(5,5,4,2) | +0.0003 |

## Discovery process

1. Built an evolutionary nm optimizer (`nm_optimizer.py`) using leaf relocation,
   SPR, pendant concentration, and edge subdivide/contract mutations.
2. Optimizer consistently converged to broom-like structures with 1 branch vertex
   but with degree sequences that didn't match standard brooms.
3. Structural analysis revealed the optimizer was finding multi-arm stars:
   trees with a high-degree hub connected to both pendant leaves AND multiple
   short paths.
4. Systematic sweep confirmed multi-arm stars beat all standard brooms.
5. Fine-grained search identified 4-arm(5,5,4,2) as the global champion for
   n ≥ 200.

## Why multi-arm stars beat brooms

Intuition: a broom concentrates all its "path structure" in a single arm.
Multi-arm stars distribute the path budget across several shorter arms, which:
- Preserves more pendant leaves (higher s for the same n)
- Creates a more balanced perturbation of the independence polynomial
- The arms interact multiplicatively with the star part, and shorter
  balanced arms produce more effective near-miss amplification

## Asymptotic behavior

Both families satisfy nm → 1 as n → ∞ with nm < 1. The scaling is
`nm ≈ 1 - C/s` where s is the star-leaf count and C depends on the arm
configuration:

| Config | C (s=1000) |
|--------|------------|
| broom(13, s) | ≈ 4.09 |
| 3-arm(6,3,2) | ≈ 4.06 |
| 4-arm(5,5,4,2) | ≈ 3.98 |

4-arm(5,5,4,2) approaches 1 approximately 3% faster than standard brooms.

## Implication for the conjecture

This does not change the fundamental picture:
- nm remains strictly below 1 at all sizes tested
- No counterexample was found
- The conjecture appears safe

But it corrects the identification of the extremal family. Papers and
computational studies claiming brooms are "closest to violation" should
be updated to note that multi-arm stars achieve higher nm.

## Code

- `nm_optimizer.py`: evolutionary optimizer that discovered this
- Results: `results/multi_arm_optimization.json`, `results/nm_optimizer_50_200.json`
