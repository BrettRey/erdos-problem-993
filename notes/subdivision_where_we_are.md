# Subdivision Lemma: Where We Are

## Summary

The subdivision lemma (edge subdivision preserves unimodality) has been empirically verified through extensive testing:
- 5.7M edge subdivisions (n≤19 exhaustive)
- 19.7M finite-core checks with no failures

The proof has been attempted through several paths:

### Path 1: Tail Ratio Monotonicity (Most Promising)

**Goal**: Show that for forest polynomials, coefficient ratios are nonincreasing in the tail.

**Challenge**: The simple lemma "product of decreasing sequences is decreasing" is FALSE.

**What works**: The specific structure of forest independence polynomials (coming from actual trees, not arbitrary sequences) seems to have special properties.

**Key observation**: In the tail region, the ratios approach 1 from below, suggesting eventual decrease.

### Path 2: First-Difference Dominance

**Goal**: Show ΔA_k <= -ΔI_k for k >= d+1

**Status**: Empirically holds. Need theoretical backing.

**Lemma proved**: A(x) <= (1+x) I(T) coefficientwise - this is a key inequality.

### Path 3: Finite Kernel Fallback

**Idea**: If we can characterize all "leaf-light" forests (which are the only potential counterexamples), we can check them directly.

**Status**: This is the most practical path. The finite-core checks already verify millions of cases.

## The Core Challenge

The subdivision lemma would be proven if we can show:

For forest polynomials F and G with tail-decreasing coefficients, their product H = F * G also has tail-decreasing coefficients.

This is true empirically but the proof requires understanding the special structure of independence polynomials.

## Recommendations

1. **Push on finite kernel**: Extend the kernel checks further - if we can get to where all potential counterexamples are covered, that's a practical proof.

2. **Look for counterexamples**: If we can find even one case where the conditions fail, we understand the proof requirements better.

3. **Alternative approach**: Maybe the lemma isn't universally true, but needs additional conditions (e.g., about the structure of the forest).

## Current Best Path

The finite-kernel approach seems most tractable:
- Show that "leaf-heavy" hubs are safe (already partially done)
- Characterize the finite set of "leaf-light" cases
- Verify directly

This would give a complete proof, even if not fully general.
