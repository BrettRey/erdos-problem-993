# Next Steps: Separation Invariant
Date: 2026-07-03

## Current Position

The computational terrain has narrowed. Log-concavity failures are no longer the right object to classify directly. Across the current stress set -- Ramos--Sun, Galvin, Bautista--Ramos, Li/Kadrawi--Levit, exhaustive `n=28` LC failures, and bounded forest products -- the useful repeated pattern is:

1. The worst post-mode tail-rise pressure occurs immediately after the mode or very near it.
2. LC defects are ratio bumps farther to the right, where the absolute ratios are already well below 1.
3. In the `d_leaf <= 1` lane, all tested stress rows remain low-mode and mean-bounded.

So the next proof-facing target should not be "almost log-concavity." It should be a separation statement about the consecutive coefficient ratios

```text
r_k = a_k / a_{k-1}.
```

## Proposed Next Move

Draft and falsify a concrete separation invariant:

> For tree/forest independence polynomials in the tested stress regimes, once `r_k < 1` near the mode, any later LC-breaking bump `r_{j+1} > r_j` occurs while `r_{j+1}` is still bounded away from 1 by a reserve that is already implied by the low-mode/mean-bound structure.

This is intentionally weaker than log-concavity and stronger than unimodality. It says LC bumps can happen, but not in a way that creates a valley recovery.

## Immediate Work Items

1. Add a ratio-profile audit script.
   - Input: a summary JSON or generated family grid.
   - Output: for each row, `max_tail_ratio`, `argmax_tail_ratio`, all LC bump locations, and `max_lc_bump_right_ratio`.
   - Purpose: make the separation invariant explicit and machine-checkable.

2. Stress larger product powers.
   - Current products stop at 3 factors and already reach near-miss `0.99238`.
   - Run powers of the top product components up to a practical bound, especially `TG_{8,8}^p`, `Galvin(21,11)^p`, and mixed products.
   - Check whether near-miss tends to 1 while LC bumps remain separated.

3. Split the invariant by lane.
   - `d_leaf <= 1`: expect low-mode and `mu < n/3` to survive.
   - outside `d_leaf <= 1`: do not require low-mode/mean-bound; instead check whether PNP transfer is the right route.

4. Only after those scans, write a proof candidate.
   - Candidate statement A: a `d_leaf <= 1` separation theorem using low-mode/mean-bound.
   - Candidate statement B: a product-stability lemma showing separated ratio bumps do not create non-unimodality under convolution.

## Stopping Rule

This is worth continuing if product powers still show:

- no non-unimodality;
- the maximum tail ratio approaches 1 only near the mode;
- LC bump right-ratios stay far below the maximum tail ratio;
- `d_leaf <= 1` products remain low-mode and mean-bounded.

If larger products produce a late LC bump with right-ratio close to 1, the separation invariant is false and we should pivot back to direct product/forest counterexample search.
