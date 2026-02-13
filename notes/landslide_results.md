# Landslide Strategy Results

## Hypothesis
The "Landslide" hypothesis proposed that the coefficient sequence of a tree's independent set polynomial could be destabilized by:
1.  Identifying a "weak layer" (a plateau or shallow descent in the coefficient sequence).
2.  "Overloading" this layer by grafting a specific "debris" tree that adds mass to the sensitive index, triggering a second peak (unimodality violation).

## Experiments

### 1. Landslide Search (`landslide_search.py`)
- **Base:** Broom trees (known for steep slopes) and Random trees.
- **Debris:** Small trees ($n \le 12$) grafted to various nodes.
- **Result:** 
    - Brooms are very stable.
    - Random trees frequently exhibit plateaus (ratio $i_{k+1}/i_k \approx 1.0$) in the tail.
    - No violations found.

### 2. Landslide Evolution (`landslide_evolution.py`)
- **Method:** Evolutionary algorithm to evolve the "debris" tree structure and the target grafting node on a fixed random base.
- **Parameters:** Base $N=60$, Debris $N=10$, 200 generations.
- **Result:**
    - The search consistently converges to a ratio of exactly `1.00000` (within float precision).
    - This corresponds to finding a "shelf" or plateau in the coefficient sequence.
    - Despite intense optimization pressure, the algorithm could not push the ratio above 1.0.

## Conclusion
The independent set polynomials of trees appear to be **locally robust** against this type of perturbation. 
- "Plateaus" ($i_k = i_{k+1}$) are easy to construct and are stable.
- "Valleys" ($i_k > i_{k+1} < i_{k+2}$) are extremely hard to create by simple grafting.
- The convolution-like operation of grafting tends to smooth out local irregularities rather than amplify them.

This suggests that if a counterexample exists, it requires a more global structural property (like the "Legs" or "Islands" found in other approaches) rather than a local "landslide" trigger on a standard tree.
