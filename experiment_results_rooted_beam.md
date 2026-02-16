# Rooted Forest Beam Search Results

## Objective
To find non-unimodal trees by iteratively constructing "rooted forest signatures" and combining them.
This implements a Beam Search over the space of all constructible trees (via root merging).
Formula: $P_{new} = \prod A_i + x \prod B_i$.

## Methodology
- **State:** Tuple $(A, B)$ where $A=I(T)$, $B=I(T-root)$.
- **Transition:** Merge $k$ subtrees (biased towards $k=1, 2$ to favor path-like structures).
- **Beam:** 1000 best signatures sorted by Log-Concavity Ratio of $A$.
- **Generations:** 20 (reached tree sizes $N > 500$).

## Results
- **Non-Unimodal Trees:** 0 found.
- **Max LC Ratio:** $\approx 0.99$ (Strictly Log-Concave).
- **Trend:** The max LC ratio starts low (0.6) and climbs *towards* 1.0 as $N$ increases, but appears to asymptote *below* 1.0 for these large "fat" trees.
- **Comparison:** The Bridge Stitching scan found ratios $\approx 1.22$ at small $N$ ($N=16$). The Beam Search, by rapidly growing $N$, might be moving *away* from the "boundary" where violations occur (likely small, sparse trees).

## Conclusion
The "Rooted Forest" construction is valid and covers the tree space. However, growing trees "randomly" or even guiding by LC ratio tends to produce stable, log-concave polynomials. The violations ($LC > 1$) seem to be specific "corner cases" (like specific brooms) that get "averaged out" when merging many subtrees.
Future work should focus on **small N** exhaustive or targeted pruning rather than massive N growth.
