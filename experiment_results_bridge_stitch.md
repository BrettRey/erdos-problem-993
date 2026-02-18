# Bridge Stitching Search Results

## Objective
To search for non-unimodal trees by constructing synthetic trees via the "Bridge Stitching" method.
This involves joining two rooted trees $T_1$ and $T_2$ by an edge between their roots $r_1$ and $r_2$.
The independence polynomial of the resulting tree $T$ is calculated using the formula:
$$ I(T) = I(T_1-r_1) \cdot I(T_2) + x \cdot I(T_1-N[r_1]) \cdot I(T_2-r_2) $$

## Methodology
- **Signature Generation:** We generated "rooted signatures" $(I(T-r), I(T-N[r]))$ for all rooted trees up to $N=15$.
    - *Correction:* A bug in `indpoly.py` previously handled forest polynomials incorrectly (returning single-component results for disconnected graphs), leading to invalid small signatures and spurious high LC ratios. This has been fixed.
- **Search Strategy:**
    - Exhaustive pairing for small component sizes.
    - Randomized sampling (100k pairs per size combination) for larger sizes where the product space exceeded 500k pairs.
- **Scope:** Covered component sizes $n_1, n_2 \in [1, 15]$, effectively probing the space of trees up to $N=30$.

## Results
- **Non-Unimodal Trees:** 0 found.
- **Max Log-Concavity Ratio:** $\approx 0.855$
    - Found at $N=29$ (from components $n_1=14, n_2=15$).
    - This ratio is $< 1.0$, meaning the polynomials are strictly log-concave.
    - The previously reported ratio of $1.22$ was due to the forest polynomial bug.
- **Validation:** The search now produces valid tree polynomials. The absence of violations up to $N=30$ reinforces the robustness of the conjecture in this size range.

## Conclusion
The bridge stitching method is a computationally efficient way to explore a vast space of tree structures. With the bug fixed, the results align with the expectation that trees up to $N=30$ are unimodal and log-concave. The "danger zone" of LC violations seems harder to reach than initially thought with random stitching.
