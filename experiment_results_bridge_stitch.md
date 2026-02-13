# Bridge Stitching Search Results

## Objective
To search for non-unimodal trees by constructing synthetic trees via the "Bridge Stitching" method.
This involves joining two rooted trees $T_1$ and $T_2$ by an edge between their roots $r_1$ and $r_2$.
The independence polynomial of the resulting tree $T$ is calculated using the formula:
$$ I(T) = I(T_1-r_1) \cdot I(T_2) + x \cdot I(T_1-N[r_1]) \cdot I(T_2-r_2) $$

## Methodology
- **Signature Generation:** We generated "rooted signatures" $(I(T-r), I(T-N[r]))$ for all rooted trees up to $N=15$.
- **Search Strategy:**
    - Exhaustive pairing for small component sizes.
    - Randomized sampling (`--random-samples` pairs per size combination) for larger sizes where the product space exceeded 500k pairs.
- **Scope:** Covered component sizes $n_1, n_2 \in [1, 15]$, effectively probing the space of trees up to $N=30$.

## Results
- **Non-Unimodal Trees:** 0 found in all runs.
- **Run A (seed=993, random-samples=100k):**
    - `checked = 6,646,354`
    - `max_lc_ratio = 1.2222222222`
- **Run B (seed=994, random-samples=300k):**
    - `checked = 12,646,354`
    - `max_lc_ratio = 1.2222222222`
- **Run C (seed=995, random-samples=500k):**
    - `checked = 18,646,354`
    - `max_lc_ratio = 1.3125`
- **Max LC witness polynomial:**
    - `[1, 4, 21, 75, 142, 146, 81, 22, 2]` (from Run C).
    - Earlier witness `[1, 3, 11, 36, 59, 46, 18, 3]` remains reproducible in Runs A/B.
- **Validation:** The search successfully found trees with LC ratios $> 1.0$, confirming that the method explores the "dangerous" territory of log-concavity violations, even if it hasn't broken unimodality yet.

## Conclusion
The bridge stitching method is a computationally efficient way to explore a vast space of tree structures ($N \approx 30$). While no counterexample was found in this sampled run, the presence of high LC ratios suggests the method is viable for deeper or more targeted searches (e.g., using evolutionary optimization on the component signatures instead of random sampling).

## Implementation Note
An early-return indentation bug in `scripts/bridge_stitch_scan.py` was fixed.
After the fix, full-size scans run to completion and produce reproducible
artifacts with explicit `seed` and `random_samples` metadata.
