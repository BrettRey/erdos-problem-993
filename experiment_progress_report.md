# Erdős Problem #993 - Experiment Progress Report

## Executive Summary
Three novel search strategies were implemented and executed to attack the unimodality conjecture for tree independence polynomials:
1.  **Evolutionary "Log-Concavity Breaker":** A genetic algorithm targeting high local LC violations.
2.  **Bridge Stitching Scan:** A constructive method combining rooted tree signatures to form larger trees ($N \le 30$).
3.  **Rooted Forest Beam Search:** An iterative "growth" algorithm merging subtrees to maximize LC ratio ($N > 500$).

**Status:** No non-unimodal counterexample found.
**Key Insight:** The search space for $N \le 30$ appears robustly unimodal. High log-concavity violations (Ratio > 1.2) exist but do not cascade into unimodality failures (Ratio > 1.0 for $i_{k-1} > i_k < i_{k+1}$ condition). Large random trees stabilize into log-concavity.

---

## 1. Evolutionary "Log-Concavity Breaker"
**Goal:** Evolve trees that maximize $R = \max_k \frac{i_{k-1} i_{k+1}}{i_k^2}$.
- **Setup (upgraded):**
  - Multi-island GA with migration.
  - Size cap enforced (`--max-nodes`) for both mutation and crossover.
  - Fitness retargeted toward near-valley behavior:
    near-miss ratio + local LC ratio around the peak window.
- **Result:**
    - The initial $N=26$ seed (LC ratio $\approx 1.145$) is no longer a strict attractor under the retargeted objective.
    - Best near-miss ratio reached `0.9271937066` at `n=60`
      (`results/best_evolutionary_tree_multiisland_seed998_g160_cap60_nmfit.json`).
    - Still no non-unimodal witness.

## 2. Bridge Stitching Scan
**Goal:** Construct $T$ by bridging $T_1, T_2$.
$$ I(T) = A_0 B_0 + x(A_0 B_1 + A_1 B_0) $$
- **Setup:** Library of rooted signatures for $n \in [1, 15]$, with seeded randomized sampling on large batches.
- **Runs:**
  - `seed=993`, `random_samples=100k`: `6,646,354` pairs checked.
  - `seed=994`, `random_samples=300k`: `12,646,354` pairs checked.
  - `seed=995`, `random_samples=500k`: `18,646,354` pairs checked.
- **Result:**
    - Sampled up to ~18.6 million pairs up to $N=30$.
    - **Max LC Ratio:** $1.3125$.
    - **Unimodality:** All constructed polynomials were unimodal.
- **Significance:** This method is extremely fast ($10^6$ checks/sec) compared to full tree generation. It successfully creates "dangerous" polynomials but hasn't broken unimodality.

## 3. Rooted Forest Beam Search
**Goal:** Grow large trees by iteratively merging $k$ subtrees to maximize LC ratio.
- **Setup:** Beam width 1000, 20 generations. Biased towards path-like growth ($k=1,2$).
- **Result:**
    - Reached $N > 500$.
    - **Max LC Ratio:** $\approx 0.99$.
    - **Observation:** The ratio climbs towards 1.0 but appears to asymptote *below* it. Large "evolved" trees tend to become strictly log-concave.
    - **Conclusion:** The "danger zone" (Ratio > 1) is likely confined to specific sparse structures at smaller $N$, which get "diluted" in large aggregations.

## Strategic Implications
- **Hardness:** The "gap" between LC failure (easy to find at small $N$) and Unimodality failure (impossible so far) is real and significant.
- **Large N Stability:** Random or guided growth to large $N$ ($>100$) seems to stabilize the polynomials into log-concavity. This suggests a counterexample, if one exists, is likely a "finely tuned" structure at moderate $N$ ($30 < N < 100$) rather than a generic property of large trees.
- **Next Steps:**
    - The "Frankenstein" approach (Step 3 in the original plan) is the most promising extension: use the bridge stitching formula but *relax* the requirement that $A, B$ come from real trees, find a numerical counterexample, and *then* try to realize it.
    - We have confirmed that generic "stitching" of valid trees doesn't easily break unimodality.
