# Erdős Problem #993 - Experiment Progress Report

## Executive Summary
Two novel search strategies were implemented and executed to attack the unimodality conjecture for tree independence polynomials:
1.  **Evolutionary "Log-Concavity Breaker":** A genetic algorithm targeting high local LC violations.
2.  **Bridge Stitching Scan:** A constructive method combining rooted tree signatures to form larger trees ($N \le 30$).

**Status:** No non-unimodal counterexample found.
**Key Insight:** The search space for $N \le 30$ appears robustly unimodal. High log-concavity violations (Ratio > 1.2) exist but do not cascade into unimodality failures (Ratio > 1.0 for $i_{k-1} > i_k < i_{k+1}$ condition).

---

## 1. Evolutionary "Log-Concavity Breaker"
**Goal:** Evolve trees that maximize $R = \max_k \frac{i_{k-1} i_{k+1}}{i_k^2}$.
- **Setup:** Pop=100, Gens=500. Seeded with known LC failures at $N=26$.
- **Result:**
    - The initial $N=26$ seed (Ratio $\approx 1.145$) is a local maximum.
    - Mutations tend to degrade this specific "sparse tail" structure.
    - Forcing growth ($N > 26$) tends to *reduce* the LC ratio, suggesting the violation might be an anomaly of small $N$ discrete combinatorics rather than a scalable property.

## 2. Bridge Stitching Scan
**Goal:** Construct $T$ by bridging $T_1, T_2$.
$$ I(T) = A_0 B_0 + x(A_0 B_1 + A_1 B_0) $$
- **Setup:** Library of rooted signatures for $n \in [1, 15]$.
- **Runs:**
  - `seed=993`, `random_samples=100k` (large batches) -> `6,646,354` pairs checked.
  - `seed=994`, `random_samples=300k` (large batches) -> `12,646,354` pairs checked.
- **Result:**
    - Sampled up to ~12.6 million pairs up to $N=30$.
    - **Max LC Ratio:** $1.222$ (at $N=16$).
    - **Unimodality:** All constructed polynomials were unimodal.
- **Significance:** This method is extremely fast ($10^6$ checks/sec) compared to full tree generation. It successfully creates "dangerous" polynomials but hasn't broken unimodality.

Implementation update:
- Fixed an early-return indentation bug in `scripts/bridge_stitch_scan.py`.
- Added explicit `--seed` and `--random-samples` controls for reproducible sampling.

## Strategic Implications
- **Hardness:** The "gap" between LC failure (easy to find) and Unimodality failure (impossible so far) is real and significant.
- **Next Steps:**
    - The "Frankenstein" approach (Step 3 in the original plan) is the most promising extension: use the bridge stitching formula but *relax* the requirement that $A, B$ come from real trees, find a numerical counterexample, and *then* try to realize it.
    - We have confirmed that generic "stitching" of valid trees doesn't easily break unimodality.
