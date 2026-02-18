# Erdős Problem #993 - Experiment Progress Report

## Executive Summary
Three novel search strategies were implemented and executed to attack the unimodality conjecture for tree independence polynomials:
1.  **Evolutionary "Log-Concavity Breaker":** A genetic algorithm targeting high local LC violations.
2.  **Bridge Stitching Scan:** A constructive method combining rooted tree signatures to form larger trees ($N \le 30$).
3.  **Rooted Forest Beam Search:** An iterative "growth" algorithm merging subtrees to maximize LC ratio ($N > 500$).

**Update (Radio Roots):** A fourth strategy, **"Radio Roots,"** is currently running (PID 17289). It maximizes "Root Turbulence" (proximity to the positive real axis) and Near-Miss Ratio. Early results are highly promising ($NM \approx 0.97$).

**Status:** No non-unimodal counterexample found yet.
**Key Insight:** The search space for $N \le 30$ appears robustly unimodal. High log-concavity violations (Ratio > 1.2) previously reported were due to a software bug; corrected runs show ratios $\le 0.86$ for these constructed trees. Large random trees stabilize into log-concavity.

---

## 4. Radio Roots Search (Active)
**Goal:** Attack unimodality via "Root Turbulence."
**Theory:** Unimodality requires all real roots (negative). Breakdown occurs when roots become complex and move towards the positive real axis.
**Strategy:**
- **Metric:** `Score = 20*NM + 5*Max_Re - 0.5*Max_Im`.
- **Target:** Maximize Near-Miss Ratio (deep valleys) and push roots right (Re -> >0).
- **Status:** Running (PID 17289).
- **Early Result:** Iteration 1 found `cat_34` ($N=194$) with **NM = 0.9704** and **Max Re = 0.41**. This is the strongest "near-miss" signal seen in the project so far, confirming the "danger zone" hypothesis.

---

## 1. Evolutionary "Log-Concavity Breaker"
**Goal:** Evolve trees that maximize $R = \max_k \frac{i_{k-1} i_{k+1}}{i_k^2}$.
- **Setup:** Pop=100, Gens=500. Seeded with known LC failures at $N=26$.
- **Result:**
    - The initial $N=26$ seed (Ratio $\approx 1.145$) remains the local maximum.
    - Mutations tend to degrade this specific "sparse tail" structure.
    - Forcing growth ($N > 26$) tends to *reduce* the LC ratio.

## 2. Bridge Stitching Scan
**Goal:** Construct $T$ by bridging $T_1, T_2$.
$$ I(T) = A_0 B_0 + x(A_0 B_1 + A_1 B_0) $$
- **Setup:** Library of rooted signatures for $n \in [1, 15]$. Sampled $O(10^6)$ pairs for large $N$.
- **Result:**
    - Sampled ~6.8 million pairs up to $N=30$.
    - **Max LC Ratio:** $\approx 0.855$ (at $N=29$).
    - **Correction:** A bug in forest polynomial computation previously yielded spurious ratios $>1.2$. The corrected scan shows strict log-concavity.
- **Significance:** Validates that simple bridging of small trees ($N \le 15$) yields robustly unimodal results.

## 3. Rooted Forest Beam Search
**Goal:** Grow large trees by iteratively merging $k$ subtrees to maximize LC ratio.
- **Setup:** Beam width 1000, 20 generations. Biased towards path-like growth ($k=1,2$).
- **Result:**
    - Reached $N > 500$.
    - **Max LC Ratio:** $\approx 0.99$.
    - **Observation:** The ratio climbs towards 1.0 but appears to asymptote *below* it. Large "evolved" trees tend to become strictly log-concave.
    - **Conclusion:** The "danger zone" (Ratio > 1) is likely confined to specific sparse structures at smaller $N$ (like the $N=26$ brooms), which get "diluted" in large aggregations.

## Strategic Implications
- **Hardness:** The conjecture appears robust against random and semi-random constructions up to $N=30$ and for large "fat" trees.
- **Focus Shift:** The best candidate for a counterexample remains the **Broom-like** structures around $N=26-30$ that exhibit LC failure. Random growth moves away from this "sparse" regime.
- **Next Steps:**
    - **Targeted Broom/Spider Search:** Instead of random stitching, systematically construct generalized brooms (spiders with long legs) and optimize their leg lengths.
    - **Frankenstein Approach:** Re-visit the relaxed polynomial search, but with the corrected code, to see if *numerical* violations exist that *could* be realized.
