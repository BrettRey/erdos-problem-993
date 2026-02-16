# ECMS Positive Shift Analysis

## Objective
Investigate the rare cases (0.05%) where edge contraction increases the mode: $	ext{mode}(I(T/e)) = 	ext{mode}(I(T)) + 1$.

## Data
- Analyzed 108,374 edges (from $n=15$ trees).
- Found 55 positive shifts (+1).
- Found 0 shifts of +2 or more.

## Structural Pattern
In all 55 cases:
1.  **Original Polynomial $I(T)$:** Has a very "flat" peak.
    *   Example: $i_5 = 908, i_6 = 907$. Mode is 5.
2.  **Contracted Polynomial $I(T/e)$:** Has a shifted peak.
    *   Example: $i_5 = 798, i_6 = 840$. Mode is 6.
3.  **Mechanism:**
    *   $I(T/e) = I(T) - \Delta$.
    *   The "removed mass" $\Delta$ is concentrated at the original mode $k$ or to its left.
    *   $\Delta_5 = 110, \Delta_6 = 67$.
    *   Removing more from 5 than 6 flips the balance from $908 > 907$ to $798 < 840$.

## Implication for Proof
- A shift of $+1$ requires a "tie-breaking" perturbation on a nearly-flat peak.
- A shift of $+2$ would require "eroding" *two* consecutive peaks.
- Since $\Delta$ is related to $I(T)$ (it's a subgraph polynomial), its shape is constrained. It cannot be arbitrarily "spiky" to suppress specific coefficients without affecting neighbors.
- This supports the conjecture that shifts $\ge 2$ are impossible.

## Conclusion
The ECMS conjecture is likely true because "positive shifts" are merely "tie-breakers" in the left direction (removing mass from the left).
