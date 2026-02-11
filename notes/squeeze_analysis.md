# Squeeze Scan Analysis (n <= 23)

**Data:** `results/squeeze_n12_networkx.json`, `results/squeeze_n23_geng.json`.

## Key findings

1) The strict squeeze condition

   `first_descent >= tail_start = ceil((2α-1)/3)`

   fails starting at n = 9.

2) Worst margin through n <= 23 is **-3**, i.e.

   `first_descent = tail_start - 3`.

3) The margin -3 first appears at n=21 and persists through n=23.

   - n=21: 15 trees total (7 double-stars, 8 other).
   - n=22: 267 trees total.
   - n=23: 3,513 trees total.

So any deterministic squeeze proof must allow a margin of at least -2.
The data now force margin at least -3.

## Implication

The “squeeze from both ends” approach needs a relaxed target, e.g.

  `first_descent >= tail_start - 3`.

This bound is tight on the current scan (n <= 23).

The remaining gap is then just two indices, which is plausibly addressable
by a sharp local inequality in a minimal-counterexample argument.

## Boundary check at t-2 and t-1 (empirical)

Let t = ceil((2α-1)/3).

1) Margin = -2 (n <= 20): when d(I)=t-2, the boundary step always satisfies
   i_{t-1} >= i_t. No counterexample found.

2) Margin = -3 (n <= 23): when d(I)=t-3, the boundary steps satisfy
   i_{t-2} >= i_{t-1} >= i_t in all cases checked.
   - n=21: 15 cases, all pass.
   - n=22: 267 cases, all pass.
   - n=23: 3,513 cases, all pass.

This suggests a two-part squeeze strategy:
1) prove the universal bound d(I) >= t-3, and
2) prove the boundary inequalities
     i_{t-2} >= i_{t-1} >= i_t
   when d(I)=t-3.

If both hold, unimodality follows from the Levit–Mandrescu decreasing tail.
