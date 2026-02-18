# Focused Mode-Tie Progress (2026-02-18)

Goal: prove

`mu(lambda_m) >= m-1`, with `lambda_m = i_{m-1}/i_m`, `m = leftmost mode(I_T)`.

## New packaged scanner

- Script: `conjecture_a_mode_tie_focused_scan.py`
- Outputs:
  - `results/whnc_mode_tie_focused_dleaf_n23.json`
  - `results/whnc_mode_tie_focused_alltrees_n22.json`

## Frontier results

1. `d_leaf<=1`, `n<=23`
- checked: `931,596`
- failures: `0`
- minimum margin: `0.36225903819956606`

2. all trees, `n<=22`
- checked: `9,114,283`
- failures: `0`
- minimum margin: `0.36225903819956606`

Global minimum witness in both runs:
- `n=21`
- `g6=T???????C?G?G?C?@??G??_?@??@???_B~o?`
- degree signature `{1:10, 2:10, 10:1}` (balanced length-2 spider `S(2^10)`)
- `m=7`, `lambda_m=0.879383429672447`, `mu(lambda_m)=6.362259038199566`

## Extremal trend on d_leaf<=1 per n

- Odd `n=21` and `n=23` minima are balanced length-2 spiders (`S(2^10)`, `S(2^11)`).
- Even `n` minima can come from near-spider mixed signatures.

## Strong candidates tested and rejected

1. Mirror pair dominance at tie point
- Candidate: `w_{m-1+j} >= w_{m-1-j}` for all `j`, where `w_k = i_k lambda_m^k`.
- False with many violations (large-gap witnesses).

2. One-sided central-mass bound
- Candidate: `sum_{k<=m-2} (m-1-k) p_k <= p_m` at tie distribution `p`.
- False on almost all trees in frontier.

These failures indicate the proof must use a more global weighted-tail mechanism,
not pointwise pairing.

## Best next proof lane

Given all-tree support through `n<=22`, try proving the focused inequality as a
**general tree theorem** (not just `d_leaf<=1`) via leaf decomposition:

- `I_T = I_{T-l} + x I_{T-{l,s}}` for leaf `l` with support `s`,
- track the focused tie functional `Phi_m(lambda) = lambda I'(lambda) - (m-1)I(lambda)`
  at `lambda=lambda_m(T)`,
- prove positivity propagates under leaf extension.

Fallback lane: prove extremality of balanced length-2 spiders for the focused
margin, then close that family analytically.

## Update: Leaf-Bridge Reduction

See `notes/mode_tie_leaf_bridge_2026-02-18.md` for a stronger staged result:

- full `d_leaf<=1` frontier (`931,596` trees through `n<=23`) has
  zero failures for existence of a leaf with both bridge terms nonnegative;
- choosing a leaf with minimum support degree has zero failures on this frontier;
- a new local target emerges: prove the bridge inequalities for support-degree-2
  leaves (tested on `4,543,370` degree-2-support leaf checks, 0 failures).

## Update: Large-k Spider Stability (staged)

To continue the focused tie lane without overflow at large `k`, we added:

- `conjecture_a_spider_mode_tie_asymptotic.py`

This script evaluates `S(2^k)` at `lambda_m(T)` with log-scaled closed forms
for `mu(T), mu(T-l), mu(T-{l,s})`, and tracks:

- `margin = mu(T,lambda_m)-(m-1)`
- `c1 = mu(T-l,lambda_m)-(m-1)`
- `c2 = mu(T-{l,s},lambda_m)-(m-2)`
- `dmu = mu(T-l,lambda_m)-mu(T-{l,s},lambda_m)`

Staged artifacts:

- `results/whnc_mode_tie_spider_asymptotic_smoke_k40.json`
- `results/whnc_mode_tie_spider_asymptotic_stage1_k41_500.json`
- `results/whnc_mode_tie_spider_asymptotic_stage2_k501_5000.json`
- `results/whnc_mode_tie_spider_asymptotic_stage3_k5001_50000.json`
- aggregate: `results/whnc_mode_tie_spider_asymptotic_staged_summary_k2_50000.json`

Aggregate (`k=2..50,000`, 49,999 cases):

- mode-formula mismatches: `0` (`m=(2k+1)//3`),
- `margin > 0`, `c1 > 0`, `c2 > 0` throughout,
- global minima:
  - `min margin = 0.33333995565772057` at `k=49906`,
  - `min c1 = 0.1666724541428266` at `k=49906`,
  - `min c2 = 0.6666799680824624` at `k=49906`,
- decomposition consistency stayed near machine precision:
  - `max |alpha+beta-1| = 2.22e-16`,
  - `max |(alpha*c1+beta*c2)-margin| = 1.10e-11`.

The smallest margins remain on the `k ≡ 1 (mod 3)` subsequence, consistent with
the observed extremal pattern.
