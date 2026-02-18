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
