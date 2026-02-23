# STRONG C2 Lane (Option 3): Determinant Split Status

Date: 2026-02-19

## Update (2026-02-19): Obstruction profile through n<=23

A dedicated obstruction scanner now isolates where rise-compensation could fail:

- script: `conjecture_a_strong_c2_obstruction_scan.py`
- merged artifact:
  `results/whnc_strong_c2_obstruction_staged_summary_n23.json`
- detailed note:
  `notes/strong_c2_obstruction_profile_2026-02-19.md`

Frontier summary (`d_leaf<=1`, `n<=23`, canonical leaf):

- checked: `931,596`
- mismatch-negative: `129`
- `combined<0`: `0`
- rise-bound failures: `0`
- `q_{m-1}<q_{m-2}` occurs only `2` times total (both mismatch-negative, both
  `mode(B)=m-1`), with large ratio slack in the normalized condition.

This narrows the strong-C2 obstruction to a very small `Q`-drop regime.

## Target

For a leaf `l` with support `s`, `A=T-l`, `B=T-{l,s}`, and `m=mode(I_T)`:

`STRONG C2` is equivalent to

`lambda_m^T >= lambda_{m-1}^B`

and (via mediant) to the cross-tree inequality

`lambda_m^A >= lambda_{m-1}^B`,

i.e.

`a_{m-1} b_{m-1} - a_m b_{m-2} >= 0`.

With `A=(1+x)B + xP` (`P = I(B-u)` for the non-leaf neighbor `u` of `s`):

`a_k = b_k + b_{k-1} + p_{k-1}`.

Substituting gives the exact split:

`a_{m-1} b_{m-1} - a_m b_{m-2}`
`= (b_{m-1}^2 - b_m b_{m-2}) + (p_{m-2} b_{m-1} - p_{m-1} b_{m-2})`
`= lc_surplus + mismatch`.

So option 3 reduces to proving:

`lc_surplus + mismatch >= 0`.

## New scanner and artifacts

- Script: `conjecture_a_strong_c2_split_scan.py`
- Frontier artifacts:
  - `results/whnc_strong_c2_split_n22.json`
  - `results/whnc_strong_c2_split_n23.json`
  - `results/whnc_strong_c2_split_staged_summary_n23.json`

Canonical leaf policy in the scan:
- choose a leaf with minimum support degree
- tie-break by smallest leaf index

## Empirical status (d_leaf<=1, n<=23)

From `results/whnc_strong_c2_split_staged_summary_n23.json`:

- checked trees: `931,596`
- negative mismatch cases: `129`
- negative combined cases (`lc_surplus + mismatch < 0`): `0`
- worst compensation ratio:
  - `max((-mismatch)/lc_surplus) = 0.0400783`
- minimum combined margin:
  - `min(lc_surplus + mismatch) = 9`

Interpretation:
- mismatch can be negative, so this is genuinely a balance inequality;
- but the negative mismatch is tiny relative to LC surplus in all checked cases;
- the determinant margin stays strictly positive on the full frontier.

## Secondary quantitative facts

In the same staged summary:

- `bound_rise_fail = 0` for the candidate bound
  `(-mismatch) <= b_{m-1}(b_{m-1}-b_{m-2})`.
- `bound_drop_fail = 1` for
  `(-mismatch) <= b_{m-2}(b_{m-1}-b_m)` (so this one is not robust).

This suggests the `rise` term is a plausible analytic control term for `-mismatch`.

## Current state of option 3

What is now done:
1. exact algebraic reduction of STRONG C2 to `lc_surplus + mismatch >= 0`,
2. full computational certificate through `n<=23` on the canonical leaf choice,
3. quantified slack showing mismatch is at most ~4.01% of LC surplus in worst case.

What remains:
- a structural inequality proving the compensation analytically (not by scan),
  e.g. an upper bound on `-mismatch` in terms of a provably nonnegative local term
  controlled by `B` around index `m-1`.
