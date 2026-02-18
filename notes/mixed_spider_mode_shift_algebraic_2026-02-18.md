# Mixed-Spider Mode Step: algebraic reduction and certificate (2026-02-18)

Target statement:

- For mixed spiders `S(2^k,1^j)`, prove
  `mode(k,j+2) = mode(k,j) + 1`.

Here `mode(k,j)` is the leftmost mode of

`I_{k,j}(x) = (1+2x)^k (1+x)^j + x(1+x)^k`.

## Algebraic reduction

If we have the closed-form mode law

`mode(k,j) = floor((4k + 3j + 3)/6)`

(on the domain of interest, e.g. `k>=6`, `j>=0`), then the parity-step claim is immediate:

`mode(k,j+2)`
`= floor((4k + 3(j+2) + 3)/6)`
`= floor((4k + 3j + 3)/6 + 1)`
`= floor((4k + 3j + 3)/6) + 1`
`= mode(k,j) + 1`.

So the step claim is an arithmetic corollary of the mode law.

## What was checked exactly

Script:

- `verify_mixed_spider_mode_shift_formula.py`

Runs:

`python3 verify_mixed_spider_mode_shift_formula.py --k-min 1 --k-max 1200 --j-max 120 --out results/whnc_mixed_spider_mode_formula_k1_1200_j120.json`

`python3 verify_mixed_spider_mode_shift_formula.py --k-min 6 --k-max 3000 --j-max 120 --out results/whnc_mixed_spider_mode_formula_k6_3000_j120.json`

Checks in both runs:

1. `mode(k,j) == floor((4k+3j+3)/6)`.
2. `mode(k,j+2) == mode(k,j)+1`.

Outcomes:

- On `k=1..1200`, `j=0..120`:
  - `mode_formula_fail_count = 1`
  - `mode_shift_fail_count = 1`
  - unique witness is the known tiny exception `(k,j)=(2,0)`.
- On target regime `k=6..3000`, `j=0..120`:
  - `mode_formula_fail_count = 0`
  - `mode_shift_fail_count = 0`.

## Status

- The parity-step statement is algebraically proved **conditional on** the closed-form mode law.
- The computational certificate supports both the mode law and the step claim on a large exact scan domain.
