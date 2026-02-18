# Sub-Claim A Harness (2026-02-18)

Target:

- For mixed spiders `S(2^k,1^j)`, prove
  `margin(k,j+2) >= margin(k,j)` for all `k>=6`, `j>=0`.

This file documents a canonical computational harness so multiple proof attempts
can share one definition and one certificate format.

## Script

- `verify_subclaim_A_parity_tail.py`

It computes, for each `(k,j)`:

- `delta_margin = margin(k,j+2) - margin(k,j)`,
- decomposition `delta_margin = delta_A + delta_E` with
  - `A = k*2λ/(1+2λ) + j*λ/(1+λ) - (m-1)`,
  - `E = r(B-A)/(1+r)`,
  - `r = λ(1+λ)^(k-j)/(1+2λ)^k`,
  - `B = 1 + k*λ/(1+λ) - (m-1)`,
- mode-shift diagnostics `m(k,j+2)-m(k,j)`.

## Main run

Command:

`python3 verify_subclaim_A_parity_tail.py --k-min 6 --k-max 8000 --j-max 80 --out results/whnc_subclaim_A_parity_tail_k6_8000_j80.json`

Outcome:

- total pairs checked: `647,595`,
- `full_fail = 0` (no failures for full claim `j>=0`),
- `even_tail_fail = 0` (`j>=4`, even),
- `odd_tail_fail = 0` (`j>=5`, odd),
- `mode_shift_hist = {1: 647,595}`.

Witnesses:

- minimum full delta:
  - `(k,j)=(8000,80)`,
  - `delta_margin = 4.5781863263982814e-05`.
- minimum odd-tail delta:
  - `(k,j)=(7999,79)`,
  - `delta_margin = 4.5806584239471704e-05`.

Observed decomposition extremes:

- minimum `delta_A` at `(k,j)=(6,2)`:
  - `delta_A = -0.004866283177902275`,
  - but `delta_E = +0.00970750149052317`,
  - so `delta_margin > 0`.
- minimum `delta_E` at `(k,j)=(6,0)`:
  - `delta_E = -0.022396362280884596`,
  - but `delta_A = +0.02888991621887449`,
  - so `delta_margin > 0`.

## Notes for parallel proof lanes

The harness suggests all mode shifts are exactly `+1` on `j -> j+2` over this
range. That is a useful structural sub-lemma candidate for analytic work:

1. prove `m(k,j+2)=m(k,j)+1` for `k>=6`, `j>=0`,
2. then prove `delta_A + delta_E >= 0` under that mode relation.
