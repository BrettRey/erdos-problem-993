# Round 21 Odd Compensation Profile (2026-03-02)

## Scope

Exhaustive through `n<=18` on support-root step-prefix `X<0` cases, boundary-correct indexing.

Fixed constant:

- `lambda0 = 0.05201381704686925`

Channels:

- `D = sum_i Lambda_old(i) P(k-i)Q(k-i)`
- `R_shift = C10 + C01 + C11`

Error split:

- `sum_err = even_err + odd_err`

## Global checks

- Total negative cases: `135,976`
- `xneg_step_not2 = 0` (all negatives at step 2)
- `even_err <= D`: 0 failures
- `sum_err <= D + lambda0*R_shift`: 0 failures

Odd-only bound status:

- `odd_err <= lambda0*R_shift` fails in `38,219` cases
- in all those cases, compensation from even slack closes the gap:
  - `odd_surplus := (odd_err - lambda0*R_shift)_+`
  - `even_slack := (D - even_err)_+`
  - observed `odd_surplus <= even_slack` in every case (equivalent to global closure + even bound)

## Extremals

Global tight witness for required lambda:

- same witness as previous Round 19 fit
- `needed_global_lambda = 0.05201381704686925`

Largest odd-only lambda demand:

- `needed_odd_lambda = odd_err / R_shift = 0.559442...`
- class `(3,13)` witness
- still has `needed_global_lambda = 0` because even slack absorbs odd surplus.

## Pair-class summary

Max needed global lambda by class:

- `(2,14): 0.0520138170` (tight)
- `(3,13): 0.0438692744`
- `(2,13): 0.0237609674`
- all other observed classes: `0`

This reproduces and explains the hard-core classes while showing that odd-only local control is too strong.

## Artifacts

- Script: `scan_round21_odd_compensation.py`
- JSON: `results/round21_odd_compensation_n18.json`
