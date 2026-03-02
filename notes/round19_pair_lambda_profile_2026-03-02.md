# Round 19 Pair-Lambda Profile (2026-03-02)

## Scope

Exhaustive scan through `n <= 18` over support-root, step-prefix instances (`k < mode(I_new)`, smallest mode index), with boundary-correct indexing.
Only `X_k < 0` cases were profiled.

Condition tested:

`sum_err <= D + lambda * (C10 + C01 + C11)`

where

- `D = sum_i Lambda_old(i) P(k-i)Q(k-i)`
- `C10 = sum_i Lambda_old(i) P(k-i-1)Q(k-i)`
- `C01 = sum_i Lambda_old(i) P(k-i)Q(k-i-1)`
- `C11 = sum_i Lambda_old(i) P(k-i-1)Q(k-i-1)`

## Global results

- Total negative cases: `135,976`
- All negatives are at step `t=2`
- Global minimal scalar through `n<=18`:
  - `lambda* = 0.05201381704686925`

Tight witness:

- `n=18`
- `g6=Q???????????????O?E??NwA_^?`
- `root=0, step=2, k=5, mode=7`
- `X=-135,037`
- `D=1,442,539`
- `sum_err=1,546,363`
- `need = sum_err - D = 103,824`
- `C10=789,600, C01=772,310, C11=434,175`
- `R_shift=1,996,085`
- `need / R_shift = 0.05201381704686925`

## Pair-class decomposition (step-2 child-size pair `(a,b)`)

Only 9 pair classes appear in all step-2 negative cases.

Classes requiring nonzero lambda:

1. `(2,14)`
   - count `68,908`
   - `lambda*_pair = 0.05201381704686925` (global tight class)
2. `(3,13)`
   - count `25,412`
   - `lambda*_pair = 0.04386927442810327`
3. `(2,13)`
   - count `28,443`
   - `lambda*_pair = 0.023760967407659456`

Classes with `lambda*_pair = 0` (already satisfy `sum_err <= D` without shifted reserve):

- `(2,12)` count `8,098`
- `(3,12)` count `2,704`
- `(2,11)` count `1,936`
- `(2,10)` count `459`
- `(3,11)` count `11`
- `(2,9)` count `5`

## Implication

This sharply reduces the proof burden:

- The odd-reserve mechanism is not uniformly hard across all step-2 pairs.
- The only genuinely hard classes through `n<=18` are `(2,14)`, `(3,13)`, `(2,13)`.
- A realistic closure strategy is:
  1. prove base dominance `sum_err <= D` for all other pair classes,
  2. derive a small shifted reserve for the three hard classes,
  3. optimize constants by pair class (instead of a single global lambda).

## Artifacts

- Script: `analyze_round19_pair_lambda.py`
- JSON: `results/round19_pair_lambda_profile_n18.json`
