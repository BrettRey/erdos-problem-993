# Round 19 Tail-Obligation Scan (2026-03-02)

Script: `scan_round19_tail_obligations.py`  
Primary artifact: `results/round19_tail_obligations_n18.json`

## Scope

- Trees exhaustive for `n=3..18` via `geng`.
- Support-step prefix regime (`k < mode(I_new)`), boundary-correct CB indexing (includes `-1`).
- Diagnostics restricted to `X_k < 0` instances.

## Core totals

- Prefix checks: 2,657,523
- `X_k < 0` cases: 135,976
- All `X_k < 0` cases occur at `step=2`.

## Diagonal structure

- Total diagonals examined over all `X_k<0` cases: 1,311,324
- Diagonals with `>1` sign changes in `D_i^(s)=W_i^(s)-W_{i+1}^(s)`: 0
- Max sign changes observed: 1

So the one-sign-change pattern remains exact through `n<=18`.

## Obligation checks

### Even-channel bound

Checked bound for even `s=2t`:

`Err_s <= Lambda_old(t) * P(k-t) * Q(k-t)`

Results:
- checks: 723,650 even diagonals
- failures: 0
- worst ratio: 0.7491

This channel bound appears robust.

### Odd residual

Odd-diagonal error `Err_s` is not negligible:
- odd diagonals with `Err_s > 0`: 38,651 / 587,674 (6.58%)
- max odd `Err_s`: 1,574,573

So the previous “odd = free” simplification is false.

### Global sum-error budget

Checked:

`sum_s Err_s <= D_k`

Results:
- failures: 85
- first failure at `n=17`, graph6 `P?????????????_?W?Dw?G~_`, `root=0`, `step=2`, `k=4`
  with `sum_err=217,793 > D_k=216,402`
- max ratio `sum_err / D_k`: 1.07197

Thus `sum Err <= D` is very close but not universally true.

## Extra decomposition stats (separate local diagnostic)

Per negative case, splitting `sum_err = even_err + odd_err`:

- `odd_err = 0` in 97,757 / 135,976 cases
- `odd_err > 0` in 38,219 / 135,976 cases
- max `odd_err / D`: 0.5489
- max `even_err / D`: 0.6157
- max `(odd_err + even_err) / D`: 1.07197

This strongly suggests the missing proof ingredient is an odd-reserve term, not an even-channel fix.

## Coarse shifted-channel reserve search (full `n<=18` negatives)

Using all 135,976 negative instances, define:

- `C10 = sum_i Lambda_old(i) P(k-i-1)Q(k-i)`
- `C01 = sum_i Lambda_old(i) P(k-i)Q(k-i-1)`
- `C11 = sum_i Lambda_old(i) P(k-i-1)Q(k-i-1)`

Grid search over `(alpha,beta,gamma) in {0,0.2,...,1.0}^3` for

`sum_err / (D + alpha*C10 + beta*C01 + gamma*C11)`.

Best coarse point:

- `alpha=1, beta=1, gamma=1`
- worst-case ratio: `0.5696` (well below 1)
- baseline worst ratio without reserve: `sum_err/D = 1.07197`

So the simple reserve

`R_shift := C10 + C01 + C11`

empirically dominates the full tail error budget on all tested negatives with strong margin.
This is a high-priority candidate for proof-level interpretation.

## Minimal scalar for shifted reserve

Using the same full `n<=18` negative corpus, solve for minimal scalar `lambda*` such that

`sum_err <= D + lambda * (C10 + C01 + C11)`

for every case.

Result:

- `lambda* = 0.05201381704686925`
- no zero-reserve pathologies (`R_shift=0` never needed where `sum_err>D`)
- tight witness:
  - `n=18`, graph6 `Q???????????????O?E??NwA_^?`, `root=0`, `step=2`, `k=5`
  - `D=1,442,539`, `sum_err=1,546,363`, `need=103,824`
  - `R_shift=1,996,085` so `need / R_shift = lambda*`.

This is substantially stronger/more interpretable than unconstrained coefficient fitting.

## Takeaway for Round 19 prompts

Most useful next target:

1. Keep the even-channel inequality as a stable base.
2. Add a minimal odd-residual reserve inequality.
3. Since all negatives are at `t=2`, focus theory and scans on step-2 closure.
