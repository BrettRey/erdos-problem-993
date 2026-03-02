# Round 22 (Instance 1): Neighbor-Transport Theorem (Replace Failed Local Lemma)

Use only your prior outputs plus this prompt. Do not claim new scans.

## Locked context

Conventions unchanged (support-root, step-prefix `k<mode(I_new)`, boundary-correct).

Known through `n<=18`:

- global closure still holds:
  `sum_err <= D + lambda0*(C10+C01+C11)` with `lambda0=0.05201381704686925`
- old per-odd local lemma fails:
  `err_{2t+1} <= lambda0*Lambda_t*(10+01+11)` fails in `38,642/587,674` odd checks
- stronger `(Lambda_t+Lambda_{t+1})` local form also fails (`1,600` failures)

New locked diagnostics:

Define
- `d_t := max(0, err_{2t+1} - lambda0*Lambda_t*(10+01+11))`
- `e_t := max(0, Lambda_t*P(k-t)Q(k-t) - err_{2t})`

Observed through `n<=18` over all `X<0` cases:

1. strict local `d_t <= e_t` fails (`393` fails over `723,650` odd slots)
2. neighbor local `d_t <= e_t + e_{t+1}` has `0` failures
3. per-case aggregate `sum_t d_t <= sum_t e_t` has `0` failures

## Task

Promote this to a theorem template replacing the failed per-odd lemma.

1. State exact theorem candidate using `{d_t,e_t}` and existing channel definitions.
2. Show how `d_t <= e_t + e_{t+1}` implies global closure.
3. Isolate one missing algebraic lemma (if needed) that is weaker than old per-odd control.
4. Explain why this naturally allows local failures of `d_t<=e_t` but still closes globally.

## Output format

1. `Revised theorem candidate`
2. `Proof skeleton`
3. `Single missing lemma`
4. `Why this matches the new diagnostics`
5. `Finite falsification recipe`
