# Round 22 (Instance 2): Derive a Compensation Identity for d_t and e_t

Use only your prior outputs plus this prompt. No new computation claims.

## Locked facts

- per-odd local bound fails; `(Lambda_t+Lambda_{t+1})` local bound also fails
- but neighbor compensation law is empirically exact through `n<=18`:
  `d_t <= e_t + e_{t+1}` (0 failures)
- and aggregate `sum_t d_t <= sum_t e_t` also has 0 failures

Definitions fixed:

- `d_t := max(0, err_{2t+1} - lambda0*Lambda_t*(10+01+11))`
- `e_t := max(0, Lambda_t*P(k-t)Q(k-t) - err_{2t})`
- `lambda0 = 0.05201381704686925`

## Task

Starting from your exact odd decomposition machinery, derive a **symbolic neighbor-compensation identity/inequality** that has this shape:

`d_t <= e_t + e_{t+1}`.

Requirements:

1. Show where the `t+1` channel enters algebraically (not heuristic).
2. Make the inequality local in `t` and in fixed arrays `(A,B,P,Q,k)`.
3. Keep free constants minimal (ideally none besides fixed `lambda0`).
4. Show summation over `t` yields global closure form.

## Output format

1. `Local compensation formula`
2. `Derivation path`
3. `Global corollary`
4. `What remains to prove`
5. `One explicit counterexample to old local lemma + why neighbor form survives`
