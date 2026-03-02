# Round 22b (Instance 2): Compensation Identity with Explicit Defect Term

Use only your prior outputs plus this prompt.

## Locked context update

Old local/neighbor claims are no longer exact globally:

- `d_t <= e_t` fails (many)
- `d_t <= e_t + e_{t+1}` was exact through `n<=18`, but at `n=19`:
  - fails in `5` odd slots
- aggregate `sum_t d_t <= sum_t e_t` was exact through `n<=18`, but at `n=19`:
  - fails in `33` cases

Definitions:
- `d_t := max(0, err_{2t+1} - lambda0*Lambda_t*(10+01+11))`
- `e_t := max(0, Lambda_t*P(k-t)Q(k-t) - err_{2t})`
- `lambda0 = 0.05201381704686925`

## Task

Derive a compensation law with an explicit nonnegative defect term:

`d_t <= e_t + e_{t+1} + b_t`

and a summable bound on `b_t` that yields a corrected global theorem.

Requirements:
1. Give explicit formula or tight upper envelope for `b_t`.
2. Show how summing gives
   `sum_err <= D + lambda0*R_shift + Defect`.
3. Provide one candidate to absorb `Defect` (extra channel or delta-lambda term).
4. Keep constants/channels minimal.

## Output format

1. `Local compensated inequality`
2. `Definition of defect term b_t`
3. `Global corollary with Defect`
4. `Single missing sublemma`
5. `How this matches the 5 and 33 failures`
