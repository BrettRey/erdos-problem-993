# Round 21 Neighbor Compensation Scan (2026-03-02)

## Scope

Exhaustive through `n<=18` in the locked regime:

- support-root processing
- step-prefix (`k < mode(I_new)`, smallest mode index)
- boundary-correct indexing
- only `X_k<0` cases

## Definitions tested

Let `lambda0 = 0.05201381704686925`.

For each odd diagonal `s=2t+1`:

- `d_t := max(0, err_{2t+1} - lambda0 * Lambda_old(t) * (10+01+11))`

For each even diagonal `s=2t`:

- `e_t := max(0, Lambda_old(t) P(k-t)Q(k-t) - err_{2t})`

with `err_s` matching the Round 19/21 script definition.

## Results

- Total negative cases: `135,976`
- All at step `t=2` (`xneg_step_not2 = 0`)

Global checks:

1. `sum_t d_t <= sum_t e_t` per case:
   - failures: `0`
2. Local check `d_t <= e_t`:
   - failures: `393` out of `723,650` odd slots
3. Local neighbor check `d_t <= e_t + e_{t+1}`:
   - failures: `0` out of `723,650` odd slots

So the data strongly support a **neighbor-compensation law**:

`d_t <= e_t + e_{t+1}`

while ruling out the stricter single-channel `d_t <= e_t`.

## Extremals

Worst ratio for failed strict check `d_t/e_t`:

- ratio `7.66616`
- witness:
  - `n=18`
  - `g6=Q???????????????O?F??f_?@~G`
  - `root=0, k=6, t=0`
  - `d_t=450,916.04`
  - `e_t=58,819`
  - `e_t+e_{t+1}=743,205` (passes neighbor form)

## Pair-class pattern

Failures of strict `d_t <= e_t` are concentrated in the same high classes:

- `(2,14)`: 241 fails
- `(2,13)`: 78 fails
- `(3,13)`: 59 fails
- smaller counts in `(2,12),(3,12),(2,11)`
- none in `(2,10),(3,11),(2,9)`

Neighbor form `d_t <= e_t + e_{t+1}` had **0 failures in every class**.

## Implication

The likely missing theorem is not per-diagonal odd control by one even channel,
but a **two-channel local transport bound** that shifts one step in `t`:

`odd_deficit_t <= even_slack_t + even_slack_{t+1}`.

This is a concrete replacement for the failed per-odd local lemma and is compatible with the verified global closure.

## Artifacts

- Script: `scan_round21_neighbor_compensation.py`
- JSON: `results/round21_neighbor_compensation_n18.json`
