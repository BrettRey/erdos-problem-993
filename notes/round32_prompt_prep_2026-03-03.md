# Round 32 Prompt Prep (2026-03-03)

Goal: adapt post-break repair strategy after n=23 worsens the global envelope break and moves both frontier witnesses into the same small-first line.

## Locked baseline at n=23

From completed mod=256 artifacts:

- `results/alpha_bookkeeping_modal_n23_n23_w256.json`
  - `min_alpha_all = 0.18243252946998884`
  - witness class `(a,b)=(2,19)` at step 2, `k=9`
- `results/lambda_frontier_modal_n23_n23_w256.json`
  - `lambda_max = 0.24039868912946966`
  - witness class `(a,b)=(2,18)` at step 2, `k=5`
- derived gap:
  - `alpha_front(23) - lambda_front(23) = -0.05796615965948082`

Drift vs n=22:

- `delta_alpha = -0.005154294480471472`
- `delta_lambda = +0.043562639089012145`
- `delta_gap = -0.04871693356948362`

## Implication for repairs

- The n=22 break could be interpreted as cross-bucket mismatch (`a<=3` vs `a>=4` extremisers).
- At n=23 both extremisers are in `a=2`, so a simple two-bucket repair is insufficient.
- Round 32 should therefore target:
  1. a finer theorem-level partition or reserve-vector repair inside small-first classes,
  2. a local correction family that resolves both `(2,18)` and `(2,19)` stress patterns,
  3. an induction v10 interface with explicit hard-step finite obligations,
  4. classifier/planner upgrade from mismatch-first to internal-deficit-first diagnostics.

## Prompt design constraints

- Keep all-diagonal transport framework (`sum_all`) and boundary-correct indexing fixed.
- Do not re-open odd-only routes.
- Do not rely on fresh scans; use locked n=23 constants and witness classes.
- Require explicit falsification/certification interfaces in each output.
