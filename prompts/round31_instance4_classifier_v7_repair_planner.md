# Round 31 (Instance 4): Classifier v7 as Repair Planner

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

- Existing score stack:
  - `Gamma_G` (global lambda demand)
  - `Gamma_L` (local neighbour stress)
  - `Gamma_M*` (alpha margin to current calibration)
  - `Gamma_D` (drift/erosion score)
- n=22 snapshot indicates global envelope break (`alpha_front < lambda_front`).

## Task

Upgrade classifier v6 into a repair planner that chooses what to fix first:

1. Define one action-space over repairs:
   - raise lower-face alpha,
   - lower upper-face lambda,
   - add/modify reserve channels,
   - class-stratify obligations.
2. Define a costed priority score per action-class pair using existing gammas.
3. Produce a queue of “highest leverage repairs” starting from n=22 witnesses.
4. Add stop criteria: when to accept a repair as sufficient for induction closure.
5. Give a minimal pipeline to compute this automatically from witness logs + channel stats.

## Output format

1. `Repair action space`
2. `Costed priority score`
3. `Ranked repair queue from n=22 baseline`
4. `Stop criteria for closure readiness`
5. `Minimal automated pipeline`

