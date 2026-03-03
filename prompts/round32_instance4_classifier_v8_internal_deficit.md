# Round 32 (Instance 4): Classifier v8 for Internal-Deficit Repair Planning

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

- Existing score stack:
  - `Gamma_G` (global lambda demand)
  - `Gamma_L` (local neighbour stress)
  - `Gamma_M*` (alpha margin to current calibration)
  - `Gamma_D` (alpha drift/erosion)
- n=23 frontier constants:
  - `alpha_front = 0.18243252946998884`
  - `lambda_front = 0.24039868912946966`
  - gap `= -0.05796615965948082`
- Frontier witnesses:
  - alpha witness class `(2,19)`
  - lambda witness class `(2,18)`
- Therefore failure is no longer primarily cross-bucket mismatch; it is an internal same-line deficit.

## Task

Upgrade classifier v7 into v8 for this new regime:

1. Add diagnostics that explicitly separate cross-class mismatch from internal within-line deficit.
2. Update action prioritization so the planner chooses correctly when alpha and lambda frontier classes are both in `a=2`.
3. Define a costed priority score per `(action, class)` using only existing gammas and frontier aggregates.
4. Produce a ranked repair queue for the n=23 baseline.
5. Define stop criteria and a minimal automated pipeline from witness logs.

## Output format

1. `Internal-deficit diagnostics`
2. `Costed priority score (v8)`
3. `Ranked repair queue from n=23 baseline`
4. `Stop criteria for closure readiness`
5. `Minimal automated pipeline`
