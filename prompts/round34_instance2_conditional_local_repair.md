# Round 34 (Instance 2): Conditional Local Repair Module from (rho,lambda) Constraints

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

- Round 33 gave the exact `a=2` repair feasibility system:
  - transfer `M = rho*Psi + lambda*Phi`
  - per witness constraints: `gap <= M <= slack`
  - feasibility is a 2D linear strip intersection in `(rho,lambda)`.
- n=23 conflict was internal to `a=2` (`(2,19)` alpha vs `(2,18)` lambda).
- n=24 global witnesses moved to different lines:
  - alpha witness `(2,20)`
  - lambda witness `(4,18)`.

## Task

Upgrade the local-repair result into a conditional module in a larger policy:

1. Define exact trigger conditions (in terms of line/class diagnostics) for when the `a=2` `(rho,lambda)` module should be activated.
2. Define a deactivation condition for when this module should be skipped (e.g., cross-line dominated regime).
3. Give a complete solver workflow:
   - witness-first calibration on stress classes,
   - full-population feasibility confirmation,
   - infeasibility certificate if no pair exists.
4. Define compatibility constraints so this local module can coexist with bucket/partition repairs without double-counting reserve.
5. Give a finite acceptance checklist and minimal telemetry fields to log each activation run.

## Output format

1. `Activation/deactivation rule for the a=2 local module`
2. `Exact feasible-system interface (reused and normalized)`
3. `Calibration and certification workflow`
4. `Composition rules with partition-level repairs`
5. `Run-time acceptance checklist`
