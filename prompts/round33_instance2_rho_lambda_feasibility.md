# Round 33 (Instance 2): Feasibility Region for a=2 Local Repair (rho, lambda)

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

- Local packet-pair backbone remains fixed.
- Round 32 repair for a=2 introduced:
  - `Psi_t` (center probe)
  - `Phi_t` (edge-band/radius-2 probe)
  - transfer mass `M_t = rho*Psi_t + lambda*Phi_t` moved from `a=2` slot to `a=3` sink.
- Target constant remains `alpha_* = 0.21034113597068071` in this local route.
- Stress classes to anchor calibration: `(2,19)` and `(2,18)`.

## Task

Convert the repair into an explicit feasible-parameter problem:

1. Derive the exact linear inequality system in `(rho, lambda)` induced by repaired `a=2` and sink `a=3` constraints (`gap <= M <= slack`).
2. State the feasible region geometry (empty/nonempty, rays, minimal corner points) in terms of witness-derived coefficients.
3. Give a calibration algorithm that returns a provably feasible pair `(rho, lambda)` if one exists.
4. Provide minimality diagnostics (when `rho` alone fails, when `lambda` alone fails, when both are necessary).
5. Give a finite verification checklist to certify the calibrated pair over the full `a=2` population.

## Output format

1. `Exact (rho, lambda) constraint system`
2. `Feasible-region geometry`
3. `Calibration algorithm from witness data`
4. `Minimality diagnostics`
5. `Finite certification checklist`
