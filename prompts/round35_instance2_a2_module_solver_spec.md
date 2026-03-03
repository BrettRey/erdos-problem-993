# Round 35 (Instance 2): A2 Local Module as Solver-Ready Optimization Spec

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

- The `a=2` local module is conditional and uses transfer `M = rho*Psi + lambda*Phi`.
- Per-instance constraints are exact strips:
  - `gap_r <= rho*Psi_r + lambda*Phi_r <= slack_r`, with `rho,lambda >= 0`.
- Round 34 added activation/deactivation rules and composition ledgers.

## Task

Produce a canonical solver specification with deterministic diagnostics:

1. Normalize the full feasible-system interface into a minimal data contract (`gap, slack, Psi, Phi`, labels).
2. Give a complete 2D feasibility algorithm that returns either:
   - feasible `(rho,lambda)` with optimality criterion, or
   - explicit infeasibility certificate (Type I/II/III).
3. Derive fast 1D special-case tests (`lambda=0` and `rho=0`) and exact conditions when both parameters are necessary.
4. Specify composition-safe integration rules with other modules using a single ledger (no double counting).
5. Define strict runtime telemetry and acceptance invariants so every activation run is auditable/replayable.

## Output format

1. `Normalized data contract`
2. `Exact feasibility/optimization algorithm`
3. `1D reductions and necessity diagnostics`
4. `Composition-safe integration contract`
5. `Auditable run ledger and acceptance invariants`
