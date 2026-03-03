# Round 37 (Instance 2): A2 Solver Implementation Contract + Ledger Commit Validator

Use only prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

- Round 36 defines typed A2 solver pseudocode, deterministic `eps` policy, and Type I/II/III infeasibility certificates.
- Ledger-safe composition and commit acceptance theorem are already specified.

## Task

Turn the contract into an implementation blueprint suitable for direct coding:

1. Define Python data classes and function signatures for records, half-planes, solver result unions, and certificates.
2. Provide executable pseudocode for the exact pipeline with deterministic ordering:
   - prechecks,
   - 1D reductions,
   - 2D vertex enumeration,
   - canonical Type III extraction.
3. Provide a deterministic numeric policy section that pins every comparison and tie-break to one `eps`.
4. Provide a commit validator algorithm that checks ledger safety invariants before append.
5. Provide a high-signal unit/integration test set, including:
   - feasible 1D-rho,
   - feasible 1D-lambda,
   - feasible 2D-vertex,
   - Type I,
   - Type II,
   - Type III.

## Output format

1. `Data classes and function signatures`
2. `Deterministic solver pipeline`
3. `Single-eps numeric policy`
4. `Ledger commit validator`
5. `Test suite blueprint`
