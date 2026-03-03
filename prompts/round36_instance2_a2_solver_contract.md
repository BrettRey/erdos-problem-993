# Round 36 (Instance 2): A2 Solver Contract with Exact Certificate Extraction

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

- Round 35 defines the normalized A2 data contract `(gap, slack, Psi, Phi)` and exact LP feasibility in `(rho,lambda)`.
- Type I/II/III infeasibility certificates are already specified.
- Composition is ledger-based and deterministic.

## Task

Convert the A2 module into a production solver contract:

1. Give strict typed pseudocode for the full solver pipeline (prechecks, 1D reductions, 2D exact solver, certificate extraction).
2. Provide worst-case complexity bounds and practical pruning rules that preserve exactness.
3. Specify deterministic numeric tolerance handling (single `eps` policy) and prove decision stability under that policy.
4. Define a canonical minimal-certificate extractor for Type III (<=3 constraints) with deterministic ordering.
5. Provide an executable acceptance theorem: if solver returns FEASIBLE and all acceptance invariants pass, the committed ledger update is composition-safe.

## Output format

1. `Typed solver pseudocode`
2. `Complexity and exactness-preserving pruning`
3. `Deterministic tolerance policy`
4. `Canonical Type III extractor`
5. `Acceptance theorem for ledger-safe commit`
