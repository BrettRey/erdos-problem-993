# Round 36 (Instance 1): Pi(n) Contract, Complexity, and Replay-Certificates

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

- Round 35 already gives canonical `Pi(n)` with deterministic refinement and certificates.
- Bucket closure is based on `(r_lower, r_upper)` and `g_B = alpha_B - lambda_B`.
- Forced witness separations at locked cutoffs are already identified (`n=23`, `n=24`).

## Task

Promote `Pi(n)` into a code-level contract with proof-grade guarantees:

1. State `Pi(n)` as strict typed pseudocode (input/output types, tie-breaks, invariants per loop step).
2. Prove time/space complexity in terms of number of records and distinct `(a,b)` classes.
3. Give a replay-certificate protocol that a second implementation can validate bit-for-bit.
4. Provide a minimal counterexample schema for `FAIL_CLASS` and `FAIL_DEGENERATE` that is both human-readable and machine-checkable.
5. State the exact frontier-minimality corollaries for n=23 and n=24 as reusable lemmas in this contract language.

## Output format

1. `Typed Pi(n) contract`
2. `Complexity and invariants`
3. `Replay-certificate protocol`
4. `Minimal failure schemas`
5. `Locked-cutoff minimality lemmas`
