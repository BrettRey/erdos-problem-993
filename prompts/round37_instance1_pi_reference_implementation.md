# Round 37 (Instance 1): `Pi(n)` Reference Implementation + Replay Verifier

Use only prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

- Round 36 provides typed `Pi(n)` contract with exact arithmetic, deterministic tie-breaks, and PASS/FAIL certificates.
- Locked forced separations at `n=23` and `n=24` are already fixed.

## Task

Produce an implementation-grade package design for `Pi(n)`:

1. Provide concrete Python module layout and exact function signatures for:
   - `compute_pi(log) -> PiResult`
   - `verify_pass_certificate(log, cert)`
   - `verify_fail_class(log, cert)`
   - `verify_fail_degenerate(log, cert)`
2. Give exact-rational comparison and canonical serialization rules that avoid float usage entirely.
3. Provide deterministic split-transcript schema and replay algorithm that reconstructs final buckets bit-for-bit.
4. Provide a minimal but complete `unittest` matrix that exercises PASS, FAIL_CLASS, FAIL_DEGENERATE, and tie-break edge cases.
5. State integration hooks so orchestrator can call `Pi(n)` as a deterministic subroutine (input contract, output contract, failure propagation).

## Output format

1. `Module layout and typed signatures`
2. `Exact arithmetic and serialization contract`
3. `Replay and certificate verification algorithm`
4. `Deterministic test matrix`
5. `Orchestrator integration hooks`
