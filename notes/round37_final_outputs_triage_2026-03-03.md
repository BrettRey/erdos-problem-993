# Round 37 Final Outputs Triage (2026-03-03)

Source: user-provided final Round 37 outputs in chat (four-instance pack).

## Scope captured

The final pack delivers implementation-grade contracts for:

1. `Pi(n)` reference implementation + replay verifier contract.
2. A2 local solver + ledger commit validator contract.
3. Obligation runtime API + deterministic escalation engine.
4. Orchestrator v13 codegen plan + artifact schemas + replay taxonomy.

No new empirical scan claims were introduced in the pack; this is architecture/spec material.

## Instance 1: `Pi(n)` implementation contract (captured)

Key deliverables extracted:

- Package layout `pi_n/` with modules:
  - `types.py`, `ratio.py`, `ordering.py`, `algorithm.py`, `transcript.py`, `serialize.py`, `verify.py`, `errors.py`.
- Deterministic public API:
  - `compute_pi(log) -> PiResult`
  - `verify_pass_certificate(...)`
  - `verify_fail_class(...)`
  - `verify_fail_degenerate(...)`.
- Exact arithmetic requirement:
  - no floats; ratio comparison via cross multiplication.
- Canonical transcript replay requirement:
  - base split `a<=3` vs `a>=4`, replayed split transcript, canonical first-failing-bucket / lambda-witness-class split checks.
- Minimal tests required:
  - PASS (no split), PASS (split), FAIL_DEGENERATE, FAIL_CLASS, tie-break determinism, input permutation determinism.

## Instance 2: A2 solver implementation contract (captured)

Key deliverables extracted:

- Typed data model for:
  - `A2Record`, half-planes, lines, solver results, Type III certificate payloads.
- Deterministic pipeline:
  - Type I/II prechecks
  - 1D reductions (`lambda=0`, `rho=0`)
  - 2D vertex enumeration
  - canonical Type III extraction (size-2 preferred, then size-3 Helly subset).
- Single-epsilon numeric policy:
  - one `eps` reused across feasibility, tie-breaks, determinant checks, clamping.
- Commit validator contract:
  - acceptance invariants checked before append-only ledger event construction.
- Minimum solver test set:
  - feasible 1D-rho, 1D-lambda, 2D-vertex, and infeasible Type I/II/III.

## Instance 3: obligation runtime + escalation contract (captured)

Key deliverables extracted:

- Runtime state schema:
  - reserve basis spec, bucket map, witness assignment map, params records (`sandwich` vs `k1`), transition history.
- Deterministic bucket validator:
  - param-shape checks, coordinatewise margin checks, LF/UF slack checks over witness x k grids, deterministic first-failure selection.
- Escalation transitions:
  - `T0` neutral rewrite,
  - `T1` singleton split of failing witness,
  - `T2` A2 threshold split,
  - `T3` reserve-basis expansion.
- Progress metrics:
  - finite potential drop under `T1/T2`; bounded-step behavior with reserve cap.
- Orchestrator handoff contract:
  - state snapshot + data bundle in, validation report + updated snapshot + work items out.

## Instance 4: orchestrator v13 codegen plan (captured)

Key deliverables extracted:

- Module map:
  - `frontier`, `partition`, `a2`, `hardstep`, `queue`, `obligations`, `replay_verify`, `run`.
- Deterministic pipeline:
  - authoritative filter -> derive -> aggregate -> candidate partitions `PI0..PI4` -> evaluate/select -> queue -> obligations.
- Strict artifact schemas:
  - `frontier_state.json`, `partition_plan.json`, `repair_queue.json`, `v13_obligations.json` with explicit `JsonReal` encoding (`number | "inf" | "-inf"`).
- Conformance harness:
  - fixture-driven tests (T0..T7), golden byte checks via canonical JSON serializer, replay-verify pass requirement.
- Authority-gating hard assertions:
  - non-authoritative runs must not output `CLOSED`; only `OPEN` or `UNKNOWN`.

## Recommended execution order (next coding phase)

1. Implement `pi_n` package first (smallest isolated deterministic core).
2. Implement A2 solver package with certificate extraction and ledger validator.
3. Implement orchestrator v13 core with candidate partition evaluation and schema-valid artifact writes.
4. Implement replay verifier and conformance harness; require green fixtures before integrating into main workflow.

## Locked baseline reminder

Authoritative locked baseline remains:

- `n=24` alpha witness `(2,20)` at `0.16161242603550297`
- `n=24` lambda witness `(4,18)` at `0.280781720999777`
- gap `-0.11916929496427403`

`n=25` remains non-authoritative until lambda dict reaches full `256/256` key completion.
