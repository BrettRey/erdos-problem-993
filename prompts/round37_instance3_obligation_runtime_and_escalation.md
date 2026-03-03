# Round 37 (Instance 3): Obligation Runtime API + Deterministic Escalation Engine

Use only prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

- Round 36 provides induction v14 in API form: `Obligation(module,bucket,witness_set,params)`.
- Module set is fixed: `ENV`, `A2_LOCAL`, `HARDSTEP`, plus `K1` branch handling.
- Deterministic escalation states/transitions (`T0..T3`) are defined conceptually.

## Task

Specify a runtime-checkable engine for obligations and escalation:

1. Define strict typed schemas for:
   - module/bucket assignment state,
   - witness sets,
   - params records (sandwich vs K1),
   - per-bucket validation outcomes.
2. Provide a deterministic `validate_state(state, data) -> report` algorithm that checks all obligations bucketwise.
3. Provide a deterministic escalation executor implementing `T0..T3` with explicit preconditions, postconditions, and tie-break ordering.
4. Provide termination/progress metrics for fixed `n` and explicit failure conditions when reserve dimension is insufficient.
5. Provide integration contract for orchestrator handoff: how state snapshots, failed witnesses, and refined buckets are serialized back into obligations artifacts.

## Output format

1. `Typed runtime state schema`
2. `Bucket validation algorithm`
3. `Escalation executor (T0..T3)`
4. `Progress and termination metrics`
5. `Orchestrator handoff contract`
