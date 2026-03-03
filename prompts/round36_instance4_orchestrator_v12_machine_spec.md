# Round 36 (Instance 4): Orchestrator v12 Machine-Executable Spec

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

- Round 35 provides orchestrator v11 with deterministic main loop, finite partition family, module semantics, and authority rules.
- Authoritative baseline remains n=24; n=25 remains partial/non-authoritative.

## Task

Upgrade to v12 as a machine-executable spec with deterministic replay:

1. Provide full end-to-end typed pseudocode including all helper procedures (`BUILD_PARTITIONS`, `EVAL_PARTITION`, `BUILD_REPAIR_QUEUE`, `BUILD_V13_OBLIGATIONS`, `CHECK_CLOSURE`).
2. Define exact objective computation and tie-break hierarchy so independent implementations produce identical `partition_plan` and `repair_queue`.
3. Add explicit partial-data gating rules as preconditions/postconditions in the algorithm (not prose-only).
4. Specify a conformance test suite: minimal synthetic logs that each force a distinct regime/module path.
5. Provide a deterministic replay verifier that checks `frontier_state`, `partition_plan`, `repair_queue`, and obligations for consistency.

## Output format

1. `Typed end-to-end pseudocode`
2. `Deterministic objective/tie-break contract`
3. `Authority-gating pre/postconditions`
4. `Conformance test matrix`
5. `Replay verifier specification`
