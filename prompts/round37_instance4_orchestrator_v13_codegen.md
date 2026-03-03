# Round 37 (Instance 4): Orchestrator v13 Codegen Plan with Replay Guarantees

Use only prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

- Round 36 provides orchestrator v12 typed pseudocode with deterministic partition selection and replay verifier.
- Authoritative baseline is locked at `n=24`; `n=25` is non-authoritative unless completed.

## Task

Produce a code-generation grade v13 plan that can be implemented directly:

1. Define concrete Python module boundaries (`frontier`, `partition`, `a2`, `hardstep`, `queue`, `obligations`, `replay_verify`) and public APIs.
2. Provide deterministic implementation pseudocode for each stage, including exact candidate partition construction (`PI0..PI4`) and key-based selection.
3. Provide strict JSON schema contracts for:
   - `frontier_state.json`,
   - `partition_plan.json`,
   - `repair_queue.json`,
   - `v13_obligations.json`.
4. Provide a conformance harness that runs the full pipeline on synthetic fixtures and checks byte-stable outputs under canonical serialization.
5. Provide authoritative/non-authoritative gating assertions as executable pre/postconditions and a failure taxonomy for verifier errors.

## Output format

1. `Module/API map`
2. `Deterministic stage pseudocode`
3. `Artifact JSON schemas`
4. `Conformance harness design`
5. `Authority-gating assertions and verifier taxonomy`
