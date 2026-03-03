# Round 37 Prompt Prep (2026-03-03)

Goal: final conversion pass from Round 36 contract text to executable, replayable implementation artifacts.

## Inputs received (Round 36, 4 outputs)

1. Typed `Pi(n)` contract with exact-rational comparisons, deterministic tie-breaks, and PASS/FAIL certificate schemas.
2. Typed A2 solver contract with deterministic `eps` policy, exact 1D/2D feasibility flow, and canonical Type III extraction.
3. Induction theorem v14 framed as a unified `Obligation(module,bucket,witness_set,params)` API with deterministic escalation.
4. Orchestrator v12 machine-executable pseudocode, strict authority gating, finite partition family, conformance tests, and replay verifier.

## Locked baseline to preserve

Authoritative complete frontier remains n=24:

- alpha witness `(2,20)`, value `0.16161242603550297`
- lambda witness `(4,18)`, value `0.280781720999777`
- gap `-0.11916929496427403`

n=25 remains partial/non-authoritative unless and until 256/256 authoritative shards complete.

## Round 37 objectives

1. Demand code-ready implementation blueprints (module boundaries, function signatures, test vectors, verifier checks), not additional abstract restatements.
2. Require deterministic replay and certificate validation as first-class outputs.
3. Force explicit integration points between `Pi(n)`, A2 solver, orchestrator, and obligation API.
4. Produce minimal unresolved-core checklist tied to current n=24 authoritative baseline and n>=25 escalation hooks.

## Prompt constraints

- Use only locked constants/objects and prior outputs.
- No fresh scan claims.
- Keep deterministic 5-section output format.
