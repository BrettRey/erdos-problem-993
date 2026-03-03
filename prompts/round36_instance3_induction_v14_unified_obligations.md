# Round 36 (Instance 3): Induction v14 with Unified Module Obligation API

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

- Round 35 provides v13 module-aligned induction (`ENV/A2_LOCAL/HARDSTEP`).
- Hard-step minimal template already separates known `n=23` and `n=24` conflict geometries.
- Bridge remains a deferred single-lemma path.

## Task

Produce v14 as a theorem/API hybrid suitable for direct orchestrator consumption:

1. Define a single obligation API signature `Obligation(module, bucket, witness_set, params)` that all branches consume.
2. State v14 theorem entirely in terms of this API and ID-based closure algebra.
3. Isolate minimal unresolved lemma obligations per module path (exactly one core unresolved item for ENV, A2_LOCAL, HARDSTEP).
4. Provide a non-circular proof skeleton that compiles from API obligations to STP2 preservation.
5. Give a deterministic escalation transition system on bucket/module states under further drift (n>=25).

## Output format

1. `Obligation API definition`
2. `Induction theorem v14 (API form)`
3. `Minimal unresolved core by module`
4. `Non-circular proof skeleton`
5. `Deterministic escalation state machine`
