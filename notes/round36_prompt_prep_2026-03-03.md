# Round 36 Prompt Prep (2026-03-03)

Goal: turn Round 35 canonical specs into implementation-ready theorem/solver/orchestrator deliverables with machine-checkable certificates.

## Inputs received (Round 35, 4 outputs)

1. Canonical `Pi(n)` algorithm with deterministic tie-breaks, termination, soundness, class-granularity completeness boundary, and certificate schemas.
2. A2 local module as exact 2D feasibility/optimization spec with Type I/II/III infeasibility certificates and composition-safe ledger rules.
3. Induction theorem v13 aligned to modules `ENV/A2_LOCAL/HARDSTEP` with minimal hard-step class template and deterministic escalation.
4. Orchestrator v11 deterministic execution spec with finite partition family, objective/tie-break policy, module semantics, authoritative-data gating, and strict output schemas.

## Locked baseline to preserve

Authoritative complete frontier remains n=24:

- alpha witness `(2,20)`, value `0.16161242603550297`
- lambda witness `(4,18)`, value `0.280781720999777`
- gap `-0.11916929496427403`

n=25 remains partial/non-authoritative.

## Round 36 objectives

1. Consolidate theorem + orchestrator interfaces into one canonical contract that can be translated directly into code and proof obligations.
2. Demand explicit pseudocode and complexity bounds for all finite algorithms (`Pi(n)`, A2 solver, hardstep list constructor).
3. Require machine-verifiable certificate formats (success + all failure modes) with deterministic replay rules.
4. Force one unresolved-core lemma statement per module path (ENV, A2_LOCAL, HARDSTEP), separating algebraic closure from empirical calibration.

## Prompt constraints

- Use only locked constants/objects and prior outputs.
- No new scan claims.
- Keep deterministic 5-section output format.
