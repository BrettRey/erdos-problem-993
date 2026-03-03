# Round 34 (Instance 4): Orchestrator v10 for Partition + Local-Module + Proof Obligations

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

- v9 planner computes per-record scores (`Gamma_G`, `Gamma_L`, `Alpha_req`, `Gamma_M*`, `Gamma_D`) and aggregates by class/line.
- Round 33 added:
  - three-bucket theorem interface,
  - exact `(rho,lambda)` local repair feasibility,
  - v11 finite hard-step theorem interface.
- n=24 complete frontiers:
  - alpha witness `(2,20)`, value `0.16161242603550297`
  - lambda witness `(4,18)`, value `0.280781720999777`
  - gap `-0.11916929496427403`.

## Task

Turn the planner into an orchestrator that outputs a complete proof-repair plan from logs:

1. Define required inputs and derived aggregates (record/class/line/bucket) with exact formulas.
2. Define deterministic regime selection that chooses among:
   - cross-line mismatch repair,
   - same-line internal repair,
   - within-class deficit repair,
   - optional activation of the `a=2` `(rho,lambda)` module.
3. Define deterministic action scoring and tie-breaks for combined actions:
   - partition updates,
   - A/B/C channel repairs,
   - local-module calibration runs,
   - theorem-interface updates (v12 obligations).
4. Produce the n=24 queue template (with only diagnostic-conditioned branches).
5. Provide an executable end-to-end pipeline with stop criteria tied to the selected partition and obligations.

## Output format

1. `Input/aggregate schema and formulas`
2. `Deterministic regime-selection state machine`
3. `Unified action scoring and tie-break rules`
4. `Ranked queue template for n=24`
5. `Executable pipeline with closure stop checks`
