# Round 33 (Instance 4): Classifier v9 Executable Repair Planner

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

- v8 added internal-deficit diagnostics (`Delta_line`, `Delta_line_mismatch`, `Delta_withinclass_max`) and shifted priority to the broken line (`a=2` at n=23).
- Locked frontier at n=23:
  - alpha witness `(2,19)`
  - lambda witness `(2,18)`
  - global gap `-0.05796615965948082`.

## Task

Turn v8 into an implementation-ready planner:

1. Define exact aggregate fields and formulas required from witness logs.
2. Specify deterministic decision logic for regime selection:
   - cross-class mismatch vs internal line deficit vs within-line mismatch.
3. Give explicit priority score formulas and tie-break rules for actions A/B/C/D.
4. Produce the n=23 ranked queue template with conditional branches only where diagnostics decide.
5. Provide a minimal end-to-end automation pipeline (inputs, aggregations, queue output, stop checks).

## Output format

1. `Required aggregates and formulas`
2. `Regime-selection decision logic`
3. `Priority score and tie-break specification`
4. `Ranked queue template for n=23`
5. `Executable pipeline and stop checks`
