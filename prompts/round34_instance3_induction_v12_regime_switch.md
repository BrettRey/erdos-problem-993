# Round 34 (Instance 3): Induction v12 with Regime-Switched Hard-Step Interfaces

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

- Branch `k=1` remains closed by the micro-lemma.
- v11 hard step (`t=2`) uses finite class partition + channel-vector sandwich.
- Known hard-step witness geometry:
  - n=23: internal conflict in `a=2` (`(2,18)` vs `(2,19)`)
  - n=24: cross-line witness split (`(2,20)` alpha vs `(4,18)` lambda).

## Task

Produce v12 as a regime-switched theorem interface that stays finite and non-circular:

1. Propose a smallest practical hard-step partition template that can represent both known regimes (same-line and cross-line).
2. State the v12 theorem with class-indexed lower/upper faces and coordinatewise margins, allowing optional neutral cancellation terms.
3. Separate exactly what must be proved immediately at hard step `t=2` from what may remain deferred at bridge `t>=3`.
4. Give a dependency graph proving non-circularity in final theorem style.
5. Add an escalation rule that is deterministic when a new internal conflict appears at n>=25.

## Output format

1. `Induction theorem v12 (regime-switched finite interface)`
2. `Hard-step partition template covering n=23 and n=24 geometries`
3. `Immediate vs deferred obligations`
4. `Dependency graph (non-circular)`
5. `Deterministic escalation rule`
