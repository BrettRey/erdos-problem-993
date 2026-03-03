# Round 32 (Instance 3): Induction Pack v10 for Internal Same-Line Break

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

- k=1 micro-lemma branch remains independent and robust.
- v9 already replaced single global envelope with repaired k>=2 branch interfaces.
- n=23 shows stronger break with both extremisers in `a=2` (`(2,19)` for alpha, `(2,18)` for lambda), so k>=2 failure is now internal to one hard-step line.

## Task

Write induction v10 with explicit post-n23 obligations:

1. Preserve three-branch structure (`k=1`, hard step `t=2`, bridge `t>=3`).
2. Give a k>=2 repaired branch schema that can handle internal same-line failure (e.g., finite hard-step class partition inside `a<=3`, or reserve-vector/channel-coordinate margins).
3. Keep unresolved interfaces minimal and non-circular.
4. State what must be proved immediately at hard step `t=2` versus what can be deferred to bridge.
5. Include escalation if the repaired hard-step interface drifts again at n>=24.

## Output format

1. `Induction theorem v10 (post-n23 break)`
2. `k>=2 repaired branch schema`
3. `Dependency graph (non-circular, minimal interfaces)`
4. `Hard-step vs bridge obligations`
5. `Escalation path under further drift`
