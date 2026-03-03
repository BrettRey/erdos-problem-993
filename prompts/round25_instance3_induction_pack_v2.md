# Round 25 (Instance 3): Near-Final Induction Pack v2

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

- k=1 micro-lemma is strongly supported (`1,407,504` checks, `0` failures).
- X<0 profile n<=19: step2 dominates; step3 has only 12 cases, all at k=1.
- Single-constant variant A is currently the clean baseline (`lambda19`).
- Split variant B must be reformulated to avoid impossible side conditions.

## Task

Write a near-final theorem package with two viable variants:

A) single-constant (`sum_err <= D + lambda19*R_shift`),
B) split-reserve using an **extended bookkeeping** form from Instance 2.

Requirements:
1. Keep explicit three-branch structure (`k=1`, `k>=2,t=2`, `k>=2,t>=3`).
2. Provide non-circular dependency graph.
3. Give exactly one unresolved lemma per variant.
4. State which variant is manuscript-default and why.

## Output format

1. `Variant A theorem`
2. `Variant B theorem (obstruction-safe)`
3. `Dependency graph`
4. `Unresolved lemma(s)`
5. `Recommendation`
