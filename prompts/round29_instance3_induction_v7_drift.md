# Round 29 (Instance 3): Induction Pack v7 with Alpha/Lambda Drift

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

- k=1 micro-lemma remains robust and closes all `k=1` updates.
- Broken bookkeeping route is disallowed.
- Alpha calibration now has drift by n:
  - n<=19: `alpha_min=0.2437206585`
  - n<=21: `alpha_*=0.21034113597068071`.
- Upper reserve side currently uses `lambda19=0.08144365672607116` (may drift at larger n).

## Task

Write a near-final induction framework that is stable under drift:

1. Use three branches:
   - `k=1` branch (micro-lemma),
   - `k>=2,t=2` branch (alpha-sandwich),
   - `k>=2,t>=3` bridge branch.
2. Parameterize theorem by `(alpha, lambda)` with explicit condition `alpha > lambda`.
3. Instantiate with current locked values (`alpha_*`, `lambda19`) to show present closure margin.
4. Provide one non-circular dependency graph and one unresolved bottleneck.

## Output format

1. `Parameterized induction theorem`
2. `Current-value instantiation (n<=21 alpha)`
3. `Dependency graph`
4. `Single unresolved bottleneck`
5. `Manuscript-ready default statement`
