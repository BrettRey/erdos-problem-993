# Round 28 (Instance 3): Induction Pack v6 (Alpha-Sandwich, No BK)

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

- `k=1` micro-lemma is robust and closes all `k=1` updates.
- Non-step2 X<0 at `n<=19` are only step3, `k=1` (already handled by micro-lemma).
- Broken bookkeeping form is forbidden.
- Calibrated constants:
  - `alpha_min = 0.2437206585182262`
  - `lambda19 = 0.08144365672607116`
  - margin `alpha_min - lambda19 > 0`.

## Task

Write a near-final induction framework with 3 branches:
1. `k=1` branch via micro-lemma,
2. `k>=2, t=2` branch via alpha-sandwich closure,
3. `k>=2, t>=3` bridge branch via same alpha-sandwich (or a strictly weaker invariant that still closes).

Requirements:
- Do not assume `Lambda >= D + R_shift - sum_err`.
- Express closure entirely via inequalities of the form
  `Lambda - D + sum_all >= alpha*R_shift`
  and
  `sum_all <= D + lambda*R_shift (+ optional tiny extras if needed)`.
- Provide non-circular dependency graph.
- State exactly one unresolved bottleneck per variant (if two variants are given).

## Output format

1. `Variant A induction theorem (single reserve)`
2. `Variant B induction theorem (if needed)`
3. `Dependency graph`
4. `Unresolved bottleneck(s)`
5. `Recommended manuscript default`
