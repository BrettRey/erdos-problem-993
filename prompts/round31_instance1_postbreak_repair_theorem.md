# Round 31 (Instance 1): Post-Break Repair Theorem (n=22 baseline)

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

- Drift envelope definitions are fixed:
  - `alpha(n) := inf (Lambda + sum_all - D)/R_shift`
  - `lambda(n) := sup (sum_all - D)/R_shift`
- n=22 modal snapshot (mod=256 complete):
  - `alpha_front(22) = 0.1875868239504603`
  - `lambda_front(22) = 0.1968360500404575`
  - gap = `alpha - lambda = -0.0092492260899972` (global sandwich break)
- Witnesses:
  - alpha witness at `(a,b)=(2,18), step=2, k=9`
  - lambda witness at `(a,b)=(4,16), step=2, k=5`

## Task

Design the smallest theorem-level repair that can restore closure despite `alpha(22) < lambda(22)`:

1. Keep all-diagonal framework (`sum_all`) and boundary-correct prefix setup.
2. Provide one repaired theorem schema that is strictly weaker than the failed global envelope but still strong enough for closure.
3. State a minimal constant/interface to certify it (single scalar or smallest finite tuple).
4. Explain why this repair avoids the exact n=22 break mechanism.
5. Give a falsification checklist for the repaired schema.

## Output format

1. `Post-break repaired theorem candidate`
2. `Minimal repair parameterization`
3. `Why this defeats the n=22 break witnesses`
4. `Single unresolved lemma interface`
5. `Finite falsification checklist`

