# Round 32 (Instance 1): Internal-Bucket Repair Theorem After n=23

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

- Same all-diagonal framework and boundary-correct prefix regime as Round 28–31.
- Drift-envelope definitions remain:
  - `alpha(n) := inf (Lambda + sum_all - D)/R_shift`
  - `lambda(n) := sup (sum_all - D)/R_shift`
- n=23 completed snapshot (mod=256):
  - `alpha_front(23) = 0.18243252946998884`
  - `lambda_front(23) = 0.24039868912946966`
  - gap = `alpha - lambda = -0.05796615965948082`
- Witnesses:
  - alpha witness: `(a,b)=(2,19)`, `step=2`, `k=9`
  - lambda witness: `(a,b)=(2,18)`, `step=2`, `k=5`
- So unlike n=22, both extremisers are in the same small-first line (`a=2`).

## Task

Design the next minimal theorem-level repair after the Round 31 two-bucket idea is no longer sufficient:

1. Keep the same base objects (`Lambda`, `D`, `R_shift`, `sum_all`) and all-diagonal transport setup.
2. Provide one repaired theorem schema that can close despite internal same-bucket break at n=23.
3. State the smallest finite parameter/interface needed to certify closure (explicit tuple, not prose).
4. Explain exactly why this repair addresses the n=23 mechanism (`(2,19)` alpha witness and `(2,18)` lambda witness in one bucket).
5. Provide a finite falsification checklist.

## Output format

1. `Internal-break repaired theorem candidate`
2. `Minimal finite parameterization`
3. `Why this fixes the n=23 same-bucket break`
4. `Single unresolved lemma interface`
5. `Finite falsification checklist`
