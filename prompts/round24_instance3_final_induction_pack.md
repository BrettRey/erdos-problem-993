# Round 24 (Instance 3): Draft a Near-Final Induction Pack

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

- Adjacent-minor bridge is dead.
- k=1 closure is strongly supported:
  - `Lambda_1 >= 0` checks: `1,407,504`
  - failures: `0`
- X<0 profile at `n<=19`:
  - total: `428,434`
  - step 2: `428,422`
  - step 3: `12` (all at `k=1`)
- For `k>=2`, reserve closure currently has two viable theorem forms:
  A) single-constant: `sum_err <= D + 0.08144365672607116 * R_shift`
  B) split-reserve: `sum_err <= D + lambda0*R_shift + c*Extra` (with explicit Extra from Round 24 Inst2)

## Task

Produce a near-final theorem package suitable for manuscript insertion, with explicit branches:

1. `k=1` branch via micro-lemma.
2. `k>=2, t=2` branch via upgraded reserve inequality.
3. `k>=2, t>=3` branch via bridge invariant.

Requirements:
1. Give two theorem variants (A and B above), and state which one you recommend.
2. Provide a dependency graph that is explicitly non-circular.
3. Isolate exactly one unresolved lemma per variant.
4. Keep notation consistent with `D`, `sum_err`, `R_shift`, `Extra`.

## Output format

1. `Variant A theorem (single-constant)`
2. `Variant B theorem (split-reserve)`
3. `Dependency graph`
4. `Unresolved lemma(s)`
5. `Recommendation and rationale`
