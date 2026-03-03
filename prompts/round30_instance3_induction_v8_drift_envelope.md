# Round 30 (Instance 3): Induction Pack v8 (Drift Envelope + Bottleneck Split)

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

- k=1 micro-lemma remains robust and independent.
- For k>=2, closure is now routed through alpha/lambda sandwich faces (not broken bookkeeping).
- Current calibration point:
  - `alpha_* = 0.21034113597068071`
  - `lambda19 = 0.08144365672607116`
  - positive margin currently remains.
- Desired architecture: stable if `alpha` and `lambda` drift with n.

## Task

Produce a near-final induction theorem (v8) with a clean unresolved interface:

1. Keep three branches:
   - `k=1` micro-lemma branch,
   - `k>=2, t=2` hard-step branch,
   - `k>=2, t>=3` bridge branch.
2. Parameterize by `(alpha(n), lambda(n))` with explicit condition for closure.
3. Separate unresolved pieces into:
   - lower-face bottleneck,
   - upper-face bottleneck,
   and show no circularity.
4. Give a manuscript-ready statement and a short proof sketch.
5. Include an explicit “if drift crosses threshold, what fails first?” clause.

## Output format

1. `Induction theorem v8 (drift-parameterized)`
2. `Three-branch closure proof sketch`
3. `Dependency graph (non-circular)`
4. `Two-bottleneck unresolved interface`
5. `Failure mode when drift threshold is crossed`

