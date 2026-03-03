# Round 27 (Instance 1): Repair Bookkeeping with a Calibrated Alpha

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

Conventions: support-root, boundary-correct indexing, prefix `k < mode(I_new)`.

The previously used bookkeeping inequality

`Lambda_k >= D + R_shift - sum_err`

is empirically false on `n<=19` X<0 cases.

New calibration scan on full `n<=19` X<0 corpus (`428,434` cases):
- For all-diagonal error (`sum_all = sum_s err_s`), the minimal feasible coefficient in
  `Lambda_k >= D + alpha*R_shift - sum_all`
  is
  `alpha_min = 0.2437206585182262`.
- Witness:
  `n=19`, `g6=R?????????????????C??o??{?DF~?`, `root=0`, `step=2`, `k=7`, `(a,b)=(2,15)`.
- For odd-only error (`sum_odd = sum_{s odd} err_s`), the same form is impossible with nonnegative alpha:
  minimum observed effective alpha is `-1.6`.

## Task

Treat this as a theorem-design problem:

1. Derive a corrected bookkeeping form with explicit calibrated coefficients.
2. Explain why odd-only replacement cannot work in the same linear form.
3. Give one clean theorem candidate using calibrated alpha (and optional extra channels).
4. Isolate one unresolved proof lemma.

## Output format

1. `Calibrated bookkeeping theorem candidate`
2. `Why odd-only linear form fails`
3. `Implication to closure inequalities`
4. `Single unresolved lemma`
5. `Finite falsification checklist`
