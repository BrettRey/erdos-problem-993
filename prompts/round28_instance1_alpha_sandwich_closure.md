# Round 28 (Instance 1): Alpha-Sandwich Closure (No Broken Bookkeeping)

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

Conventions: support-root, boundary-correct indexing, prefix `k < mode(I_new)`.

1. The old bookkeeping form
`Lambda_k >= D + R_shift - sum_err`
is invalid on `n<=19` X<0 corpus.
- all-diagonal version (`sum_all = sum_s err_s`): 9,223 failures
- odd-only version (`sum_odd`): 428,434 failures

2. Calibrated lower-bound scan on the same `n<=19` X<0 corpus gives:
- minimal feasible `alpha` in
  `Lambda_k >= D + alpha*R_shift - sum_all`
  is
  `alpha_min = 0.2437206585182262`.
- odd-only linear analogue is impossible with nonnegative alpha.

3. Locked global reserve constant from prior rounds:
`lambda19 = 0.08144365672607116`.

Hence numerical margin:
`alpha_min - lambda19 = 0.162276... > 0`.

## Task

Build a theorem candidate that closes via an **alpha-sandwich**:
- lower side uses calibrated-alpha inequality with `sum_all`,
- upper side uses global `sum_all` reserve inequality,
- closure is obtained from positive margin `alpha - lambda`.

Do **not** use the invalid bookkeeping form.

## Output format

1. `Alpha-sandwich theorem candidate`
2. `Derivation of closure from alpha-lambda margin`
3. `Why odd-only cannot replace sum_all here`
4. `Single unresolved lemma (or minimal pair)`
5. `Finite falsification checklist`
