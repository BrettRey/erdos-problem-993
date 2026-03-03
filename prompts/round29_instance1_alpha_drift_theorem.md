# Round 29 (Instance 1): Alpha-Sandwich with Drift (n<=21)

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

Conventions: support-root, boundary-correct indexing, prefix `k < mode(I_new)`.

From the Modal artifact `results/alpha_bookkeeping_modal_n3_n21.json`:
- `min_alpha_all(n<=19) = 0.2437206585182262`
- `min_alpha_all(n<=21) = 0.21034113597068071`
- odd-only analogue remains negative (`min_alpha_odd(n<=21) = -2.088235294117647`)
- witness for new alpha minimum at `n=21`, `(a,b)=(2,17)`, `step=2`, `k=8`.

Previously locked reserve constant on n<=19:
- `lambda19 = 0.08144365672607116`.

## Task

Produce a drift-aware theorem statement:

1. Replace fixed `alpha_min` by an explicit `alpha_*` that is valid through `n<=21`.
2. State closure in a way that is robust to future lambda drift (parameterized threshold form).
3. Quantify the current certified margin using `alpha_* = 0.21034113597068071`.
4. Give one explicit failure criterion for the alpha-sandwich route.

Do not reuse invalid bookkeeping (`Lambda >= D + R_shift - sum_err`).

## Output format

1. `Drift-aware alpha-sandwich theorem candidate`
2. `Current certified margin and closure implication`
3. `Failure criterion / danger threshold`
4. `Single unresolved lemma (or minimal pair)`
5. `Finite falsification checklist`
