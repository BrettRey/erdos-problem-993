# Round 29 Prompt Prep (2026-03-03)

Goal: move from fixed n<=19 alpha-sandwich to drift-aware theorem/proof plan after n<=21 calibration.

## New locked facts

From `results/alpha_bookkeeping_modal_n3_n21.json`:
- `min_alpha_all (n<=21) = 0.21034113597068071`
- `min_alpha_all (n<=19) = 0.2437206585182262`
- `min_alpha_odd (n<=21) = -2.088235294117647`
- new alpha witness: `n=21`, `step=2`, `k=8`, `(a,b)=(2,17)`.

Carry-forward constants:
- `lambda19 = 0.08144365672607116`
- current gap at n<=21-alpha level: `0.21034113597068071 - 0.08144365672607116 = 0.12889747924460955`.

## Prompt strategy

1. Keep alpha-sandwich architecture, but parameterize by `(alpha, lambda)`.
2. Push on lower-face mechanism (`G_k >= alpha R_shift`) via exact identity route.
3. Use drift-aware induction statement and classifier risk thresholds.
4. Continue forbidding broken bookkeeping (`Lambda >= D + R_shift - sum_err`).
