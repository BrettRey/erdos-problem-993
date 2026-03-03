# Round 28 Prompt Prep (2026-03-03)

Goal: steer next prompts away from broken bookkeeping and toward alpha-sandwich closure.

## Hard constraints to carry forward

- Broken bookkeeping is invalid on n<=19 X<0 corpus:
  - `Lambda >= D + R_shift - sum_all`: 9,223 failures
  - `Lambda >= D + R_shift - sum_odd`: 428,434 failures

- Calibrated lower bound (same corpus):
  - minimal feasible `alpha` in
    `Lambda >= D + alpha*R_shift - sum_all`
    is `0.2437206585182262`.

- Odd-only linear alpha route is impossible (`min effective alpha = -1.6`).

- Current global reserve constant from prior rounds:
  - `lambda19 = 0.08144365672607116`
  - positive margin `alpha_min - lambda19 ≈ 0.16228`.

- Split-domination obstruction remains:
  - global side condition `(1-lambda0)R_shift >= c*Extra` fails due cases with `R_shift=0`, `Extra>0`.

## Prompt direction

- Use all-diagonal `sum_all` for calibrated lower/upper sandwich.
- Ban reuse of broken bookkeeping form.
- Keep k=1 micro-lemma branch unchanged.
- Update classifier with a proof-margin score against alpha-sandwich failure.
