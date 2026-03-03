# Round 28 (Instance 4): Classifier v4 with Proof-Margin Score

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

Existing scores:
- `Gamma_G`: global reserve hardness (`(sum_err-D)_+ / R_shift`)
- `Gamma_L`: local neighbour-compensation stress

New calibrated signal:
- `alpha_min = 0.2437206585182262` from
  `Lambda >= D + alpha*R_shift - sum_all` (all-diagonal form).
- Current sharp global reserve level: `lambda19 = 0.08144365672607116`.

## Task

Extend classifier with a **proof-margin score** that measures distance to alpha-sandwich failure.

Requirements:
1. Define explicit formula for a new score (`Gamma_M` or equivalent), instance-level and class-level.
2. Explain interaction among `Gamma_G`, `Gamma_L`, and the new score.
3. Use it to refine next-jump ranking beyond prior `(3,15),(2,16),...` ordering.
4. Give a concrete “danger threshold” criterion for when a future lambda jump would break the alpha-sandwich route.
5. Separate proved algebraic statements vs conjectural bridge assumptions.

## Output format

1. `Three-score (or four-score) definitions`
2. `Interpretation and coupling`
3. `Refined next-jump ranking with danger thresholds`
4. `Proved vs conjectural`
5. `Minimal validation pipeline`
