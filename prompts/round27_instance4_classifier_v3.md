# Round 27 (Instance 4): Classifier v3 with Calibrated-Bookkeeping Signals

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

Existing two-score framework:
- `Gamma_G`: global reserve hardness,
- `Gamma_L`: local neighbour-compensation stress.

New signal from calibrated bookkeeping:
- minimum feasible alpha in `Lambda >= D + alpha*R_shift - sum_all` is about `0.2437206585` on n<=19 X<0.

## Task

Extend the classifier with one additional score that captures bookkeeping pressure.

Requirements:
1. Define explicit formula for a third score (`Gamma_B` or equivalent).
2. Explain how it interacts with `Gamma_G` and `Gamma_L`.
3. Use it to refine next-jump ranking beyond the prior `(3,15),(2,16),...` list.
4. Clearly mark proved algebraic parts vs conjectural bridge assumptions.

## Output format

1. `Three-score classifier definitions`
2. `Interpretation of each score`
3. `Refined next-jump ranking`
4. `Proved vs conjectural`
5. `Minimal validation pipeline`
