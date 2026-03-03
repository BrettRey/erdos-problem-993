# Round 24 (Instance 4): Quantitative Classifier Update After n=19

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

- Through `n<=18`: hard classes were `(2,14),(3,13),(2,13)`.
- At `n=19`: new hard classes `(2,15),(3,14)` appear.
- Pair-class maxima on full `n<=19` X<0 corpus:
  - `(3,14): max_lambda = 0.08144365672607116`
  - `(2,15): max_lambda = 0.07972646169405961`
  - `(2,14): 0.05201381704686925`
  - `(3,13): 0.04386927442810327`
  - `(2,13): 0.023760967407659456`
- `(4,13)` has X<0 cases but `max_lambda = 0` in this regime.
- Coarse geometry predicate is saturated and non-discriminative.

## Task

Upgrade the classifier into a **two-score** system that separates:

1. Global reserve hardness (`lambda` requirement), and
2. Local neighbor-compensation stress.

Requirements:
1. Scores must be explicit formulas in local coefficients/minors/channels.
2. Explain how `(2,15),(3,14)` become hard while `(4,13)` stays globally easy.
3. Provide a predictive criterion for the *next* class likely to activate at larger n.
4. Clearly label what is proved algebraically vs conjectural bridge assumptions.

## Output format

1. `Two-score classifier definitions`
2. `Why new hard classes activate`
3. `Why (4,13) can remain lambda=0`
4. `Next-jump predictive criterion`
5. `Proved vs conjectural`
