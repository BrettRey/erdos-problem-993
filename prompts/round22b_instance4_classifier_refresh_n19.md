# Round 22b (Instance 4): Refresh Quantitative Classifier with n=19 Data

Use only your prior outputs plus this prompt.

## Locked context update

Binary geometry is still saturated and non-discriminative.

At `n=19`, global reserve hardness changed:
- fixed `lambda0` fails in `33` cases
- new max needed global lambda is `0.08144365672607116`
- new hard classes appear (nonzero global lambda):
  - `(2,15)`, `(3,14)`, plus prior `(2,14)`, `(3,13)`, `(2,13)`
- class `(4,13)` shows local compensation stress (neighbor failures) but zero global lambda in current scan.

## Task

Update the quantitative classifier so it predicts both:
1. global reserve hardness (needed lambda), and
2. local compensation stress (neighbor-law defect risk).

You must define two explicit scores, e.g.:
- `Gamma_global := ((sum_err-D)_+) / R_shift`
- `Gamma_local := d_t / (e_t + e_{t+1})` or an analytic surrogate.

Requirements:
- explicit coefficient/minor formulas,
- clear link from each score to theorem terms,
- finite inequality families for `a=2,3,4`.

## Output format

1. `Two-score classifier definition`
2. `Analytic link to reserve/defect terms`
3. `Finite inequality families (a=2,3,4)`
4. `Proved vs conjectural parts`
5. `Why this resolves the n<=18 to n=19 shift`
