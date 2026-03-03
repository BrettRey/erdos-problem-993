# Round 26 (Instance 4): Constants, Tradeoffs, and Next-Jump Forecast v2

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

Pair-class maxima on n<=19 X<0 corpus:
- (3,14): 0.08144365672607116
- (2,15): 0.07972646169405961
- (2,14): 0.05201381704686925
- (3,13): 0.04386927442810327
- (2,13): 0.023760967407659456

Tradeoff artifact exists (`lambda` vs `c` for mixed reserves).

## Task

Refine the quantitative picture into a decision table for manuscript strategy:

1. pure single-constant theorem,
2. split theorem with moderate Extra coefficient,
3. defect+tiny-residual theorem.

Requirements:
1. Provide explicit constants and interpretability tradeoffs.
2. Rank options by (a) simplicity, (b) plausibility of proof, (c) empirical sharpness.
3. Give a concrete next-jump forecast rule for classes at larger n.
4. State exactly which claims are algebraic vs conjectural.

## Output format

1. `Decision table`
2. `Recommended default theorem form`
3. `Forecast rule for next hard class`
4. `Algebraic vs conjectural split`
5. `Minimal validation plan`
