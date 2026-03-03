# Round 25 (Instance 4): Two-Score Classifier with Next-Jump Forecast

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

Class maxima on n<=19 X<0 corpus:
- (3,14): 0.08144365672607116
- (2,15): 0.07972646169405961
- (2,14): 0.05201381704686925
- (3,13): 0.04386927442810327
- (2,13): 0.023760967407659456
- (4,13): max_lambda = 0 in this regime.

Coarse geometry is saturated and non-discriminative.

## Task

Refine the two-score classifier (`Gamma_G`, `Gamma_L`) into a practical forecast tool.

Requirements:
1. Give explicit formulas in minors/channels.
2. Explain `(2,15),(3,14)` activation and `(4,13)` global-easy/local-stress split.
3. Produce a concrete ranked forecast list for next candidate classes at larger n.
4. Provide a minimal finite-check pipeline for validating the forecast.

## Output format

1. `Final two-score definitions`
2. `Activation mechanism`
3. `Ranked next-jump forecast`
4. `Minimal check pipeline`
5. `Proved vs conjectural`
