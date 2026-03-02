# Round 22 (Instance 4): Quantitative Classifier via Compensation Load

Use only your prior outputs plus this prompt.

## Locked context

Binary classifier `W<=3 + boundary-touch` is saturated (true for all negative cases in all 9 classes), so it does not separate hard/easy.

New quantitative diagnostics suggest using compensation load:

- per case: needed global lambda
  `Gamma_case := ((sum_err - D)_+) / (C10+C01+C11)`
- per odd slot: local deficit/slack
  `d_t`, `e_t`, `e_t+e_{t+1}`

Observed through `n<=18`:

- hard classes with nonzero lambda: `(2,14),(3,13),(2,13)`
- zero-lambda classes: the other six observed classes
- `d_t<=e_t` sometimes fails, but `d_t<=e_t+e_{t+1}` never fails

## Task

Construct an explicit quantitative classifier score that predicts reserve hardness and is compatible with the above.

Requirements:

1. Define score(s) explicitly from coefficient-level quantities.
2. Link score analytically to required reserve constant (global or class-wise).
3. Provide finite inequality family for `a=2` and `a=3`.
4. Explain why binary overlap criterion saturates but your score still separates classes.

## Output format

1. `Score definition(s)`
2. `Analytic connection to lambda`
3. `Finite family (a=2,a=3)`
4. `Proved vs conjectural parts`

No binary-only classifier claims.
