# Round 21 (Instance 4): Refine Classifier (Current Hard Criterion is Too Coarse)

You proposed:
`Hard(instance) := exists k<prefix with W(k)<=3 and boundary touch`.

Use only your prior output + this prompt.

## Locked falsification of coarse classifier

On all negative step-2 cases through `n<=18` (all 9 pair classes), evaluated at the actual negative `k`:

- `W(k)<=3` is true in 100% of cases
- boundary touch is true in 100% of cases

So this binary criterion does **not** separate nonzero-lambda classes from zero-lambda classes.

## Locked target structure

Nonzero lambda classes: `(2,14),(3,13),(2,13)`.
Other six classes have lambda 0 under global reserve fit.

## Task

Upgrade to a quantitative classifier score `H_score(instance,k)` that can separate:

- high-reserve-needed classes,
- zero-reserve classes.

Constraints:

1. Must be explicit in coefficients / local ratios / minor magnitudes (not black-box).
2. Must reduce to finite inequality families for `a=2` and `a=3`.
3. Must explain why binary `W<=3 + boundary` saturates but magnitude score still separates.

Suggested direction: normalized obstruction magnitude, e.g.

`Gamma = (sum_err - D)_+ / (C10+C01+C11)`

or a local surrogate built from your `Omega(k)` expansion coefficients.

## Output format

1. `Quantitative classifier definition`
2. `Analytic link to reserve constant`
3. `Finite inequality family (a=2,a=3)`
4. `What needs proof vs what is empirical`

No binary classifier claims.
