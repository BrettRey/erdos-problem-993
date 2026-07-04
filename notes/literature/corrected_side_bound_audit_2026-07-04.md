# Corrected Side-Bound Audit
Date: 2026-07-04

## Purpose

This follows the corrected signed conditional reduction after the missing
X-side `beta_X` boundary term was found. The goal was to make the corrected
X-side and reflected Y-side lower bounds computable and then stress-test the
naive side-selection target

```text
V * max(X_side_bound, Y_side_bound) >= 1/4.
```

This is a diagnostic for issue #5. It is not a proof of the signed reserve
lemma.

## Implementation

`scripts/analyze_signed_conditionals.py` now computes:

```text
x_reduction_bound
y_reduction_bound
best_side_reduction_bound = max(x_reduction_bound, y_reduction_bound)
```

using the corrected X-side reciprocal boundary

```text
beta_X = a_{n_X} b_{n_X+1-D}/c_D.
```

It also checks the reciprocal identities

```text
c_{D-1}/c_D = E_pi[h_X] + beta_X,
c_{D-1}/c_D = E_pi[k_Y] + beta_Y,
```

and the product identities for `R_+/R_-`.

A regression test was added for the exact audit counterexample

```text
X = Bernoulli(1/2),
Y = Binomial(10,1/2).
```

The test verifies:

```text
D = -3,
Delta = 3/10,
old omitted-boundary expression = 4/11,
beta_X = 42/55,
corrected X-side bound = -1/55.
```

## Commands

Corrected top-row audit:

```bash
python3 scripts/analyze_signed_conditionals.py \
  results/signed_ratio_drop_breaker_extended_2026-07-04.json \
  --top 25 \
  --out results/signed_conditional_reduction_audit_2026-07-04.json
```

Smoke breaker with `side_reduction` as an additional optimizer target:

```bash
python3 scripts/signed_ratio_drop_breaker.py \
  --random-samples 100 \
  --max-side-n 400 \
  --de-maxiter 4 \
  --de-popsize 3 \
  --de-restarts 1 \
  --top 10 \
  --out results/signed_ratio_drop_breaker_smoke_2026-07-04.json
```

Joint Newton-or-reserve smoke breaker:

```bash
python3 scripts/signed_ratio_drop_breaker.py \
  --random-samples 100 \
  --max-side-n 400 \
  --de-maxiter 4 \
  --de-popsize 3 \
  --de-restarts 1 \
  --fallback-constant 0.75 \
  --top 10 \
  --out results/signed_ratio_drop_breaker_joint_smoke_2026-07-04.json
```

Half-heavy+dust grid after adding the explicit family:

```bash
python3 scripts/signed_ratio_drop_breaker.py \
  --random-samples 0 \
  --max-side-n 400 \
  --de-restarts 0 \
  --fallback-constant 0.75 \
  --top 10 \
  --out results/signed_ratio_drop_breaker_joint_grid_2026-07-04.json
```

Joint optimizer smoke:

```bash
python3 scripts/signed_ratio_drop_breaker.py \
  --random-samples 0 \
  --max-side-n 400 \
  --de-maxiter 2 \
  --de-popsize 2 \
  --de-restarts 1 \
  --fallback-constant 0.75 \
  --top 10 \
  --out results/signed_ratio_drop_breaker_joint_optimizer_smoke_2026-07-04.json
```

Tests:

```bash
python3 -m unittest test_all.py -v
```

## Results

The corrected top-row audit processed `68` unique rows and found maximum
identity error

```text
6.66e-16.
```

On those saved near-sharp rows, the smallest observed corrected side bound was

```text
V * best_side_reduction_bound ~= 0.333876.
```

The smoke breaker processed:

```text
attempted = 38039
analyzed = 36485
feasible with V >= 1 = 35689
effective-drop failures below 1/4 = 0
raw-reserve failures below 1/4 = 0
rows with V * best_side_reduction_bound below 1/4 = 146
```

The joint smoke breaker used the same small run with fallback threshold `0.75`
and found:

```text
rows with V * best_side_reduction_bound below 1/4 = 146
rows also having V * reserve below 0.75 = 0
rows also having V * Delta_eff below 0.75 = 0
```

After adding the explicit half-heavy+dust grid, the no-optimizer grid run
processed:

```text
attempted = 63607
analyzed = 62266
feasible with V >= 1 = 61470
rows with V * best_side_reduction_bound below 1/4 = 10591
rows also having V * reserve below 0.75 = 0
rows also having V * Delta_eff below 0.75 = 0
```

The joint optimizer smoke processed:

```text
attempted = 64375
analyzed = 62746
feasible with V >= 1 = 61950
rows with V * best_side_reduction_bound below 1/4 = 10623
rows also having V * reserve below 0.75 = 0
rows also having V * Delta_eff below 0.75 = 0
```

The best effective-drop row remained near the sparse one-sided boundary:

```text
V * Delta_eff ~= 0.333877,
V * best_side_reduction_bound ~= 0.333876.
```

The low side-bound rows are not reserve threats. In the smoke output, among
rows with `V * best_side_reduction_bound < 1/4`, the smallest stored actual
slack was

```text
min V * reserve ~= 0.827112,
min V * Delta_eff ~= 0.799800.
```

The low side-bound rows are strongly half-heavy. In the joint smoke output,
the lowest-reserve low-side-bound examples had half-heavy variance fraction
between about `0.86` and `0.97`; the smallest-side-bound optimizer examples
had half-heavy variance fraction essentially `1`.

The explicit half-heavy+dust grid found the closest fallback rows so far:

```text
X = Bernoulli(0.25) + Binomial(6, 1/2),
Y = Binomial(4, 1/2),
V * best_side_reduction_bound ~= 0.234159,
V * Delta_eff ~= 0.785269,
V * reserve ~= 0.950680.
```

and

```text
X = Binomial(5, 0.004) + Binomial(7, 1/2),
Y = Bernoulli(0.02) + Binomial(4, 1/2),
V * best_side_reduction_bound ~= 0.243064,
V * Delta_eff ~= 0.789679,
V * reserve ~= 0.789747.
```

Thus `0.75` remains a plausible fallback search threshold in these runs, but
`0.8` is already too aggressive as a conjectural universal fallback constant.

The worst side-bound row was the balanced half-heavy case

```text
X = Binomial(500, 1/2),
Y = Binomial(500, 1/2),
V = 250,
D = 1,
V * best_side_reduction_bound ~= 0.005962,
V * Delta_eff ~= 0.997012,
V * reserve ~= 1.49402.
```

## Interpretation

The naive side-selection lemma

```text
V * max(X_side_bound, Y_side_bound) >= 1/4
```

is false as a proof target. The corrected conditional Newton bound is too
weak in high-variance, half-heavy, normal-like cases. This does not threaten
the signed reserve route because those rows already have large raw and
effective reserve.

The better target is a two-clause lemma:

> At the signed first descent, either the corrected side-bound gives
> `V * best_side_reduction_bound >= c`, or the raw/effective reserve is already
> at least `c'/V` by a direct smoothing or half-heavy argument.

Empirically, the side-bound clause carries the near-sharp sparse boundary, and
the fallback clause carries the half-heavy balanced regime.

The most concrete fallback target suggested by this audit is:

> If at least a fixed fraction of `V` comes from Bernoulli parameters near
> `1/2`, and the first descent is not carried by the corrected conditional
> side-bound, then `V * Delta_eff` is bounded below by an absolute constant
> well above `1/4`.

The `0.75` value survived the smoke/grid/optimizer passes above, but it should
be treated only as a search calibration, not a proposed theorem constant.

The exact fair-binomial model is proved separately in

```text
notes/literature/fair_binomial_signed_fallback_2026-07-04.md
```

It gives the theorem-level base case

```text
V * Delta_eff >= 5/8,
V * reserve >= 3/4
```

for `X~Binomial(m,1/2)` and `Y~Binomial(n,1/2)`. Both constants are sharp at
total fair count `m+n=4`; constants near `0.8` are disfavored by the
half-heavy+dust rows above and cannot be inferred from the exact model.

## Overclaim Guard

This audit does not prove the signed reserve theorem, the effective-drop
diagnostic, the hub-bouquet reserve, or Erdos 993. It also disproves one
over-strong proof target: the corrected conditional side-bound alone cannot be
the whole argument.
