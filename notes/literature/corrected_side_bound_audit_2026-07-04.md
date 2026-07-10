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

A hardened rerun after adding relative-drop descent tolerance and non-finite
diagnostic counters gave the same top-row result:

```text
attempted = 71
duplicates = 3
analyzed = 68
terminal descents = 0
non-finite expectation terms = 0
non-finite identity errors = 0
max identity error = 6.661e-16
```

This identity error certifies the displayed identities at the computed descent
index. It does not, by itself, certify descent placement. Non-dyadic
near-plateau rows should be rechecked with a relative-drop tolerance or exact
arithmetic before their descent index is treated as a certificate.

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

These are finite smoke runs over the listed families and parameter ranges; the
zero failure counts are not universal statements. Rows dropped between
`attempted` and `analyzed` include terminal descents with `Delta=1`, which
trivially satisfy the reserve targets but are not useful for minimum searches.

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

In these runs, the low side-bound rows were not reserve threats. In the smoke
output, among rows with `V * best_side_reduction_bound < 1/4`, the smallest
stored actual slack was

```text
min V * reserve ~= 0.827112,
min V * Delta_eff ~= 0.799800.
```

In these runs, the low side-bound rows were strongly half-heavy. In the joint
smoke output, the lowest-reserve low-side-bound examples had half-heavy
variance fraction between about `0.86` and `0.97`; the smallest-side-bound
optimizer examples had half-heavy variance fraction essentially `1`.

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
`0.8` is already too aggressive as a conjectural universal fallback constant;
the second displayed half-heavy+dust row is an exact rational counterexample.

## Exact Side-Bound Disproof

The smaller displayed half-heavy row also gives a compact exact disproof of
the side-bound-only target. Take

```text
X = Bernoulli(1/4) + Binomial(6,1/2),
Y = Binomial(4,1/2).
```

Exact rational evaluation of the corrected analyzer definitions gives

```text
D = 2,
V = 43/16,
X_side_bound = 399961/4594590,
Y_side_bound = 8069/92610,
V * max(X_side_bound,Y_side_bound)
  = 346967/1481760
  < 1/4.
```

The exact gap below `1/4` is

```text
23473/1481760.
```

The row is not a reserve threat:

```text
V * Delta_eff = 19393/24696,
V * reserve   = 559/588.
```

This certificate removes the former dependence on a large floating-point row
for the theorem-level retirement of the side-bound-only target.

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

The `Binomial(500,1/2)` row has also been exact-certified. With
`C=binom(1000,501)`, its two corrected side bounds are

```text
B_X = 1001/(501*502^2) + 1/(502*C),
B_Y = 1001*1504/(500*501*502^2) - 1/(501*502*C),
```

with `B_Y>B_X`. Hence

```text
250 * max(B_X,B_Y)
  = 188188/31563501 - 250/(501*502*C)
  < 1/4.
```

The saved float differs from the exact value by less than `4.8e-14`. This
large row explains how the conditional side bound can collapse in a balanced
normal-like regime, while the small rational row above is the preferred
disproof certificate.

## Interpretation

The naive side-selection lemma

```text
V * max(X_side_bound, Y_side_bound) >= 1/4
```

should be retired as a proof target. The corrected conditional Newton bound is
too weak in the high-variance, half-heavy, normal-like cases found here; the
`Binomial(500,1/2)` pair shows this weakness can be severe. This does not
threaten the signed reserve route because those rows already have large raw
and effective reserve.

The better target is a two-clause lemma:

> At the signed first descent, either the corrected side-bound gives
> `V * best_side_reduction_bound >= c`, or the raw/effective reserve is already
> at least `c'/V` by a direct smoothing or half-heavy argument.

Empirically in these runs, the side-bound clause carries the near-sharp sparse
boundary, and the fallback clause carries the half-heavy balanced regime.

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

The source note now records the corrected constants from the 2026-07-04 audit.
It gives the theorem-level base case

```text
V * Delta_eff >= 5/8,
V * reserve >= 3/4
```

for `X~Binomial(m,1/2)` and `Y~Binomial(n,1/2)` with
`V=(m+n)/4 >= 1`. Both constants are sharp at total fair count `m+n=4`.
Constants near `0.8` are ruled out for the displayed half-heavy+dust
perturbation by an exact rational regression test, and cannot be inferred from
the exact fair-binomial model.

## Overclaim Guard

This audit does not prove the signed reserve theorem, the effective-drop
diagnostic, the hub-bouquet reserve, or Erdos 993. It also disproves one
over-strong proof target: the corrected conditional side-bound alone cannot be
the whole argument. Issue #5 remains open.
