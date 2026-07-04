# Local-Mode Mean Breaker Pass
Date: 2026-07-04

## Target

This pass tried to falsify the local-mode mean lemma:

```text
0 <= w_i <= 1, e_2(w) > 0, e_3(w) >= e_2(w)
    =>  mu = sum_i w_i/(1+w_i) >= 5/2.
```

This note records computational pressure only. It is not a proof of the lemma.

## Harness

The new breaker harness is:

```bash
python3 scripts/local_mode_mean_breaker.py
```

It combines four search modes:

1. grouped constrained SLSQP over integer multiplicity patterns;
2. ungrouped constrained SLSQP restarts;
3. structured families around five weights equal to one;
4. coarse two-group grids and grouped differential-evolution penalty runs.

The grouped pass is intended to put pressure on KKT-style finite-support extremals, where an optimizer would be expected to have few distinct nonzero weights if the constraint is active.

## Default Run

Command:

```bash
python3 scripts/local_mode_mean_breaker.py \
  --out results/local_mode_mean_breaker_2026-07-04.json
```

Result:

```text
processed = 57904
feasible = 37070
counterexamples = 0
near_counterexamples_strict_tolerance = 0
best_mean = 2.4999999999892917
```

The best rows are numerical copies of the expected boundary:

```text
w_1 = ... = w_5 = 1, with the remaining weights zero.
```

Tiny negative values of `e_3-e_2` in some best rows are floating-point SLSQP tolerance effects, not strict feasible violations.

## Extended Run

Command:

```bash
python3 scripts/local_mode_mean_breaker.py \
  --max-support 22 --max-groups 5 --group-restarts 6 \
  --ungrouped-max-n 14 --ungrouped-restarts 6 \
  --de-restarts 2 --de-maxiter 180 \
  --out results/local_mode_mean_breaker_extended_2026-07-04.json
```

Result:

```text
processed = 117395
feasible = 89009
counterexamples = 0
near_counterexamples_strict_tolerance = 0
best_mean = 2.4999999999041704
```

Again, the best row is the five-one boundary plus zeros, up to numerical tolerance.

## Status

This breaker pass did not find a counterexample to the local-mode mean lemma. It strengthens the evidence that the boundary is five fair weights, but it does not prove the lemma. The next proof-side target remains an extremal argument showing that any minimizer under `e_3 >= e_2` can be reduced to the five-one boundary or a higher-mean case.
