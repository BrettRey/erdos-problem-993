# Poisson Boundary For The Effective Signed Ratio Drop
Date: 2026-07-03

## Purpose

The signed ratio-drop optimizer found rows with

```text
V * Delta_eff = 0.336426...
```

where

```text
Delta_eff = 1 - (c_{D+1}/c_D)/(c_D/c_{D-1}).
```

This note explains that value. It is not a new theorem about all Poisson-binomial laws. It is an exact calculation in the one-sided Poisson limit, used to calibrate the possible constant for this sufficient proof route.

## One-Sided Poisson Calculation

Let

```text
X ~ Pois(lambda),
p_k = P(X = k).
```

Then

```text
p_k / p_{k-1} = lambda/k.
```

The first strict descent index is

```text
D = floor(lambda) + 1
```

for all `lambda > 0`; at integer `lambda = m`, the plateau at `m-1,m` means the first strict descent is `D = m+1`.

At that `D`,

```text
p_D/p_{D-1} = lambda/D,
p_{D+1}/p_D = lambda/(D+1).
```

Therefore

```text
Delta_eff
  = 1 - (lambda/(D+1))/(lambda/D)
  = 1/(D+1).
```

Since `Var X = lambda`,

```text
V * Delta_eff = lambda/(D+1).
```

For `V = lambda >= 1`, the minimum over the Poisson family is attained at the left endpoint `lambda = 1`, where `D = 2`, giving

```text
V * Delta_eff = 1/3.
```

## Interpretation

This explains the optimizer row near

```text
V = 1,
D = 2,
c_D/c_{D-1} ~= 1/2,
c_{D+1}/c_D ~= 1/3.
```

The finite optimizer value `0.336426...` is slightly above `1/3`, consistent with a finite binomial approximation to `Pois(1)` plus a very small reflected component.

This gives only a limiting calibration for the effective-drop route:

1. A universal effective-drop theorem with constant greater than `1/3` faces this limiting obstruction. Turning that into a finite Poisson-binomial upper-bound proposition still requires writing the standard approximation step.
2. The working constant `1/4` remains below this boundary.
3. This boundary concerns the sufficient diagnostic `Delta_eff`, not the raw reserve. The raw reserve has a different sparse boundary, discussed in `notes/literature/skellam_sparse_limit_reserve_2026-07-03.md`.

## Consequence For Issue #5

The effective-drop route is still plausible at the nonsharp `1/4` scale, but it is now clear why it has less slack than the raw-reserve route. The proof should aim for

```text
Delta_eff >= 1/(4V)
```

or another comparably weak constant. It should not be framed as evidence for any constant above `1/3`.
