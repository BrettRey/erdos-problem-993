# One-Sided Effective-Drop Reduction
Date: 2026-07-03

## Purpose

This note isolates the one-sided version of the effective signed ratio-drop route. It does not prove the signed lemma. It records what Newton's inequalities already give and what extra localization would be needed for the `1/4` target.

## Setup

Let

```text
S = sum_i Bernoulli(p_i),
0 <= p_i <= 1/2,
a_k = P(S = k),
V = Var S.
```

Let

```text
D = min { k >= 1 : a_k < a_{k-1} }
```

and write

```text
rho_k = a_{k+1}/a_k.
```

The one-sided effective drop at first descent is

```text
Delta_eff = 1 - rho_D/rho_{D-1}.
```

## Newton Reduction

Write

```text
a_k = C e_k(w_1,...,w_n),
w_i = p_i/(1-p_i) <= 1.
```

Newton's inequalities give

```text
rho_D/rho_{D-1}
  <= D(n-D)/((D+1)(n-D+1))
  <= D/(D+1).
```

Therefore

```text
Delta_eff >= 1/(D+1).
```

Consequently, the desired nonsharp bound

```text
Delta_eff >= 1/(4V)
```

would follow from the mode-variance localization

```text
D + 1 <= 4V.
```

This is a sufficient reduction, not an equivalence.

## What Is Already Proved

The existing low-probability reserve note proves the weaker chain

```text
D+1 <= mu+3 <= 2V+3.
```

It follows immediately that

```text
Delta_eff >= 1/(5V)
```

for `V >= 1`, and that the `1/(4V)` target follows from this crude chain whenever `V >= 3/2`.

The remaining one-sided gap for this route is therefore the low-variance band

```text
1 <= V < 3/2.
```

In that band, the natural sharpened target is exactly

```text
D + 1 <= 4V.
```

## Probe

I added:

```bash
python3 scripts/probe_one_sided_effective_drop.py \
  --out results/one_sided_effective_drop_probe_2026-07-03.json
```

The probe scans random grouped low-probability laws, a binomial grid, and finite Poisson approximants. It checks:

```text
Delta_eff >= 1/(4V),
D + 1 <= 4V.
```

This is only a falsification probe. It is not a proof of either statement.

The default run processed `8,074` one-sided rows, including `6,652` with `V >= 1`, and found:

```text
quarter_effective_failures = 0
mode_variance_4_failures  = 0
```

The smallest observed effective-drop values by variance cutoff were:

| Variance cutoff | Source | `V` | First descent | `V * Delta_eff` | `(D+1)/V` |
|---:|---|---:|---:|---:|---:|
| 1 | finite Poisson, `lambda=1.01` | `1.009798` | `2` | `0.336733` | `2.970891` |
| 1.25 | random grouped | `1.290047` | `2` | `0.434337` | `2.325497` |
| 1.5 | finite Poisson, `lambda=2` | `1.999600` | `3` | `0.500050` | `2.000400` |
| 2 | random grouped | `2.019659` | `3` | `0.512845` | `1.980533` |
| 5 | binomial grid | `5.046597` | `6` | `0.743125` | `1.387073` |

The largest observed values of `(D+1)/V` were instead near small high-probability binomial rows:

| Source | `V` | First descent | `(D+1)/V` | `V * Delta_eff` |
|---|---:|---:|---:|---:|
| `Bin(5,1/2)` | `1.25` | `4` | `4.000000` | `0.750000` |
| random grouped, `n=5` | `1.000105` | `3` | `3.999580` | `0.624193` |

Thus the two empirical boundaries differ:

1. The effective-drop constant is stressed by sparse Poisson-like laws near `Pois(1)`, where the limiting value is `1/3`.
2. The sufficient localization `D+1 <= 4V` is stressed by small nearly-fair binomial laws.

This distinction matters because proving `D+1 <= 4V` may be easier or harder than proving the effective-drop bound directly; it is only one sufficient route.

## Proof Target

A useful next lemma would be:

> If `0 <= p_i <= 1/2`, `V >= 1`, and `D` is the first strict descent, then `D+1 <= 4V`.

Together with Newton, this would prove the one-sided effective-drop bound at the `1/4` scale. It would still not prove the signed case, but it would remove one uncertainty from the route.
