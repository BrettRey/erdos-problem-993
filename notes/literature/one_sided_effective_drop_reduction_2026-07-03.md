# One-Sided Effective-Drop Reduction
Date: 2026-07-03

**2026-07-10 status update.** The one-sided theorem and its plateau-safe weak
endpoint remain valid. A later Hillion--Johnson cubic-curvature proof now
establishes the general finite signed theorem in
`notes/literature/universal_poisson_binomial_effective_drop_2026-07-10.md`.

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

### Support domain

Delete every zero parameter and let `m` be the number of positive parameters.
The variance assumption `V>=1` and the bound
`p_i(1-p_i)<=1/4` imply `m>=4`.

Let `M` be the largest mode. If `m=4`, then `V>=1` forces all four
parameters to equal `1/2`, so `M=2`. If `m>=5`, Darroch localization and
`p_i<=1/2` give

```text
M < mu+1 <= m/2+1,
```

which implies `M<=m-2` after integer rounding. Since a Poisson-binomial mass
sequence is log-concave, `D=M+1<=m-1`. Therefore
`a_{D-1}`, `a_D`, and `a_{D+1}` are positive, and both ratios in
`Delta_eff` are defined.

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

With `m` denoting positive support size as above, the exact first inequality
also gives

```text
Delta_eff
  >= (m+1)/((D+1)(m-D+1))
  > 1/(D+1).
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

The remaining localization target is narrowed further in:

```text
notes/literature/one_sided_localization_reduction_2026-07-04.md
```

That follow-up note reduces the only one-sided gap to the local-mode mean lemma. The lemma is proved in:

```text
notes/literature/local_mode_mean_bound_proof_2026-07-04.md
```

Consequently, the sufficient localization `D+1 <= 4V` is now proved for one-sided low-probability PB laws with `V >= 1`, and hence the one-sided effective-drop bound `Delta_eff >= 1/(4V)` follows from Newton. This argument alone does not prove the signed case; the 2026-07-10 universal theorem cited at the top now does so by a different route.

## Plateau-Safe Endpoint

For later signed perturbation, define the first weak descent

```text
D_weak = min { k >= 1 : a_k <= a_{k-1} }.
```

The support and localization argument above applies with
`D_weak<=D`. Newton therefore gives

```text
1 - a_{D_weak+1}a_{D_weak-1}/a_{D_weak}^2
  >= 1/(D_weak+1)
  >= 1/(4V).
```

This weak-descent formulation is essential at plateaux. For
`X~Binomial(5,1/2)`, the first strict descent is `4` and its effective drop is
`3/5`, whereas the first weak descent is `3` and its effective drop is `1/2`.
If `Y~Bernoulli(q)` with `q>0`, the first strict descent of `X-Y` is `3` and

```text
Delta_eff(X-Y)
  = 1 - 10(5-4q)/(10-5q)^2
  -> 1/2
```

as `q->0+`, even though `Var(Y)->0`. Thus strict-first-descent reserve is not
uniformly continuous in a vanishing reflected variance. A perturbative signed
proof must anchor at `D_weak` and handle descent-index selection explicitly.
