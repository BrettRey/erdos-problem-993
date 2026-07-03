# Signed Low-Probability Poisson-Binomial Reserve Probe
Date: 2026-07-03

## Purpose

The one-sided low-probability reserve lemma handles sums

```text
X = sum_i Bernoulli(p_i), p_i <= 1/2.
```

After splitting variables with `p_i > 1/2`, a general Poisson-binomial law is a deterministic shift plus

```text
X - Y,
```

where `X` and `Y` are independent low-probability Poisson-binomial sums. The one-sided Newton proof does not automatically apply to this signed law: if one shifts the Laurent polynomial for `X-Y` into an ordinary polynomial, Newton's ratio-drop index sees the artificial shift. The invariant quantity is the first descent in the original signed support coordinates.

This probe stress-tests the signed-support reserve directly.

## Script

```bash
python3 scripts/probe_signed_pb_reserve.py \
  --random-samples 2000 \
  --binomial-grid-size 31 \
  --max-groups 6 \
  --out results/signed_pb_reserve_probe_2026-07-03.json
```

The script forms

```text
c_z = P[X - Y = z],
V = Var(X) + Var(Y),
D = min { z : c_z < c_{z-1} },
reserve = 1 - c_{D+1}/c_D.
```

It records `V * reserve` and checks the candidate inequalities

```text
reserve >= 1/(4V),
reserve >= 1/(5V)
```

for rows with `V >= 1`.

## Coverage

The run processed `11,820` signed laws:

- `2,000` random grouped low-probability pairs,
- a `31 x 31` grid of two-binomial pairs over several `(n_x,n_y)` choices,
- finite Skellam approximants to `Pois(lambda) - Pois(eta)`.

There were `11,683` rows with variance at least `1`.

The run found:

```text
quarter_failures = 0
fifth_failures   = 0
```

## Best Rows By Variance Cutoff

| Variance cutoff | Source | `V` | Pressure | `V * reserve` |
|---:|---|---:|---:|---:|
| 1 | random grouped | `1.0410653829` | `0.4606454539` | `0.5615033470` |
| 2 | two binomial | `2.0161400000` | `0.6375323963` | `0.7307854345` |
| 5 | random grouped | `5.0628662021` | `0.8255957773` | `0.8829852444` |
| 10 | random grouped | `10.2764525264` | `0.9100315277` | `0.9245567346` |
| 20 | random grouped | `20.0139966764` | `0.9514693263` | `0.9712927419` |
| 50 | two binomial | `57.8717500000` | `0.9829115306` | `0.9889396274` |

The smallest observed `V * reserve` above variance `1` was about `0.5615`, still above the working constants `1/4` and `1/5`.

## Interpretation

This is evidence, not a theorem. Its useful information is negative: the obvious signed-law stress families did not create a new cancellation mechanism. The hardest rows still look close to the one-sided sparse boundary or to a deterministic shift plus a low-variance low-probability component.

The next proof target is therefore precise:

> Prove a signed low-probability reserve lemma for `h + X - Y`, with `X` and `Y` independent low-probability Poisson-binomial sums, using the first descent in signed support coordinates.

One possible algebraic form is the coefficient sequence of

```text
x^h prod_i (1 + u_i x) prod_j (x + v_j),
```

but any proof must be invariant under the arbitrary shift `h`. Equivalently, one can try to prove a direct local-ratio inequality for the convolution of a low-probability law with a reflected low-probability law.

## Repository Artifacts

```text
scripts/probe_signed_pb_reserve.py
results/signed_pb_reserve_probe_2026-07-03.json
```
