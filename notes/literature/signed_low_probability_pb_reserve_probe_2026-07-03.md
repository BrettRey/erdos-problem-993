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

## Follow-Up Optimizer

I added an adversarial grouped optimizer:

```bash
python3 scripts/optimize_signed_pb_reserve.py \
  --out results/signed_pb_reserve_optimizer_2026-07-03.json
```

It searches laws of the form

```text
X = sum_g Bin(x_g, p_g),   Y = sum_h Bin(y_h, q_h),
p_g, q_h <= 1/2,
```

and minimizes `V * reserve` over grouped block counts and probabilities. The default run covered six `(x_n,y_n)` pairs and six variance cutoffs, for `31` feasible optimizer runs.

The best rows by variance cutoff were:

| Variance cutoff | `(x_n,y_n)` | `V` | Pressure | `V * reserve` |
|---:|---:|---:|---:|---:|
| 1 | `(200,50)` | `1.0000007257` | `0.4985807522` | `0.5014196117` |
| 2 | `(100,100)` | `2.0054842731` | `0.6643453157` | `0.6731501905` |
| 5 | `(50,50)` | `5.0595761830` | `0.8336119443` | `0.8418530436` |
| 10 | `(100,100)` | `10.7669509552` | `0.9146638892` | `0.9188097196` |
| 20 | `(200,200)` | `21.9441263825` | `0.9562838118` | `0.9593135594` |
| 50 | `(200,50)` | `53.5762478932` | `0.9816026958` | `0.9856585309` |

The optimizer improves the broad probe's minimum from about `0.5615` to about `0.5014`. That is essentially the same sparse one-sided boundary seen in the unsigned grouped optimizer: the best row has `Y` almost deterministic at zero and `X` close to a sparse low-mean law.

This strengthens the empirical message but not the theorem. The working constant `c = 1/4` still has a large buffer in the adversarial signed search, while a sharp constant above `1/2` remains implausible.

## Balanced Side-Variance Runs

The optimizer now also records `x_variance` and `y_variance` and can enforce

```text
min(Var X, Var Y) >= rho * (Var X + Var Y).
```

I ran three constrained optimizations:

```bash
python3 scripts/optimize_signed_pb_reserve.py \
  --min-side-variance-fraction 0.05 \
  --out results/signed_pb_reserve_optimizer_balanced_f05_2026-07-03.json

python3 scripts/optimize_signed_pb_reserve.py \
  --min-side-variance-fraction 0.10 \
  --out results/signed_pb_reserve_optimizer_balanced_f10_2026-07-03.json

python3 scripts/optimize_signed_pb_reserve.py \
  --min-side-variance-fraction 0.25 \
  --out results/signed_pb_reserve_optimizer_balanced_f25_2026-07-03.json
```

Best observed `V * reserve` values by variance cutoff:

| Side variance floor | `V >= 1` | `V >= 2` | `V >= 5` | `V >= 10` | `V >= 20` | `V >= 50` |
|---:|---:|---:|---:|---:|---:|---:|
| none | `0.5014196117` | `0.6731501905` | `0.8418530436` | `0.9188097196` | `0.9593135594` | `0.9856585309` |
| `0.05 V` | `0.5121224485` | `0.6832518560` | `0.8730412077` | `0.9292849474` | `0.9626594641` | `0.9842424147` |
| `0.10 V` | `0.5264484058` | `0.7077536943` | `0.8816073736` | `0.9314281715` | `0.9626594641` | `0.9847391662` |
| `0.25 V` | `0.5774912085` | `0.8103035475` | `0.8865137720` | `0.9335001870` | `0.9626594641` | `0.9894354882` |

Interpretation: the optimizer tries to sit exactly on the side-variance constraint at low variance. As the required two-sided variance share increases, the worst observed reserve moves away from the `1/2` sparse boundary. This suggests a proof split:

1. If one side has very small variance, reduce to the one-sided low-probability lemma plus a perturbation by a low-variance reflected law.
2. If both sides carry a fixed variance fraction, prove a stronger genuinely two-sided reserve, likely by a local-ratio or smoothing argument.

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
scripts/optimize_signed_pb_reserve.py
results/signed_pb_reserve_probe_2026-07-03.json
results/signed_pb_reserve_optimizer_2026-07-03.json
results/signed_pb_reserve_optimizer_balanced_f05_2026-07-03.json
results/signed_pb_reserve_optimizer_balanced_f10_2026-07-03.json
results/signed_pb_reserve_optimizer_balanced_f25_2026-07-03.json
```
