# Two-Sided Sparse Poisson Limit For The Variance Reserve
Date: 2026-07-03

## Purpose

The grouped Poisson-binomial optimizer found hardest rows with many probabilities near `0` and some probabilities near `1`. This note identifies the limiting local model for that shape and calibrates a limiting upper bound on any universal constant in the variance-reserve lemma.

## Limiting Model

Take independent Bernoulli variables in two blocks:

```text
N variables with p = lambda/N,
M variables with p = 1 - eta/M.
```

If the near-deterministic block is factored as a shift by `M`, then

```text
sum Bernoulli(1 - eta/M) = M - Bin(M, eta/M),
```

so the shifted local law converges to

```text
Z = Pois(lambda) - Pois(eta).
```

The variance tends

```text
V = lambda + eta.
```

This is the two-sided sparse boundary of the Poisson-binomial problem. It matches the apparent shape seen in the grouped optimizer: a large deterministic shift, a low-probability success component, and a small low-probability failure component from variables with `p_i` near `1`.

## Probe

I added:

```bash
python3 scripts/probe_skellam_reserve.py \
  --variance-values 1,2,5,10,20,50,100 \
  --grid-size 801 \
  --boundary-exponents 8 \
  --out results/skellam_reserve_probe_2026-07-03.json
```

The scan computes the first strict descent of the distribution of `Pois(lambda)-Pois(eta)` and the first post-descent pressure.

Best rows by fixed variance:

| `V` | `lambda` | `eta` | First descent | Pressure | `V * reserve` |
|---:|---:|---:|---:|---:|---:|
| 1 | `0.99999999` | `1e-8` | 1 | `0.4999999942` | `0.5000000058` |
| 2 | `1.99999998` | `2e-8` | 2 | `0.6666666578` | `0.6666666844` |
| 5 | `4.99999995` | `5e-8` | 5 | `0.8333333200` | `0.8333333998` |
| 10 | `9.9999999` | `1e-7` | 10 | `0.9090908931` | `0.9090910689` |
| 20 | `19.9999998` | `2e-7` | 20 | `0.9523809346` | `0.9523813078` |
| 50 | `49.9999995` | `5e-7` | 50 | `0.9803921378` | `0.9803931092` |
| 100 | `99.999999` | `1e-6` | 100 | `0.9900989904` | `0.9901009611` |

The near-sharp `V=1` boundary is:

| `eta` | `lambda` | Pressure | `V * reserve` |
|---:|---:|---:|---:|
| `1e-1` | `0.9` | `0.4433983167` | `0.5566016833` |
| `1e-2` | `0.99` | `0.4941852661` | `0.5058147339` |
| `1e-3` | `0.999` | `0.4994168540` | `0.5005831460` |
| `1e-4` | `0.9999` | `0.4999416685` | `0.5000583315` |
| `1e-5` | `0.99999` | `0.4999941667` | `0.5000058333` |
| `1e-6` | `0.999999` | `0.4999994167` | `0.5000005833` |
| `1e-7` | `0.9999999` | `0.4999999417` | `0.5000000583` |
| `1e-8` | `0.99999999` | `0.4999999942` | `0.5000000058` |

## Consequence

This gives a limiting obstruction to any universal constant larger than `1/2` in the variance-reserve lemma. In the limit

```text
lambda = 1 - eta,
eta down to 0,
V = lambda + eta = 1,
```

the first strict descent occurs at `1`, the post-descent pressure tends to `1/2`, and therefore

```text
V * reserve -> 1/2.
```

The homogeneous Poisson endpoint `lambda=1, eta=0` has a plateau at `0,1`, so its first strict descent moves to `2` and gives `V * reserve = 2/3`. The near-deterministic failure component breaks that plateau and exposes the limiting boundary approached by the grouped optimizer. This is consistent with why the finite grouped optimizer found best values just above `0.5`.

To turn this from calibration into a formal upper-bound proposition, one still has to write the standard finite approximation step: take `N` variables with `p=lambda/N` and `M` variables with `p=1-eta/M`, then pass the fixed local probabilities and adjacent ratios to the limit.

## Proof Implication

The right theorem target is now:

> For a Poisson-binomial law with variance `V >= 1`, the first post-descent reserve satisfies
>
> ```text
> 1 - a_{D+1}/a_D >= c/V
> ```
>
> for some absolute `c <= 1/2`.

The conservative working target `c = 1/4` remains appropriate. If the full theorem is true, this boundary suggests that any sharp constant is at most `1/2`; proving sharpness would require a separate argument.

For the hub-bouquet lane, sharpness is unnecessary. It is enough to prove any explicit constant below the sparse-boundary value.
