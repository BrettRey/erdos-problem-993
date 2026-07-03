# Binomial Variance Reserve Lemma
Date: 2026-07-03

## Purpose

This proves the binomial base case for the variance-reserve program in issue #5. It is the model calculation behind the proposed Poisson-binomial lemma.

## Statement

Let

```text
b_k = binom(n,k) p^k (1-p)^{n-k},
0 < p < 1,
q = 1-p,
V = npq.
```

Let `D` be the first strict descent:

```text
D = min { k : b_k < b_{k-1} }.
```

Assume `D < n`, so the first post-descent ratio `b_{D+1}/b_D` is defined. Then

```text
V * (1 - b_{D+1}/b_D) >= V / (V+3).
```

In particular, if `V >= 1`, then

```text
b_{D+1}/b_D <= 1 - 1/(4V).
```

So the binomial law satisfies the candidate variance reserve with constant `c = 1/4` outside the finite/small-variance regime.

## Exact Formula

The consecutive coefficient ratios are

```text
r_k = b_k / b_{k-1}
    = ((n-k+1)/k) * (p/q).
```

Thus `r_k < 1` exactly when

```text
(n+1)p < k.
```

Set

```text
a = (n+1)p,
D = floor(a) + 1,
theta = D - a.
```

Then

```text
0 < theta <= 1,
```

where `theta = 1` in the integer case.

The first post-descent ratio is

```text
b_{D+1}/b_D
  = ((n-D)/(D+1)) * (p/q).
```

Therefore

```text
1 - b_{D+1}/b_D
  = ((D+1)q - (n-D)p) / ((D+1)q)
  = (D+1 - (n+1)p) / ((D+1)q)
  = (theta+1) / ((D+1)q).
```

Multiplying by `V = npq` gives the exact expression

```text
V * (1 - b_{D+1}/b_D)
  = np(theta+1)/(D+1).
```

## Lower Bound

Since `theta > 0`,

```text
theta + 1 >= 1.
```

Also

```text
D+1 = (n+1)p + theta + 1
    = np + p + theta + 1
    <= np + 3,
```

because `p <= 1` and `theta <= 1`.

Hence

```text
V * (1 - b_{D+1}/b_D)
  = np(theta+1)/(D+1)
  >= np/(np+3).
```

Since `V = npq <= np` and the function `x/(x+3)` is increasing on `x >= 0`,

```text
np/(np+3) >= V/(V+3).
```

This proves

```text
V * (1 - b_{D+1}/b_D) >= V/(V+3).
```

If `V >= 1`, then `V+3 <= 4V`, so

```text
V/(V+3) >= 1/4.
```

Equivalently,

```text
1 - b_{D+1}/b_D >= 1/(4V).
```

## Sharpness And Interpretation

The proof is intentionally conservative. For large, central binomial laws, `V * reserve` approaches `1`, not `1/4`. The constant `1/4` is chosen because it is easy to prove and strong enough for the hub-bouquet program.

There can be no uniform positive lower bound on `V * reserve` without a lower-variance hypothesis: as `p -> 0` with fixed `n`, `V -> 0` and the product can go to `0`.

## Consequence For Issue #5

The full Poisson-binomial target should now be framed as:

```text
For a Poisson-binomial law with variance V >= 1,
the first post-descent reserve is at least c/V.
```

The binomial calculation proves this with `c = 1/4` in the homogeneous case. It also suggests that the heterogeneous proof should use:

1. mode/Darroch localization to place `D` just to the right of the mean;
2. log-concavity to reduce to the adjacent first-tail ratio;
3. a conditional expectation or local-ratio inequality to replace the exact binomial formula.
