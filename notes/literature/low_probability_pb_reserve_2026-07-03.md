# Low-Probability Poisson-Binomial Reserve Lemma
Date: 2026-07-03

## Purpose

This records a proved partial case of the variance-reserve lemma from issue #5. It covers the one-sided sparse regime

```text
0 <= p_i <= 1/2.
```

The full grouped-optimizer obstruction is two-sided, because variables with `p_i` near `1` become small failure probabilities after a shift. Still, this lemma gives a real theorem-level foothold: the low-probability half already has the desired `1/V` reserve.

## Statement

Let

```text
S = B_1 + ... + B_n,
B_i ~ Bernoulli(p_i),
0 <= p_i <= 1/2,
a_k = P(S = k),
mu = E S,
V = Var S.
```

Let

```text
D = min { k >= 1 : a_k < a_{k-1} }
```

and assume `D+1` is in the support. If `V >= 1`, then

```text
1 - a_{D+1}/a_D > 1/(5V).
```

So the variance-reserve lemma is true in this one-sided low-probability regime, with the explicit nonsharp constant `c = 1/5`.

## Proof

Write the probability generating polynomial as

```text
f(x) = prod_i (1-p_i+p_i x)
     = C prod_i (1+w_i x),
w_i = p_i/(1-p_i).
```

Thus, up to the positive constant `C`,

```text
a_k = e_k(w_1, ..., w_n),
```

where `e_k` is the elementary symmetric polynomial.

Newton's inequalities imply that

```text
e_k^2 / binom(n,k)^2 >=
  (e_{k-1}/binom(n,k-1)) (e_{k+1}/binom(n,k+1)).
```

Equivalently,

```text
e_{k+1}/e_k <= [k(n-k)/((k+1)(n-k+1))] e_k/e_{k-1}.
```

In particular,

```text
e_{k+1}/e_k <= [k/(k+1)] e_k/e_{k-1}.
```

Let

```text
rho_k = a_{k+1}/a_k = e_{k+1}/e_k.
```

Since `D` is the first strict descent,

```text
rho_{D-1} = a_D/a_{D-1} < 1.
```

The Newton inequality gives

```text
rho_D <= [D/(D+1)] rho_{D-1} < D/(D+1).
```

Therefore

```text
1 - a_{D+1}/a_D
  = 1 - rho_D
  > 1/(D+1).
```

Now let `M` be the largest mode. The first strict descent is one step after the last mode, so

```text
D = M + 1.
```

Use Darroch's mode theorem in the weaker form

```text
|M - mu| <= 1.
```

Then

```text
D+1 = M+2 <= mu+3.
```

Because all `p_i <= 1/2`,

```text
p_i(1-p_i) >= p_i/2,
```

so

```text
V >= mu/2,
mu <= 2V.
```

Combining these estimates,

```text
1 - a_{D+1}/a_D
  > 1/(D+1)
  >= 1/(mu+3)
  >= 1/(2V+3).
```

If `V >= 1`, then `2V+3 <= 5V`, and therefore

```text
1 - a_{D+1}/a_D > 1/(5V).
```

This proves the claim.

## Sharpness

The constant `1/5` is not sharp. The proof loses in two places:

1. replacing the stronger Newton factor

```text
k(n-k)/((k+1)(n-k+1))
```

by `k/(k+1)`;

2. replacing the actual mode location by the coarse Darroch bound.

The Poisson boundary at `V = 1` has `V * reserve = 2/3` in the one-sided endpoint, so the theorem is deliberately conservative.

## Consequence For Issue #5

This proves the variance-reserve lemma in the low-probability PB case. The remaining hard case is not arbitrary heterogeneity; it is specifically the two-sided shifted law obtained when some variables have `p_i > 1/2`:

```text
S = h + X - Y,
```

where `X` and `Y` are low-probability Bernoulli sums. The sparse limiting obstruction

```text
Pois(lambda) - Pois(eta)
```

belongs to this two-sided case. The next proof target should therefore be a signed/shifted version of the Newton-ratio drop argument, not another general PB optimizer.

## Source Note

The only external input is Darroch's mode localization theorem for Poisson-binomial laws. Tang and Tang's survey, *The Poisson Binomial Distribution -- Old & New*, records this as Darroch's rule and states that a mode differs from the mean by at most `1`:

```text
https://arxiv.org/abs/1908.10024
```
