# Fixed-`r` Fugacity-Shift Sublemma

## Purpose

This note formalizes the least explicit link in the fixed-`r` certificate
lemma: why the hub-on term changes the Route-2 fugacity and hub-off mean by at
most the quantity used in

```text
gpt_attack/fixed_r_hubon_route2_perturbation.py
```

The result is deliberately conservative.  It is meant to match the current
certificate code, not to optimize constants.

## Setup

For a fixed lane, write

```text
I(T;x) = F(x) + G(x),
```

and let `m` be the certified mode.  Define

```text
lambda_0 = F_{m-1}/F_m,
lambda   = (F_{m-1}+G_{m-1})/(F_m+G_m).
```

Assume:

```text
3/4 <= lambda_0 <= 2,
T >= max(G_{m-1}, G_m) / F_m,
T <= 1/2.
```

In the current scripts, `T` is bounded by `best_witness_term(...)`, which in
fact bounds `max_k G_k/F_m`, so it also bounds the two coefficients needed
here.

## Lemma 1: Fugacity Shift

Under the assumptions above,

```text
|lambda - lambda_0| <= 4T
```

and

```text
lambda >= 1/2.
```

### Proof

Compute directly:

```text
lambda - lambda_0
 = (F_m G_{m-1} - F_{m-1}G_m) / (F_m(F_m+G_m)).
```

Since all coefficients are nonnegative,

```text
|lambda - lambda_0|
 <= (F_m G_{m-1} + F_{m-1}G_m) / F_m^2
 <= T + lambda_0 T
 <= 3T
 <= 4T.
```

Also,

```text
lambda
 = (F_{m-1}+G_{m-1})/(F_m+G_m)
 >= F_{m-1}/(F_m+G_m)
 >= lambda_0/(1+T)
 >= (3/4)/(3/2)
 = 1/2.
```

This proves the lemma.

## Lemma 2: Hub-Off Mean Shift

Let `H(x)` be any polynomial of degree at most `N` with nonnegative
coefficients.  If `lambda_0, lambda >= 1/2`, then

```text
|mu_H(lambda) - mu_H(lambda_0)|
  <= (N^2/2) |lambda - lambda_0|.
```

Consequently, under Lemma 1,

```text
|mu_H(lambda) - mu_H(lambda_0)| <= 2N^2 T.
```

### Proof

For the Gibbs distribution with weights proportional to `h_k lambda^k`,

```text
d mu_H(lambda) / d lambda = Var_lambda(X) / lambda.
```

Since `0 <= X <= N`,

```text
Var_lambda(X) <= N^2/4.
```

On the interval between `lambda_0` and `lambda`, Lemma 1 gives
`lambda >= 1/2` after using the conservative certificate condition `T <= 1/2`.
Therefore

```text
|d mu_H / d lambda| <= (N^2/4)/(1/2) = N^2/2.
```

Integrating from `lambda_0` to `lambda` and using
`|lambda-lambda_0| <= 4T` gives

```text
|mu_H(lambda) - mu_H(lambda_0)|
  <= (N^2/2) 4T
  = 2N^2T.
```

This is exactly the `lambda_shift_bound` form used in the perturbation script.

## Lemma 3: Mixture Mean Shift

Let

```text
B(x) = H(x) + J(x)
```

with nonnegative coefficients and degree at most `N`.  For any `lambda > 0`,
if

```text
J(lambda)/H(lambda) <= U,
```

then

```text
|mu_B(lambda) - mu_H(lambda)| <= N U.
```

### Proof

The `B`-weighted distribution is a mixture of the `H`-weighted and
`J`-weighted distributions.  The mixture weight on the `J` part is

```text
w = J(lambda)/(H(lambda)+J(lambda)) <= J(lambda)/H(lambda) <= U.
```

Both component means lie in `[0,N]`, so replacing the `H` component by the
mixture can move the mean by at most `Nw <= NU`.

## Fixed-`r` Specialization

For the Route-2 bridge polynomial

```text
B_{a,r}(x) = F^-_a(x) + G^-_a(x),
F^-_a(x) = P_r(x)(1+2x)^(a-1),
G^-_a(x) = xP_{r-1}(x)(1+x)^(a-1),
```

use

```text
N <= a+r.
```

The current script proves:

```text
T(a,r,q) <= best_witness_term(...)
```

and, for `lambda in [1/2,2]`,

```text
G^-_a(lambda)/F^-_a(lambda) <= C_r (3/4)^(a-1),
```

where

```text
C_r = 2 P_{r-1}(2) / P_r(1/2).
```

Combining Lemmas 2 and 3 gives the perturbation bound used by the code:

```text
|mu_{B_{a,r}}(lambda) - mu_{F^-_a}(lambda_0)|
 <= 2(a+r)^2 T(a,r,q)
    + (a+r) C_r (3/4)^(a-1).
```

When the comparison target is proportional to `1/a`, the certificate must
verify that `a` times each perturbation term decreases as `a -> a+3` in each
residue class.  The corrected scripts check exactly this weighted monotonicity
for the fugacity-shift term and the mixture term.

## What This Closes

This note supplies the missing human-readable justification for the
`lambda_shift_bound` and `mixture_bound` terms in
`fixed_r_hubon_route2_perturbation.py`.

Remaining proof-writing obligations:

```text
1. State the coefficient-positivity principle for shifted rational functions.
2. State the path-polynomial evaluation recurrence for hub-off reserve.
3. Prove or certify the global F-mode margin from checked local margins.
```
