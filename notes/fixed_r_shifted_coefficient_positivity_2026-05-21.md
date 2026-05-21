# Shifted-Coefficient Positivity Principle

## Purpose

The fixed-`r` certificates repeatedly prove inequalities of rational functions
in a residue-class parameter.  This note states the elementary positivity
principle used by the scripts.

This is the justification for checks such as:

```text
shift t by the threshold;
check every numerator and denominator coefficient is positive.
```

## Principle

Let

```text
R(t) = N(t)/D(t)
```

be a rational function with rational coefficients.  Fix an integer `T >= 0`.
If both shifted polynomials

```text
N(u+T)
D(u+T)
```

have strictly positive coefficients, then

```text
R(t) > 0
```

for every integer `t >= T`.

The same conclusion holds if the coefficients are nonnegative and each shifted
polynomial has at least one strictly positive coefficient.  The scripts use
the stricter all-positive condition.

## Proof

For every integer `t >= T`, write

```text
t = u + T,    u >= 0.
```

If every coefficient of `N(u+T)` is positive, then `N(t)>0` for every
`u>=0`.  Likewise `D(t)>0`.  Therefore

```text
R(t)=N(t)/D(t)>0.
```

This proves the principle.

## Residue-Class Use

For fixed `r`, each asymptotic certificate is split by

```text
a = 3t + q,    q in {0,1,2}.
```

Given a threshold `A`, the relevant lower bound on `t` is

```text
T_q = ceil((A-q)/3).
```

To prove an inequality for all

```text
a >= A,    a == q mod 3,
```

the scripts substitute `a=3t+q`, form the rational expression in `t`, then
apply the shifted-coefficient principle at `T_q`.

## Script Mapping

The Python symbolic checker uses:

```text
shifted_poly(poly, threshold_t)
shifted_poly_positive(poly, threshold_t)
shifted_positive(expr, threshold_t)
```

in:

```text
gpt_attack/fixed_r_huboff_certificate.py
gpt_attack/fixed_r_hubon_mode_certificate.py
gpt_attack/fixed_r_hubon_route2_perturbation.py
```

The C++/GMP reserve helper receives already-shifted integer coefficient lists
from:

```text
gpt_attack/fixed_r_huboff_cpp_certificate.py
```

and checks that the final numerator and denominator coefficient lists are
strictly positive.

## What This Certifies

The principle turns a finite coefficient check into a proof of the inequality
for every integer point in the infinite tail of one residue class.  It is not a
sampling argument.

For the fixed-`r` Route-2 certificate, it is used to prove:

```text
left and right F-boundary margins;
lambda_0 bounds;
hub-off reserve positivity;
hub-on monotonicity/decrease conditions.
```
