# Fixed-`r` Hub-Off Reserve Recurrence Certificate

## Purpose

This note states the exact recurrence used to certify the hub-off reserve in
the fixed-`r` Route-2 proof scaffold.  It explains what the Python and
C++/GMP scripts are proving, independently of implementation details.

## Hub-Off Object

For

```text
T_{a,r}=S(2^a,r),
F_a(x)=P_r(x)(1+2x)^a,
```

let

```text
m = (2a+D_q)/3,
lambda_0 = F_{m-1}/F_m.
```

After removing one length-2 arm, the hub-off bridge approximation is

```text
F^-_a(x)=P_r(x)(1+2x)^(a-1).
```

The hub-off reserve certificate proves

```text
mu_{F^-_a}(lambda_0) >= m - 4/3 + 1/(R a)
```

for all `a >= A` in each residue class.

Equivalently,

```text
mu_{F^-_a}(lambda_0) - m + 4/3 - 1/(R a) >= 0.
```

## Mean Formula

For any `lambda>0`,

```text
mu_{F^-_a}(lambda)
 = lambda P'_r(lambda)/P_r(lambda)
   + 2(a-1)lambda/(1+2lambda).
```

The second term is the mean contribution of `(1+2x)^(a-1)`.  The first term is
the path-polynomial contribution.

Thus the certificate must prove:

```text
lambda_0 P'_r(lambda_0)/P_r(lambda_0)
+ 2(a-1)lambda_0/(1+2lambda_0)
- m + 4/3
- 1/(R a)
>= 0.
```

## Cleared Path Recurrences

Write

```text
lambda_0 = L/Z.
```

Let

```text
d_n = deg P_n = ceil(n/2),
E_n = Z^d_n P_n(L/Z),
N_n = Z^d_n (L/Z) P'_n(L/Z).
```

The path polynomials satisfy

```text
P_0=1,
P_1=1+x,
P_n=P_{n-1}+xP_{n-2}.
```

The cleared recurrences are:

```text
E_0 = 1,          N_0 = 0,
E_1 = Z+L,        N_1 = L.
```

For `n >= 2`:

```text
if n is even:
  E_n = E_{n-1} + L E_{n-2},
  N_n = N_{n-1} + L(E_{n-2}+N_{n-2});

if n is odd:
  E_n = Z E_{n-1} + L E_{n-2},
  N_n = Z N_{n-1} + L(E_{n-2}+N_{n-2}).
```

These formulas are exactly what the scripts implement.  They avoid expanding
`P_r(L/Z)` and `(L/Z)P'_r(L/Z)` as nested rational expressions.

## Derivation

The recurrence for `E_n` follows from

```text
P_n=P_{n-1}+xP_{n-2}
```

after multiplying by `Z^d_n`.  Since

```text
d_n=d_{n-1}=d_{n-2}+1      when n is even,
d_n=d_{n-1}+1=d_{n-2}+1    when n is odd,
```

the even and odd cases above follow.

For `N_n`, differentiate:

```text
lambda P'_n(lambda)
 = lambda P'_{n-1}(lambda)
   + lambda P_{n-2}(lambda)
   + lambda^2 P'_{n-2}(lambda).
```

Multiplying by `Z^d_n` and using the same degree identities gives the stated
recurrences.

## Reserve Numerator

After the recurrences compute

```text
E_r = Z^d P_r(L/Z),
N_r = Z^d (L/Z)P'_r(L/Z),
```

the path mean is

```text
N_r/E_r.
```

The arm mean is

```text
2(a-1)L/(Z+2L).
```

Therefore the hub-off reserve inequality is:

```text
N_r/E_r
+ 2(a-1)L/(Z+2L)
- m + 4/3
- 1/(R a)
>= 0.
```

Equivalently, after multiplying by the positive denominator

```text
E_r (Z+2L) (3Ra),
```

one checks nonnegativity of the integer polynomial numerator:

```text
3Ra N_r (Z+2L)
+ 3Ra E_r 2(a-1)L
+ E_r (Z+2L) (Ra(4-3m)-3)
>= 0.
```

The scripts check a slightly rearranged but algebraically identical expression.

## Residue-Class Certification

For each residue class:

```text
a = 3t+q,
m = (2a+D_q)/3,
lambda_0 = L(t)/Z(t).
```

The scripts shift

```text
t = u + T_q,    T_q = ceil((A-q)/3),
```

convert `L` and `Z` to a common integer coefficient scale, run the recurrence,
and check that the final numerator and denominator polynomials have strictly
positive coefficients in `u`.

By the shifted-coefficient positivity principle, this proves the reserve for
every integer `t >= T_q`, hence every `a >= A` in that residue class.

## Common Scaling

The Python driver may replace `L,Z` by a common positive integer multiple:

```text
L -> cL,
Z -> cZ.
```

This does not change `lambda_0=L/Z`.  The cleared `E_n,N_n` are scaled by
positive powers of `c`, so coefficient positivity and the sign of the reserve
inequality are unchanged.

## Script Mapping

Python implementation:

```text
gpt_attack/fixed_r_huboff_certificate.py
```

Large-lane GMP implementation:

```text
gpt_attack/fixed_r_huboff_cpp_certificate.py
gpt_attack/fixed_r_huboff_reserve_cpp.cpp
```

The C++ helper performs the same recurrence and final numerator/denominator
positivity check using exact integer arithmetic.  Threading affects only the
order of exact polynomial multiplication, not the certificate statement.

## What This Closes

Together with:

```text
notes/fixed_r_shifted_coefficient_positivity_2026-05-21.md
notes/fixed_r_fugacity_shift_sublemma_2026-05-21.md
notes/fixed_r_global_f_mode_margin_sublemma_2026-05-21.md
```

this supplies the local proof justifications needed by:

```text
notes/fixed_r_certificate_lemma_2026-05-21.md.
```

The remaining arbitrary fixed-`r` problem is now the production of certificate
data for every fixed `r`, or a direct analytic lower bound replacing the
per-lane reserve positivity certificate.
