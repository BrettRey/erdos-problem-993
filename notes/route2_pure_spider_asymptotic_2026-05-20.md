# Route-2 Pure Spider Asymptotic (2026-05-20)

## Purpose

The wider spider-lane scans suggest that the large-`a` obstruction is not a
new counterexample family.  It is the binomial spine of the pure spider
`S(2^a)`, with Route-2 slack tending to `1/6` from above.

The closed-form checker is:

```bash
python3 gpt_attack/pure_spider_route2_asymptotic.py \
  --exact-max 250 --float-max 5000 \
  --out results/pure_spider_route2_asymptotic.json
```

Output:

- exact `a = 2..250`,
- mode formula failures: `0`,
- minimum exact Route-2 slack in that window: `0.169322709163` at `a=250`,
- minimum exact margin above `1/6`: `0.00265604249668` at `a=250`,
- worst exact error versus the binomial-spine approximation:
  `-0.00472836578512` at `a=12`.

Large samples:

```text
a=  100  m=   67  slack=0.173267326733  spine=0.173267326733
a=  199  m=  133  slack=0.17            spine=0.17
a=  500  m=  333  slack=0.170658682635  spine=0.170658682635
a= 1000  m=  667  slack=0.167332667333  spine=0.167332667333
a= 5000  m= 3333  slack=0.167066586683  spine=0.167066586683
```

An attempted exact rational run to `a=1000` was killed as not worth the
runtime; exact rational arithmetic grows too quickly and does not add much
once the asymptotic formula is isolated.

## Formulas

For the pure spider with `a` arms of length `2`,

```text
I_a(x) = (1 + 2x)^a + x(1 + x)^a.
```

The coefficient is:

```text
c_k(a) = 2^k binom(a,k) + binom(a,k-1).
```

The leftmost mode formula, proved below and checked exactly for
`2 <= a <= 250`, is:

```text
m_a = 2,                    a = 2,
m_a = floor((2a+1)/3),      a >= 3.
```

For Route-2, remove one length-2 arm, so `B = S(2^(a-1))`, and set:

```text
lambda_a = c_{m_a-1}(a) / c_{m_a}(a).
Delta_a  = mu_{S(2^(a-1))}(lambda_a) - (m_a - 3/2).
```

The mean of `S(2^b)` has the exact two-branch form:

```text
rho_b(lambda) = lambda ((1+lambda)/(1+2lambda))^b,

mu_b(lambda) =
[
  2b lambda/(1+2lambda)
  + rho_b(lambda) (1 + b lambda/(1+lambda))
]
/
[
  1 + rho_b(lambda)
].
```

Dropping the hub-on branch gives the binomial-spine fugacity:

```text
lambda_a^0 = m_a / (2(a-m_a+1)).
```

At this fugacity, the binomial-spine Route-2 slack is:

```text
Delta_a^0
  = 2(a-1)lambda_a^0/(1+2lambda_a^0) - (m_a - 3/2)
  = 3/2 - 2m_a/(a+1).
```

By residue class:

```text
a = 3t:     Delta_a^0 = (t+3)/(2(3t+1)),
a = 3t+1:   Delta_a^0 = (t+2)/(2(3t+2)),
a = 3t+2:   Delta_a^0 = (t+5)/(6t+6).
```

Each expression is strictly larger than `1/6` and tends to `1/6`.

## Mode Formula Proof

This section proves only the pure-spider mode formula.  It does not by itself
prove the Route-2 inequality; that still needs the asymptotic error bound in
the next section.

The mode formula is not just empirical.  Let

```text
d_k(a) = c_k(a) - c_{k-1}(a).
```

For `1 <= k <= a+1`, using

```text
binom(a,k) = binom(a,k-1) (a-k+1)/k,
binom(a,k-2) = binom(a,k-1) (k-1)/(a-k+2),
```

we get

```text
d_k(a)
= binom(a,k-1)
  [
    2^(k-1) (2a - 3k + 2)/k
    + (a - 2k + 3)/(a - k + 2)
  ].
```

The binomial prefactor is positive for this range, so only the bracket
matters.

Set

```text
A_k = 2a - 3k + 2.
```

For `a >= 3`, put `m = floor((2a+1)/3)`.

### Increase Through `m`

If `k <= m`, then `A_k >= 1`.  The first bracket term is positive.  The
second bracket term is greater than `-1`, since

```text
(a - 2k + 3)/(a - k + 2) > -1
iff
2a - 3k + 5 > 0,
```

and `2a - 3k + 5 = A_k + 3 >= 4`.

For `k = 1`, the first term is already positive enough.  For `k >= 2`,

```text
2^(k-1) A_k/k >= 2^(k-1)/k >= 1,
```

with equality only at `k=2`, while the second term is still strictly greater
than `-1`.  Hence `d_k(a) > 0` for every `1 <= k <= m`.

### Descent After `m`

For `k >= m+2`, the residue classes give `A_k <= -3`.  The first bracket term
is at most

```text
-3 * 2^(k-1)/k < -1,
```

while the second term is always less than `1` for `k > 1`.  Hence
`d_k(a) < 0` for every `k >= m+2`.

It remains to check `k = m+1`.  Write `a` by residue class:

```text
a = 3t:     m=2t,     A_{m+1}=-1,
a = 3t+1:   m=2t+1,   A_{m+1}=-2,
a = 3t+2:   m=2t+1,   A_{m+1}= 0.
```

For the first two classes the first term is negative and the second term is
nonpositive, so `d_{m+1}(a) < 0`.  In the third class, the first term vanishes
and the second term is

```text
(1-t)/(t+2).
```

Thus `d_{m+1}(a) = 0` only when `a=5`; it is negative for `a=8,11,...`.

Therefore the coefficient sequence strictly increases through `m`, then
decreases after `m`, except for the harmless tie `c_3(5)=c_4(5)`.  The
leftmost mode is

```text
m_a = floor((2a+1)/3),      a >= 3,
```

with the separate small case `m_2=2`.

## Route-2 Slack Proof

This lane now has a conservative proof of the Route-2 inequality
`Delta_a > 0`.

The exact fugacity can be written as:

```text
lambda_a =
  lambda_a^0
  (1 + epsilon_{m_a-1})/(1 + epsilon_{m_a}),

epsilon_k = k / (2^k (a-k+1)).
```

For `a >= 30`, the mode formula gives

```text
m_a - 1 >= (2a-5)/3.
```

Also `lambda_a` stays in a fixed compact interval inside `(0,infinity)`.
Indeed, for `a >= 30`, the residue formulas for `lambda_a^0` give
`lambda_a^0 >= 7/8`, while `lambda_a^0 <= 1`.  The same lower bound on
`m_a-1` gives

```text
epsilon_{m_a}, epsilon_{m_a-1}
  <= a / 2^((2a-5)/3)
  < 1/8,
```

where the right-hand side is already true at `a=30` and improves thereafter.
Hence

```text
7/9 <= (7/8)/(1+1/8) <= lambda_a <= 1+1/8 < 2.
```

For the crude bound below it is enough that `1/2 <= lambda_a <= 2`.

Let

```text
f_b(lambda) = 2b lambda/(1+2lambda),
h_b(lambda) = 1 + b lambda/(1+lambda),
rho_b(lambda) = lambda ((1+lambda)/(1+2lambda))^b,
```

so that

```text
mu_b(lambda) = [f_b(lambda) + rho_b(lambda) h_b(lambda)]/[1+rho_b(lambda)].
```

For `b=a-1`,

```text
|f_b(lambda_a) - f_b(lambda_a^0)|
  <= 2a |lambda_a - lambda_a^0|
  <= 4a^2 / 2^((2a-5)/3).
```

The hub-on branch contributes at most

```text
rho_b(lambda_a) |h_b(lambda_a)-f_b(lambda_a)|
  <= 4a (3/4)^(a-1),
```

because `lambda_a >= 1/2` gives
`(1+lambda_a)/(1+2lambda_a) <= 3/4`, while `lambda_a <= 2` and
`|h_b-f_b| <= h_b+f_b <= 2a`.

Therefore

```text
|Delta_a - Delta_a^0|
  <= 4a^2 / 2^((2a-5)/3) + 4a (3/4)^(a-1).       (E)
```

Both terms on the right are decreasing for `a >= 30`, since their successive
ratios are bounded by

```text
((31/30)^2) 2^(-2/3) < 1,
(31/30)(3/4) < 1.
```

At `a=30`, the right side of `(E)` is

```text
0.0394729524415 < 1/6.
```

For an exact certificate avoiding the decimal, use the weaker threshold bound

```text
4(30)^2/2^18 + 4(30)(3/4)^29
  = 1524235892972445 / 36028797018963968
  < 1/6.
```

Since `Delta_a^0 > 1/6` in every residue class, this proves `Delta_a > 0` for
all `a >= 30`.

The finite range is checked exactly by:

```bash
python3 gpt_attack/pure_spider_route2_bound.py
```

Output:

```text
finite exact range: a=2..29
finite failures: 0
minimum finite slack: a=28 slack=0.189596603771

asymptotic threshold: a >= 30
simple error bound at threshold: 0.0394729524415
rational upper at threshold: 1524235892972445/36028797018963968 ~= 0.0423060445835
reserve 1/6 - bound: 0.127193714225
rational reserve: 13441690830564649/108086391056891904 ~= 0.124360622083
term1 ratio bound at threshold: 0.672657849416
term2 ratio bound at threshold: 0.775
```

Thus Route-2 holds for the pure spider lane `S(2^a)`, removing a length-2 arm.

## Interpretation

Pure `S(2^a)` is not a threat.  It explains the apparent `1/6` asymptote in
the wider spider scans and gives a template for fixed-`r` lanes
`S(2^a,r)`: the finite path factor should shift the mode by `O(1)`, while the
hub-on branch remains exponentially small.
