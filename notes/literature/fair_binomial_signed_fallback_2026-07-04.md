# Fair-Binomial Signed Fallback
Date: 2026-07-04

## Purpose

This records the exact half-heavy model calculation for the fallback branch of
issue #5. It is a proof-track note for the Newton-or-reserve split; it does
not prove the full signed reserve theorem.

## Lemma

Let

```text
X ~ Binomial(m, 1/2),
Y ~ Binomial(n, 1/2),
N = m+n,
V = Var(X-Y) = N/4.
```

Assume `V >= 1`, equivalently `N >= 4`. Let `D` be the first strict descent of
the signed law `Z=X-Y`, and define

```text
R_- = P(Z=D)/P(Z=D-1),
R_+ = P(Z=D+1)/P(Z=D),
Delta_eff = 1 - R_+/R_-.
```

Then

```text
V * Delta_eff >= 5/8,
V * (1 - R_+) >= 3/4.
```

The constants are sharp for this exact model: both minima occur at `N=4`.

## Proof

Since `n-Y ~ Binomial(n,1/2)` and is independent of `X`,

```text
Z+n = X + (n-Y) ~ Binomial(N, 1/2).
```

Thus the signed sequence is just a shifted fair-binomial coefficient sequence.
Write `k = D+n`. The first strict descent occurs at the first `k` with

```text
binom(N,k) < binom(N,k-1),
```

or equivalently

```text
k > (N+1)/2.
```

Therefore

```text
k = floor((N+1)/2) + 1.
```

For the fair binomial sequence,

```text
R_- = binom(N,k)/binom(N,k-1) = (N-k+1)/k,
R_+ = binom(N,k+1)/binom(N,k) = (N-k)/(k+1).
```

Hence

```text
Delta_eff
  = 1 - R_+/R_-
  = 1 - k(N-k)/((k+1)(N-k+1))
  = (N+1)/((k+1)(N-k+1)).
```

There are two parity cases.

If `N=2r`, then `k=r+1`, `r>=2`, and

```text
V * Delta_eff
  = (r/2) * (2r+1)/(r(r+2))
  = (2r+1)/(2(r+2))
  >= 5/8.
```

The last inequality is minimized at `r=2`.

Also

```text
1 - R_+ = 1 - (r-1)/(r+2) = 3/(r+2),
V * (1-R_+) = 3r/(2(r+2)) >= 3/4.
```

If `N=2r+1`, then `k=r+2`, `r>=2`, and

```text
V * Delta_eff
  = ((2r+1)/4) * ((2r+2)/(r(r+3)))
  = (2r+1)(r+1)/(2r(r+3))
  >= 3/4.
```

For the last inequality, cross-multiplication gives
`2(r-1)(r-2) >= 0`.

Also

```text
1 - R_+ = 1 - (r-1)/(r+3) = 4/(r+3),
V * (1-R_+) = (2r+1)/(r+3) >= 1.
```

Combining the two parity cases gives the stated bounds.

At the sharp case `N=4`, the post-descent ratio `R_+` uses the last support
point. Thus the endpoint case must not be excluded by an interior-support
assumption.

## Calibration

An earlier version of this note incorrectly stated that the sharp
effective-drop constant was `1/2`, attained at `N=5`. The correct values at
`N=5` are

```text
R_- = 1/2,
R_+ = 1/5,
Delta_eff = 3/5,
V * Delta_eff = 3/4.
```

The global minimum under `V>=1` is instead attained at `N=4`, where

```text
R_- = 2/3,
R_+ = 1/4,
Delta_eff = 5/8,
V * Delta_eff = 5/8,
V * (1-R_+) = 3/4.
```

This exact model supports a conservative half-heavy fallback branch, but it
does not by itself prove stability under sparse dust or arbitrary
near-half-heavy perturbations.

## Overclaim Guard

This lemma covers only the exact fair-binomial model. It does not prove the
thresholded half-heavy stability lemma, the signed reserve theorem, the
hub-bouquet reserve, or Erdos 993.
