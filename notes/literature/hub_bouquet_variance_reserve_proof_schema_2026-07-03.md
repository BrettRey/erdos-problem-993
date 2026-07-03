# Hub-Bouquet Variance Reserve Proof Schema
Date: 2026-07-03

## Purpose

This note turns issue #5 into a proof checklist. The computational evidence says that the hardest crossing-pressure examples are high-leaf hub/broom families. The right proof object is not the full tree recurrence. It is the coefficient distribution of

```text
A(x) = (1+x)^s Q(x),
```

where `Q` is a product of path independence polynomials. The hub-included term

```text
x R(x)
```

is a perturbation to handle after the product term is understood.

## Exact Setup

For a hub with `s` pendant leaves and path arms `a_1, ..., a_m`,

```text
I(x) = (1+x)^s Q(x) + x R(x),
Q(x) = prod_j I(P_{a_j}; x),
R(x) = prod_j I(P_{a_j-1}; x).
```

For the product term

```text
A(x) = (1+x)^s Q(x) = sum_k A_k x^k,
Q(x) = sum_y q_y x^y,
```

the coefficients are

```text
A_k = sum_y q_y binom(s, k-y).
```

Normalize by `A(1)`. Then `A_k/A(1)` is the law of

```text
X + Y,
X ~ Bin(s, 1/2),
Y distributed by q_y / Q(1),
```

and if `Q` is a path-product polynomial, `Y` is itself a Poisson-binomial random variable. Hence `X+Y` is Poisson-binomial.

## Proven Reduction 1: Worst Tail Ratio Is Adjacent To The Mode

Path independence polynomials are real-rooted with nonpositive roots. Therefore `Q` and `(1+x)^s Q` are real-rooted with nonnegative coefficients. Their coefficient sequences are log-concave.

Let

```text
r_k = A_k / A_{k-1}
D = min { k : r_k < 1 }.
```

Since log-concavity makes `r_k` nonincreasing, for every `j >= D`,

```text
A_{j+1}/A_j = r_{j+1} <= r_{D+1}.
```

Thus for the product term, crossing-only separation reduces to the single adjacent ratio

```text
A_{D+1}/A_D < 1.
```

This is already automatic from log-concavity. The quantitative task is to bound its reserve away from zero at the correct scale.

## Exact Ratio Identity

The adjacent ratio has an exact conditional expectation form. For every `k` with `A_k > 0`, define a probability distribution on arm sizes `Y` by

```text
P_k(Y=y) = q_y binom(s, k-y) / A_k.
```

Then

```text
A_{k+1}/A_k
  = E_k [ binom(s, k+1-Y) / binom(s, k-Y) ]
  = E_k [ (s-k+Y) / (k+1-Y) ].
```

Equivalently,

```text
1 - A_{k+1}/A_k
  = E_k [ (2(k-Y)+1-s) / (k+1-Y) ].
```

This identity is the most concrete algebraic handle in the problem. It says the reserve is governed by the conditional leaf-count excess over the middle of the binomial layer.

## Asymptotic Form Of The Desired Bound

Let

```text
mu = E[X+Y],
V  = Var(X+Y).
```

For a Poisson-binomial law, Darroch localizes the mode within distance `< 1` of `mu`. If `D` is the first strict descent index, then `D - mu` is bounded. A local-normal expansion gives

```text
log(A_{D+1}/A_D)
  = - (D + 1/2 - mu) / V + lower-order terms,
```

so

```text
1 - A_{D+1}/A_D
  = (D + 1/2 - mu) / V + lower-order terms.
```

The computations match this: for the dominant broom product term at `s = 5000`, `V * reserve` stays order `1` as the broom handle grows, while `s * reserve` does not.

## Candidate Lemma

The proof target should be:

> **Variance Reserve Lemma.** Let `A(x)` be a Poisson-binomial generating polynomial with variance `V` and coefficient sequence `(A_k)`. Let `D` be the first strict descent index. Under a mild nondegeneracy condition, there is an absolute constant `c > 0` such that
>
> ```text
> A_{D+1}/A_D <= 1 - c/V.
> ```

For fixed or controlled arms in `(1+x)^s Q(x)`, `V = s/4 + O(s)` or `s/4 + O(1)`, yielding the earlier `1 - c'/s` reserve as a corollary. For very long broom handles, `V` is the correct denominator.

The natural first constant to try is not large. The data suggest `V * reserve` often lies between `1` and `2` in the hard range, so even `c = 1/4` would be a useful theorem.

## Perturbation By The Hub-Included Term

The full hub-bouquet polynomial is

```text
I(x) = A(x) + xR(x).
```

For fixed arms, `R` has bounded degree while the first descent of `A` is near `s/2 + O(1)`. Thus `xR` contributes nothing at the relevant coefficients once `s` is large enough. This gives an immediate fixed-arm theorem after the Variance Reserve Lemma.

For growing arms, `xR` is not automatically absent. The next perturbation statement should be:

> If `R` is the one-step-shortened path product associated with `Q`, then near the first descent of `A`, the coefficient ratios of `A+xR` differ from those of `A` by `o(1/V)` or by a controlled fraction of the reserve.

This is the real broom case. It should be attacked after the product-term lemma, not before.

## Immediate Proof Tasks

1. Prove or source the Variance Reserve Lemma for Poisson-binomial laws.
2. If no off-the-shelf theorem exists, prove the weaker version needed here using the exact ratio identity and Darroch mode localization.
3. Certify the fixed-arm hub-bouquet corollary: for fixed path-product `Q`, the term `xR` is absent near first descent for all large `s`.
4. Extend to broom handles `l = l(s)` using the perturbation estimate between `(1+x)^s I(P_l)` and `(1+x)^s I(P_l) + xI(P_{l-1})`.

## Why This Is Promising

The earlier LC-failure corpora were useful for falsifying proof mechanisms, but they were not optimizing the crossing obstruction. The hub-ridge computations are optimizing the obstruction directly. They point to a standard probabilistic object: local ratios of sums of independent Bernoulli variables near the mode. That is a much smaller problem than tree independence polynomial unimodality.
