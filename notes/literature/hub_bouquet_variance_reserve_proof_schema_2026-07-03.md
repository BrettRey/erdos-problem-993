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

## Binomial Base Case: Proved

The homogeneous case is now closed in

```text
notes/literature/binomial_variance_reserve_lemma_2026-07-03.md
```

For `b_k = binom(n,k)p^k(1-p)^{n-k}`, with variance `V=np(1-p)`, let `D` be the first strict descent and assume `D<n`. Then

```text
V * (1 - b_{D+1}/b_D) >= V/(V+3).
```

In particular, if `V >= 1`,

```text
b_{D+1}/b_D <= 1 - 1/(4V).
```

The proof is exact. Writing `a=(n+1)p`, `D=floor(a)+1`, and `theta=D-a`, one has

```text
1 - b_{D+1}/b_D = (theta+1)/((D+1)(1-p)),
V * (1 - b_{D+1}/b_D) = np(theta+1)/(D+1).
```

The lower bound follows from `theta+1 >= 1`, `D+1 <= np+3`, and `V <= np`.

I also added an exact rational checker:

```bash
python3 scripts/verify_binomial_variance_reserve.py \
  --max-n 200 \
  --max-den 100 \
  --out results/binomial_variance_reserve_check_2026-07-03.json
```

It checked 949,721 rational parameter rows with zero failures. The smallest sampled `V * reserve` among rows with `V >= 1` was about `0.561`, so the theorem's `1/4` constant is conservative.

## General Poisson-Binomial Falsification Probe

I added a separate falsification probe for the candidate lemma:

```bash
python3 scripts/probe_pb_variance_reserve.py \
  --out results/pb_variance_reserve_probe_2026-07-03.json
```

It scanned 7,645 binomial-grid rows and 2,440 deterministic random Poisson-binomial rows. It did not find evidence that `V * reserve` collapses toward zero. The smallest observed values by variance cutoff were:

| Variance cutoff | Source | `V` | Pressure | `V * reserve` |
|---:|---|---:|---:|---:|
| 1 | random PB, `beta_0.5_5` | `1.3000146658` | `0.5349365442` | `0.6045893131` |
| 2 | binomial `n=20, p=0.142` | `2.43672` | `0.7033799534` | `0.7227800000` |
| 5 | binomial `n=50, p=0.117` | `5.16555` | `0.8328749393` | `0.8632928571` |
| 10 | binomial `n=50, p=0.294` | `10.3782` | `0.9109419263` | `0.9242625000` |
| 20 | binomial `n=100, p=0.297` | `20.8791` | `0.9539760474` | `0.9609387097` |
| 50 | binomial `n=500, p=0.499` | `124.9995` | `0.9920398247` | `0.9950179283` |
| 500 | binomial `n=5000, p=0.887` | `501.155` | `0.9977801313` | `1.1124983097` |

Interpretation: the candidate lemma still looks plausible, and the extremal behavior may already be visible in ordinary binomial laws. The probe also warns against trying to prove a sharp constant from the tree data. A conservative constant such as `c = 1/4` or `c = 1/2`, with an explicit lower-variance exception range, would be enough for the hub lane.

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
