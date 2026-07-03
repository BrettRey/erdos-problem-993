# Signed Poisson-Binomial Conditional Ratio Identity
Date: 2026-07-03

## Purpose

This records the algebraic handle for the signed reserve lemma in issue #5. The current computational picture says the signed law

```text
Z = X - Y
```

has its worst reserve near the one-sided sparse boundary. The next proof should not work with an arbitrarily shifted ordinary polynomial. It should work directly with signed-support coefficients.

## Setup

Let

```text
a_i = P(X = i),
b_j = P(Y = j),
c_z = P(X - Y = z) = sum_j a_{z+j} b_j.
```

Assume `X` and `Y` are independent low-probability Poisson-binomial sums, so all Bernoulli parameters on both sides are at most `1/2`.

Let

```text
r_i = a_{i+1}/a_i,
s_j = b_{j+1}/b_j
```

where the ratios are defined inside the support.

## Exact X-Side Identity

For each signed support point `z` with `c_z > 0`, define the conditional distribution

```text
pi_z(j) = a_{z+j} b_j / c_z.
```

This is the law of `Y` conditional on `X - Y = z`. Then

```text
c_{z+1}/c_z
  = sum_j a_{z+1+j} b_j / c_z
  = E_{pi_z} [ a_{z+1+Y}/a_{z+Y} ]
  = E_{pi_z} [ r_{z+Y} ].
```

Similarly,

```text
c_z/c_{z-1}
  = E_{pi_{z-1}} [ r_{z-1+Y} ].
```

Thus the post-descent pressure is an average of the one-sided local ratios of `X`, but under the conditional law at the signed level `z`.

## Exact Y-Side Identity

The same ratio can be written using the `Y` side. Since

```text
c_z = sum_i a_i b_{i-z},
```

let

```text
rho_z(i) = a_i b_{i-z} / c_z,
```

the law of `X` conditional on `X - Y = z`. Then

```text
c_{z+1}/c_z
  = sum_i a_i b_{i-z-1} / c_z
  = E_{rho_z} [ b_{X-z-1}/b_{X-z} ]
  = E_{rho_z} [ 1/s_{X-z-1} ].
```

This is the reflected local-ratio form. It should be useful when the `Y` side, rather than the `X` side, carries the variance controlling the signed descent.

## Newton Input On Each Side

For a low-probability Poisson-binomial law, after removing the positive constant,

```text
a_i = e_i(u_1,...,u_m),
u_l = p_l/(1-p_l) <= 1.
```

Newton's inequalities give

```text
r_i <= [i/(i+1)] r_{i-1}.
```

Equivalently,

```text
1 - r_i/r_{i-1} >= 1/(i+1)
```

whenever the ratios are positive. The same statement holds for the `s_j` ratios of `Y`.

The signed problem is therefore not missing log-concavity. The missing step is to convert conditional averages of these one-sided ratio drops into a shift-invariant bound in terms of

```text
V = Var X + Var Y.
```

## Proof Split Suggested By The Optimizer

The balance-constrained optimizer suggests the following two-regime target.

### Near-One-Sided Regime

If

```text
min(Var X, Var Y) <= epsilon V,
```

then `Z = X - Y` should be treated as a low-variance perturbation of the dominant one-sided law. The one-sided theorem already gives

```text
reserve >= 1/(5 Var X)
```

or the reflected analogue. The remaining task is a stability estimate: convolving with a reflected law of variance `epsilon V` should not reduce the first-descent reserve by more than a controlled fraction. For issue #5, a very crude stability loss is enough, because the target constant is only `1/4` in the product-term lemma and `1/5` in the proved one-sided lemma.

### Genuinely Two-Sided Regime

If

```text
Var X >= epsilon V
and
Var Y >= epsilon V,
```

then the optimizer sees a larger reserve. The proof should use the two conditional identities above plus smoothing: both conditional laws have enough spread that the averaged one-sided Newton drops cannot concentrate entirely at the sparse boundary.

The computational side-variance runs give the right calibration:

| `epsilon` | best observed `V * reserve` at `V >= 1` |
|---:|---:|
| `0` | `0.5014196117` |
| `0.05` | `0.5121224485` |
| `0.10` | `0.5264484058` |
| `0.25` | `0.5774912085` |

These are empirical guideposts, not constants in a theorem.

## Immediate Lemma Targets

1. **Conditional-index localization.** At the first signed descent `D`, bound the relevant conditional index, for example `E_{\pi_D}[D+Y]`, by `O(Var X + Var Y)`.
2. **Conditional Newton drop.** Turn

```text
r_i <= [i/(i+1)] r_{i-1}
```

into a drop for `E_{\pi_D}[r_{D+Y}]` relative to the preceding signed ratio.
3. **Perturbative side.** Prove that if `Var Y <= epsilon V`, then the signed first-descent reserve differs from the one-sided `X` reserve by `O(epsilon/V)` or another loss small enough to preserve a constant.
4. **Balanced side.** Prove a direct reserve bound under `Var X, Var Y >= epsilon V`.

The first target is the most concrete next calculation: it replaces the artificial shifted polynomial index by the actual conditional success count on one side.
