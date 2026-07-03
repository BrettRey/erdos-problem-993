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

## X-Side Identity With Boundary Term

For each signed support point `z` with `c_z > 0`, define the conditional distribution

```text
pi_z(j) = a_{z+j} b_j / c_z.
```

This is the law of `Y` conditional on `X - Y = z`. If no lower boundary term enters, then

```text
c_{z+1}/c_z
  = sum_j a_{z+1+j} b_j / c_z
  = E_{pi_z} [ a_{z+1+Y}/a_{z+Y} ]
  = E_{pi_z} [ r_{z+Y} ].
```

In general there is a lower-boundary correction. With `a_i = 0` and `b_j = 0` outside support,

```text
c_{z+1}/c_z
  = E_{pi_z} [ r_{z+Y} ] + a_0 b_{-z-1}/c_z.
```

The extra term is zero unless `-z-1` is in the support of `Y`. This matters for rows whose first descent is negative.

Similarly,

```text
c_z/c_{z-1}
  = E_{pi_{z-1}} [ r_{z-1+Y} ] + a_0 b_{-z}/c_{z-1}.
```

Thus the post-descent pressure is an average of the one-sided local ratios of `X`, but under the conditional law at the signed level `z`.

## Reflected Y-Side Identity With Boundary Term

The same ratio can be written using the `Y` side. Define

```text
t_j = b_{j-1}/b_j,   t_0 = 0.
```

Then

```text
c_{z+1}/c_z
  = E_{pi_z} [ t_Y ] + a_{z+1+n_Y} b_{n_Y}/c_z,
```

where `n_Y` is the maximum support point of `Y`. The upper-boundary correction is zero unless `z+1+n_Y` is in the support of `X`.

Equivalently, using the law of `X` conditional on `X-Y=z`,

```text
rho_z(i) = a_i b_{i-z} / c_z,
```

the interior part is

```text
E_{rho_z} [ b_{X-z-1}/b_{X-z} ].
```

This reflected local-ratio form should be useful when the `Y` side, rather than the `X` side, carries the variance controlling the signed descent.

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

## Conditional Analyzer

I added a numerical checker for the identities and the conditional index statistics:

```bash
python3 scripts/analyze_signed_conditionals.py \
  --out results/signed_pb_conditional_analysis_2026-07-03.json
```

It reads the signed probe and optimizer certificates, extracts the best rows, reconstructs the PMFs, and evaluates the boundary-corrected conditional identities. The run processed `69` unique signed rows and found maximum identity error

```text
7.77e-16
```

up to floating-point precision.

The worst-reserve rows have the expected one-sided shape:

- `V` is essentially `1`;
- the first descent is at `D = 1`;
- `min(Var X, Var Y)/V` is tiny;
- `E[X | X-Y=D]` is about `1.00` to `1.03`;
- `V * E[1/(X+1) | X-Y=D]` is about `0.49` to `0.50`.

This numerically explains why the sparse boundary gives `V * reserve` near `1/2`: the conditional one-sided index is sitting at `X = 1`, where the Newton factor has its largest drop.

The same analysis gives a concrete target for localization. Among the analyzed rows,

```text
(E[X | X-Y=D] + 1) / V
```

was at most about `2.89`. This is only empirical, but it makes the first lemma target sharper: prove a shift-invariant bound of this form, with boundary terms handled explicitly.

## Conditional-Index Stress Probe

The analyzer above only looks at the best rows from existing certificates. I added a broader generated-corpus probe:

```bash
python3 scripts/probe_signed_conditional_index.py \
  --out results/signed_pb_conditional_index_probe_2026-07-03.json
```

This regenerated the signed random/grid/finite-Skellam corpus, analyzed all `11,820` rows, and tested the tentative bound

```text
(E[X | X-Y=D] + 1)/V <= 3.
```

That bound is false for this corpus: there were `2` X-side failures, with maximum observed ratio about `3.6131`. The worst row had

```text
V = 1.1108317260,
D = 3,
V * reserve = 0.8568451527,
(E[X | X-Y=D] + 1)/V = 3.6130670245.
```

The failures do not threaten the reserve conjecture; their reserve is large. They only show that the localization lemma should be stated with a looser constant.

A rerun with candidate bound `4` is recorded in:

```bash
python3 scripts/probe_signed_conditional_index.py \
  --candidate-bound 4 \
  --out results/signed_pb_conditional_index_probe_bound4_2026-07-03.json
```

It processed the same `11,820` rows and found zero X-side or Y-side failures. Thus `4V` is the current empirical localization target, not `3V`.

## Immediate Lemma Targets

1. **Conditional-index localization.** At the first signed descent `D`, bound the relevant conditional index, for example `E_{\pi_D}[D+Y]`, by `O(Var X + Var Y)`. The current empirical target is a constant around `4`, not `3`.
2. **Conditional Newton drop.** Turn

```text
r_i <= [i/(i+1)] r_{i-1}
```

into a drop for `E_{\pi_D}[r_{D+Y}]` relative to the preceding signed ratio.
3. **Perturbative side.** Prove that if `Var Y <= epsilon V`, then the signed first-descent reserve differs from the one-sided `X` reserve by `O(epsilon/V)` or another loss small enough to preserve a constant.
4. **Balanced side.** Prove a direct reserve bound under `Var X, Var Y >= epsilon V`.

The first target is the most concrete next calculation: it replaces the artificial shifted polynomial index by the actual conditional success count on one side.
