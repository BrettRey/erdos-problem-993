# Signed Conditional Reserve Reduction
Date: 2026-07-04

## Purpose

This note advances the signed part of issue #5 after the one-sided local-mode
mean lemma was proved. It does not prove the signed reserve theorem. Its role
is to isolate an exact conditional reduction for

```text
Z = X - Y,
```

where `X` and `Y` are independent low-probability Poisson-binomial sums, and
to separate theorem-level inputs from the remaining signed obstacles.

The main output is a boundary-aware inequality of the form

```text
signed effective drop
  >= inverse-index Newton term
     - conditional-dispersion penalty
     - boundary penalty.
```

This gives a smaller proof target than "prove the signed reserve lemma" and
explains why the one-sided theorem is useful but not by itself sufficient.

## Proved Inputs Available Now

The following inputs are now theorem-level within the project notes.

1. For a one-sided low-probability Poisson-binomial law with variance `V >= 1`,
   the first-descent effective ratio drop satisfies

```text
Delta_eff >= 1/(4V).
```

This follows from Newton's inequalities, the localization `D+1 <= 4V`, and
the local-mode mean bound proved in:

```text
notes/literature/local_mode_mean_bound_proof_2026-07-04.md
notes/literature/one_sided_effective_drop_reduction_2026-07-03.md
```

2. Pointwise Newton ratio drops hold on each low-probability side. If

```text
a_i = P(X=i),       r_i = a_{i+1}/a_i,
```

then for `i >= 1`,

```text
r_i/r_{i-1} <= i/(i+1) = 1 - 1/(i+1).
```

The same statement holds for the ratios of `Y`.

3. The boundary-corrected conditional identities for signed coefficients are
recorded in:

```text
notes/literature/signed_pb_conditional_ratio_identity_2026-07-03.md
```

The identities are algebraic and have been numerically checked on the current
signed corpus. The numerical checks are not proof of the reserve theorem.

## Signed Target

Let

```text
c_z = P(X-Y=z) = sum_j a_{z+j} b_j,
V = Var X + Var Y.
```

Let `D` be the first strict descent in signed coordinates:

```text
c_D < c_{D-1}.
```

The conditional reductions below are for positive-support descents, so assume
`c_D>0`. A terminal descent to zero is a signed support-edge case rather than
the interior first-descent obstruction targeted here. Also assume the
non-deterministic Bernoulli parameters on each side are strictly below `1`,
so the side supports are contiguous intervals after any deterministic shifts
are stripped off.

Write

```text
R_- = c_D/c_{D-1},
R_+ = c_{D+1}/c_D,
Delta = 1 - R_+/R_-.
```

Since `R_- < 1`, a lower bound on `Delta` implies at least the same lower
bound for the raw reserve:

```text
1 - R_+ > Delta.
```

Thus a signed effective-drop theorem

```text
Delta >= c/V
```

would be sufficient for the signed reserve route. The optimizer data suggest
that `c=1/4` is plausible but not sharp. The one-sided Poisson boundary shows
that this effective-drop route should not be expected to support any universal
constant above `1/3`.

## Corrected X-Side Conditional Reduction

Let `pi` be the conditional law of `Y` given `X-Y=D`:

```text
pi(y) = a_{D+y} b_y / c_D.
```

Under `pi`, write

```text
N = D + Y.
```

Thus `N` is the corresponding conditional value of `X`.

Define the X-side ratios

```text
r_n = a_{n+1}/a_n
```

inside the support, with `r_n=0` if `n+1` is outside the support. Define

```text
h_X = 1_{N>=1} / r_{N-1}.
```

The reciprocal denominator ratio has an upper-boundary correction. If `n_X`
is the maximum support point of `X`, then

```text
c_{D-1}/c_D = E_pi[h_X] + beta_X,
beta_X = a_{n_X} b_{n_X+1-D}/c_D,
```

where the boundary term is interpreted as zero when `n_X+1-D` is outside
the support of `Y`. Each term with `1 <= N <= n_X` maps from `a_N` to
`a_{N-1}` and contributes to `E_pi[h_X]`; the extra term corresponds to
`N=n_X+1`, which appears in `c_{D-1}` but not in the conditional law at
`D`.

The numerator ratio has one lower-boundary correction:

```text
c_{D+1}/c_D = E_pi[r_N] + nu_X,
nu_X = a_0 b_{-D-1}/c_D,
```

where the boundary term is interpreted as zero when `-D-1` is outside the
support of `Y`.

Therefore

```text
R_+/R_- = (E_pi[r_N] + nu_X) (E_pi[h_X] + beta_X).
```

Define

```text
I_X = E_pi[1/(N+1)],
H_X = E_pi[r_N] E_pi[h_X] - E_pi[r_N h_X].
```

The pointwise Newton drop gives, for every conditional atom,

```text
r_N h_X <= 1 - 1/(N+1).
```

For `N=0`, both sides are `0`; for `N>=1`, this is exactly
`r_N/r_{N-1} <= N/(N+1)`.

Consequently,

```text
Delta
  = 1 - (E_pi[r_N] + nu_X) (E_pi[h_X] + beta_X)
  >= I_X - H_X - nu_X E_pi[h_X]
     - (E_pi[r_N] + nu_X) beta_X.
```

This inequality is proved by algebra plus Newton. It is not a signed reserve
theorem because the last three penalty terms may, a priori, consume the
inverse-index gain.

Interpretation:

```text
I_X                  inverse-index Newton gain
H_X                  conditional-dispersion penalty from changing measures
nu_X E_pi[h_X]       lower-boundary penalty for the X-side identity
(E_pi[r_N]+nu_X)beta_X
                     upper-boundary penalty for the reciprocal denominator
```

If `Y` is deterministic at `0`, then `pi` is deterministic, `H_X=0`,
`nu_X=0`, and `beta_X=0` away from the top support edge. In that case this
reduces to the one-sided pointwise Newton drop at the signed descent. The
one-sided localization theorem then supplies the `1/(4V)` scale. This is the
precise sense in which the new one-sided result is a black-box endpoint for
the signed route.

### Why The `beta_X` Term Matters

An earlier version of this reduction omitted `beta_X`. That omission makes
the X-side bound false, not merely incomplete.

For example, take

```text
X = Bernoulli(1/2),
Y = Binomial(10,1/2).
```

The first strict descent is `D=-3`, and exact arithmetic gives

```text
Delta = 3/10.
```

The omitted-boundary X-side expression gives

```text
I_X - H_X - nu_X E[h_X] = 4/11,
```

which is larger than the true `Delta`. The missing term is

```text
beta_X = 42/55.
```

A low-probability example shows the same mechanism:

```text
X = Bernoulli(0.25),
Y = Binomial(20,0.25).
```

Here the first strict descent is `D=-4`; numerically,

```text
Delta ~= 0.229176,
old omitted-boundary X-side expression ~= 0.265599.
```

Thus the corrected X-side reduction must include the upper-boundary penalty.
In regimes where `Y` has much wider support than `X`, this penalty can be
large, so the reflected Y-side reduction may be the useful side.

## Reflected Y-Side Analogue

There is an analogous but slightly less clean reflected reduction. Let

```text
n_Y = max support point of Y,
t_y = b_{y-1}/b_y,       t_0 = 0.
```

Under the same conditional law `pi` at `D`, set

```text
k_Y = 1_{Y<n_Y} b_{Y+1}/b_Y = 1_{Y<n_Y}/t_{Y+1}.
```

The reflected numerator identity is

```text
R_+ = E_pi[t_Y] + nu_Y,
nu_Y = a_{D+1+n_Y} b_{n_Y}/c_D,
```

again with the boundary term interpreted as zero when the indicated `X`
index is outside support.

The reciprocal denominator has a lower reflected boundary:

```text
c_{D-1}/c_D = E_pi[k_Y] + beta_Y,
beta_Y = a_{D-1} b_0 / c_D.
```

Newton on the `Y` side gives

```text
t_Y k_Y <= 1 - 1/(Y+1).
```

Thus, with

```text
I_Y = E_pi[1/(Y+1)],
H_Y = E_pi[t_Y] E_pi[k_Y] - E_pi[t_Y k_Y],
```

one obtains

```text
Delta
  >= I_Y - H_Y
     - E_pi[t_Y] beta_Y
     - nu_Y (E_pi[k_Y] + beta_Y).
```

When the reflected boundary terms vanish, this has the same structure as the
X-side reduction.

The two reductions should be used together. Boundary terms that are awkward
on one side often correspond to an edge of signed support where the reflected
side is the more natural coordinate.

## What Remains Conjectural

The exact reductions above show that the signed problem is not missing
one-sided Newton drops. It is missing control of three signed effects.

### 1. Conditional Inverse-Index Localization

One needs a shift-invariant lower bound such as

```text
I_X = E[1/(X+1) | X-Y=D] >= c_0/V
```

or the reflected analogue for `Y`, at least on a side where the dispersion and
boundary penalties are controlled.

The existing conditional-index probes suggest constants in this range. For
example, the empirical bound

```text
(E[X | X-Y=D] + 1)/V <= 4
```

survived the generated signed corpus, while `3` failed in rows with large
reserve. By Jensen, an expectation bound of the form

```text
E[X | X-Y=D] + 1 <= C V
```

would imply

```text
I_X >= 1/(C V).
```

This is still a conjectural signed localization statement.

### 2. Conditional-Dispersion Control

The penalty

```text
H_X = E[r_N] E[h_X] - E[r_N h_X]
```

measures the loss from comparing consecutive signed ratios under different
conditional laws. This term is absent in the deterministic one-sided endpoint
and should be small in a perturbative near-one-sided regime.

A useful near-one-sided lemma would be:

> If `Var Y <= epsilon (Var X + Var Y)` and `V>=1`, then at the signed first
> descent either the X-side boundary terms are zero/negligible and
>
> ```text
> H_X + nu_X E[h_X] + (E[r_N]+nu_X) beta_X <= (1/2) I_X,
> I_X >= c_0/V,
> ```
>
> or the reflected Y-side gives a direct reserve bound.

This would give

```text
Delta >= c_0/(2V)
```

in the near-one-sided regime. The constant need not be sharp.

### 3. Boundary Absorption

The X-side boundary penalties are

```text
nu_X E[h_X],
(E[r_N]+nu_X) beta_X.
```

The reflected Y-side has two boundary penalties:

```text
E[t_Y] beta_Y,
nu_Y (E[k_Y] + beta_Y).
```

These terms only occur at signed support edges, but a proof cannot ignore
them. A clean route would prove a side-selection lemma:

> At the signed first descent, at least one of the X-side or reflected Y-side
> reductions has total boundary penalty at most a fixed fraction of its
> inverse-index gain, unless the raw reserve is already at least `c/V`.

This is narrower than the full signed reserve theorem and directly targets
the boundary issue exposed by the corrected identities.

## Two-Regime Proof Route

The current proof route should split signed laws by side variance.

### Near-One-Sided Regime

Assume, after possibly swapping/reflection,

```text
Var Y <= epsilon V.
```

Use the one-sided theorem as the endpoint. The new task is a stability lemma:
small reflected variance should not create enough conditional-dispersion or
boundary loss to erase the one-sided `1/(4V)` effective drop.

A conservative target is:

> For some absolute `epsilon>0` and `c>0`, if `Var Y <= epsilon V` and `V>=1`,
> then
>
> ```text
> Delta >= c/V.
> ```

The conditional reduction suggests proving this by bounding
`H_X + nu_X E[h_X] + (E[r_N]+nu_X) beta_X` relative to `I_X`.

### Genuinely Two-Sided Regime

Assume

```text
Var X >= epsilon V,
Var Y >= epsilon V.
```

The optimizer data show a larger reserve in this regime, but the proof should
not rely on the data. The likely lemma is a two-sided smoothing statement:

> At the signed first descent, the X-side and reflected Y-side reductions
> cannot both have large dispersion-plus-boundary losses relative to their
> inverse-index gains.

Equivalently, prove that at least one side satisfies

```text
I_side - H_side - Boundary_side >= c(epsilon)/V.
```

This would turn the algebraic reduction directly into the signed effective
drop.

## Best Next Lemma Target

The most focused next lemma is the following conditional Newton-drop lemma.

> **Signed conditional Newton-drop lemma.** Let `X` and `Y` be independent
> low-probability Poisson-binomial sums, `V=Var X+Var Y>=1`, and let `D` be
> the first strict descent of `c_z=P(X-Y=z)`. With the X-side quantities
> `I_X,H_X,nu_X,beta_X` and the reflected Y-side quantities
> `I_Y,H_Y,beta_Y,nu_Y` defined above, prove that at least one side has
>
> ```text
> I_side - H_side - Boundary_side >= c/V
> ```
>
> for an absolute constant `c>0`.

This lemma is still essentially the signed reserve problem, but it has a
smaller attack surface: all terms are explicit conditional expectations and
boundary masses. For a first proof attempt, split it into the two smaller
lemmas:

1. **Perturbative near-one-sided lemma.** Prove the displayed inequality when
   `min(Var X, Var Y) <= epsilon V`, using the one-sided `1/(4V)` theorem as
   the endpoint.
2. **Balanced smoothing lemma.** Prove the displayed inequality when both side
   variances are at least `epsilon V`.

The perturbative lemma is the better immediate target because it is anchored
to the proved one-sided result and to the observed sparse boundary. The
balanced lemma is probably cleaner after the dispersion term `H_side` is
better understood.

## Overclaim Guard

This note proves only the conditional reductions and the conditional
implications stated explicitly above. It does not prove:

```text
Delta >= 1/(4V)
```

for signed laws, the raw signed reserve theorem, the hub-bouquet perturbation
step, or Erdos 993. Issue #5 should remain open.
