# Fixed-`r` Global `F`-Mode Margin Sublemma

## Purpose

The certificate lemma needs a global hub-off margin:

```text
F_m - F_k >= F_m/(M a)    for every k != m.
```

The hub-off scripts directly check only the adjacent margins:

```text
F_m - F_{m-1} >= F_m/(M a),
F_m - F_{m+1} >= F_m/(M a).
```

This note records why those adjacent checks are enough for

```text
F(x)=P_r(x)(1+2x)^a.
```

## Lemma

Let

```text
F(x)=sum_k F_k x^k
```

have nonnegative coefficients and be unimodal.  Suppose `m` satisfies

```text
F_m - F_{m-1} >= eta F_m,
F_m - F_{m+1} >= eta F_m
```

with the convention that missing coefficients are `0`.  If `m` is a mode
position for the unimodal sequence, then for every `k != m`,

```text
F_m - F_k >= eta F_m.
```

### Proof

Since the coefficient sequence is unimodal with mode at `m`,

```text
k < m  =>  F_k <= F_{m-1},
k > m  =>  F_k <= F_{m+1}.
```

Therefore:

```text
k < m  =>  F_m - F_k >= F_m - F_{m-1} >= eta F_m,
k > m  =>  F_m - F_k >= F_m - F_{m+1} >= eta F_m.
```

This proves the global margin.

## Why `F=P_r(1+2x)^a` Is Unimodal

The path independence polynomial `P_r(x)` is real-rooted with nonpositive
roots.  This follows, for example, from the Chudnovsky-Seymour theorem for
independence polynomials of claw-free graphs, since paths are claw-free; it can
also be proved directly from the path recurrence

```text
P_r(x)=P_{r-1}(x)+xP_{r-2}(x).
```

The factor `(1+2x)^a` is also real-rooted.  Hence

```text
F(x)=P_r(x)(1+2x)^a
```

is real-rooted with nonnegative coefficients.  Newton's inequalities imply its
coefficient sequence is log-concave with no internal zeros, hence unimodal.

Thus the adjacent margin checks performed in

```text
gpt_attack/fixed_r_huboff_certificate.py
gpt_attack/fixed_r_huboff_cpp_certificate.py
```

are sufficient to supply the global `F`-margin required by the mode-domination
certificate.

## Caveat

The lemma assumes the adjacent inequalities are checked at the intended mode
index `m=(2a+D_q)/3`.  If either adjacent check fails, then this lemma does not
identify the mode.  The scripts avoid this by proving both shifted adjacent
margin rational functions have positive numerator and denominator coefficients
for all `a >= A` in each residue class.

## What This Closes

Together with

```text
notes/fixed_r_fugacity_shift_sublemma_2026-05-21.md
```

this closes the two previously informal links in

```text
notes/fixed_r_certificate_lemma_2026-05-21.md.
```

The remaining theorem-level gap is no longer the composition of the
certificates.  It is the production of certificate data for arbitrary fixed
`r`, or an analytic replacement for the per-lane hub-off reserve positivity
checks.
