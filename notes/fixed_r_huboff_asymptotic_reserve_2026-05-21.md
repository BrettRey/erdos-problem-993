# Fixed-`r` Hub-Off Reserve Asymptotic

## Purpose

The hub-off reserve certificate is the current hard part of the arbitrary
fixed-`r` theorem.  The GMP checker proves the reserve exactly for individual
lanes, but an arbitrary fixed-`r` theorem needs a conceptual lower bound.

This note records the first asymptotic extraction.

## Object

For

```text
F_a(x)=P_r(x)(1+2x)^a,
m=2a/3 + delta,
lambda_0=F_{m-1}/F_m,
```

the hub-off reserve is

```text
R_{a,r,q}
 = mu_{P_r(1+2x)^(a-1)}(lambda_0) - m + 4/3.
```

The certificate needs:

```text
R_{a,r,q} >= 1/(1000a)
```

for all sufficiently large `a` in each residue class.

## Expansion Target

The expansion script computes:

```text
R_{a,r,q} = C_{r,q}/a + O_r(a^-2),
```

where

```text
C_{r,q} = lim_{a -> infinity, a=q mod 3} a R_{a,r,q}.
```

The useful theorem target is now:

```text
For every fixed r and every residue q, the stabilized shift D_{r,q} gives
C_{r,q} > 0.
```

Once this is known, the hub-off reserve follows for all sufficiently large
`a`; the remaining problem is making the threshold effective.

## Script

Script:

```bash
python3 gpt_attack/fixed_r_huboff_asymptotic.py \
  --r-values 8,20,80,120,160,200,240,280,320,400 \
  --threshold 650 \
  --out results/fixed_r_huboff_asymptotic_r400.json
```

The script works with `h=1/a` and expands the fixed-path convolution term by
term.  It avoids the dense rational functions in the residue parameter `t`.

For each path coefficient index `s`, it expands

```text
2^-u binom(a,m-u)/binom(a,m)
```

as a truncated series in `h`, where `u=s` for `F_m` and `u=s+1` for
`F_{m-1}`.  This gives `lambda_0` to second order in `h`, enough to extract
`C_{r,q}`.

## Closed Moment Formula

Write

```text
P_r(x)=sum_s p_s x^s,
Pr(S=s)=p_s/P_r(1),
M_j=E[S^j],
delta=D_{r,q}/3.
```

The script now checks, exactly, that the extracted coefficient equals:

```text
C_{r,q}
 = (-24 delta + 54 M_1^3 - 9 M_1^2 - 81 M_1 M_2
    + 24 M_1 + 9 M_2 + 27 M_3 + 16)/12.
```

Equivalently, if

```text
V = Var(S),
K_3 = E[(S-M_1)^3],
```

then

```text
C_{r,q} = (9V + 27K_3 + 24(M_1-delta) + 16)/12.
```

This is the first real compression of the hub-off reserve problem: the leading
constant depends only on the first three moments of the fixed path polynomial
and the stabilized mode shift.

The same script now also records the first-order mode-shift candidates from
`M_1`; see `notes/fixed_r_shift_rule_from_path_moments_2026-05-21.md`.

The positivity problem is proved further in
`notes/fixed_r_c_positivity_reduction_2026-05-21.md`: with
`alpha=M_1-D_{r,q}/3`,

```text
C_{r,q} = (9(V+3K_3) + 24alpha + 16)/12.
```

The shift interval gives `alpha>-2/3` for `r>=4`, so it is enough to prove
`V+3K_3>=0` for the path-size distribution.  The note proves this inequality
by an exact Fibonacci/Lucas calculation.

## Internal Cancellation Check

For every tested lane, the script also verifies the expected cancellations:

```text
coefficient of a in R_{a,r,q}:      0
constant coefficient in R_{a,r,q}:  0
```

Thus the first nonzero term is indeed the `C_{r,q}/a` term.

## Results

Minimum and maximum `C_{r,q}` across the three residue classes:

```text
r=8:    min 0.366240921613, max 1.699574254946
r=20:   min 0.666288254007, max 1.999621587341
r=80:   min 2.776776058863, max 4.110109392197
r=120:  min 3.961544868065, max 5.294878201398
r=160:  min 4.479647010600, max 5.812980343933
r=200:  min 5.664415819801, max 6.997749153135
r=240:  min 6.849184629003, max 8.182517962336
r=280:  min 8.033953438204, max 9.367286771538
r=320:  min 9.218722247406, max 10.552055580739
r=400:  min 10.921593199143, max 12.254926532476
```

All tested constants are comfortably positive.  The minimum grows with `r` in
the tested range except for small residue-pattern oscillations.

## First-Order Fugacity Term

Let

```text
mu_r = P'_r(1)/P_r(1),
delta = D_{r,q}/3.
```

The expansion gives:

```text
lambda_0 = 1 + L_{r,q}/a + O_r(a^-2),
L_{r,q} = (9/2)(delta - mu_r) - 3.
```

This first-order term explains the cancellation of the constant part of the
reserve:

```text
mu_r
+ (a-1) * 2lambda_0/(1+2lambda_0)
- m
+ 4/3
```

has no `a^0` term after substituting `L_{r,q}`.

The positivity question therefore lives at the next coefficient, `C_{r,q}`.

## Theorem Direction

The next analytic lemma should prove positivity of the closed moment formula
for the stabilized shifts.  A practical statement would be:

```text
For each fixed r and residue q, if D_{r,q} is the stabilized shift produced by
the hub-off mode equation, then C_{r,q}>0.
```

A stronger useful statement would control it uniformly enough to give an
effective threshold:

```text
R_{a,r,q} >= C_{r,q}/(2a)
```

for all `a >= A_reserve(r)`.

## Consequence for the Current Proof Plan

This reframes the hub-off reserve problem.  The GMP positivity checker is
proving a finite-threshold version of a phenomenon whose leading asymptotic
constant is already positive and growing in the tested lanes.

The remaining arbitrary fixed-`r` gap is:

```text
1. characterize the stabilized shifts D_{r,q};
2. prove the closed moment formula gives C_{r,q}>0 for those shifts;
3. make the O_r(a^-2) remainder effective.
```

The first item has been reduced to the first-order interval
`3M_1-1 <= D <= 3M_1+2`, plus boundary handling for `r=2,3`.
The second item is now reduced to, and proved by, the path moment inequality
`V+3K_3>=0`.

Lane-by-lane GMP certificates remain useful as stress tests, but the proof
target is now the effective remainder bound, plus the final mode-shift
remainder needed to justify the stabilized shifts from the first-order rule.
