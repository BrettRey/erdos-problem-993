# Fixed-`r` Shift Rule from Path Moments

## Purpose

This records the first-order rule for the stabilized mode shifts in

```text
F_{a,r}(x)=P_r(x)(1+2x)^a.
```

It is the mode-shift half of the fixed-`r` theorem gap.  The hub-on term still
has to be dominated separately.

## Setup

Let

```text
m=(2a+D)/3,
delta=D/3,
q=a mod 3,
D == q mod 3.
```

Let `S` be the independent-set size on the fixed path `P_r`, with probability
proportional to the coefficients of `P_r` at fugacity `1`:

```text
P_r(x)=sum_s p_s x^s,
Pr(S=s)=p_s/P_r(1),
M_1=E[S].
```

## First-Order Ratio

For

```text
lambda_0=F_{m-1}/F_m,
```

the expansion gives

```text
lambda_0 = 1 + L(D)/(a) + O_r(a^-2),
L(D) = (9/2)(delta-M_1) - 3.
```

The other adjacent ratio is the same expansion with `m` replaced by `m+1`,
equivalently `delta` replaced by `delta+1`:

```text
F_m/F_{m+1} = 1 + L(D+3)/(a) + O_r(a^-2).
```

Thus the first-order adjacent mode inequalities

```text
F_{m-1} <= F_m,
F_{m+1} <= F_m
```

are:

```text
L(D) <= 0,
L(D+3) >= 0.
```

Equivalently,

```text
M_1 - 1/3 <= delta <= M_1 + 2/3
```

or

```text
3M_1 - 1 <= D <= 3M_1 + 2.
```

For a fixed residue `q`, the first-order candidate shifts are therefore the
integers

```text
D == q mod 3
3M_1 - 1 <= D <= 3M_1 + 2.
```

The interval has length `3`, so each residue has one candidate except possible
boundary cases where one residue has two.

## Computed Boundary Check

Implemented in:

```text
gpt_attack/fixed_r_huboff_asymptotic.py
```

The script now records:

```text
first_order_shift_candidates
shift_matches_first_order
```

in `results/fixed_r_huboff_asymptotic_r400.json`.

The checked lanes through `r=400` all match the first-order rule exactly, and
each tested residue has a single candidate.

Additional boundary check run:

```bash
PYTHONPATH=gpt_attack python3 - <<'PY'
from fixed_r_huboff_asymptotic import first_order_shift_candidates, raw_moments_at_one
from route2_spider_lane_scan import path_polys
paths = path_polys(2000)
wide=[]
boundary=[]
for r in range(2,2001):
    candidates = first_order_shift_candidates(paths[r])
    if any(len(v) != 1 for v in candidates.values()):
        wide.append((r,candidates,raw_moments_at_one(paths[r])[0]))
    m1=raw_moments_at_one(paths[r])[0]
    if (3*m1).denominator == 1:
        boundary.append((r,m1))
print('wide_count', len(wide))
print('first_wide', wide[:10])
print('boundary_count', len(boundary))
print('first_boundary', boundary[:10])
PY
```

Output:

```text
wide_count 2
first_wide [(2, {0: [3], 1: [1, 4], 2: [2]}, 2/3),
            (3, {0: [3], 1: [4], 2: [2, 5]}, 1)]
boundary_count 2
first_boundary [(2, 2/3), (3, 1)]
```

So, through `r=2000`, the only first-order boundary cases are `r=2` and
`r=3`.  The exact stabilized shifts choose the lower boundary candidate in
those cases:

```text
r=2: shifts {0:3, 1:1, 2:2}
r=3: shifts {0:3, 1:4, 2:2}
```

## Proof Target

The mode-shift theorem can now be split into two cleaner claims:

```text
1. For r>=4, each residue q has a unique integer D == q mod 3 in
   [3M_1-1, 3M_1+2], and this D is the eventual F-mode shift.

2. For r=2,3, the first-order boundary cases are settled by the second-order
   ratio term, choosing the lower candidate listed above.
```

Once the hub-on domination lemma is applied, this gives the full-polynomial
stabilized mode shifts for `S(2^a,r)`.

## Divisibility Lemma for Boundary Cases

The boundary condition is exactly:

```text
3M_1 is an integer.
```

For the path polynomial,

```text
P_r(1)=F_{r+2}.
```

If

```text
A_r=P'_r(1),
```

then differentiating `P_r=P_{r-1}+xP_{r-2}` at `x=1` gives

```text
A_r=A_{r-1}+A_{r-2}+F_r,
A_0=0,
A_1=1.
```

Solving this recurrence gives

```text
A_r = (r L_{r+1}+2F_r)/5.
```

Since

```text
L_{r+1}=F_r+F_{r+2},
```

we have

```text
5A_r = (r+2)F_r + rF_{r+2}.
```

If `3M_1=3A_r/F_{r+2}` is an integer, then `F_{r+2}` divides `3A_r`.
Multiplying by `5` and reducing modulo `F_{r+2}` gives

```text
F_{r+2} | 3(r+2)F_r.
```

But `gcd(F_r,F_{r+2})=gcd(F_r,F_2)=1`, so

```text
F_{r+2} | 3(r+2).
```

For `r>=7`, `F_{r+2}>3(r+2)`, so no boundary is possible.  The remaining
direct checks are:

```text
r=2:  M_1=2/3,  boundary
r=3:  M_1=1,    boundary
r=4:  M_1=5/4,  not boundary
r=5:  M_1=20/13, not boundary
r=6:  M_1=38/21, not boundary
```

Thus, for `r>=4`, each residue has a unique first-order shift candidate.  The
only first-order boundary cases are `r=2` and `r=3`.

## Remaining Gap

The uniqueness of the first-order shift is now proved for `r>=4`.  What
remains for a theorem-ready mode-shift lemma is:

```text
1. write the adjacent-ratio expansion with an explicit O_r(a^-2) remainder,
   so the unique first-order candidate is eventually the true F-mode shift;
2. settle the r=2,3 boundary cases by the second-order term or direct exact
   formulas;
3. apply the existing hub-on domination lemma to transfer the F-mode shift to
   the full polynomial.
```
