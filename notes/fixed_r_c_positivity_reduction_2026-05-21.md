# Fixed-`r` Reduction for Positivity of `C_{r,q}`

## Purpose

The hub-off reserve asymptotic has leading constant

```text
R_{a,r,q}=C_{r,q}/a+O_r(a^-2).
```

This note records a sharper positivity target for `C_{r,q}`.  It reduces the
problem to a path-distribution moment inequality.

## Moment Formula

Let

```text
P_r(x)=sum_s p_s x^s,
Pr(S=s)=p_s/P_r(1),
M=E[S],
V=Var(S),
K_3=E[(S-M)^3],
delta=D_{r,q}/3,
alpha=M-delta.
```

The extracted formula is:

```text
C_{r,q} = (9V + 27K_3 + 24alpha + 16)/12
        = (9(V+3K_3) + 24alpha + 16)/12.
```

## Shift-Interval Consequence

The first-order shift rule gives

```text
M - 1/3 <= delta <= M + 2/3,
```

so

```text
-2/3 <= alpha <= 1/3.
```

For `r>=4`, the first-order candidate is unique, so the lower endpoint is not
attained:

```text
alpha > -2/3.
```

For the two boundary cases:

```text
r=2: actual shifts {0:3, 1:1, 2:2}
r=3: actual shifts {0:3, 1:4, 2:2}
```

the exact shift chooses the lower boundary candidate, giving `alpha=1/3` in
the ambiguous residue rather than `alpha=-2/3`.

Therefore the following path moment inequality is enough to prove
`C_{r,q}>0` for all stabilized shifts:

```text
V + 3K_3 >= 0.
```

For `r>=4`, strict positivity of `C_{r,q}` follows from
`24alpha+16>0`.  For `r=2,3`, it follows directly from the listed exact
shifts.

## Verification

The reduction was checked through `r=5000` with the actual boundary choices:

```bash
PYTHONPATH=gpt_attack python3 - <<'PY'
from fixed_r_huboff_asymptotic import central_moments_at_one, first_order_shift_candidates, reserve_constant_moment_formula
from route2_spider_lane_scan import path_polys
paths=path_polys(5000)
min_b=None
min_c=None
for r in range(2,5001):
    _,v,k=central_moments_at_one(paths[r])
    b=v+3*k
    if min_b is None or b < min_b[0]:
        min_b=(b,r)
    if r == 2:
        shifts={0:3,1:1,2:2}
    elif r == 3:
        shifts={0:3,1:4,2:2}
    else:
        shifts={q: ds[0] for q, ds in first_order_shift_candidates(paths[r]).items()}
    for q,d in shifts.items():
        c=reserve_constant_moment_formula(paths[r],d)
        if min_c is None or c < min_c[0]:
            min_c=(c,r,q,d)
print('min_V_plus_3K3', min_b[0], 'at r', min_b[1], 'float', float(min_b[0]))
print('min_C_actual_shift_rule', min_c[0], 'at r,q,D', min_c[1:], 'float', float(min_c[0]))
PY
```

Output:

```text
min_V_plus_3K3 0 at r 2 float 0.0
min_C_actual_shift_rule 7312/19965 at r,q,D (8, 0, 9) float 0.36624092161282246
```

## Real-Rooted Reformulation

Since `P_r` has only real negative roots, the size distribution at fugacity
`1` is a sum of independent Bernoulli variables.  The path roots are

```text
x_j = -1/(4 cos^2(j pi/(r+2))),
j=1,...,floor((r+1)/2).
```

Set

```text
y_j = 4 cos^2(j pi/(r+2)),
p_j = y_j/(1+y_j).
```

For a Bernoulli variable with parameter `p`,

```text
Var + 3 third_central = p(1-p)(4-6p).
```

Thus the needed path inequality becomes the finite trigonometric sum

```text
sum_j 2 y_j(2-y_j)/(1+y_j)^3 >= 0.
```

This is now the cleanest analytic target for positivity of the hub-off reserve
constant.

## Proof of the Path Moment Inequality

The inequality `V+3K_3>=0` also has a direct Fibonacci proof.

Let

```text
T_r=P_r(1),
S_j(r)=sum_s s^j p_s.
```

Thus

```text
M=S_1/T_r,
V+3K_3
 = S_2/T_r - (S_1/T_r)^2
   + 3(S_3/T_r - 3S_1S_2/T_r^2 + 2S_1^3/T_r^3).
```

The path recurrence gives

```text
T_r=T_{r-1}+T_{r-2},
S_1(r)=S_1(r-1)+S_1(r-2)+T_{r-2},
S_2(r)=S_2(r-1)+S_2(r-2)+2S_1(r-2)+T_{r-2},
S_3(r)=S_3(r-1)+S_3(r-2)+3S_2(r-2)+3S_1(r-2)+T_{r-2}.
```

Solving these recurrences gives:

```text
T_r = F_{r+2},

S_1(r) = (r/2+2/5)F_r + (r/10)L_r,

S_2(r) = (r^2/5+r/2+8/25)F_r - (r/50)L_r,

S_3(r) = (r^3/10+21r^2/50+23r/50+4/25)F_r
         -(r^3/50+3r^2/50+3r/50)L_r.
```

Substituting and using

```text
F_{r+2}=(3F_r+L_r)/2,
L_r^2-5F_r^2=4(-1)^r
```

gives

```text
250 F_{r+2}^3 (V+3K_3)
 = (80r+304)F_r^3 + (36r+136)F_r^2L_r
   - (-1)^r[
       (15r^3+51r^2-22r-200)F_r
       +(9r^3+41r^2+50r)L_r
     ].
```

For odd `r>=3`, the bracket is added rather than subtracted, and both bracket
polynomials are positive.  Hence the expression is positive.

For even `r>=6`, use `0<L_r<3F_r` and ignore the positive
`(36r+136)F_r^2L_r` term:

```text
250 F_{r+2}^3 (V+3K_3)
 >= F_r[(80r+304)F_r^2 - (42r^3+174r^2+128r-200)].
```

Since `F_r>=r` for `r>=6`,

```text
(80r+304)F_r^2 - (42r^3+174r^2+128r-200)
 >= 38r^3+130r^2-128r+200 > 0.
```

The remaining small cases are:

```text
r=0: V+3K_3=0
r=1: V+3K_3=1/4
r=2: V+3K_3=0
r=4: V+3K_3=5/32
```

Therefore

```text
V+3K_3 >= 0
```

for every path length `r`, with equality only at `r=0,2`.

## Consequence

The leading hub-off reserve constant is now positive for the first-order
stabilized-shift candidates:

```text
C_{r,q} = (9(V+3K_3) + 24alpha + 16)/12 > 0.
```

For `r>=4`, the shift interval is strict at the lower end, so
`24alpha+16>0`.  For `r=2,3`, the exact boundary choices listed above give
positive `C_{r,q}` directly.

The remaining hub-off asymptotic work is no longer positivity of the leading
constant.  It is making the `O_r(a^-2)` remainder effective.
