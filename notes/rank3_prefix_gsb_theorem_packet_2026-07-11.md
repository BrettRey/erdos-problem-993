# Rank-three prefix-GSB theorem packet

## Frozen statement

Let `T` be a finite tree with independence number `alpha>=7`, and let `i_k`
denote the number of independent `k`-sets. Then

```text
5 i_3 i_5 <= 4 i_4^2 + i_3 i_4.                 (R3-GSB)
```

Equivalently, for `mu_r=(r+1)i_(r+1)/i_r`,

```text
mu_4 <= mu_3+1.
```

This is exactly prefix GSB at rank `r=3`. Indeed,
`3<=ceil((2alpha-1)/3)-2` is equivalent to `alpha>=7`.

The replayable exact certificate is
`scratch_rank3_gsb_certificate_20260711.py`. It uses only integer/rational
arithmetic and symbolic polynomial identities. A full run prints
`symbolic_certificates: passed` and the finite-base table below.

## Exact coefficient formulas

Put `t=n-1`, `d_v=deg(v)`, and

```text
S_j  = sum_v d_v^j,
P    = sum_(uv in E) d_u d_v,
P21  = sum_(uv in E) d_u d_v(d_u+d_v),
L    = sum_v (sum_(u~v) d_u)^2.
```

Inclusion--exclusion over edge subsets gives

```text
i_3 = (3S_2+t^3-6t^2-t)/6,

i_4 = (-24P+12tS_2-4S_3+t^4-14t^3+23t^2-2t)/24,

i_5 = (60L-120tP+60P21
       +30t^2S_2-90tS_2+55S_2-20tS_3-30S_3+5S_4
       +t^5-25t^4+125t^3-115t^2-6t)/120.
```

The certificate replays all three formulas against the repository's exact
tree DP for every nonisomorphic tree through order 14.

Write

```text
Delta = 4i_4^2+i_3i_4-5i_3i_5.
```

The proof divides according to `h`, the number of nonleaf vertices.

## Cores with at most three nonleaf vertices

### `h=1`

The tree is `K_(1,t)`. For `k>=2`, `i_k=C(t,k)`, and directly

```text
Delta = t^2(t-3)(t-2)^2(t-1)^2/72 > 0
```

in the prefix range `t=alpha>=7`.

### `h=2`

The nonleaf core is one edge. Let its endpoint excess degrees be positive
integers `a,b`, put `m=a+b` and `p=ab`. The tree has matching number two, so
`alpha=m`; hence the prefix hypothesis is `m>=7`. Also

```text
m-1 <= p <= m^2/4.
```

As a cubic polynomial in `p`, `72Delta` has the following degree-three
Bernstein coefficients on that interval:

```text
B0=(m-3)(m-2)^2(m-1)^2(m+2)(m+6),
B1=(m-2)^2(m-1)(13m^4+18m^3-158m^2-109m+114)/12,
B2=m(m-2)^2(51m^4-72m^3-355m^2+120m+172)/48,
B3=m^2(m-2)^2(30m^3-73m^2-94m-16)/32.
```

Each is positive for `m>=7`, proving `Delta>0`.

### `h=3`

The nonleaf core is a path. Let its endpoint excess degrees be `a,b>=1`,
let the central excess degree be `c>=1`, and put `m=a+b`, `p=ab`. Again
`m-1<=p<=m^2/4`.

- If `c=1`, the matching number is two and `alpha=m+1`, so the prefix starts
  at `m>=6`.
- If `c>=2`, the matching number is three and `alpha=m+c-1`, so the prefix
  condition is `m+c>=8`.

The exact cubic `72Delta` in `p` has nonnegative Bernstein coefficients on
the whole interval in each of the following shifted nonnegative orthants:

```text
(m,c)>=(6,1), (5,2), (4,4), (3,5), (2,6).
```

Their union contains every integer prefix case just listed. The script
constructs these four Bernstein coefficients symbolically and verifies that,
after each shift, every coefficient as a polynomial in the two nonnegative
shift variables is nonnegative. This step proves slightly more than the
needed prefix statement; in particular `(5,2)` safely includes one extra
non-prefix boundary row.

## At least four nonleaf vertices, large order

Assume `h>=4` and `t>=14`. Put

```text
x_v=d_v-1,                s=sum_v x_v=t-1,
Q=sum_v x_v^2,            R=sum_v x_v^3,
U=sum_v x_v^4,            C=sum_(uv in E) x_u x_v.
```

The positive `x_v` are precisely the weights on the nonleaf core.

### Joint radius-two bound

Let

```text
a_v=sum_(u~v)x_u,
D=sum_(uv in E)x_u x_v(x_u+x_v),
Z=sum_v sum_({u,w} subset N(v)) x_u x_w.
```

Exact expansion gives

```text
L+P21 = 2R+8Q+8C+D+2Z+10t-6.
```

Since a tree gives each distance-two pair a unique middle vertex,

```text
D+2Z <= sum_(u<v)x_u x_v(x_u+x_v)=sQ-R.
```

The inequality uses `2x_ux_v<=x_ux_v(x_u+x_v)` for positive weights at
distance two; adjacent pairs contribute exactly `D`, and all omitted pairs
are nonnegative. Therefore

```text
L+P21 <= R+(t+7)Q+8C+10t-6.                    (J)
```

### Remaining moment bounds

All directions and hypotheses used below are:

```text
Q <= s^2,                                      (sum x_v=s)
U <= (t-4)R,                                   (h>=4)
C >= t-2,                                      (h>=2)
R >= Q^2/(t-1).                                (Cauchy--Schwarz)
```

For the second inequality, every positive weight is at most
`s-(h-1)=t-h<=t-4`, so `x_v^4<=(t-4)x_v^3`.

For the third, let `K` be the nonleaf core. On every core edge,
`xy>=x+y-1`; hence

```text
C >= sum_v deg_K(v)x_v-(h-1)
  = s+sum_v(deg_K(v)-1)x_v-(h-1)
  >= s+(h-2)-(h-1)=t-2.
```

The last inequality is `(sum x_v^2)^2<=(sum x_v)(sum x_v^3)`.

### Directions of substitution

After the degree formulas are rewritten in `Q,R,U,C,L+P21`, the coefficient
of `L+P21` in `144Delta` is

```text
-60A,
```

and the coefficient of `U` is

```text
-5A,
```

where

```text
A=3Q+t^3-6t^2+8t-3=6i_3>0.
```

Thus the upper bounds `(J)` and `U<=(t-4)R` give a lower bound `E_t(C,Q,R)`
for `144Delta`; the substitution directions are correct because both
coefficients are negative.

The exact derivative identities are

```text
dE/dC = 24[(t-1)(3t^3-29t^2+46t-8)
           +48C+8R+9(t-1)((t-1)^2-Q)],

dE/dR = (t-1)(7t^3-56t^2+84t-32)
        +192C+32R+(51t-186)((t-1)^2-Q).
```

Every displayed summand is nonnegative for `t>=14`, and the first terms are
strictly positive. Consequently `E_t` is increasing in both `C` and `R` on
the relaxed moment domain. Using the lower bounds on `C,R` therefore gives

```text
144Delta >= E_t(t-2,Q,Q^2/(t-1)) = N_t(Q)/(t-1)^2,
```

where

```text
N_t(Q)=16Q^4+(-51t^2+237t-186)Q^3
 +(7t^5+35t^4-442t^3+1040t^2-1121t+481)Q^2
 +(-9t^7+90t^6-439t^5+1262t^4-1740t^3+877t^2+136t-177)Q
 +4t^9-87t^8+805t^7-3998t^6+11458t^5-18977t^4
  +16209t^3-3412t^2-4172t+2170.
```

### Positivity of the residual quartic

The global minimum of `N_t''(Q)` is

```text
m_t=(t-1)(896t^4-2427t^3+13519t^2-78784t+42220)/64 > 0.
```

At `t=u+14`, the quartic factor is

```text
896u^4+47749u^3+965281u^2+8707168u+29350016,
```

so the curvature bound is strict. Set

```text
Q0=(t-1)(9t-50)/14.
```

Strong convexity gives, for every real `Q`,

```text
N_t(Q) >= N_t(Q0)+N_t'(Q0)(Q-Q0)+m_t(Q-Q0)^2/2
       >= N_t(Q0)-N_t'(Q0)^2/(2m_t).
```

The last expression has denominator

```text
941192(896t^4-2427t^3+13519t^2-78784t+42220)
```

and its numerator, after `t=u+14`, is

```text
933662464u^13+145897280186u^12+10411810306803u^11
+448767361382575u^10+13013050452178914u^9
+267481862079265246u^8+3996891304693378827u^7
+43754278751461351515u^6+348347452443392936972u^5
+1967759034918318346702u^4+7491802131968761085156u^3
+17295916918434354058336u^2+18486494986684683921664u
+675585773681130314496.
```

It is strictly positive for `u>=0`. Hence `Delta>0` whenever `h>=4` and
`t>=14`.

## Exact finite base

It remains to check `h>=4`, `alpha>=7`, and `t<14`, i.e. orders below 15.
There are no such cases at orders 8 or 9: `alpha>=7` there forces matching
number at most two, already covered by `h<=3`. Exact nonisomorphic-tree
enumeration gives:

| `n` | rows checked | minimum `Delta` | graph6 attaining minimum |
|---:|---:|---:|---|
| 10 | 22 | 5415 | `I???F?[c_` |
| 11 | 138 | 16318 | `J???EA_U@K?` |
| 12 | 476 | 42161 | `K??????w?{QK` |
| 13 | 1270 | 94984 | `L???????F_@waK` |
| 14 | 3122 | 197692 | `M?????????]?BoaM?` |

All comparisons are exact integers. Together with the symbolic large-order
certificate and the complete `h<=3` classifications, this finite base is
exhaustive rather than a bounded-order heuristic.

## Audit boundary

This packet proves only rank-three prefix GSB. Combined with the already
proved rank-one/rank-two bases, it moves any failure of prefix GSB to rank at
least four. It does not prove prefix GSB at arbitrary rank and does not by
itself resolve Erdős Problem 993.
