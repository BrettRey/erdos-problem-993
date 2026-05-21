# Fixed-`r` Route-2 Proof Note Draft

## Purpose and Scope

This is a new proof-note draft.  It is not a manuscript edit.

The aim is to state the fixed-`r` spider-lane result in a form that separates:

```text
1. what is already fully certificate-backed for sampled lanes;
2. what is analytic and applies to arbitrary fixed r;
3. what remains a certificate schema rather than a manuscript theorem.
```

No claim is made here for the full `d_leaf<=1` problem.

## Notation

Let `P_j(x)` denote the independence polynomial of the path on `j` vertices.
For fixed `r>=2`, define

```text
T_{a,r}=S(2^a,r),
```

the spider with `a` arms of length `2` and one arm of length `r`.

After removing one length-`2` arm, define

```text
B_{a,r}=S(2^(a-1),r).
```

Write the independence polynomial of `T_{a,r}` as

```text
I(T_{a,r};x)=F_{a,r}(x)+G_{a,r}(x),
F_{a,r}(x)=P_r(x)(1+2x)^a,
G_{a,r}(x)=xP_{r-1}(x)(1+x)^a.
```

Similarly,

```text
I(B_{a,r};x)=F^-_{a,r}(x)+G^-_{a,r}(x),
F^-_{a,r}(x)=P_r(x)(1+2x)^(a-1),
G^-_{a,r}(x)=xP_{r-1}(x)(1+x)^(a-1).
```

For a polynomial `H(x)=sum_k h_k x^k`, write

```text
mu_H(lambda)=lambda H'(lambda)/H(lambda).
```

Let `m=m(a,r)` be the leftmost mode of `I(T_{a,r};x)` and set

```text
lambda=i_{m-1}(T_{a,r})/i_m(T_{a,r}).
```

The fixed-`r` Route-2 target is

```text
mu_{B_{a,r}}(lambda) >= m - 3/2.
```

## Theorem A: Sampled Fixed-`r` Certificate Proposition

For

```text
r in {4,8,12,16,20,24,32,40,60,80},
```

the Route-2 target holds for every `a>=1`:

```text
mu_{B_{a,r}}(lambda) >= m - 3/2.
```

### Proof Status

This is a certificate-backed proposition.

The proof is split at `a=200`:

```text
finite range:       1 <= a <= 199;
asymptotic range:   a >= 200.
```

For the finite range, exact integer arithmetic verifies the Route-2 inequality
after clearing denominators.

For the asymptotic range, the certificate supplies:

```text
1. a stabilized residue-class mode formula;
2. hub-off mode margins and hub-off reserve;
3. hub-on mode domination;
4. hub-on Route-2 perturbation domination.
```

After correcting the hub-on monotonicity condition to check decrease of
`a * perturbation`, the sampled suite was rerun and still reports:

```text
all_ok=True.
```

This theorem is the cleanest fully supported output at present.

## Certificate Lemma

Fix `r` and a threshold `A`.  Suppose that for each residue class
`q in {0,1,2}` the certificate supplies an integer `D_q` such that

```text
a == q mod 3,
m=(2a+D_q)/3.
```

Assume the following four conditions.

### C0. Finite Verification

For every `1<=a<A`, with `m` the actual leftmost mode of `I(T_{a,r};x)`,

```text
mu_{B_{a,r}}(lambda) >= m - 3/2.
```

### C1. Mode Verification

For every `a>=A`, the certified value

```text
m=(2a+D_q)/3
```

is the leftmost mode of `I(T_{a,r};x)`.

A sufficient certificate has two parts:

```text
F_{a,r,m}-F_{a,r,k} >= F_{a,r,m} eta_mode/a     for every k != m,
2 max_k G_{a,r,k} < F_{a,r,m} eta_mode/a.
```

The first line follows from adjacent `F`-mode margins plus the
real-rooted/log-concave coefficient sequence of `P_r(x)(1+2x)^a`.

### C2. Hub-Off Reserve

Let

```text
lambda_0=F_{a,r,m-1}/F_{a,r,m}.
```

For every `a>=A`,

```text
3/4 <= lambda_0 <= 2
```

and

```text
mu_{F^-_{a,r}}(lambda_0) >= m - 4/3 + eta_reserve/a.
```

### C3. Hub-On Perturbation

For every `a>=A`,

```text
|mu_{F^-_{a,r}+G^-_{a,r}}(lambda)
  - mu_{F^-_{a,r}}(lambda_0)|
<= eta_reserve/a.
```

### Conclusion

Under C0-C3, the Route-2 target holds for every `a>=1`.

For `a<A`, this is C0.  For `a>=A`, C1 makes `m` the actual leftmost mode and
therefore makes `lambda` the correct tie fugacity.  Then C2 and C3 give

```text
mu_{B_{a,r}}(lambda)
 >= m - 4/3
 >  m - 3/2.
```

## Analytic Lemma 1: First-Order Shift Rule

Let

```text
m=(2a+D)/3,
delta=D/3,
q=a mod 3,
D == q mod 3.
```

Let `S` be the independent-set size on the fixed path `P_r`, distributed by
the coefficients of `P_r` at fugacity `1`, and let

```text
M_1=E[S]=P'_r(1)/P_r(1).
```

For

```text
lambda_0=F_{m-1}/F_m,
```

we have the expansion

```text
lambda_0 = 1 + L(D)/a + O_r(a^-2),
L(D)=(9/2)(D/3-M_1)-3.
```

Similarly,

```text
F_m/F_{m+1}=1+L(D+3)/a+O_r(a^-2).
```

The first-order adjacent-mode conditions are therefore

```text
L(D)<=0,
L(D+3)>=0,
```

equivalently

```text
3M_1-1 <= D <= 3M_1+2,
D == q mod 3.
```

For `r>=4`, each residue class has a unique integer `D` in this interval.

Indeed,

```text
P_r(1)=F_{r+2},
P'_r(1)=(rL_{r+1}+2F_r)/5.
```

If `3M_1` were an integer, then `F_{r+2}` would divide `3(r+2)`.  This is
impossible for `r>=7`, and `r=4,5,6` are checked directly.  Thus the only
first-order boundary cases are `r=2,3`.

Boundary cases:

```text
r=2: shifts {0:3, 1:1, 2:2}
r=3: shifts {0:3, 1:4, 2:2}
```

For a manuscript proof, `r=2,3` should be handled as explicit small cases or
by separate shifted-positive certificates.

## Analytic Lemma 2: Hub-Off Reserve Constant

Define the hub-off reserve

```text
R_{a,r,q}=mu_{F^-_{a,r}}(lambda_0)-m+4/3.
```

For fixed `r` and residue `q`,

```text
R_{a,r,q}=C_{r,q}/a+O_r(a^-2).
```

Let

```text
M=E[S],
V=Var(S),
K_3=E[(S-M)^3],
alpha=M-D/3.
```

Then

```text
C_{r,q}=(9(V+3K_3)+24alpha+16)/12.
```

By the first-order shift interval,

```text
alpha >= -2/3,
```

and for `r>=4` the lower endpoint is not attained.  The `r=2,3` boundary
choices give positive constants directly.

Thus positivity of `C_{r,q}` reduces to the path moment inequality

```text
V+3K_3 >= 0.
```

## Analytic Lemma 3: Path Moment Positivity

Let

```text
P_r(x)=sum_s p_s x^s,
T_r=P_r(1),
S_j(r)=sum_s s^j p_s.
```

The path recurrence gives

```text
T_r=T_{r-1}+T_{r-2},
S_1(r)=S_1(r-1)+S_1(r-2)+T_{r-2},
S_2(r)=S_2(r-1)+S_2(r-2)+2S_1(r-2)+T_{r-2},
S_3(r)=S_3(r-1)+S_3(r-2)+3S_2(r-2)+3S_1(r-2)+T_{r-2}.
```

Solving these recurrences gives

```text
T_r = F_{r+2},

S_1(r) = (r/2+2/5)F_r + (r/10)L_r,

S_2(r) = (r^2/5+r/2+8/25)F_r - (r/50)L_r,

S_3(r) = (r^3/10+21r^2/50+23r/50+4/25)F_r
         -(r^3/50+3r^2/50+3r/50)L_r.
```

Substituting into `V+3K_3` and using

```text
L_r^2-5F_r^2=4(-1)^r
```

yields

```text
250F_{r+2}^3(V+3K_3)
 = (80r+304)F_r^3 + (36r+136)F_r^2L_r
   - (-1)^r[
       (15r^3+51r^2-22r-200)F_r
       +(9r^3+41r^2+50r)L_r
     ].
```

For odd `r>=3`, the bracket is added and all terms are positive.

For even `r>=6`, use `0<L_r<3F_r`:

```text
250F_{r+2}^3(V+3K_3)
 >= F_r[(80r+304)F_r^2-(42r^3+174r^2+128r-200)].
```

Since `F_r>=r` for `r>=6`,

```text
(80r+304)F_r^2-(42r^3+174r^2+128r-200)
 >= 38r^3+130r^2-128r+200 > 0.
```

The remaining small cases are direct:

```text
r=0: V+3K_3=0
r=1: V+3K_3=1/4
r=2: V+3K_3=0
r=4: V+3K_3=5/32
```

Therefore

```text
V+3K_3>=0
```

for every path length, with equality only at `r=0,2`.

Consequently, the hub-off leading constant `C_{r,q}` is positive for the
stabilized fixed-`r` shifts.

## Analytic Lemma 4: Shifted-Positive Threshold Termination

Let `Q(t)` be a nonzero polynomial with positive leading coefficient.  Then
there exists `T` such that all coefficients of

```text
Q(u+T)
```

are positive.

Reason: the coefficient of `u^j` in `Q(u+T)` is `Q^{(j)}(T)/j!`.  Every
nonzero derivative of `Q` has positive leading coefficient after finitely many
shifts.

For a rational function `N(t)/D(t)` with positive eventual sign, apply the
same argument to numerator and denominator after normalizing the leading
coefficient of the denominator to be positive.

This gives a terminating certificate procedure for eventual positivity.

## Analytic Lemma 5: Hub-On Perturbation Bound

Let

```text
lambda=(F_{m-1}+G_{m-1})/(F_m+G_m).
```

Assume

```text
3/4 <= lambda_0 <= 2,
T >= max(G_{m-1},G_m)/F_m,
T <= 1/2.
```

Then

```text
|lambda-lambda_0| <= 4T,
lambda >= 1/2.
```

For a nonnegative polynomial `H` of degree at most `N`,

```text
|mu_H(lambda)-mu_H(lambda_0)| <= 2N^2T.
```

For the bridge mixture,

```text
|mu_{F^-+G^-}(lambda)-mu_{F^-}(lambda)|
 <= (a+r)C_r(3/4)^(a-1).
```

Therefore

```text
|mu_{B_{a,r}}(lambda)-mu_{F^-_{a,r}}(lambda_0)|
 <= 2(a+r)^2T + (a+r)C_r(3/4)^(a-1).
```

When this bound is compared with a reserve of the form `eta/a`, the
monotonicity check must be applied to

```text
a * [2(a+r)^2T + (a+r)C_r(3/4)^(a-1)],
```

not merely to the unweighted perturbation.  This is the corrected
monotonicity condition.

## Theorem B: Fixed-`r` Computable-Threshold Criterion

Fix `r>=4`.  Suppose the following objects are produced:

```text
1. the unique first-order shift D_q in each residue class;
2. a shifted-positive threshold for adjacent F-mode margins;
3. a shifted-positive threshold for the hub-off reserve
   R_{a,r,q} >= eta_{r,q}/a;
4. a shifted-positive or analytic threshold for hub-on mode domination;
5. a shifted-positive or analytic threshold for hub-on perturbation domination.
```

Then there is a computable threshold `A(r)` such that the Route-2 target holds
for all `a>=A(r)`.

If the finite range `1<=a<A(r)` is checked exactly, then the Route-2 target
holds for all `a>=1` in that fixed lane.

### Current Effective Threshold Data

For the hub-off/mode/lambda pieces, shifted-coefficient probing with
`eta=C_{r,q}/2` gives:

```text
r=4:   A=12
r=8:   A=24
r=20:  A=94
r=24:  A=22
r=40:  A=26
r=80:  A=44
r=120: A=62
```

For effective hub-on thresholds, using half the actual mode-margin constants
and half `C_{r,q}`, step-5 grids give:

```text
r=80:  A=145
r=120: A=200
r=160: A=255
r=200: A=305
r=240: A=360
r=280: A=415
r=320: A=470
r=400: A=575
```

Thus the practical threshold driver is currently hub-on perturbation, not the
hub-off reserve.

## Boundary Cases

The arbitrary fixed-`r` criterion above is cleanest for `r>=4`.  The cases
`r=2,3` should be handled separately, because they are first-order boundary
cases for the shift rule.

Acceptable treatments:

```text
1. direct exact formulas for their adjacent ratios and reserve;
2. separate shifted-positive certificates;
3. reduction of r=2 to the pure-spider lane plus a direct r=3 certificate.
```

Until this is written, avoid stating an unqualified arbitrary fixed-`r`
theorem for all `r>=2`.

## Limitations

This proof note does not prove the full `d_leaf<=1` case.

It handles structured spider lanes

```text
S(2^a,r)
```

and gives a certificate framework for fixed `r`.  The general tree problem
still requires a bridge from arbitrary `d_leaf<=1` trees to the Route-2
endpoint inequality.

## Next Use

This note is suitable as the base for a future manuscript section only after:

```text
1. r=2,3 are disposed of cleanly;
2. the global F-mode margin lemma is written in final form;
3. the certificate procedure is stated without implementation-specific
   Python language;
4. Theorem A and Theorem B are kept visibly separate.
```
