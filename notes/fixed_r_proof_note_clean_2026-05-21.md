# Fixed-`r` Route-2 Proof Note

## Purpose and Scope

This is a clean proof-note draft for the fixed-`r` spider lanes.  It is not a
manuscript edit.

It separates three levels of claim:

```text
1. a fully certificate-backed sampled-lane proposition;
2. analytic lemmas that apply to arbitrary fixed r;
3. a computable-threshold criterion for fixed r, stated as a certificate
   criterion rather than as an unconditional manuscript theorem.
```

This note does not prove the full `d_leaf<=1` case.

## Notation

Let `P_j(x)` be the independence polynomial of the path on `j` vertices.  For
fixed `r>=2`, define

```text
T_{a,r}=S(2^a,r),
```

the spider with `a` arms of length `2` and one arm of length `r`.  Removing
one length-`2` arm gives

```text
B_{a,r}=S(2^(a-1),r).
```

The spider product formula gives

```text
I(T_{a,r};x)=F_{a,r}(x)+G_{a,r}(x),
F_{a,r}(x)=P_r(x)(1+2x)^a,
G_{a,r}(x)=xP_{r-1}(x)(1+x)^a.
```

For the bridge tree,

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

### Certificate Artifacts

The finite range is verified by exact integer arithmetic in:

```text
results/fixed_r_finite_route2_selected_a199.json
```

The asymptotic certificate suite, rerun after correcting the hub-on
monotonicity condition, is:

```text
results/fixed_r_sampled_certificate_suite_corrected_monotonicity.json
```

with result:

```text
all_ok=true.
```

The relevant scripts are:

```text
gpt_attack/fixed_r_finite_route2_check.py
gpt_attack/fixed_r_huboff_certificate.py
gpt_attack/fixed_r_hubon_mode_certificate.py
gpt_attack/fixed_r_hubon_route2_perturbation.py
gpt_attack/fixed_r_sampled_certificate_suite.py
```

### Proof

The proof is split at `a=200`.

For

```text
1 <= a <= 199,
```

the finite checker verifies the Route-2 inequality directly after clearing
denominators.

For

```text
a >= 200,
```

the certificate supplies:

```text
1. a stabilized residue-class mode formula;
2. hub-off mode margins and hub-off reserve;
3. hub-on mode domination;
4. hub-on Route-2 perturbation domination.
```

The general certificate lemma below turns those four checks into the Route-2
conclusion.

## Certificate Lemma

Fix `r` and a threshold `A`.  For each residue class `q in {0,1,2}`, suppose
the certificate supplies an integer `D_q` and sets

```text
a == q mod 3,
m=(2a+D_q)/3.
```

Assume C0-C3.

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

A sufficient certificate is:

```text
F_{a,r,m}-F_{a,r,k} >= F_{a,r,m} eta_mode/a     for every k != m,
2 max_k G_{a,r,k} < F_{a,r,m} eta_mode/a.
```

Then, for every `k != m`,

```text
(F_m+G_m)-(F_k+G_k)
 >= (F_m-F_k)-|G_m-G_k|
 > 0.
```

Thus `m` is the unique mode of `I(T_{a,r};x)`, hence the leftmost mode.

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
therefore makes `lambda` the correct tie fugacity.  By C2 and C3,

```text
mu_{B_{a,r}}(lambda)
 >= m - 4/3
 >  m - 3/2.
```

## Lemma 1: Adjacent `F` Margins Give Global `F` Margins

Let `(c_k)` be a nonnegative real-rooted coefficient sequence, hence
log-concave and unimodal.  Suppose `m` is a mode and

```text
c_m-c_{m-1} >= eps c_m,
c_m-c_{m+1} >= eps c_m
```

for some `eps>0`.  Then

```text
c_m-c_k >= eps c_m
```

for every `k != m`.

Proof.  Log-concavity implies that the ratios

```text
c_k/c_{k-1}
```

are nonincreasing where defined.  Moving left from `m`, the coefficients
cannot exceed `c_{m-1}`; moving right from `m`, they cannot exceed `c_{m+1}`.
Thus the adjacent deficits dominate all farther deficits.

Apply this to

```text
F_{a,r}(x)=P_r(x)(1+2x)^a.
```

The path polynomial `P_r(x)` is real-rooted with nonpositive roots, and
`(1+2x)^a` is real-rooted, so `F_{a,r}` is real-rooted.  Its coefficient
sequence is therefore log-concave.

## Lemma 2: First-Order Shift Rule

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

the adjacent-ratio expansion is

```text
lambda_0 = 1 + L(D)/a + O_r(a^-2),
L(D)=(9/2)(D/3-M_1)-3.
```

Similarly,

```text
F_m/F_{m+1}=1+L(D+3)/a+O_r(a^-2).
```

The first-order adjacent-mode conditions are

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

These boundary cases should be handled separately before stating an
unqualified theorem for all `r>=2`.

## Lemma 3: Hub-Off Reserve Constant

Define

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
alpha in [-2/3,1/3].
```

For `r>=4`, the lower endpoint is not attained.  For `r=2,3`, the listed
boundary shifts give positive constants directly.  Therefore positivity of
`C_{r,q}` reduces to

```text
V+3K_3 >= 0.
```

## Lemma 4: Path Moment Positivity

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

The following closed forms are verified by induction from the recurrences and
initial values:

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

## Lemma 5: Shifted-Positive Threshold Termination

Let `Q(t)` be a nonzero polynomial with positive leading coefficient.  Then
there exists `T` such that every coefficient of

```text
Q(u+T)
```

is positive.

Proof.  The coefficient of `u^j` in `Q(u+T)` is `Q^{(j)}(T)/j!`.  Every
nonzero derivative of `Q` has positive leading coefficient, so after a
sufficiently large shift all such values are positive.

For a rational function `N(t)/D(t)`, first cancel common factors and normalize
so that the denominator has positive leading coefficient.  Apply the same
shifted-positive test to numerator and denominator.  Positivity of the shifted
denominator coefficients guarantees that the denominator has no zeros past
the chosen threshold.

This gives a terminating certificate procedure for eventual positivity of
each rational inequality in a fixed residue class.

## Lemma 6: Hub-On Perturbation Bound

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

For the bridge mixture, using `lambda in [1/2,2]`,

```text
G^-_{a,r}(lambda)/F^-_{a,r}(lambda)
 <= C_r(3/4)^(a-1),
```

where

```text
C_r = 2 P_{r-1}(2)/P_r(1/2).
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

This is a criterion: it assumes production of the listed threshold
certificates.  It should not be restated as an unconditional arbitrary
fixed-`r` theorem until those certificate constructions are presented as a
complete terminating procedure.

## Current Effective Threshold Data

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

The criterion above is stated for `r>=4`.  The cases `r=2,3` should be handled
separately, because they are first-order boundary cases for the shift rule.

Acceptable treatments:

```text
1. direct exact formulas for their adjacent ratios and reserve;
2. separate shifted-positive certificates;
3. reduction of r=2 to the pure-spider lane plus a direct r=3 certificate.
```

Until this is written, avoid stating an unqualified arbitrary fixed-`r`
theorem for all `r>=2`.

## Limitations

This note handles structured spider lanes

```text
S(2^a,r).
```

It does not prove the full `d_leaf<=1` case.  The general tree problem still
requires a bridge from arbitrary `d_leaf<=1` trees to the Route-2 endpoint
inequality.

## Next Use

This note is suitable as the base for future manuscript work only after:

```text
1. r=2,3 are disposed of cleanly;
2. the certificate procedure is stated without implementation-specific Python
   language;
3. Theorem A and Theorem B remain visibly separate.
```
