# Fixed-`r` Route-2 Proof Assembly Scratch

## Status

This is a new scratch document.  It is not a manuscript edit.

Purpose: assemble the current fixed-`r` Route-2 argument in theorem order,
without touching `paper/main_v2.tex` or the existing notes.

## Target Theorem

For a fixed integer `r>=2`, let

```text
T_{a,r}=S(2^a,r)
```

be the spider with `a` arms of length `2` and one arm of length `r`.  Remove
one length-`2` arm to obtain

```text
B_{a,r}=S(2^(a-1),r).
```

Let `m=m(a,r)` be the leftmost mode of `I(T_{a,r};x)`, and set

```text
lambda = i_{m-1}(T_{a,r}) / i_m(T_{a,r}).
```

The fixed-`r` Route-2 target is:

```text
mu_{B_{a,r}}(lambda) >= m - 3/2
```

for every `a>=1`.

The theorem shape now available is:

```text
For every fixed r>=2, there is a computable threshold A(r) such that
the Route-2 target holds for all a>=A(r).  The finite range
1<=a<A(r) is then checked by exact integer arithmetic.
```

This is still not a uniform-in-`r` theorem.

## Polynomial Split

Let `P_j(x)` be the independence polynomial of the path on `j` vertices.  Then

```text
I(T_{a,r};x)=F_{a,r}(x)+G_{a,r}(x),
F_{a,r}(x)=P_r(x)(1+2x)^a,
G_{a,r}(x)=xP_{r-1}(x)(1+x)^a.
```

After removing a length-`2` arm,

```text
I(B_{a,r};x)=F^-_{a,r}(x)+G^-_{a,r}(x),
F^-_{a,r}(x)=P_r(x)(1+2x)^(a-1),
G^-_{a,r}(x)=xP_{r-1}(x)(1+x)^(a-1).
```

The proof strategy is:

```text
1. localize the mode using F, then show G cannot move it;
2. prove a hub-off reserve for F^- at lambda_0=F_{m-1}/F_m;
3. show replacing lambda_0 by lambda and F^- by F^-+G^- costs less than
   the hub-off reserve.
```

## Lemma 1: First-Order Shift Rule

Let

```text
m=(2a+D)/3,
delta=D/3,
q=a mod 3,
D == q mod 3.
```

Let `S` be the independent-set size on the fixed path `P_r`, distributed by
the coefficients of `P_r` at fugacity `1`, and put

```text
M_1=E[S]=P'_r(1)/P_r(1).
```

For

```text
lambda_0=F_{m-1}/F_m,
```

the adjacent-ratio expansion is:

```text
lambda_0 = 1 + L(D)/a + O_r(a^-2),
L(D) = (9/2)(delta-M_1)-3.
```

Similarly,

```text
F_m/F_{m+1}=1+L(D+3)/a+O_r(a^-2).
```

Thus the first-order adjacent-mode conditions are:

```text
L(D)<=0,
L(D+3)>=0,
```

equivalently

```text
3M_1-1 <= D <= 3M_1+2,
D == q mod 3.
```

For `r>=4`, each residue has a unique candidate in this interval.

Reason: `P_r(1)=F_{r+2}` and

```text
P'_r(1)=(rL_{r+1}+2F_r)/5.
```

If `3M_1` is integral, then `F_{r+2}` divides `3(r+2)`.  This is impossible
for `r>=7`, and `r=4,5,6` are checked directly.  The only first-order
boundary cases are `r=2,3`, where direct exact shifts select the lower
candidate:

```text
r=2: {0:3, 1:1, 2:2}
r=3: {0:3, 1:4, 2:2}
```

Remaining manuscript obligation: write the `O_r(a^-2)` remainder as an
explicit rational-function threshold procedure, or cite the shifted-positive
certificate procedure for each fixed `r`.

## Lemma 2: Hub-Off Reserve Leading Constant

Set

```text
R_{a,r,q}
 = mu_{F^-_{a,r}}(lambda_0) - m + 4/3.
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
C_{r,q} = (9(V+3K_3)+24alpha+16)/12.
```

The shift interval gives `alpha>-2/3` for `r>=4`.  For `r=2,3`, the direct
boundary shifts give positive constants.

It remains only to prove:

```text
V+3K_3>=0
```

for the path-size distribution.

## Lemma 3: Path Moment Inequality

Let

```text
T_r=P_r(1),
S_j(r)=sum_s s^j p_s
```

where `P_r(x)=sum_s p_s x^s`.

The path recurrence gives:

```text
T_r=T_{r-1}+T_{r-2},
S_1(r)=S_1(r-1)+S_1(r-2)+T_{r-2},
S_2(r)=S_2(r-1)+S_2(r-2)+2S_1(r-2)+T_{r-2},
S_3(r)=S_3(r-1)+S_3(r-2)+3S_2(r-2)+3S_1(r-2)+T_{r-2}.
```

Solving:

```text
T_r = F_{r+2},

S_1(r) = (r/2+2/5)F_r + (r/10)L_r,

S_2(r) = (r^2/5+r/2+8/25)F_r - (r/50)L_r,

S_3(r) = (r^3/10+21r^2/50+23r/50+4/25)F_r
         -(r^3/50+3r^2/50+3r/50)L_r.
```

Substitution and `L_r^2-5F_r^2=4(-1)^r` give:

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

Small cases:

```text
r=0: V+3K_3=0
r=1: V+3K_3=1/4
r=2: V+3K_3=0
r=4: V+3K_3=5/32
```

Therefore `V+3K_3>=0` for all paths, with equality only at `r=0,2`.

Conclusion:

```text
C_{r,q}>0
```

for all fixed-`r` stabilized shifts.

## Lemma 4: Effective Hub-Off Threshold

Choose a rational

```text
0<eta_{r,q}<C_{r,q}.
```

For each fixed residue class, substitute `a=3t+q`.  Then

```text
R_{a,r,q}-eta_{r,q}/a
```

is an exact rational function of `t`.  Its leading term is positive, so an
eventual threshold is computable.

Practical certificate method:

```text
shifted-coefficient positivity after t -> u+T.
```

The current effective-threshold probe uses `eta=C_{r,q}/2` and the cleared
path recurrence for the reserve numerator.  Verified shifted-positive
thresholds for the hub-off/mode/lambda pieces:

```text
r=4:   A=12
r=8:   A=24
r=20:  A=94
r=24:  A=22
r=40:  A=26
r=80:  A=44
r=120: A=62
```

Stress note: `r=160` shifted probing was stopped before the first residue
completed.  This is a practical scaling issue, not a sign issue.

## Lemma 5: Hub-On Mode Domination

The full mode is controlled by showing that `G` is too small to move the
certified `F` mode.

A sufficient condition is:

```text
2 max_k G_k < F_m * eta_mode/a,
```

where `eta_mode` is a positive lower bound for the first-order adjacent
`F`-mode margin.

The code bounds

```text
max_k G_k/F_m
```

by choosing a witness term in `F_m` and comparing the binomial spines.
Monotonicity in `a -> a+3` is checked by a shifted polynomial.

This part is not the bottleneck in tested lanes.

## Lemma 6: Hub-On Route-2 Perturbation

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

For any nonnegative polynomial `H` of degree at most `N`,

```text
|mu_H(lambda)-mu_H(lambda_0)| <= 2N^2T.
```

For the bridge mixture,

```text
|mu_{F^-+G^-}(lambda)-mu_{F^-}(lambda)|
 <= (a+r) C_r (3/4)^(a-1).
```

Thus

```text
|mu_{B_{a,r}}(lambda)-mu_{F^-_{a,r}}(lambda_0)|
 <= 2(a+r)^2T + (a+r)C_r(3/4)^(a-1).
```

Correction to carry into any manuscript proof:

Since the reserve comparison target is `eta/a`, it is not enough to show that
the perturbation terms decrease.  The certificate must show that

```text
a * perturbation
```

decreases in each residue class.  The scripts have been corrected to check
this weighted monotonicity.

Effective hub-on thresholds using `eta_reserve=C_{r,q}/2` and step-5 grids:

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

Current diagnosis:

```text
hub-on perturbation is the practical threshold driver.
```

## Certificate Proposition for Sampled Lanes

For

```text
r in {4,8,12,16,20,24,32,40,60,80},
```

the existing certificate suite proves the Route-2 target for all `a>=1`.

Split:

```text
finite exact check: a<=199;
asymptotic certificates: a>=200.
```

After correcting the hub-on monotonicity condition, the sampled suite was
rerun and still reports:

```text
all_ok=True.
```

This supports a clean computed proposition:

```text
For each r in {4,8,12,16,20,24,32,40,60,80}, Route-2 holds for S(2^a,r)
after removing a length-2 arm, for every a>=1.
```

## Remaining Gaps Before Manuscript Text

The proof components are now coherent, but a manuscript-ready theorem still
needs careful scope control.

Open items:

```text
1. Decide whether to state the arbitrary fixed-r result as a computable
   threshold theorem or keep it as a theorem schema plus certified lanes.

2. Write the effective threshold procedure abstractly enough that it does not
   depend on Python implementation details.

3. Either prove a sharper analytic hub-on threshold bound, or state hub-on
   domination as a shifted-polynomial certificate condition.

4. Keep the sampled-lane proposition separate from the arbitrary fixed-r
   theorem schema.

5. Do not claim the full d_leaf<=1 theorem from these spider lanes alone.
```

## Current Call

The next best writing target is not another computation.  It is a polished
proof note with:

```text
Theorem A: sampled fixed-r certificate proposition.
Theorem B: arbitrary fixed-r computable-threshold schema.
Lemma package: shift rule, path moment positivity, hub-off reserve,
hub-on perturbation.
```

The proof should explicitly say where a condition is analytic and where it is
certificate-backed.
