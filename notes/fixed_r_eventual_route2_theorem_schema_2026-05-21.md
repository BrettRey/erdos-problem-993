# Fixed-`r` Eventual Route-2 Theorem Schema

## Purpose

This note records the upgraded theorem shape after the path-moment reduction.
It does not replace the lane certificates.  It explains why, for every fixed
`r`, the current proof template should have a computable eventual threshold.

The finite range below that threshold still has to be checked exactly.

## Statement Shape

For every fixed `r>=2`, there should be a computable threshold `A(r)` such
that the Route-2 inequality holds for

```text
T_{a,r}=S(2^a,r)
```

for all `a>=A(r)`.  Together with exact finite checking for `a<A(r)`, this
proves the fixed-`r` lane.

The proof has four eventual pieces:

```text
1. eventual F-mode shift;
2. eventual hub-off reserve;
3. eventual hub-on mode domination;
4. eventual hub-on Route-2 perturbation domination.
```

## 1. Eventual `F`-Mode Shift

For

```text
F_{a,r}(x)=P_r(x)(1+2x)^a,
m=(2a+D)/3,
q=a mod 3,
```

the adjacent ratios have expansions:

```text
F_{m-1}/F_m = 1 + L(D)/a + O_r(a^-2),
F_m/F_{m+1} = 1 + L(D+3)/a + O_r(a^-2),
```

where

```text
L(D) = (9/2)(D/3-M_1)-3,
M_1=P'_r(1)/P_r(1).
```

Thus the first-order candidate interval is:

```text
3M_1-1 <= D <= 3M_1+2,
D == q mod 3.
```

The divisibility lemma in
`notes/fixed_r_shift_rule_from_path_moments_2026-05-21.md` proves that, for
`r>=4`, each residue has a unique candidate.  Therefore both adjacent
inequalities have nonzero first-order margins for the chosen `D`, and the
`O_r(a^-2)` term is eventually smaller than the `1/a` margin.

For `r=2,3`, the only first-order boundary cases are handled separately by
the second-order term or by direct exact formulas.

Because the adjacent ratios are rational functions of `a` on each residue
class, this eventual threshold is computable by exact rational arithmetic and
root isolation.

## 2. Eventual Hub-Off Reserve

Set

```text
lambda_0=F_{m-1}/F_m,
F^-_{a,r}(x)=P_r(x)(1+2x)^(a-1).
```

The reserve is

```text
R_{a,r,q}
 = mu_{F^-_{a,r}}(lambda_0)-m+4/3.
```

The asymptotic expansion is:

```text
R_{a,r,q}=C_{r,q}/a+O_r(a^-2).
```

The leading constant has the closed form

```text
C_{r,q}=(9(V+3K_3)+24alpha+16)/12,
alpha=M_1-D/3.
```

The note `notes/fixed_r_c_positivity_reduction_2026-05-21.md` proves
`V+3K_3>=0` for path polynomials.  Since the shift interval gives
`alpha>-2/3` for `r>=4`, and the `r=2,3` boundary shifts are checked
directly, we have:

```text
C_{r,q}>0
```

for the stabilized candidates.

Choose any rational

```text
0 < eta_{r,q} < C_{r,q}.
```

Then

```text
R_{a,r,q} >= eta_{r,q}/a
```

for all sufficiently large `a` in residue class `q`.

This is also computable: after substituting `a=3t+q`, the difference

```text
R_{a,r,q} - eta_{r,q}/a
```

is an exact rational function of `t`.  Its numerator and denominator have
known positive leading behavior, so exact root isolation gives a threshold
after which the inequality holds.

The compact fugacity interval

```text
3/4 <= lambda_0 <= 2
```

is eventual for the same reason, since `lambda_0 -> 1`.

## 3. Eventual Hub-On Mode Domination

The full polynomial is

```text
I(T_{a,r};x)=F_{a,r}(x)+G_{a,r}(x),
G_{a,r}(x)=xP_{r-1}(x)(1+x)^a.
```

Near the `F` saddle, `G` is exponentially smaller than `F`.  More precisely,
for fixed `r`, the path polynomial contributes only a bounded convolution,
while the binomial spines compare as

```text
max_k [x^k](1+x)^a
```

against

```text
2^m binom(a,m)
```

with `m=2a/3+O_r(1)`.  Stirling or standard binomial large-deviation bounds
therefore give constants

```text
K_r>0, B_r>=0, 0<rho<1
```

such that

```text
max_k G_{a,r,k}/F_{a,r,m} <= K_r a^{B_r} rho^a.
```

One fully explicit conservative route is:

```text
max_k G_{a,r,k} <= P_{r-1}(1) max_j binom(a,j) <= P_{r-1}(1) 2^a,
F_{a,r,m} >= 2^m binom(a,m).
```

Once `a` is large enough that `1/2 <= m/a <= 3/4`,

```text
binom(a,m) >= 2^{a H(m/a)}/(a+1) >= 2^{a H(3/4)}/(a+1),
2^m >= 2^{a/2}.
```

Hence

```text
max_k G_{a,r,k}/F_{a,r,m}
 <= P_{r-1}(1)(a+1) 2^{a(1/2-H(3/4))}.
```

Since `H(3/4)>1/2`, this is exponentially small.  The current code uses
sharper constants, but this bound is already enough for the theorem schema.

Consequently this ratio is eventually smaller than any prescribed rational
margin of the form `eta/a`, including the eventual `F`-mode margin.

This transfers the `F`-mode shift to the full polynomial.

## 4. Eventual Hub-On Route-2 Perturbation

For the bridge polynomial

```text
B_{a,r}(x)=F^-_{a,r}(x)+G^-_{a,r}(x),
G^-_{a,r}(x)=xP_{r-1}(x)(1+x)^(a-1),
```

the already-separated perturbation lemma gives:

```text
|mu_{B_{a,r}}(lambda)-mu_{F^-_{a,r}}(lambda_0)|
 <= 2(a+r)^2T(a,r,q)
    +(a+r)C_r(3/4)^(a-1),
```

where `T(a,r,q)` is controlled by the hub-on/full-mode domination ratio.

By the previous exponential domination, the right-hand side is eventually
smaller than

```text
eta_{r,q}/a.
```

Thus the hub-off reserve survives the hub-on perturbation.

## Consequence

For every fixed `r`, the proof template now gives an eventual theorem:

```text
There is a computable A(r) such that Route-2 holds for S(2^a,r)
for all a>=A(r).
```

The exact finite checker can then verify `1<=a<A(r)`.

## Effective-Threshold Probe

Implemented in:

```text
gpt_attack/fixed_r_effective_threshold_probe.py
```

The script constructs exact rational functions for:

```text
left F-margin  >= eta_L/a,
right F-margin >= eta_R/a,
3/4 <= lambda_0 <= 2,
hub-off reserve >= eta_Route2/a.
```

For the margin and reserve targets it uses half of the proved leading
constants.  The hub-off reserve is built using the cleared path recurrence,
not by expanding `P_r(lambda_0)` directly.

The script now has three threshold methods:

```text
isolate: exact real-root isolation;
shift:   shifted-coefficient positivity with exponential/binary search;
cauchy:  crude Cauchy positive-root bound.
```

The shifted-coefficient method is the practical one.  It gives:

```text
r=4:   A=12
r=8:   A=24
r=20:  A=94
r=24:  A=22
r=40:  A=26
r=80:  A=44
r=120: A=62
```

The `r=120` per-residue output is:

```text
q=0, D=99,  A=57; limiting check: hub-off reserve, degree 2745/2746
q=1, D=100, A=58; limiting check: left F-margin, degree 44/45
q=2, D=101, A=62; limiting check: left F-margin, degree 44/45
```

This validates the rational-function threshold route in nontrivial cases and
shows that the hub-off/mode/lambda pieces are not the current asymptotic
threshold bottleneck through `r=120`.

Scaling result:

```text
r=160 shifted-coefficient probing was stopped before q=0 output.
```

So larger practical hub-off threshold certificates still need either more
engineering around the recurrence arithmetic or the existing GMP path.

## Effective Hub-On Probe

Implemented in:

```text
gpt_attack/fixed_r_hubon_effective_threshold_probe.py
```

This probe uses:

```text
eta_mode    = half of the first-order F-margin constant,
eta_reserve = C_{r,q}/2.
```

It also corrects the monotonicity condition: since the comparison target is
`const/a`, the certificate checks decrease of `a` times each perturbation
term, not just decrease of the perturbation term itself.

Step-5 grid results:

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

Thus the hub-on perturbation remains the dominant practical threshold after
the hub-off leading-constant improvement.

## What Remains for a Manuscript-Ready Proof

The main missing work is no longer the sign of the leading hub-off constant.
The remaining proof-writing tasks are:

```text
1. write the adjacent-ratio expansion with an explicit rational-function
   threshold procedure;
2. write the hub-off reserve threshold procedure using root isolation rather
   than shifted coefficient positivity;
3. state a clean binomial large-deviation bound for the hub-on term;
4. decide whether the manuscript presents this as an existence/computable
   threshold theorem or as certified lanes up to a stated r frontier.
```

The current certificates are stronger for the sampled lanes because they
produce concrete thresholds and exact positivity files.  This theorem schema
explains why the same template is not merely empirical for fixed `r`.
