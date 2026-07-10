# Universal Poisson-Binomial Effective-Drop Theorem

Date: 2026-07-10
Status: proved and independently audited; the compact range has independent
exact Sturm and positive-Bernstein certificates

## Theorem

Let

\[
W=\sum_{i=1}^n\operatorname{Bernoulli}(p_i)
\]

be a finite nondegenerate Poisson-binomial law. Deterministic parameters may
first be stripped and absorbed into a shift. Write

\[
f_k=\Pr(W=k),\qquad V=\operatorname{Var}(W)\ge1,
\]

and suppose the first strict descent is nonterminal:

\[
D=\min\{k\ge1:f_k<f_{k-1}\}\le n.
\]

With \(f_{n+1}=0\),

\[
\boxed{
V\left(1-\frac{f_{D-1}f_{D+1}}{f_D^2}\right)\ge\frac14.
}
\]

If the largest mode is the upper support endpoint, the next descent has
\(f_D=0\), so the quotient is undefined. For the raw-reserve application this
terminal case is immediate: the next ratio is zero and the reserve is one.

## Signed Corollary

Every signed low-probability law

\[
Z=X-Y
\]

becomes an ordinary Poisson-binomial law after shifting by the number of
Bernoulli factors in \(Y\):

\[
Z+m_Y=X+\sum_{j=1}^{m_Y}(1-Y_j).
\]

The shift preserves the variance, the first-descent local ratios, and the
effective drop. Therefore the theorem proves the full nonterminal finite
signed effective-drop target, not only the one- and two-reflected-factor
cases. A terminal signed descent has immediate raw reserve, but its effective
quotient is undefined.

At a nonterminal first descent, write

\[
R_-=\frac{f_D}{f_{D-1}}<1,\qquad
R_+=\frac{f_{D+1}}{f_D}.
\]

Then

\[
1-R_+
\ge
1-\frac{R_+}{R_-},
\]

so the same theorem also yields the raw reserve

\[
V(1-R_+)\ge\frac14.
\]

This closes the finite signed raw-reserve lemma, with the terminal case
handled as above. It does not by itself prove the
remaining hub-bouquet comparison, the full Case-B tree bound, or Erdős
Problem 993.

## Cubic Curvature Propagation

Strip deterministic parameters so that the positive support is
\(\{0,\ldots,n\}\), and set

\[
f_{-1}=f_{n+1}=0.
\]

Define the normalized Turán deficit

\[
\delta_k
=1-\frac{f_{k-1}f_{k+1}}{f_k^2},
\qquad 0\le k\le n.
\]

In particular,

\[
\delta_0=\delta_n=1.
\]

Hillion and Johnson's cubic inequalities for Bernoulli sums imply

\[
(f_{k-1}^2-f_{k-2}f_k)f_{k+1}
\le
f_{k-1}(f_k^2-f_{k-1}f_{k+1}),
\]

\[
(f_{k+1}^2-f_kf_{k+2})f_{k-1}
\le
f_{k+1}(f_k^2-f_{k-1}f_{k+1}).
\]

After division by positive central masses, these become

\[
\boxed{
\delta_{k-1}(1-\delta_k)\le\delta_k,
\qquad
\delta_{k+1}(1-\delta_k)\le\delta_k.
}
\tag{1}
\]

The endpoint convention is essential. The inequalities remain valid adjacent
to the support boundary; for example,

\[
\delta_n(1-\delta_{n-1})\le\delta_{n-1}.
\]

A recurrence-only statement without \(\delta_0=\delta_n=1\) is false. For
example, the normalized weights proportional to
\((1000,1000,900,729)\) have variance greater than one and local deficit
\(1/10\), but fail the terminal cubic inequality.

The primary source is Hillion--Johnson,
[A proof of the Shepp--Olkin entropy concavity
conjecture](https://arxiv.org/abs/1503.01570), Appendix A.2, equations
(A.7)--(A.8).

## Explicit Mass Windows

Let

\[
\delta=\delta_D.
\]

If \(\delta\ge\tfrac14\), the theorem follows immediately from \(V\ge1\).
Assume henceforth that

\[
0<\delta<\frac14.
\]

Iterating \(x\mapsto x/(1-x)\) in (1) gives

\[
\delta_{D\pm r}\le\frac{\delta}{1-r\delta}
\qquad(r\delta<1).
\tag{2}
\]

The strict positivity of \(\delta\) also follows from (2) and the endpoint
values: zero curvature would propagate all the way to an endpoint and
contradict \(\delta_0=\delta_n=1\).

If an endpoint occurred at distance \(r\) with
\((r+1)\delta<1\), then (2) would give

\[
1=\delta_{\mathrm{endpoint}}
\le\frac{\delta}{1-r\delta}<1,
\]

a contradiction. Thus all masses in the window below genuinely lie inside
the support.

Translate the largest mode \(D-1\) to zero, and put

\[
p=f_{D-1}=\max_k f_k,
\qquad
q_k=\frac{f_k}{f_{k-1}}.
\]

Because \(q_{D-1}\ge1>q_D\), (1) gives

\[
q_D
=q_{D-1}(1-\delta_{D-1})
\ge\frac{1-2\delta}{1-\delta}
=:a.
\tag{3}
\]

Let

\[
K=\max\{r\ge1:(r+1)\delta<1\}.
\]

Telescoping (2)--(3) yields, for \(1\le r\le K\),

\[
\frac{f_{D-1+r}}p
\ge
R_r
:=
a^r\prod_{j=1}^{r-1}(1-j\delta),
\tag{4}
\]

\[
\frac{f_{D-1-r}}p
\ge
L_r
:=
(1-\delta)^{-r}\prod_{j=2}^{r+1}(1-j\delta).
\tag{5}
\]

## Two Variance Bounds

Set

\[
w_0=1,\qquad w_r=R_r,\qquad w_{-r}=L_r,
\]

and define

\[
A(\delta)
=\frac12\sum_{i,j=-K}^K w_iw_j(i-j)^2.
\]

The pairwise variance identity gives

\[
V
=\frac12\sum_{i,j}f_if_j(i-j)^2
\ge p^2A(\delta).
\tag{6}
\]

There is also a sharp elementary max-atom bound:

\[
V\ge\frac{p^{-2}-1}{12}.
\tag{7}
\]

To prove it, let \(U\sim\operatorname{Unif}[-\tfrac12,\tfrac12]\) be
independent. The density of \(W+U\) is bounded by \(p\), and

\[
\operatorname{Var}(W+U)=V+\frac1{12}.
\]

Among real densities bounded by \(p\), the centered uniform density on an
interval of length \(1/p\) minimizes variance. Hence
\(\operatorname{Var}(W+U)\ge1/(12p^2)\), proving (7).

It now suffices to establish the scalar inequality

\[
\boxed{
A(\delta)\ge\frac{3+\delta}{4\delta^2}.
}
\tag{8}
\]

Indeed, if \(V<1/(4\delta)\), then (7) implies

\[
p^2>\frac{\delta}{3+\delta},
\]

while (6) and (8) imply the incompatible inequality

\[
p^2<\frac1{4\delta A(\delta)}
\le\frac{\delta}{3+\delta}.
\]

## Scalar Reparameterization

Put

\[
H=\frac{1-\delta}{\delta}=\frac1\delta-1>3.
\]

The left weights simplify to

\[
L_r=b_r
:=
\prod_{s=1}^r\left(1-\frac{s}{H}\right).
\tag{9}
\]

The right weights are no smaller. If \(C_r=R_r/L_r\), then \(C_1=1\) and

\[
\frac{C_{r+1}}{C_r}
=
\frac{(H-1)(H+1-r)}
     {(H+1)(H-r-1)}
\ge1,
\]

because the numerator exceeds the denominator by \(2r\). Thus

\[
R_r\ge L_r.
\tag{10}
\]

For symmetric weights \(b_{\pm r}=b_r\), write

\[
S=1+2\sum_r b_r,
\qquad
T=2\sum_r r^2b_r.
\]

Their first moment is zero, so their pairwise expression is \(ST\). Since
the pairwise expression is coordinatewise increasing in all weights,
(10) gives

\[
A(\delta)\ge ST.
\]

The right side of (8) becomes

\[
Q(H)
:=
\frac{(3H+4)(H+1)}4.
\tag{11}
\]

The remaining proof is \(ST\ge Q(H)\), except on the smallest cell where the
exact asymmetric weights are needed.

## Exact Cell \(3<H<4\)

Here \(K=3\). Use the exact six weights
\(L_1,L_2,L_3,R_1,R_2,R_3\), rather than the symmetric relaxation.
After subtracting \(Q(H)\), the denominator is the positive quantity

\[
4H^5(H+1)^3,
\]

and the numerator is

\[
\begin{aligned}
P(H)={}&-3H^{10}-16H^9+750H^8-3676H^7+6613H^6\\
&-5460H^5+800H^4+4696H^3-7176H^2+4272H-912.
\end{aligned}
\]

The exact Sturm count on \((3,4)\) is zero, while

\[
P(3)=22272,\qquad P(4)=3486896.
\]

Therefore \(P(H)>0\) throughout this cell.

## Exact Compact Cells \(4\le H\le16\)

On \(H\in[m,m+1]\), for \(m=4,\ldots,15\), use

\[
S_m=1+2\sum_{r=1}^m b_r,
\qquad
T_m=2\sum_{r=1}^m r^2b_r.
\]

The extra \(b_m\) vanishes at \(H=m\), so this convention is exact at both
integer boundaries. Define

\[
P_m(H)=4H^{2m}\bigl(S_mT_m-Q(H)\bigr).
\]

Each \(P_m\) is an integer polynomial of degree \(2m+2\). Exact Sturm
sequences give zero roots on every interval \([m,m+1]\), and both endpoint
values are positive. The replay script constructs and checks all twelve
polynomials, root counts, and endpoint integers; the JSON certificate records
their degrees, root counts, and exact positivity outcomes.

Thus

\[
S_mT_m\ge Q(H)
\qquad(4\le H\le16).
\]

These compact Sturm checks, together with the preceding \(3<H<4\) Sturm
check, are the only finite computer-assisted parts of the proof.

### Lean-friendly Bernstein replacement

There is also an exact Sturm-free certificate for the same thirteen finite
cells. After writing \(H=m+t\), \(0\le t\le1\), expand each numerator in its
full-degree Bernstein basis. Every coefficient is strictly positive: 11 for
the asymmetric cell and 264 across the twelve compact symmetric cells, for
275/275 positive coefficients in total. Exact reconstruction verifies every
identity.

Thus the finite range can equivalently be discharged by coefficient
positivity and nonnegativity of

\[
\binom di t^i(1-t)^{d-i},
\]

which is substantially easier to replay in Lean than Sturm root counting.
The independent certificate is documented in
`notes/literature/universal_pb_finite_bernstein_replacement_2026-07-10.md`
and replayed by `scripts/verify_universal_pb_finite_bernstein.py`.

## Analytic Range \(H\ge16\)

Choose \(R\) so that

\[
\frac{R(R+1)}2\le H\le\frac{(R+1)(R+2)}2.
\]

The Weierstrass product bound gives, for \(1\le r\le R\),

\[
b_r
\ge
c_r
:=
1-\frac{r(r+1)}{2H}
\ge0.
\]

Put

\[
S_0=1+2\sum_{r=1}^R c_r
=2R+1-\frac{R(R+1)(R+2)}{3H},
\]

\[
T_0=2\sum_{r=1}^Rr^2c_r
=
\frac{R(R+1)(2R+1)}3
-\frac{Q_4+Q_3}{H},
\]

where

\[
Q_3=\frac{R^2(R+1)^2}{4},
\]

\[
Q_4=
\frac{R(R+1)(2R+1)(3R^2+3R-1)}{30}.
\]

It is enough to show \(S_0T_0\ge Q(H)\). Let

\[
N_R(H)=H^2\bigl(S_0T_0-Q(H)\bigr).
\]

On the triangular interval, its degree-four Bernstein coefficients are
positive prefactors divided by \(2880\), multiplied respectively by

\[
\begin{aligned}
P_0(R)&=121R^4-430R^3-965R^2-510R-736,\\
P_1(R)&=121R^5-188R^4-1081R^3-1889R^2-1757R-696,\\
P_2(R)&=121R^6+54R^5-1097R^4-3068R^3-3778R^2-2492R-480,\\
P_3(R)&=121R^5+54R^4-1121R^3-2885R^2-3189R-1350,\\
P_4(R)&=121R^4+54R^3-1529R^2-3246R-2520.
\end{aligned}
\]

The prefactors are, in order,

\[
R^2(R+1)^2,\quad
R(R+1)^2,\quad
(R+1)^2,\quad
(R+1)^2(R+2),\quad
(R+1)^2(R+2)^2.
\]

For \(R=u+6\), the five polynomials become

\[
\begin{aligned}
&121u^4+2474u^3+17431u^2+46014u+25400,\\
&121u^5+3442u^4+37967u^3+199405u^2+480475u+384510,\\
&121u^6+4410u^5+65863u^4+512764u^3\\
&\hspace{2cm}+2172926u^2+4668316u+3829440,\\
&121u^5+3684u^4+43735u^3+249961u^2+671859u+644400,\\
&121u^4+2958u^3+25579u^2+88782u+91440.
\end{aligned}
\]

Every coefficient is positive, proving the claim for \(R\ge6\).

The remaining \(R=5\) range needed is \(16\le H\le21\). Its degree-four
Bernstein coefficients are

\[
2360,\quad 7500,\quad \frac{25055}{2},
\quad \frac{254205}{16},\quad \frac{31115}{2},
\]

all positive. This completes the analytic range and proves (8).

## Proof Boundary

The theorem proves the universal finite Poisson-binomial effective-drop and
raw-reserve inequalities at the quarter scale. Its only external
mathematical input is the Hillion--Johnson cubic inequality. The compact
scalar range is discharged by exact symbolic certificates, not floating-point
sampling.

The result does not determine the sharp constant. Sparse Poisson limits still
suggest \(1/3\) as the true boundary, but the present proof claims only
\(1/4\).

## Lean Formalization Status

The first formal layer is complete in
`formalization/pb_effective_drop_aristotle/PBReserve/Core.lean`. It verifies in
Lean 4.28:

- iteration of the normalized cubic-curvature recurrence over the window used
  in this proof;
- exclusion of a curvature-one endpoint from that window;
- the first-crossing ratio lower bound;
- transfer from effective drop to raw reserve.

The file has no `sorry` or additional axioms. This layer assumes the normalized
Hillion--Johnson recurrence rather than deriving it from a Poisson-binomial
law, so it is not yet an end-to-end formalization of the theorem. The audited
Aristotle run and exact scope boundary are recorded in
`formalization/pb_effective_drop_aristotle/ARISTOTLE_RESULT.md`.

## Replay

Run:

    python3 scripts/verify_universal_pb_effective_drop.py
    python3 scripts/verify_universal_pb_finite_bernstein.py

The machine-readable certificate is:

    results/universal_pb_effective_drop_certificate_2026-07-10.json
    results/universal_pb_finite_bernstein_certificate_2026-07-10.json
