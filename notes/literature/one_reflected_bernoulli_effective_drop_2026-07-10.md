# One Reflected Bernoulli Effective-Drop Theorem

Date: 2026-07-10
Status: proved, conditional only on the established one-sided localization
\(S+1\le4V_X\)

## Theorem

Let

\[
X=\sum_{i=1}^m \operatorname{Bernoulli}(p_i),
\qquad 0<p_i\le\frac12,
\qquad V_X=\operatorname{Var}(X)\ge1,
\]

and let \(Y\sim\operatorname{Bernoulli}(q)\) be independent, with
\(0\le q\le\tfrac12\). Write

\[
c_k=\Pr(X-Y=k),\qquad
V=V_X+q(1-q).
\]

Define the first strict descent in signed coordinates by

\[
D=\min\{k\ge0:c_k<c_{k-1}\}.
\]

Then \(c_D>0\) and

\[
\boxed{
V\left(1-\frac{c_{D+1}c_{D-1}}{c_D^2}\right)\ge\frac14.
}
\]

This is the first nontrivial finite signed case of the proposed effective-drop
theorem.

## Setup

Let

\[
a_k=\Pr(X=k),\qquad r_k=\frac{a_{k+1}}{a_k},
\]

with \(a_k=0\) outside \(0\le k\le m\). Delete any zero parameters before
defining \(m\). The ratios \(r_k\) are strictly decreasing.

Let \(S\) be the first strict descent of \(a\):

\[
a_S<a_{S-1}.
\]

Thus

\[
r_{S-2}\ge1,\qquad r_{S-1}<1.
\]

Moreover, \(S\ge2\). Indeed, if \(S=1\), then with
\(w_i=p_i/(1-p_i)\),

\[
\sum_i w_i=\frac{a_1}{a_0}<1,
\]

whereas

\[
V_X=\sum_i\frac{w_i}{(1+w_i)^2}<\sum_iw_i<1,
\]

contrary to the hypothesis.

The one-sided support and localization results give

\[
S\le m-1,\qquad S+1\le4V_X.
\]

Put

\[
t=\frac q{1-q}\in[0,1].
\]

The signed coefficients and their consecutive ratios are

\[
c_k=(1-q)a_k+qa_{k+1}
=(1-q)a_k(1+t r_k),
\]

\[
R_k:=\frac{c_k}{c_{k-1}}
=r_{k-1}\frac{1+t r_k}{1+t r_{k-1}}.
\]

Since \(r_k\le r_{k-1}\),

\[
r_k\le R_k\le r_{k-1}.
\]

## Descent Selection

The signed sequence cannot descend at its new left endpoint:

\[
c_0-c_{-1}=(1-2q)a_0+qa_1\ge0.
\]

For every \(1\le k\le S-2\), \(r_k\ge1\), so \(R_k\ge1\). At \(k=S\),
\(R_S\le r_{S-1}<1\). Therefore

\[
\boxed{D\in\{S-1,S\}.}
\]

There is also an exact selection threshold. If

\[
u=r_{S-2}\ge1,\qquad v=r_{S-1}<1,
\]

then

\[
R_{S-1}<1
\quad\Longleftrightarrow\quad
t>\frac{u-1}{u(1-v)}.
\]

Hence \(D=S-1\) above the threshold and \(D=S\) below it. Equality creates
a plateau at \(S-2,S-1\), so the first *strict* descent is again \(S\).

This explicit threshold is the plateau-safe replacement for the false
strict-descent continuity heuristic.

## A Curvature Bound for the Bernoulli Transform

For every interior \(k\),

\[
\begin{aligned}
c_k^2-c_{k-1}c_{k+1}
={}&(1-q)^2(a_k^2-a_{k-1}a_{k+1})\\
&+q^2(a_{k+1}^2-a_ka_{k+2})\\
&+q(1-q)(a_ka_{k+1}-a_{k-1}a_{k+2}).
\end{aligned}
\]

Newton's inequalities imply the three lower bounds

\[
a_k^2-a_{k-1}a_{k+1}\ge\frac{a_k^2}{k+1},
\]

\[
a_{k+1}^2-a_ka_{k+2}\ge\frac{a_{k+1}^2}{k+2},
\]

and, by multiplying two consecutive ratio drops,

\[
a_ka_{k+1}-a_{k-1}a_{k+2}
\ge\frac{2a_ka_{k+1}}{k+2}.
\]

Divide the resulting inequality by
\(c_k^2=(1-q)^2a_k^2(1+t r_k)^2\). A direct simplification gives

\[
\boxed{
\Delta_k(c):=
1-\frac{c_{k+1}c_{k-1}}{c_k^2}
\ge
\frac1{k+2}
+\frac1{(k+1)(k+2)(1+t r_k)^2}.
}
\]

This estimate is useful precisely in the case where the strict descent does
not move.

## Case 1: The Descent Moves Left

Suppose \(D=S-1\). Define

\[
d_j=c_{j-1}.
\]

Then \(d\) is the probability mass of

\[
X+(1-Y),
\]

an ordinary Poisson-binomial sum. Its effective drop at index \(S\) equals
the signed effective drop at \(D\). Newton's inequality for this
\((m+1)\)-factor polynomial gives

\[
\Delta_D(c)=\Delta_S(d)\ge\frac1{S+1}.
\]

Using \(S+1\le4V_X\),

\[
\bigl(V_X+q(1-q)\bigr)\Delta_D(c)
\ge\frac{V_X}{S+1}\ge\frac14.
\]

No low-probability hypothesis is needed for the added factor \(1-q\);
ordinary Newton inequalities apply to every Poisson-binomial law.

## Case 2: The Descent Stays Put

Suppose \(D=S\), and set

\[
\ell=S+1.
\]

Because \(r_S<1\), the curvature bound at \(k=S\) gives

\[
\Delta_D(c)
\ge
\frac1{\ell+1}
\left(1+\frac1{\ell(1+t)^2}\right).
\]

Since

\[
V_X\ge\frac\ell4,\qquad
q(1-q)=\frac{t}{(1+t)^2},
\]

we obtain

\[
4V\Delta_D(c)
\ge
\frac{\ell+4t/(1+t)^2}{\ell+1}
\left(1+\frac1{\ell(1+t)^2}\right).
\]

Subtracting \(1\) from that product yields exactly

\[
\frac1{\ell+1}
\left(
\frac{t(2-t)}{(1+t)^2}
+\frac{4t}{\ell(1+t)^4}
\right)\ge0,
\]

because \(0\le t\le1\). Hence \(V\Delta_D(c)\ge\tfrac14\).

## Exact Calibration Examples

### Plateau jump

For \(X\sim\operatorname{Bin}(5,\tfrac12)\), \(S=4\). At \(q=0\), the
strict descent is \(D=4\) and
\[
\Delta_D=\frac35.
\]
For every \(q>0\), the threshold is zero, \(D=3\), and
\[
\Delta_D\longrightarrow\frac12
\qquad(q\downarrow0).
\]
The theorem handles this discontinuity through Case 1.

### Nonzero threshold

For
\[
X\sim\operatorname{Bin}\left(5,\frac{11}{31}\right),
\]
we have
\[
S=3,\qquad
u=\frac{11}{10},\qquad
v=\frac{11}{20},\qquad
t_*=\frac{20}{99},\qquad
q_*=\frac{20}{119}.
\]

At \(q=q_*\), the signed law has a plateau at indices \(1,2\) and its first
strict descent remains \(D=3\). At \(q=\tfrac12\), \(D=2\); at
\(q=\tfrac1{10}\), \(D=3\).

## Exact Falsification Harness

`scripts/verify_one_reflected_bernoulli.py` checks the selection rule, the
three-term curvature identity, each Newton lower bound, both proof branches,
and the final quarter inequality using `fractions.Fraction` throughout. Its
default grid covers 47,850 rows from 4,785 eligible unique profiles, including
47 exact threshold equalities and the two calibration families above. The
2026-07-10 run had zero failures; its certificate is
`results/one_reflected_bernoulli_verify_2026-07-10.json`.

This finite exact scan is a falsification harness, not part of the proof.

## Scope and Next Obstruction

The proof uses more than a small-variance perturbation. It exploits the exact
two-term coefficient transform, which forces the descent to move by at most
one and supplies a mixed Newton curvature term.

For two reflected factors, a genuinely multi-step curvature identity is
required; naive iteration fails. The later universal theorem in
`notes/literature/universal_poisson_binomial_effective_drop_2026-07-10.md`
subsumes every finite number of reflected factors through Hillion--Johnson
cubic curvature propagation. This one-factor proof remains useful for its
explicit plateau-selection mechanism.
