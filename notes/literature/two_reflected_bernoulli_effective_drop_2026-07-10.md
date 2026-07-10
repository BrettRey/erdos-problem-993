# Two Reflected Bernoullis: Effective-Drop Theorem

Date: 2026-07-10
Status: proved and independently audited, conditional only on the established
one-sided localization \(S+1\le4V_X\)

## Theorem

Let

\[
X=\sum_{i=1}^m\operatorname{Bernoulli}(p_i),
\qquad 0<p_i\le\frac12,
\qquad V_X=\operatorname{Var}(X)\ge1,
\]

and let

\[
Y=\operatorname{Bernoulli}(q_1)+\operatorname{Bernoulli}(q_2),
\qquad 0\le q_1,q_2\le\frac12,
\]

be independent of \(X\). Write

\[
c_k=\Pr(X-Y=k),\qquad
V=V_X+\sum_{i=1}^2q_i(1-q_i).
\]

If \(D\) is the first strict descent of the signed mass sequence, then
\(c_D>0\) and

\[
\boxed{
V\left(1-\frac{c_{D+1}c_{D-1}}{c_D^2}\right)\ge\frac14.
}
\]

The proof is not an iteration of the one-reflected-Bernoulli theorem. It uses
a new three-term curvature identity.

## Descent Localization

Let

\[
a_k=\Pr(X=k),\qquad r_k=\frac{a_{k+1}}{a_k},
\]

and let \(S\) be the first strict descent of \(a\). The ratios are decreasing,
and the established one-sided results give

\[
S\ge2,\qquad S\le m-1,\qquad S+1\le4V_X,
\]

\[
r_{S-2}\ge1,\qquad r_{S-1}<1.
\]

Put

\[
t_i=\frac{q_i}{1-q_i}\in[0,1],\qquad
e=t_1+t_2,\qquad f=t_1t_2.
\]

After deleting the positive normalizing factor
\((1-q_1)(1-q_2)\), the signed coefficients are

\[
\widetilde c_k=a_k+e\,a_{k+1}+f\,a_{k+2}.
\]

Extend \(a_k\) by zero outside its support. The exact first-difference identity
is

\[
\widetilde c_k-\widetilde c_{k-1}
=(a_k-a_{k-1})+e(a_{k+1}-a_k)+f(a_{k+2}-a_{k+1}).
\]

Every term is nonnegative for \(k\le S-3\). At \(k=S\), all three
parenthesized differences are negative and the first has coefficient one.
Consequently,

\[
\boxed{D\in\{S-2,S-1,S\}.}
\]

In particular, \(D\ge0\), and

\[
c_D\ge(1-q_1)(1-q_2)a_D>0.
\]

The selection rule is explicit. Define

\[
\begin{aligned}
A&=a_{S-2}-a_{S-3},&
B&=a_{S-1}-a_{S-2},\\
C&=a_{S-1}-a_S,&
E&=a_S-a_{S+1}.
\end{aligned}
\]

Here \(A,B\ge0\) and \(C,E>0\). Put

\[
\delta_0=A+eB-fC,
\qquad
\delta_1=B-eC-fE.
\]

Then

\[
D=
\begin{cases}
S-2,&\delta_0<0,\\
S-1,&\delta_0\ge0\text{ and }\delta_1<0,\\
S,&\delta_0\ge0\text{ and }\delta_1\ge0.
\end{cases}
\]

Equality in either shoulder expression is a plateau, so the first *strict*
descent waits for the next candidate.

## The Leftmost Branch

If \(D=S-2\), shift by two:

\[
d_j=c_{j-2}.
\]

This is the mass sequence of the ordinary Poisson-binomial sum

\[
X+(1-Y_1)+(1-Y_2).
\]

Its descent index is \(S\). The full Newton factor for its \(m+2\) Bernoulli
terms gives

\[
\Delta_D(c)=\Delta_S(d)
\ge\frac{m+3}{(S+1)(m-S+3)}
\ge\frac1{S+1}.
\]

Since \(S+1\le4V_X\le4V\), the quarter bound follows.

The other two branches require a shift-invariant curvature estimate.

## Three-Term Curvature Identity

Fix \(k\ge0\), and abbreviate

\[
N_j=a_{k+j}^2-a_{k+j-1}a_{k+j+1}
\qquad(j=0,1,2),
\]

\[
C_{01}=a_ka_{k+1}-a_{k-1}a_{k+2},
\]

\[
C_{12}=a_{k+1}a_{k+2}-a_ka_{k+3},
\]

\[
L=a_{k+1}^2-a_{k-1}a_{k+3}.
\]

Direct expansion and regrouping give

\[
\boxed{
\widetilde c_k^2-\widetilde c_{k-1}\widetilde c_{k+1}
=
N_0+(t_1^2+t_2^2)N_1+fL+f^2N_2+eC_{01}+efC_{12}.
}
\]

Every term on the right is nonnegative. More precisely, put \(K=k+1\).
Consecutive Newton ratio drops give

\[
N_0\ge\frac{a_k^2}{K},\qquad
N_1\ge\frac{a_{k+1}^2}{K+1},\qquad
N_2\ge\frac{a_{k+2}^2}{K+2},
\]

\[
C_{01}\ge\frac{2a_ka_{k+1}}{K+1},\qquad
C_{12}\ge\frac{2a_{k+1}a_{k+2}}{K+2}.
\]

For the four-step gap,

\[
\frac{a_{k-1}a_{k+3}}{a_{k+1}^2}
\le
\frac{k(k+1)}{(k+2)(k+3)},
\]

so

\[
L\ge
\frac{2(2K+1)}{(K+1)(K+2)}a_{k+1}^2.
\]

## Normalized Curvature Bound

Set

\[
x=\frac{a_{k+1}}{a_k},\qquad
y=\frac{a_{k+2}}{a_{k+1}}.
\]

At either of the two remaining descent indices, \(k\ge S-1\), so

\[
0\le y\le x<1.
\]

After dividing the preceding lower bound by

\[
\widetilde c_k^2
=a_k^2(1+ex+fxy)^2,
\]

we obtain \(\Delta_k(\widetilde c)\ge G_K(t_1,t_2,x,y)\), where

\[
\begin{aligned}
G_K(u,v,x,y)
=\frac{1}{(1+(u+v)x+uvxy)^2}\Bigg[
&\frac1K+\frac{2(u+v)x}{K+1}
+\frac{(u^2+v^2)x^2}{K+1}\\
&+\frac{2uv(2K+1)x^2}{(K+1)(K+2)}
+\frac{2uv(u+v)x^2y}{K+2}
+\frac{u^2v^2x^2y^2}{K+2}
\Bigg].
\end{aligned}
\]

The next lemma is the exact variance payment:

\[
\boxed{
\left(
K+4\frac{u}{(1+u)^2}+4\frac{v}{(1+v)^2}
\right)G_K(u,v,x,y)\ge1
}
\]

for \(K\ge1\) and \(0\le y\le x\le1\), \(0\le u,v\le1\).

## Monotonicity Reduction

The derivative in \(y\) is

\[
\frac{\partial G_K}{\partial y}
=
-\frac{2uvx\,H}
{K(K+1)(K+2)(1+(u+v)x+uvxy)^3},
\]

where

\[
\begin{aligned}
H={}&2K^2uvx^2-K^2uvxy+K^2ux+K^2vx+K^2\\
&+Ku^2x^2-Kuvxy+3Kux+Kv^2x^2+3Kvx+3K+2.
\end{aligned}
\]

Because \(y\le x\),

\[
2K^2uvx^2-K^2uvxy\ge K^2uvx^2,
\]

and

\[
Ku^2x^2+Kv^2x^2-Kuvxy
\ge K(u^2-uv+v^2)x^2\ge0.
\]

Thus \(H\ge0\), and \(G_K\) decreases in \(y\).

On the diagonal \(y=x\),

\[
\frac{d}{dx}G_K(u,v,x,x)
=-\frac{2J}
{K(K+1)(K+2)(1+ux)^3(1+vx)^3},
\]

where

\[
\begin{aligned}
J={}&Kuv(u-v)^2x^3+3Kuv(u+v)x^2\\
&+(8K+4)uvx+(K+2)(u+v)\ge0.
\end{aligned}
\]

Therefore

\[
G_K(u,v,x,y)\ge G_K(u,v,x,x)\ge G_K(u,v,1,1).
\]

## Exact Bernstein Certificate

Let

\[
\beta_i(t)=\binom4i t^i(1-t)^{4-i}.
\]

Exact simplification gives

\[
\left(
K+4\frac{u}{(1+u)^2}+4\frac{v}{(1+v)^2}
\right)G_K(u,v,1,1)-1
=
\frac{K^2C_2(u,v)+KC_1(u,v)+C_0(u,v)}
{K(K+1)(K+2)(1+u)^4(1+v)^4}.
\]

Each \(C_r\) has a nonnegative degree-\((4,4)\) Bernstein expansion

\[
C_r(u,v)=\sum_{i,j=0}^4 M^{(r)}_{ij}\beta_i(u)\beta_j(v).
\]

With rows and columns indexed by \(0,\ldots,4\), the exact coefficient
matrices are

\[
M^{(2)}=
\begin{pmatrix}
0&1/2&3/2&3&4\\
1/2&15/8&55/12&9&14\\
3/2&55/12&191/18&62/3&100/3\\
3&9&62/3&40&64\\
4&14&100/3&64&96
\end{pmatrix},
\]

\[
M^{(1)}=
\begin{pmatrix}
0&2&5&9&12\\
2&15/2&16&29&46\\
5&16&581/18&57&272/3\\
9&29&57&98&152\\
12&46&272/3&152&224
\end{pmatrix},
\]

\[
M^{(0)}=
\begin{pmatrix}
0&2&4&6&8\\
2&6&31/3&15&20\\
4&31/3&52/3&25&100/3\\
6&15&25&36&48\\
8&20&100/3&48&64
\end{pmatrix}.
\]

All entries are nonnegative, proving the scalar inequality on the full square
\([0,1]^2\).

## Completion of the Proof

If \(D=S\), take \(k=D\) and \(K=S+1\). The curvature lemma and one-sided
localization give

\[
\Delta_D(c)\ge
\frac1{K+4q_1(1-q_1)+4q_2(1-q_2)}
\ge\frac1{4V}.
\]

If \(D=S-1\), then \(K=S\), while

\[
4V\ge S+1+4q_1(1-q_1)+4q_2(1-q_2)
>K+4q_1(1-q_1)+4q_2(1-q_2).
\]

The same curvature lemma again gives \(V\Delta_D(c)\ge\tfrac14\).
Together with the shifted-Newton argument for \(D=S-2\), this proves the
theorem.

## Why Naive Iteration Fails

The theorem cannot be obtained by simply applying the one-factor result
twice. For example, start with

\[
X\sim\operatorname{Bin}(4,\tfrac12)
\]

and add one reflected \(\operatorname{Bernoulli}(q)\), with
\(0<q<\tfrac12\). The first transform can leave the relevant descent in
place, but the shifted intermediate law has descent \(T=4\) and

\[
T+1=5>4\bigl(1+q(1-q)\bigr).
\]

Thus the one-sided localization hypothesis needed for a second application
has already failed. The direct multi-term curvature identity is essential.

The subsequent universal theorem in
`notes/literature/universal_poisson_binomial_effective_drop_2026-07-10.md`
uses a different cubic-curvature propagation argument and subsumes this
two-factor result. The present proof remains an independently audited finite
model of how long-gap Newton minors pay for reflected variance.

## Exact Falsification Harness

`scripts/verify_two_reflected_bernoulli.py` checks:

1. the three-place descent localization;
2. the six-term curvature identity;
3. every Newton lower bound, including the four-step gap;
4. the two derivative identities;
5. reconstruction of all three Bernstein matrices;
6. both the proof lower bound and the actual quarter inequality on exhaustive
   rational grids.

The default 2026-07-10 run checked 97,488 exact rows from 2,708 eligible
profiles, including both observed descent branches and 71 exact shoulder
plateaus. It performed 584,928 individual Newton/long-gap checks, reconstructed
all 75 Bernstein coefficients, and found zero failures. The minimum observed
scaled drop was

\[
\frac{1302716}{2970375}\approx0.43856954.
\]

No \(D=S-2\) row appeared in this grid or a separate 200,000-profile random
probe. That branch may be structurally impossible, but the theorem does not
claim this sharpening. The result certificate is
`results/two_reflected_bernoulli_verify_2026-07-10.json`.

The verifier is a falsification artifact and is not part of the proof.
