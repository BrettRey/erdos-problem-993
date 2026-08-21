# A universal Pascal-smoothing lemma for erasure shadows

Date: 2026-08-20. Companion exact checks:
`verify_pascal_smoothing_20260820.py`.

## Statement

Let \(\Delta\) be any nonempty downward-closed family on at most \(M\) labelled
coordinates, and let \(a_j\) be its number of \(j\)-faces. Coordinates
that occur in no face are allowed. Define

\[
E_d=\sum_j a_j\binom{M-j}{d}. \tag{1}
\]

Equivalently, if

\[
F(t)=(1+t)^M A_\Delta\!\left(\frac{t}{1+t}\right),
\qquad A_\Delta(u)=\sum_j a_ju^j,
\]

then \(E_d=[t^{M-d}]F(t)\). The following holds without any graph,
poset, or forest hypothesis.

> **Pascal-smoothing lemma.** Put \(m=M-d\). For
> \(1\le d\le M-1\),
> \[
> \frac{E_d^2}{E_{d-1}E_{d+1}}
> \ge
> \frac89\frac{(d+1)(m+1)}{dm}. \tag{2}
> \]

Consequently:

1. every upper-tail log-concavity inequality through defect depth eight
   (\(d\le8\)) holds;
2. the entire coefficient sequence of \(F\) is log-concave whenever
   \(M\le33\);
3. at the depth needed in the current Erdős 993 frontier,
   \[
   E_3^2\ge
   \frac{32(M-2)}{27(M-3)}E_2E_4>E_2E_4. \tag{3}
   \]

The last inequality gives a quantitative reserve, not merely weak
log-concavity.

## Proof

Normalize the face numbers by the Boolean layers:

\[
b_j=\frac{a_j}{\binom Mj}.
\]

Double-count incidences between \(j\)-faces and \((j+1)\)-faces. Each
\((j+1)\)-face has \(j+1\) lower faces, while a \(j\)-face has at most
\(M-j\) one-element extensions. Hence

\[
(j+1)a_{j+1}\le(M-j)a_j,
\qquad b_{j+1}\le b_j. \tag{4}
\]

Set

\[
S_m=\sum_{j=0}^m\binom mj b_j.
\]

The absorption identity

\[
\binom Mj\binom{M-j}d=\binom Md\binom{M-d}j
\]

turns (1) into

\[
E_d=\binom Md S_{M-d}. \tag{5}
\]

It remains to control three adjacent Pascal transforms of a decreasing
sequence. Write

\[
\begin{aligned}
A&=\sum_j\binom{m-1}j b_j,\\
B&=\sum_j\binom{m-1}j b_{j+1},\\
C&=\sum_j\binom{m-1}j b_{j+2}.
\end{aligned}
\]

Then Pascal's identity gives

\[
S_{m-1}=A,\qquad S_m=A+B,\qquad S_{m+1}=A+2B+C. \tag{6}
\]

Monotonicity (4) gives \(0\le C\le B\le A\). With \(x=B/A\),

\[
\frac{S_m^2}{S_{m-1}S_{m+1}}
\ge\frac{(1+x)^2}{1+3x}\ge\frac89, \tag{7}
\]

because

\[
9(1+x)^2-8(1+3x)=(3x-1)^2.
\]

Combining (5) and (7) proves (2). The sufficient condition for the
\(d\)-th inequality is

\[
d(M-d)\le8(M+1). \tag{8}
\]

It holds automatically for \(d\le8\). It also holds for every \(d\)
when \(M\le33\), since

\[
d(M-d)\le\left\lfloor\frac{M^2}{4}\right\rfloor\le8(M+1).
\]

Taking \(d=3\) yields (3).

## Application to the maximum layer

For the order-ideal code of the free forest poset \(P\), the antichains
of \(P\) form a downward-closed family. The constant coordinates from
forced matching bags and singleton bags simply increase \(M\). Thus
the extendable erasure profile \((e_d)\) satisfies the whole lemma. In
the remaining rung \(M=\alpha\le19\), every adjacent inequality of the
extendable profile is therefore proved, not just the depth-3 one.

This removes the proposed forest-poset shadow lemma entirely. The only
remaining issue in the matching-bag model is the nonextendable correction

\[
s_d=e_d+b_d,
\]

where \(b_d\) counts blocked alternating-path assignments. Formula
(3) supplies explicit slack against that correction.

## Why this does not solve the original conjecture by itself

The face family used above is the antichain family associated with the
maximum-assignment code. Its transform counts precisely the partial
assignments that extend to a maximum assignment. The full independence
complex of a tree need not be pure: a partial independent set can be
maximal without extending to a maximum independent set. Those are
exactly the blocked terms \(b_d\). Applying (4) directly with
\(M=\alpha(T)\) would therefore be invalid; the ambient tree has
\(n\), not \(\alpha\), vertices.
