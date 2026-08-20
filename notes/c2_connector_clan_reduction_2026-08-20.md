# Clan cancellation proves the two-branch-vertex theorem

Date: 2026-08-20. Status: proof draft complete; exact core verification
passing. The even-connector binomial lemma left open in the first draft
is proved below. The 14,400-case exact audit of its Laurent core is in
`verify_two_hub_B_logconcavity_20260820.py`.

This note starts from
`notes/two_hub_B_logconcavity_proof_2026-08-20.md` and asks whether the
same proof can bypass LAW-V/LAW-W and establish log-concavity directly
for every tree with two branch vertices joined by a path.

## Imbalance across a connector

Let the connector have \(\ell\) edges, color the first hub positive,
and let \(p,q\) be the numbers of odd positive arm-prefixes at its two
ends.  A unit connector component containing both hubs has signed
bipartition difference

\[
D=\sum_{j=0}^{\ell}(-1)^j-p+(-1)^{\ell+1}q.
\]

Hence

\[
|D|=\begin{cases}
|p-q|,&\ell\text{ odd},\\
|1-p-q|,&\ell\text{ even}.
\end{cases}\tag{1}
\]

Apply the local Li transforms at the two ends as in the adjacent proof.
Write \(r=p-1\), \(s=q-1\), and let the two partner constants be
\(c,d\ge1\).  Put \(X=z+z^{-1}\) and \(H_j=z^j+z^{-j}\).

## Connector-state partition

Equation (1) is needed only for clan maps whose positive component
contains both hubs.  In every nonzero two-row contribution, that forces
every internal connector multiplicity to be one.  If a connector
multiplicity is at least two next to a positive vertex, the clan
component contains a triangle and is 2-Schur-zero.  If a zero (or an
isolated cloned vertex) breaks the connector, every remaining component
contains at most one hub.  It is a spider, a path, or an isolated
\(K_2\).

Partition first at the connector breaks.  On each separated side apply
the published local spider transform, allowing the positive connector
prefix to serve as an ordinary arm.  The two sides have disjoint support,
so their two-state Laurent blocks tensor; convolution preserves central
unimodality.  Consequently all broken-connector maps are already
2-Schur-positive in disjoint local blocks.  The sections below address
the sole new stratum: the connector is unit-positive from end to end.
Here the local transforms may be chosen on pendant odd prefixes, and the
connector contributions are exactly those in (1).

## Odd connector lengths: closed

For odd \(\ell\), the four hub-state maps have the same imbalance core
as the adjacent case.  If \(\ell=1\), their sum is

\[
(cX^r+H_r)(dX^s+H_s)-H_{r+s}.
\]

If \(\ell\ge3\), the both-hubs-absent state also contains a nonempty
balanced internal connector path.  It contributes one additional
positive multiple of \(X^{r+s}\).  The adjacent core is centrally
unimodal by the four-map theorem, and convolution/addition of the
balanced path term preserves central unimodality.  Thus the clan proof
establishes log-concavity for every odd connector length.  The
one-active blocks are unchanged up to the same harmless balanced-path
factor: an unmatched present endpoint has zero or one odd arm-prefix,
so its two-state sum is again a binomial row with its two endpoints
raised.

## Even connector lengths: one algebraic lemma

For even \(\ell\), the internal connector left after deleting both hubs
has imbalance one.  The active four-map core becomes

\[
F_{r,s;c,d}(z)
=cdX^{r+s+1}+dX^sH_{r+1}+cX^rH_{s+1}+H_{r+s+1}.
\tag{2}
\]

All dependence on actual connector length and arm lengths has vanished;
only \(r,s,c,d\) remain.  The full C2 theorem therefore follows from:

**Even-connector lemma.** For all integers \(r,s,c,d\ge1\), the
Laurent coefficient sequence of (2) is centrally unimodal.

There is no additional endpoint case hidden here.  With only one
active local block, deleting the active hub leaves either the odd
internal connector path or the unmatched endpoint component.  For
\(q=0,1\) the resulting two-map sum is again
\(cX^m+H_m\) for the appropriate \(m\), hence centrally unimodal.  If
neither endpoint is active, (1) gives imbalance at most one.  Thus (2)
is the sole even-parity obligation.

The constants reduce cleanly.  With
\(A=X^{r+s+1}\), \(B=X^sH_{r+1}\),
\(C=X^rH_{s+1}\), \(D=H_{r+s+1}\),

\[
cdA+dB+cC+D
=(c-1)(d-1)A+(c-1)(A+C)+(d-1)(A+B)+(A+B+C+D).
\]

The first three summands are centrally unimodal: for example
\(A+C=X^r(X^{s+1}+H_{s+1})\), a binomial row convolved with a binomial
row whose two endpoints have each been raised by one.  Thus it is
enough to prove the unit-constant case

\[
X^{r+s+1}+X^sH_{r+1}+X^rH_{s+1}+H_{r+s+1}.
\tag{3}
\]

In ordinary rank variable \(y\), (3) has coefficient polynomial

\[
(1+y)^N+(1+y)^s(1+y^{r+1})
 +(1+y)^r(1+y^{s+1})+1+y^N,
\quad N=r+s+1.\tag{4}
\]

### Proof of the even-connector lemma

It remains, by the constants reduction, to prove that (4) is unimodal.
Assume \(r\ge s\) by symmetry and write

\[
\Delta_n(j)=\binom n{j+1}-\binom n j.
\]

We use two elementary estimates.  The quantity \(\Delta_n(j)\) is
increasing in \(n\) once \(n\ge2j+1\), directly from

\[
\Delta_n(j)=\binom n j\frac{n-2j-1}{j+1}.
\tag{5}
\]

If \(r\ge s\ge j+1\), put \(N=r+s+1\).  Since
\(N\ge r+j+2\), Vandermonde and (5) give

\[
\begin{aligned}
\frac{\Delta_N(j)}{\binom rj}
&\ge \frac{r-j+1}{j+1}
 \left(1+(j+2)\frac{\binom r{j-1}}{\binom rj}\right)\\
&=\frac{r+j^2+j+1}{j+1}\ge2.
\end{aligned}
\tag{6}
\]

We also need the boundary version: if \(r\ge j+1\), then

\[
\Delta_{r+j+1}(j)\ge\binom rj.\tag{7}
\]

The same calculation, now using
\(\binom{r+j+1}j\ge\binom rj+(j+1)\binom r{j-1}\), reduces (7), after
writing \(r=j+1+t\), to

\[
j^2t+j^2-j+t^2+2t\ge0.
\]

Let \(f_j\) be the coefficient of \(y^j\) in (4), and consider
\(f_{j+1}-f_j\) before the middle.

- At \(j=0\), direct substitution gives
  \(f_1-f_0=2r+2s-3>0\).
- For \(1\le j<s\), none of the shifted copies has begun, so
  \[
  f_{j+1}-f_j=\Delta_N(j)+\Delta_r(j)+\Delta_s(j).\tag{8}
  \]
  If at least one of the last two terms is nonnegative, group the other
  with \(\Delta_N\).  The relevant grouped polynomial is
  \((1+y)^s((1+y)^{r+1}+1+y^{r+1})\), or its symmetric counterpart;
  it is a convolution of symmetric unimodal sequences, so that grouped
  difference is nonnegative.  If both terms are negative, their total
  magnitude is at most
  \(\binom rj+\binom sj\le2\binom rj\), which (6) says is at most
  \(\Delta_N(j)\).  Thus (8) is nonnegative.
- At \(j=s<r\), the disappearing last coefficient of \((1+y)^s\)
  cancels the first coefficient of the shifted copy of \((1+y)^r\).
  Hence the difference is \(\Delta_N(s)+\Delta_r(s)\).  It is trivial
  if the second term is nonnegative; otherwise its magnitude is at most
  \(\binom rs\), and (7) applies because \(N=r+s+1\).
- For \(s<j<N/2\), the two terms involving the \(s\)-row vanish at both
  adjacent indices.  The remaining coefficients agree locally with
  \[
  (1+y)^r\bigl((1+y)^{s+1}+1+y^{s+1}\bigr),
  \]
  a convolution of two symmetric unimodal sequences.  Their difference
  is therefore nonnegative.
- If \(r=s\), the same argument covers \(j<s\), while the two central
  coefficients \(f_s,f_{s+1}\) are equal by symmetry.

Thus (4) increases to its center.  It is symmetric, so it is centrally
unimodal.  The constants decomposition proves (2) for all
\(c,d\ge1\). \(\square\)

Combining the odd-connector argument, the even-connector lemma, and the
unmatched endpoint cases proves: **the independence polynomial of every
tree with at most two vertices of degree at least three is
log-concave.** This is the C2 theorem sought in
`notes/c2_bounded_pendant_core_2026-08-08.md`, obtained without LAW-V,
LAW-W, or the connector recurrence.

## Verification and scope

`verify_two_hub_B_logconcavity_20260820.py` checks (2) for all
\(1\le r,s\le30\) and \(1\le c,d\le4\), totaling 14,400 exact integer
cores. The same verifier audits the adjacent four-map block and 3,330
bounded arm pairs. These computations check the displayed algebra; they
do not replace the requested independent review of the clan-map partition
and normalization.
