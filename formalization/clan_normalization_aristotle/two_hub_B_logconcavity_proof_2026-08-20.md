# Log-concavity of the adjacent two-hub independence polynomial

Date: 2026-08-20. Status: proof draft, exact core audit passing, with the
normalization kernel formally verified but the global partition still open.
Companion verifier: `verify_two_hub_B_logconcavity_20260820.py`. This repairs
the map-level collision isolated in
`notes/two_hub_clan_cancellation_attack_2026-08-20.md`.
The local transform/inverse, including the new \(p=2\) use, has a
separate finite audit in `verify_local_li_block_20260820.py`.

The independent Aristotle audit in
`formalization/clan_normalization_aristotle/RESULT.md` proves the clan
normalization, normalized imbalance formula, binary product law, local
`p=2` injection, four-map identity, and the needed block inequality for
derived scalars `c,d>=1`. It refutes the unnecessarily broad version for
arbitrary positive rationals and does not formalize the arbitrary-even-arm
weight derivation or the exhaustive global block partition asserted in
Sections 3--5. Thus the theorem statement below remains the target of this
proof draft rather than a fully independently closed theorem.

## Theorem

Let \(D(a,b)\) be a tree consisting of adjacent vertices \(u,v\), with
any finite multisets of pendant paths attached at \(u\) and \(v\).  Then
the independence polynomial of \(D(a,b)\) is log-concave.

Equivalently, in the notation of the C2 program,

\[
B=Q_aQ_b+x(R_aQ_b+Q_aR_b)
\]

is log-concave for every pair of arm multisets.  This proves Fact 1 of
`notes/c2_bounded_pendant_core_2026-08-08.md`.

The proof extends the 2-Schur cancellation used by Li, Li, Yang, and
Zhang for spiders (arXiv:2501.04245, Theorem 3.1; source text on disk at
`notes/literature/li_yang_zhang_zhang_2025_symmetric_function_log_concavity.txt`).

## 1. Two-row Schur coefficients as an imbalance sequence

We use the notation and results of Li--Li--Yang--Zhang.  For
\(\alpha:V(G)\to\mathbb N\), let \(G^\alpha\) be the clan graph and
\(X_G^\alpha\) its normalized chromatic symmetric function.  They prove

\[
Y_G=\sum_\alpha X_G^\alpha,
\]

and that \(Y_G\) is 2-Schur-positive if and only if \(I(G;x)\) is
log-concave (their Corollary 2.2).

Here is a convenient reformulation of their Proposition 2.4.  Suppose
every component of \(G^\alpha\) is bipartite.  For each component choose
one orientation of its bipartition and let \(\delta_j\ge0\) be the
difference between its color-class sizes.  Up to the positive
normalization factors from cloned isolated vertices, the Laurent
polynomial

\[
L_\alpha(z)=\prod_j(z^{\delta_j}+z^{-\delta_j})
\tag{1}
\]

counts its semi-ordered stable bipartitions by signed size difference.
If \(A_e=[z^e]L_\alpha(z)\), Proposition 2.4 says that the coefficient
at the two-row partition whose parts differ by \(e\) is

\[
A_e-A_{e+2}.\tag{2}
\]

Consequently the two-row part is Schur-positive exactly when the
symmetric coefficient sequence of \(L_\alpha\) is nondecreasing from
either extreme toward zero.  We call this *central unimodality*.

If a clan component is nonbipartite, its two-row Schur truncation is
zero (Li--Li--Yang--Zhang, Corollary 2.5), so it may be discarded.
Otherwise every nontrivial component involving a positive-multiplicity
edge has multiplicity one at all its vertices: multiplicity at least
two next to a positive vertex creates a triangle.  Components away from
the two hubs are therefore paths or isolated \(K_2\)'s.

We shall use twice the elementary fact that the convolution of two
nonnegative symmetric centrally unimodal sequences is again symmetric
and centrally unimodal.  One proof decomposes each sequence into a
nonnegative linear combination of centered interval indicators; the
convolution of two centered intervals is a symmetric trapezoid.

## 2. The local two-state block at one hub

Fix all multiplicities on the pendant paths at one hub \(w\), and
suppose \(\alpha(w)=1\).  Its positive component contains a unit prefix
of each incident arm.  Let \(p\) be the number of odd-length prefixes.
The connected spider component has bipartition imbalance

\[
r=p-1.\tag{3}
\]

When \(p\ge2\), apply the local map from the published spider proof:
delete the hub, replace the first shortest odd prefix by
\(2,0,2,0,\ldots,2\), and, except in their length-one case, use the
shorter \(2,0,\ldots,2,0\) marker on a second odd prefix.  The transform
preserves total multiplicity and is injective; the unequal full and
marker runs recover the original prefixes and hub.  The published proof
states the map for its bad range \(p\ge3\), but the same definition and
inverse work unchanged for \(p=2\).

Call the hub-present state \(H_p\) and its hub-absent partner \(P_p\).
Their imbalance Laurent polynomials, after factoring all common outside
path components, are

\[
H_p(z)=z^r+z^{-r},\qquad
P_p(z)=c(z+z^{-1})^r,\qquad r=p-1,\ c\ge1.\tag{4}
\]

Indeed, the transform removes one odd component and leaves \(p-1\)
odd path components.  Even paths contribute only a positive scalar;
the two orientations of every cloned isolated \(K_2\) cancel its clan
normalization factorial, exactly as in the published calculation.  Thus
\(c\) is a positive integer (a power of two from balanced ordinary path
components).

The two-state sum

\[
A_{r,c}(z)=c(z+z^{-1})^r+z^r+z^{-r}\tag{5}
\]

is centrally unimodal.  In coefficient notation it is the binomial row
\(c\binom rj\), with one added at each endpoint.  If \(r=1\) it is
constant on its support.  If \(r\ge2\), its first inward coefficient is
\(cr\ge c+1\), and the remaining binomial coefficients increase toward
the center.

Because the local transform is injective and changes the hub state from
one to zero, the local state space decomposes into disjoint two-state
blocks \(\{H_p,P_p\}\) for \(p\ge2\), plus unmatched states.  Every
unmatched hub-present state has \(p\le1\) and hence imbalance at most
one; unmatched hub-absent states are disjoint unions of paths.  They are
already 2-Schur-positive.

## 3. Two active blocks: the four-map identity

Take an active block on each side, with

\[
r=p-1\ge1,\qquad s=q-1\ge1,
\]

and positive constants \(c,d\).  There are four global clan maps,
according as each hub is present (\(H\)) or replaced by its local
partner (\(P\)).  In the first three cases containing an absent hub the
two sides are disconnected, so their Laurent polynomials multiply.
When both hubs are present, the central edge forces their colorings to
be opposite.  Its connected bipartition imbalance is \(|p-q|=|r-s|\).
The sum of the four Laurent polynomials is therefore

\[
\begin{aligned}
F(z)={}&P_pP_q+H_pP_q+P_pH_q
       +(z^{r-s}+z^{s-r})\\
={}&A_{r,c}(z)A_{s,d}(z)
       -(z^{r+s}+z^{-(r+s)}).\tag{6}
\end{aligned}
\]

The identity uses

\[
(z^r+z^{-r})(z^s+z^{-s})
=z^{r+s}+z^{-(r+s)}+z^{r-s}+z^{s-r}.
\]

Both factors in the product on the second line of (6) are centrally
unimodal by (5), so their product is.  The subtracted monomials affect
only the two extreme coefficients, reducing each by one.  This can only
strengthen the inward inequalities.  Hence every four-map block is
2-Schur-positive.

This is the capacity repair: a map-level partner can be shared by the
two single-hub sources and the double-hub source, but the *sum of all
four maps* is positive.  No target map is counted twice.

## 4. One active block and one unmatched side

If the other hub is absent, the two-state polynomial (5) is merely
multiplied by the centrally unimodal Laurent polynomial of a path
forest, so positivity is preserved.

Suppose instead that the other hub is present but unmatched.  It has
\(q\le1\) odd prefixes.

- If \(q=0\), its component has imbalance one.  Joining the hubs and
  adding the partner map gives, up to positive common factors,
  \[
  c(z+z^{-1})^{r+1}+z^{r+1}+z^{-(r+1)},
  \]
  which is (5) with exponent \(r+1\).
- If \(q=1\), its component is balanced.  The corresponding sum is a
  positive scalar multiple of \((z+z^{-1})^r\), with one added at each
  endpoint, again of form (5).

Thus every two-map block in this case is 2-Schur-positive.

## 5. No active block

If neither side is active, every present hub has at most one odd arm
prefix.  A component containing both hubs then has imbalance
\(|p-q|\le1\).  All other components are paths or isolated \(K_2\)'s.
Every remaining clan map is individually 2-Schur-positive (or
2-Schur-zero if it has a nonbipartite clan component).

Sections 3--5 partition all clan maps into disjoint positive blocks.
Therefore \(Y_{D(a,b)}\) is 2-Schur-positive.  Corollary 2.2 of
Li--Li--Yang--Zhang now implies that \(I(D(a,b);x)=B\) is log-concave.
\(\square\)

## Verification and scope

The companion verifier checks the symmetric-unimodality argument in
(5)--(6) for \(1\le r,s\le60\) and \(1\le c,d\le6\), checks the two
one-active cases, and independently recomputes \(B\)'s Turan gaps over
3,330 arm pairs of bounded total weight using the exact C2 engine.  The
finite checks audit the algebra; the proof above is uniform in the arm
lengths and counts.

This theorem supplies the adjacent base fact. The later proof draft
`c2_connector_clan_reduction_2026-08-20.md` extends the same partition to
every connector length and bypasses LAW-V, LAW-W, and the LAW-U head strip.
The independent audit has closed the normalization and algebraic parts, but
the arbitrary-even-arm weight derivation and exhaustive global block
partition remain before either theorem is treated as final.
