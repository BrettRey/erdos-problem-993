# A clan-graph cancellation attack on the adjacent two-hub polynomial

Date: 2026-08-20. Status: superseded later the same session by the
capacity repair in `notes/two_hub_B_logconcavity_proof_2026-08-20.md`.
This note preserves the local construction and the collision that forced
the correct Laurent-polynomial grouping. Companion arithmetic audit:
`verify_two_hub_clan_attack_20260820.py`.

## Target

Let \(D(a,b)\) be the tree formed from adjacent hubs \(u,v\), with
pendant paths of lengths in \(a\) at \(u\) and in \(b\) at \(v\).  Its
independence polynomial is

\[
B=Q_aQ_b+x(R_aQ_b+Q_aR_b).
\]

The missing Fact 1 in the C2 program is that \(B\) is log-concave for
all arm multisets.  This note extends the cancellation geometry in the
published spider proof of Li--Li--Yang--Zhang (arXiv:2501.04245,
Theorem 3.1) to the only genuinely new clan component: one containing
both adjacent hubs.

The local construction works with a strictly positive binomial reserve.
It is not yet a proof of Fact 1 because the most direct global map is
many-to-one.  The collision is exact and recorded below; the remaining
target is a capacitated injection into *stable partitions* of the
partner clan graph, rather than an injection into partner maps alone.

## Published framework being extended

For a map \(\alpha:V(G)\to\mathbb N\), let \(G^\alpha\) be its clan
graph and \(X_G^\alpha\) the normalized chromatic symmetric function.
Li--Li--Yang--Zhang use the identity

\[
Y_G=\sum_\alpha X_G^\alpha
\]

and prove that \(Y_G\) is 2-Schur-positive exactly when \(I(G;x)\) is
log-concave.  A connected bipartite clan component with bipartition
sizes \(K\ge L\) has one negative length-two Schur coefficient exactly
when \(K-L\ge2\): the coefficient of \(s_{(K-1,L+1)}\) is \(-1\).
Their spider proof pairs each such map with a map whose clan graph has
enough stable bipartitions of type \((K-1,L+1)\) to cancel it.

Any clan component with an edge is either 2-Schur-zero (a positive
multiplicity at least two meets another positive vertex and creates a
triangle), or every multiplicity in the component is one.  Away from
the hubs, every surviving component is therefore a path or an isolated
\(K_2\), and is 2-Schur-positive.  The only component not already
covered by the published spider argument is a unit-multiplicity induced
subtree containing both hubs.

## Exact imbalance law for a two-hub component

Color \(u\) positive and \(v\) negative.  In the component containing
both hubs, let

- \(p\) be the number of odd-length positive arm-prefixes at \(v\);
- \(q\) be the number of odd-length positive arm-prefixes at \(u\).

Even prefixes contribute zero to the signed bipartition difference.
An odd prefix at \(v\) contributes \(+1\), one at \(u\) contributes
\(-1\), and the adjacent hubs cancel.  Therefore

\[
K-L=|p-q|.
\]

This is the complete badness test: assume by symmetry that
\(p\ge q+2\).

## The local Li transform

On one hub side having at least two odd prefixes, apply the transform
from the published spider proof.

1. Choose the first shortest odd prefix \(A\).
2. Set the hub multiplicity to zero.
3. If \(|A|=1\), replace its sole vertex by an isolated \(K_2\)
   (multiplicity two).
4. If \(|A|=2t+1\ge3\), replace \(A\) by the pattern
   \(2,0,2,0,\ldots,2\), and on the first \(2t\) vertices of a second
   odd prefix of length at least \(|A|\), use
   \(2,0,\ldots,2,0\), leaving the next unit vertex in place.

The full transformed prefix gains one unit of total multiplicity, which
exactly replaces the deleted hub.  All new clan components are paths or
isolated \(K_2\)'s.  The asymmetric full/partial pair is a prefix code:
the two runs contain respectively \(t+1\) and \(t\) copies of \(K_2\),
so the original unit prefixes and the hub can be recovered.  This is
the injectivity mechanism in Li--Li--Yang--Zhang's proof.

## Local cancellation for every two-hub imbalance

### Case 1: \(q\le2\)

Apply the local Li transform only at \(v\).  After deleting \(v\), the
component at \(u\) is a spider with \(q\) odd prefixes.  Its two color
classes differ by \(|1-q|\le1\), so its individual chromatic symmetric
function is 2-Schur-positive.  Everything else is a path or a \(K_2\).

Counting orientations of the odd components gives the following
surplus of stable bipartitions of difference \(p-q-2\) over those of
difference \(p-q\):

\[
\begin{array}{c|c}
q&\text{surplus}\\ \hline
0&\binom p1-\binom p0=p-1,\\
1&\binom{p-1}1-\binom{p-1}0=p-2,\\
2&\binom p2-\binom p1=p(p-3)/2.
\end{array}
\]

The assumptions \(p\ge q+2\) make every quantity strictly positive.
Factors from balanced path components are positive and factors from
the orientations of cloned \(K_2\)'s cancel their normalization
factorials, exactly as in the published proof.

### Case 2: \(q\ge3\)

Now both sides have at least three odd prefixes.  Apply the local Li
transform independently at both hubs.  Both hubs disappear; one odd
path component is replaced by balanced \(K_2\)'s on each side.  There
remain

\[
r=(p-1)+(q-1)=p+q-2
\]

odd path components.  A stable bipartition of difference \(p-q\)
chooses \(q-1\) of them against the majority orientation, whereas one
of difference \(p-q-2\) chooses \(q\).  The exact surplus is

\[
\binom{p+q-2}{q}-\binom{p+q-2}{q-1}>0,
\]

because \(p\ge q+2\) places \(q\) on the increasing side of the
binomial row.  Equivalently, the ratio of the first binomial to the
second is \((p-1)/q>1\).

Thus every bad two-hub component has a locally constructed partner with
strictly more positive capacity at its unique negative Schur index.

## Why this is not yet a proof: the global collision

The local map is injective when only one spider center is in scope, but
its direct extension to the whole double-spider map space is not.
Consider five unit arms at \(v\) and three at \(u\).

- A double-hub bad map puts multiplicity one on both hubs and all eight
  leaves.  The double transform (the length-one case on each side)
  sends it to: both hubs zero, one leaf of multiplicity two on each
  side, and the other six leaves of multiplicity one.
- A single-hub bad map puts multiplicity one on \(v\) and its five
  leaves, zero on \(u\), and multiplicities \((2,1,1)\) on the three
  \(u\)-leaves.  The single transform at \(v\) has exactly the same
  image.

Both sources and the image have total multiplicity ten.  Therefore a
map-level injection cannot simply combine the published one-center map
with the new two-center map.

This collision does **not** kill the method.  The common image is a
disjoint union of two \(K_2\)'s and six isolated vertices, with many
stable bipartitions.  The local calculations above measured precisely
this multiplicity but then discarded it by mapping only to the clan
map.  The correct target is an enriched pair

\[
(\beta,\ \text{a positive stable partition counted in }
[s_{(K-1,L+1)}]X_G^\beta).
\]

The next proof obligation is a sign-reversing, capacity-respecting
injection from negative stable partitions of all bad maps into these
enriched targets.  The binomial surpluses prove that enough *local*
capacity exists.  What remains is to assign the capacities consistently
when single- and double-hub sources collide.

## Consequences and priority

1. The hard part of the two-hub extension is not finding a positive
   clan partner; one exists with an explicit strict binomial reserve for
   every imbalance.
2. The sole identified gap is global bookkeeping across overlapping
   image strata.  This is finite combinatorial injection design, not an
   analytic inequality over arm lengths.
3. A successful capacitated injection would prove the open
   log-concavity of \(B\), one of the three remaining C2 base facts,
   without LAW-V or LAW-W.
4. The imbalance law \(K-L=|p-q|\) and the two transforms are independent
   of arm lengths except for the published prefix code.  This makes the
   route naturally formalizable after the capacity assignment is fixed.
