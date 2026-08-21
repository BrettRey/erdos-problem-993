# Maximum-matching bags: a small-defect model for the depth-3 frontier

Date: 2026-08-20. Companion exact verifier:
`verify_matching_bag_csp_20260820.py`.

## Representation theorem

Let \(T\) be a tree, let \(M\) be any maximum matching, and contract
each edge of \(M\).  Leave every unmatched vertex as a singleton bag.
The quotient \(S\) is again a tree.  Its bags have size two (matched
edges) or one (unmatched vertices), and

\[
|V(S)|=n-|M|=\alpha(T)
\]

by Konig's theorem.

At each two-vertex bag choose one endpoint or choose the empty state;
at a singleton bag choose its vertex or the empty state.  Every quotient
edge remembers the original endpoint (the *port*) at each end.  The
only constraint is that the two incident port states may not be chosen
simultaneously.  Then

\[
I(T;x)=\sum_{\sigma\text{ feasible}}x^{\#\text{nonempty bags}}.
\tag{1}
\]

This is a bijection: a tree independent set selects at most one endpoint
of each matching edge, and every feasible bag assignment is an
independent set.  Thus the reversed coefficient

\[
i_{\alpha-d}(T)
\]

counts feasible assignments with exactly \(d\) empty bags.

The number of singleton bags is the matching deficiency

\[
\delta=n-2|M|=2\alpha-n.\tag{2}
\]

## Consequence for the distance-3 window

The remaining depth-3 rung has \(\alpha\in\{17,18,19\}\) and
\(33\le n\le38\).  In the bag model it is therefore a port-labeled tree
on only 17--19 bags, with at most five singleton bags.  A depth-3 LC
break is exactly

\[
s_3^2<s_2s_4,
\qquad s_d=\#\{\text{feasible assignments with }d\text{ empty bags}\}.
\tag{3}
\]

This replaces a 531-billion-tree question by a bounded structural
question about the first five defect layers of a binary tree CSP.

## Extendable and blocked defects

Fix an empty bag set \(D\).  The occupied bags outside \(D\) impose
boundary restrictions on the one- or two-state variables in \(D\).
The partial independent set extends to a maximum independent set if and
only if this residual port CSP is satisfiable.

This separates

\[
s_d=e_d+b_d,
\]

where \(e_d\) counts extendable erasures of maximum assignments and
\(b_d\) counts blocked partial assignments.  The extendable part is an
erasure-shadow profile of the binary code of maximum independent sets.
The blocked part is the genuinely tree-specific correction.

For \(d=1\), a matched bag is blocked precisely when selected boundary
neighbors hit both of its ports (a singleton is blocked when its sole
port is hit).  This gives a structural explanation for the observed
depth-1 LC failures: they are created by two-sided boundary saturation,
not by the erasure shadow itself.  For \(d\le4\), every obstruction is
an unsatisfiable port CSP on a forest of at most four bags.  There are
only finitely many obstruction shapes.

### Unit-propagation lemma

The obstruction shapes are simpler than a general binary CSP suggests.
Each skeleton edge forbids only one ordered pair of endpoint values.  If
one endpoint still has both values available, it always supplies a
compatible value to either value at the other endpoint.  A constraint
can therefore propagate only from a vertex already forced to a single
value, and it can delete at most one value at its neighbor.

On a tree, ordinary leaf elimination (equivalently, arc consistency) is
complete: after processing a leaf, retain at its neighbor exactly the
values having a compatible leaf value.  Induction on the number of
vertices proves that the residual CSP is satisfiable iff no domain is
emptied.  Consequently every minimal blocked residual component is an
implication path whose endpoint boundary forces propagate incompatible
values to the same vertex or across the same final edge.  A length-zero
path is the familiar depth-1 event in which boundary neighbors hit both
ports of one empty bag.

For \(|D|\le4\), the only minimal blocked cores are thus forced paths on
one, two, three, or four empty bags, with finitely many endpoint-port and
edge-port parities.  Branching residual trees introduce no new minimal
core: prune every branch not lying on the two colliding implication
chains.  This converts step 1 below from a generic CSP enumeration into
a four-row path-signature table.

### Alternating-port characterization

The path-signature table collapses further.  At a matched bag, label its
two endpoint states 0 and 1.  A selected boundary neighbor incident at
port \(a\) forbids state \(a\), hence forces state \(1-a\).  If an
internal edge uses ports \((a,b)\), a force propagates across it only
when the current forced state is \(a\); it then forces state \(1-b\) at
the next bag.

It follows that a force traverses an empty matched bag exactly when the
two path incidences use **opposite ports** at that bag.  A second selected
boundary neighbor creates a contradiction exactly when its incidence is
opposite to the incoming path incidence.  Therefore every minimal
blocked core containing only matched bags is precisely:

\[
\text{selected boundary port}
\;--\;\underbrace{\text{empty bags using opposite ports}}_{m\text{
bags}}\;--\;\text{selected boundary port}.
\tag{4}
\]

Both boundary occupied bags must be selected at the displayed incident
ports.  There is no additional parity label: local endpoint swaps gauge
all such transmitting paths to the same signature.

A singleton bag has only its retained state.  If it is empty, any force
arriving at its sole port is already a contradiction.  Thus the second
kind of minimal core is an alternating-port path with one selected
boundary source and a singleton terminal.  For defect depth at most
four, the complete obstruction list is consequently:

- two-boundary alternating paths through 1, 2, 3, or 4 empty matched
  bags;
- one-boundary alternating paths ending at a singleton, with total empty
  length 1, 2, 3, or 4.

This is an exact classification, not a finite-search conjecture.  The
remaining difficulty is overlap bookkeeping: one partial assignment may
contain several such cores, so \(b_d\) is the size of their union rather
than the raw number of paths.

## Poset gauge (stronger reduction)

The bipartition of the original tree gauges away all port parities.  After
complementing maximum independent sets to minimum vertex covers, orient a
matched-edge variable by which bipartition endpoint lies in the cover.  Every
nonmatching edge then becomes an inequality \(x_j\le x_i\).  The constraint
graph is a forest poset, while unmatched vertices impose consistent unary
values.  Consequently the code of maximum assignments is an order-ideal code
of a forest poset, times constant coordinates.

The full proof and its graph translation are in
`matching_bag_poset_reduction_2026-08-20.md`.  In particular, the extendable
erasure polynomial is a binomial factor times the independence polynomial of
the canonical bipartite graph \(B(P)\) with edge \(i^1j^0\) whenever
\(i\le_Pj\).  Alternating-port obstruction paths are simply violated
transitive comparisons in this gauge.  This supersedes the need for a
port-parity signature table.

## Updated proof attack

The extendable half is now closed by the universal Pascal-smoothing lemma
in `pascal_smoothing_shadow_lemma_2026-08-20.md`. For every rank, its
upper-tail inequalities through defect depth eight hold; for the present
rank at most 19, the
whole extendable profile is log-concave. At depth three it has the explicit
reserve

\[
e_3^2\ge \frac{32(\alpha-2)}{27(\alpha-3)}e_2e_4.
\]

The remaining work is therefore confined to the blocked term:

1. Express \(b_2,b_3,b_4\) as counts of the classified local signatures in the
   matching skeleton.
2. Use (3) to isolate the exact linear/quadratic combination of signature
   counts that could be negative.
3. Charge that combination against the displayed extendable reserve.
4. Prove that any remaining negative combination requires either at least 20 bags or
   at least six singleton bags.

Either conclusion closes the entire current rung: it has at most 19
bags and at most five singletons.  Unlike full enumeration, the attack
scales with defect depth, not with the number of original trees.

## Verification

The companion script constructs the bags and port constraints from an
independently computed maximum matching, evaluates (1) by a three-state
tree DP, and checks it coefficient-for-coefficient against the ordinary
tree independence-polynomial DP on every nonisomorphic tree through
order 12.  It also verifies (2) on every instance.
