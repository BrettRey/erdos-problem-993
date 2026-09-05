# Formalization request: fivefold bounds for blocked sets in finite forests

## Required result and exact scope

Produce a standalone Lean 4 / Mathlib proof of the following statement:
for every finite forest G with independence number alpha >= 5, if every independent
set of size alpha-2 is contained in a maximum independent set, then
5*b_3 <= e_3 and 5*b_4 <= e_4. Here e_d counts independent sets of size alpha-d
contained in a maximum independent set and b_d counts those not so contained.
Disconnected forests and isolated vertices are included. There is NO bound on
order, maximum degree, or matching deficiency.

This is an auxiliary result in the Erdős #993 project, not the full problem.
A complete written proof is supplied below; independently audit it and return
an exact refutation or missing inference if it fails.

Frozen source: notes/joint_density_analytic_sprint_2026-09-04.md
SHA-256: 6b554778d602baa14bcb1db59650d940efb9b0689988286ca3caae566566ff19
Repository HEAD: cfd88a5d704599644f930c7ad34c9169433968d0 (source uncommitted).

## Lean definitions: locally compiled against Mathlib v4.28.0

The following definitions were actually compiled with Lean 4.28.0 before
submission. Equivalent representations are allowed only with proved bridges
to these meanings. The guard d <= G.indepNum is essential: negative-size
layers must be zero, not the empty-set layer produced by truncated Nat subtraction.

```lean
import Mathlib.Combinatorics.SimpleGraph.Clique
import Mathlib.Combinatorics.SimpleGraph.Acyclic

namespace FivefoldForest

variable {V : Type*} [Fintype V]

/-- Containment in a maximum independent set, not merely a maximal one. -/
def Extendable (G : SimpleGraph V) (S : Finset V) : Prop :=
  ∃ M : Finset V, G.IsMaximumIndepSet M ∧ S ⊆ M

noncomputable def extendableCount (G : SimpleGraph V) (d : ℕ) : ℕ := by
  classical
  exact if d ≤ G.indepNum then
    (Finset.univ.filter fun S : Finset V =>
      S.card = G.indepNum - d ∧ G.IsIndepSet (S : Set V) ∧ Extendable G S).card
  else 0

noncomputable def blockedCount (G : SimpleGraph V) (d : ℕ) : ℕ := by
  classical
  exact if d ≤ G.indepNum then
    (Finset.univ.filter fun S : Finset V =>
      S.card = G.indepNum - d ∧ G.IsIndepSet (S : Set V) ∧ ¬ Extendable G S).card
  else 0

/-- Specification only; this definition is not a proof of the assertion. -/
def FivefoldConclusion (G : SimpleGraph V) : Prop :=
  5 * blockedCount G 3 ≤ extendableCount G 3 ∧
  5 * blockedCount G 4 ≤ extendableCount G 4

#check @FivefoldForest.extendableCount
#check @FivefoldForest.blockedCount
#check @FivefoldForest.FivefoldConclusion

end FivefoldForest
```

Required top-level theorem shape (NOT a theorem already proved by the definition):

```text
theorem fivefold_of_forest
    {V : Type*} [Fintype V] (G : SimpleGraph V)
    (hforest : G.IsAcyclic)
    (halpha : 5 <= G.indepNum)
    (hb2 : blockedCount G 2 = 0) :
    FivefoldConclusion G
```

A finite labelled-vertex representation is acceptable if a general finite-type
transport theorem is proved. Fixed-size computations are not the target.

## Deliverable and honest grading

G0 (required for COMPLETE): the actual universal forest theorem above.
G1: expose the supporting lemmas (b2=0 forces b1=0 for bipartite alpha>=4;
component factorization; star absorption; pendant-pair and pendant-star
recurrences; essential-root cone identity; small exceptional star cases)
as individually traceable declarations.
G2: a README mapping the written proof steps to Lean declarations, and
#print axioms for the final theorem and the named supporting lemmas.

No sorry/admit, added axioms, opaque unsupported hypotheses, implemented_by,
native_decide, bounded enumeration substituting for the theorem, or hidden
assumption that all smaller trees or all recurrence profiles already satisfy
the desired conclusion. Proving only numerical profile propagation from
assumed graph recurrences is PARTIAL, not G0. If a support/decomposition theorem
or a graph-to-profile bridge remains open, list it explicitly.

Use the proof route below or a sound simpler one. Strong induction must be
on vertex count; the conditional component lemma is an independently proved
algebraic result, not an induction assumption on a same-size forest.
A complete-looking conditional theorem is not a substitute for G0.
An honest PARTIAL with real proved steps is preferable to silent weakening.

Return a buildable pinned Lean project. Include the main theorem in the default
build targets. The local verifier will rebuild it, inspect its statement, search
for escape hatches, and inspect its axiom dependencies before accepting any verdict.

## Supporting pendant-pair definitions and recurrence

## 6. Exact pendant-\(P_2\) recurrence

Attach a new path \(w-p-l\) to a vertex \(w\in T\), with new inner vertex
\(p\) and leaf \(l\), and call the result \(T^+\). Then

\[
I_{T^+}(x)=(1+x)I_T(x)+xI_{T-w}(x),
\qquad
\alpha(T^+)=\alpha(T)+1. \tag{5}
\]

Let \(E_T\) be the face enumerator of the subcomplex generated by the maximum
independent sets of \(T\). Let \(Q\) be the corresponding enumerator generated
only by maximum sets avoiding \(w\), and put

\[
H=I_{T-w}-Q.
\]

The maximum facets of \(T^+\) are \(\{l\}\cup M\) for every maximum set
\(M\) of \(T\), together with \(\{p\}\cup M_0\) for maximum sets \(M_0\) of
\(T\) which avoid \(w\). Therefore

\[
E_{T^+}=(1+x)E_T+xQ,
\qquad
B_{T^+}=(1+x)B_T+xH. \tag{6}
\]

If \(q_d=[x^{\alpha-d}]Q\) and \(h_d=[x^{\alpha-d}]H\), the reversed profiles
satisfy

\[
e_d^+=e_{d-1}+e_d+q_d,
\qquad
b_d^+=b_{d-1}+b_d+h_d, \tag{7}
\]

with negative subscripts interpreted as zero.

If \(b_2^+=0\), nonnegativity in (7) forces \(b_1=b_2=h_2=0\). Assuming
(3)--(4) for the base tree, they propagate to \(T^+\) from only the two
root-conditioned boundary inequalities

\[
5h_3\le e_2+q_3,
\qquad
5h_4\le q_4. \tag{8}
\]

Indeed,

\[
5b_3^+=5(b_3+h_3)\le e_3+e_2+q_3=e_3^+
\]

and

\[
5b_4^+=5(b_3+b_4+h_4)\le e_3+e_4+q_4=e_4^+.
\]

The exact rooted audit in
`../scripts/audit_pendant_p2_boundary_20260904.py` checks all
nonisomorphic base trees through order 14 and every root. Among 5,446 trees,
196 satisfy \(b_1=b_2=0\); their 2,469 roots yield 838 eligible extensions
with \(h_2=0\). There are zero failures of either inequality in (8) or either
propagated fivefold inequality. The \(h_4\) inequality is sharp: its largest
observed ratio is one, on rooted graph6 `FqPA?` at root 1. The \(h_3\)
ratio is at most \(5/26\) in this corpus. Runtime was 8.23 seconds.

Exact certificate: `../results/pendant_p2_boundary_n14_20260904.json`.

Contracting a fixed maximum matching turns \(T\) into a tree of matched and
singleton bags. If the matching deficiency is zero or one, that quotient has
at least two leaf bags and at most one singleton bag, so at least one leaf bag
is matched. In the original tree it is exactly a pendant \(P_2\). Thus
(6)--(8) are especially well aligned with the previously uncovered
deficiency-zero and deficiency-one cells. This does not yet prove that
\(D\ge0\) in those cells; it supplies a recursive mechanism for the
\(b_2=0\) subcase.

## Written proof and supporting elementary lemmas

## 9. Bounded forest-product pass (2026-09-04)

Brett requested short reasoning passes with regular checkpoints and low file
churn; substantial local compute is acceptable. This pass stops at the
following elementary lemmas. It does not close the rooted pair or #993.

**Root deletion identifies the forest profiles.** For \(\alpha\ge3\),
\(h_2=0\) implies that a maximum independent set of \(T\) avoids \(w\).
Otherwise \(Q=0\), whereas \(T-w\) has independence number at least
\(\alpha-1\), so it has an independent set of size \(\alpha-2\) and
\(h_2>0\). Consequently \(\alpha(T-w)=\alpha(T)\), and \(Q,H\) are exactly
the extendable and blocked enumerators of the forest \(F=T-w\).

**Lemma 3 (bipartite graphs).** If \(G\) is bipartite, \(\alpha(G)\ge4\),
and \(b_2(G)=0\), then \(b_1(G)=0\).

**Proof.** Suppose a blocked independent set \(S\) has size \(\alpha-1\).
It is maximal: adding one vertex would give a maximum set containing it.
For each \(v\in S\), the set \(S-\{v\}\) extends to a maximum set by
\(b_2=0\). This maximum set excludes \(v\), since otherwise it would
contain \(S\), and therefore supplies two vertices outside \(S\).
Each supplied vertex has a neighbour in \(S\) by maximality, and its only
such neighbour is \(v\). The supplied pairs for different \(v\) are thus
disjoint. Hence \(|V(G)|\ge3(\alpha-1)\). Bipartiteness gives
\(|V(G)|\le2\alpha\), forcing \(\alpha\le3\), a contradiction. ∎

In particular, the eligible rooted forests in Section 7 have
\(b_1(F)=b_2(F)=0\). The order restriction matters: the six-vertex path
has \(\alpha=3\), \(b_1>0\), and \(b_2=0\).

**Lemma 4 (conditional product closure).** Let a forest have components
\(F_i\), each with \(b_{i,1}=b_{i,2}=0\). If every component satisfies
\(5b_{i,3}\le e_{i,3}\) and \(5b_{i,4}\le e_{i,4}\), then the forest
satisfies both fivefold bounds.

**Proof.** Use reversed enumerators
\(\mathcal E_i(z)=\sum_d e_{i,d}z^d\) and
\(\mathcal B_i(z)=\sum_d b_{i,d}z^d\). Extendability is componentwise, so

\[
\mathcal E_F=\prod_i\mathcal E_i,\qquad
\mathcal B_F=\prod_i(\mathcal E_i+\mathcal B_i)-\prod_i\mathcal E_i.
\]

Each \(\mathcal B_i\) starts at degree at least three. Through degree four,
only terms with one blocked factor contribute. At degree three their
component-defect patterns are \((3,0,\ldots)\); at degree four they are
\((4,0,\ldots)\) and \((3,1,0,\ldots)\). Applying the component bounds
replaces these terms by distinct terms of \(\mathcal E_F\). The remaining
terms of \(\mathcal E_F\) are nonnegative. ∎

The component hypotheses \(b_{i,1}=b_{i,2}=0\) follow from
\(b_1(F)=b_2(F)=0\): the degree-one product coefficient forces all
\(b_{i,1}=0\), and then degree two forces all \(b_{i,2}=0\), since every
\(e_{i,0}>0\).

**Small-component exception.** Among connected tree components with
\(\alpha_i\le4\) and \(b_{i,1}=b_{i,2}=0\), the only possible nonzero
blocked profile is that of \(K_{1,4}\). For \(\alpha_i\le3\), the remaining
layers are empty sets or negative sizes and have no blocked faces. For
\(\alpha_i=4\), a blocked singleton cannot belong to an independent pair,
since all pairs extend by \(b_{i,2}=0\); its vertex is universal, forcing
the tree to be \(K_{1,4}\). This star has
\(e=(1,4,6,4,1)\), \(b=(0,0,0,1,0)\), and fails \(5b_3\le e_3\).
Thus product closure alone does not finish the rooted argument. The connected
component bounds for \(\alpha_i\ge5\), compensation for the star exception,
and a decomposition supporting a noncircular induction remain to be supplied.

**Exact checks in this pass.** A direct independent-set scan of all 201 trees
of orders 1 through 10 found the small star exception and the order-six
\(b_1>0,b_2=0\) exception at \(\alpha=3\). A second scan reused
`independent_masks`, `profile`, and `rooted_profile` from
`scripts/audit_pendant_p2_boundary_20260904.py`, over all 5,446 trees of
orders 2 through 14 and every eligible root. It recovered 838 eligible roots,
zero with \(h_1>0\), zero components with \(b_1>0\) or \(b_2>0\), and 60
occurrences of small blocked components, all four-leaf stars. Wall time was
8.45 seconds. These were read-only console checks; the existing certificates
were not overwritten. The lemmas above rest on their displayed proofs.

## 10. Absorbing the four-leaf star (2026-09-04)

This short follow-up resolves the small-component exception left in Section 9.
The connected-tree fivefold bounds remain assumptions, not proved results.

**Lemma 5 (star absorption).** Let \(G\) be a nonempty forest with
\(b_1=b_2=0\), \(5b_3\le e_3\), and \(5b_4\le e_4\).
Then \(G\sqcup K_{1,4}\) satisfies the same conditions.

**Proof.** Put \(m_d=e_d-5b_d\). The reversed extendable and blocked
enumerators of the star are \((1+z)^4\) and \(z^3\). Multiplication of
the enumerators, with \(b_0=b_1=b_2=0\), gives exactly

\[
m_3'=m_3+4e_2+6e_1-e_0,
\qquad
m_4'=m_4+4m_3+6e_2-e_1+e_0. \tag{9}
\]

Write \(a=\alpha(G)\ge1\) and \(n=|V(G)|\le2a\). Counting incidences
between maximum sets and their one-vertex deletions gives
\(ae_0\le(n-a+1)e_1\), hence \(e_0\le2e_1\).
If \(a\ge2\), the same count between extendable sets of sizes
\(a-1\) and \(a-2\) gives
\((a-1)e_1\le(n-a+2)e_2\), hence \(e_1\le4e_2\).
These inequalities make both right-hand sides of (9) nonnegative.
If \(a=1\), then \(e_2=0\), \(e_1=1\), and \(e_0\ge1\), which
also suffice. The zero blocked coefficients at defects one and two are
preserved by multiplication. ∎

The nonempty hypothesis is essential: the lone star fails the defect-three
bound. Two stars provide the base when there is no other component:

\[
(e_3,e_4)=(56,70),\quad (b_3,b_4)=(2,8),\quad (m_3,m_4)=(46,30).
\]

The defect-four estimate can be sharp: \(K_{1,4}\sqcup K_1\) has
\(e_4=5\) and \(b_4=1\).

**Corollary (conditional reduction to connected components).** Let \(F\)
be a forest with \(\alpha(F)\ge5\) and \(b_2(F)=0\). If every component
with independence number at least five satisfies both fivefold bounds, then
\(F\) satisfies both bounds.

Indeed, Lemma 3 gives \(b_1(F)=0\), and Section 9 passes both zero
conditions to its components. Every smaller component satisfies the bounds
unless it is \(K_{1,4}\). If any non-star components exist, combine them
using Lemma 4, then add the stars using Lemma 5. Otherwise there are at least
two stars; start with the displayed two-star base and apply Lemma 5.

For the rooted problem at \(\alpha\ge5\), this proves
\(5h_3\le q_3\) and \(5h_4\le q_4\), **conditional on the corresponding
connected-component bounds for \(T-w\)**. The first is stronger than the
required \(5h_3\le e_2+q_3\). This removes forest multiplication as a
separate obstruction; it does not supply the missing connected-tree theorem
or a complete induction. In particular, a proof restricted by matching
deficiency would need a scope check: root deletion here raises total
deficiency by one, so the original deficiency-five census does not itself
cover every possible residual component.

**Exact checks.** Read-only Python checks reused `independent_masks` and
`profile` from `scripts/audit_pendant_p2_boundary_20260904.py`. All 987
nonisomorphic trees of orders 1 through 12 passed the two incidence bounds
used in (9). Seventy bases satisfied the absorption hypotheses; adding twenty
successive stars to each gave 1,400 exact truncated-polynomial checks of (9)
and both fivefold bounds, with zero failures. Pure star unions were checked
for 2 through 100 stars. The combined run took 0.40 seconds. Separately,
direct independent-set enumeration of each base and its star union replayed
(9) for all 12 bases with \(b_1=b_2=0\) through order eight, including
the exceptional lone-star base. All identities and applicable inequalities
passed. These console checks support the displayed proof; no existing
certificate was overwritten.

Checkpoint: star compensation is settled. Stop this pass here; the next
question is a connected-tree reduction that can support induction.

## 11. Connected-tree reduction (2026-09-04; induction assembly pending)

The next bounded pass supplies the missing structural case split. A complete
strong-induction proof, with every small exception checked together, has not
yet been assembled; this section does not promote the fivefold conjecture
to a theorem.

Every non-star tree has a vertex \(p\) adjacent to \(m\ge1\) leaves and
exactly one other vertex \(w\): root at a non-leaf vertex and choose a
deepest non-leaf vertex. Removing \(p\) and its leaf neighbours leaves a
nonempty connected tree \(G\). The case \(m=1\) is the existing pendant-pair
recurrence. For \(m\ge2\), every maximum independent set excludes \(p\),
since replacing \(p\) by its leaves strictly increases size. Writing
\(a=\alpha(G)\), the resulting tree \(T\) satisfies

\[
\alpha(T)=a+m,\quad
E_T(x)=(1+x)^m E_G(x),\quad
B_T(x)=(1+x)^m B_G(x)+xI_{G-w}(x). \tag{10}
\]

In particular \(b_2(T)=0\) forces \(b_1(G)=b_2(G)=0\).
For \(m=2\), the last term contributes an independent set of size
\(a-1\) in \(G-w\) to \(b_2(T)\); such a set always exists because
\(\alpha(G-w)\ge a-1\). Thus **the two-leaf case is impossible** under
\(b_2(T)=0\).

For \(m\ge4\), \(T\) is obtained from \(G\sqcup K_{1,m}\) by adding
the edge \(wp\). Every maximum set of the disjoint union excludes \(p\),
so this edge leaves the maximum sets, and therefore all extendable counts,
unchanged. It can only remove blocked sets. Consequently any fivefold
bounds for the disjoint union pass to \(T\). Section 10 handles the
four-leaf star; for \(m\ge5\), the isolated star itself satisfies both
bounds (only \(m=5\) has a nonzero blocked coefficient at defects three
or four, with \(b_4=1\) and \(e_4=5\)).

**Three leaves.** When \(m=3\), equation (10) and \(b_2(T)=0\) force
\(\alpha(G-w)=a-1\), so every maximum set of \(G\) contains \(w\).
Put \(F=G-w\), and let \(f_d\) be its extendable profile relative to
\(\alpha(F)=a-1\). Because \(b_1(G)=0\), every maximum set of \(F\)
extends by \(w\); because \(b_2(G)=0\), every independent set of size
\(a-2\) in \(F\) also extends to a maximum set of \(G\), and hence to
one of \(F\). Thus \(b_1(F)=0\) and

\[
E_G(x)=(1+x)E_F(x),\quad
e_{G,0}=f_0,\quad e_{G,1}=f_0+f_1,\quad e_{G,2}=f_1+f_2.
\]

With \(M_d=e_{G,d}-5b_{G,d}\), equation (10) gives

\[
e_{T,3}-5b_{T,3}=M_3+3f_2+6f_1-f_0,
\]
\[
e_{T,4}-5b_{T,4}=M_4+3M_3+3f_2-f_1+f_0. \tag{11}
\]

Both added expressions are nonnegative for nonempty \(F\). Section 10's
incidence bound gives \(f_0\le2f_1\). For \(r=\alpha(F)\ge3\), its
second incidence bound gives
\(f_1\le(r+2)f_2/(r-1)\le(5/2)f_2\), so
\(3f_2-f_1+f_0\ge0\). For \(r=2\), \(f_2=1\), \(f_1\le|V(F)|\le4\),
and \(f_0\ge1\) suffice. For \(r=1\), \(f_2=0\), \(f_1=1\), and
\(f_0\ge1\) suffice. In the target range \(\alpha(T)\ge5\),
\(r=a-1\ge1\). Therefore the three-leaf operation preserves both bounds
whenever \(G\) satisfies them.

The exceptional base \(G=K_{1,4}\) also survives this operation: its
\((M_3,M_4)=(-1,1)\), the eligible root is a leaf, and
\((f_0,f_1,f_2)=(1,3,3)\). Formula (11) yields margins \((25,5)\).
This is exactly the nine-vertex local-injection obstruction of Section 4,
whose \((e_3,e_4)=(35,35)\) and \((b_3,b_4)=(2,6)\).

**Exact check.** A read-only independent-set scan of all 5,446 trees of
orders 2 through 14 found 193 with \(\alpha\ge5\) and \(b_2=0\), including
47 without a pendant pair. Every tree had a removable support as described
above. There were no two-leaf supports in these trees. All 28 occurrences
of the three-leaf case passed the essential-root condition, the cone
identities for \(E_G\), and both exact margin identities (11).
The first three-leaf example was the nine-vertex obstruction. Runtime was
7.73 seconds, using `independent_masks` and `profile` from
`scripts/audit_pendant_p2_boundary_20260904.py`; no certificates were
overwritten.

Checkpoint: the structural obstruction to a pendant-pair-only induction is
resolved by (10)--(11). The next short pass should assemble the strong
induction, auditing the \(\alpha=4\) exceptions in the pendant-pair case
and the interaction of those bases with the forest reduction. No larger
census is needed for that audit.

## 12. Fivefold forest theorem: assembled proof (2026-09-04)

**Theorem.** Let \(F\) be a finite forest with independence number
\(\alpha\ge5\). Write \(e_d\) for the number of independent sets of size
\(\alpha-d\) contained in a maximum independent set, and \(b_d\) for the
number not so contained. If \(b_2=0\), then

\[
5b_3\le e_3,\qquad 5b_4\le e_4. \tag{12}
\]

No matching-deficiency restriction is imposed. Coefficients at negative set
sizes are zero. The proof uses the elementary lemmas and identities displayed
in Sections 9--11, rather than the bounded computational evidence.

**Proof for trees by strong induction on vertex count.** Let \(T\) satisfy
the hypotheses, and assume the assertion for every smaller tree. Lemma 3
also gives \(b_1(T)=0\).

For a star \(K_{1,m}\), the hypotheses give \(m\ge5\). The only blocked
set is the centre singleton. When \(m=5\), \(b_3=0\), \(b_4=1\), and
\(e_4=5\); when \(m\ge6\), both blocked counts vanish. Thus stars satisfy
(12). For a non-star, use Section 11's decomposition into a smaller tree
\(G\), an attached support \(p\), its \(m\) leaf neighbours, and the
attachment vertex \(w\in G\).

**One leaf (\(m=1\)).** Put \(a=\alpha(G)=\alpha(T)-1\ge4\).
The pendant-pair identity gives
\(b_1(G)=b_2(G)=h_2=0\), with \(H=I_{G-w}-Q\) as in Section 6.
Section 9 then identifies \(Q=E_{G-w}\), \(H=B_{G-w}\), and
\(\alpha(G-w)=a\).

First suppose \(G\ne K_{1,4}\). If \(a\ge5\), the induction hypothesis
gives (12) for \(G\). If \(a=4\), Section 9's small-component
classification gives \(b_3(G)=b_4(G)=0\), so the same bounds hold.

For the residual forest \(R=G-w\), Lemma 3 and \(h_2=0\) give
\(b_1(R)=b_2(R)=0\), and these zero conditions pass to each component.
If \(a\ge5\), each component with independence number at least five is
a smaller tree satisfying the induction hypothesis. The *conditional*
forest corollary of Section 10 therefore proves (12) for \(R\), yielding
\(5h_3\le q_3\) and \(5h_4\le q_4\).
If \(a=4\) and \(R\ne K_{1,4}\), every component satisfies the bounds by
the small-component classification, so Lemma 4 gives the same conclusion.
The only remaining residual possibility is \(R=K_{1,4}\), when
\((h_3,h_4,q_3,q_4)=(1,0,4,1)\). Here the weaker required inequality
\(5h_3\le e_{G,2}+q_3\) still holds, since a maximum set of \(G\)
alone supplies \(\binom42=6\) extendable pairs; the other rooted bound
is immediate. In all these cases Section 6's propagation proves (12) for
\(T\).

It remains to handle the exceptional base \(G=K_{1,4}\). A leaf root is
impossible: it belongs to every maximum set, so \(Q=0\) and \(h_2>0\).
Thus \(w\) is the centre. Directly,

\[
E_T(x)=(1+2x)(1+x)^4,\qquad B_T(x)=x+x^2.
\]

Consequently \(\alpha(T)=5\),
\((e_3,e_4)=(14,6)\), and \((b_3,b_4)=(1,1)\), proving (12).

**At least two leaves (\(m\ge2\)).** Equation (10) forces
\(b_1(G)=b_2(G)=0\). As shown in Section 11, \(m=2\) contradicts
\(b_2(T)=0\).

For \(m=3\), put \(a=\alpha(G)=\alpha(T)-3\ge2\).
If \(G\ne K_{1,4}\), it satisfies (12) either by induction
(\(a\ge5\)) or by the small-component classification (\(a\le4\)).
The two margin identities (11) and the nonnegative added expressions proved
there give (12) for \(T\). If \(G=K_{1,4}\), Section 11's explicit
calculation gives margins \((25,5)\), so this exception also passes.

For \(m\ge4\), form \(U=G\sqcup K_{1,m}\). It has
\(\alpha(U)=\alpha(T)\ge5\) and \(b_1(U)=b_2(U)=0\).
Each of its two components is strictly smaller than \(T\). Each component
with independence number at least five therefore satisfies (12) by
induction. Apply Section 10's conditional forest corollary to obtain (12)
for \(U\). Adding the edge \(wp\) preserves all maximum sets and all
extendable counts, while only removing independent sets. Hence
\(e_d(T)=e_d(U)\) and \(b_d(T)\le b_d(U)\), proving (12) for \(T\).

This exhausts the support cases and completes the tree induction.

**Passage to forests.** Given any forest in the theorem, Lemma 3 gives
\(b_1=0\), and the zero conditions pass to its components. The tree theorem
just established applies to each component with independence number at least
five. Section 10's conditional corollary now proves (12) for the forest. ∎

**Dependency audit.** The tree induction only invokes the theorem for
strictly smaller *connected* trees. Lemmas 3--5, the small-component
classification, and the conditional forest corollary were proved separately
without the tree theorem. The comparison forest \(U\) in the \(m\ge4\)
case has the same order as \(T\), but its components are smaller; that case
uses the independently established component corollary, not an induction
hypothesis for \(U\). Root deletion can raise deficiency, which causes no
scope problem because (12) has no deficiency hypothesis.

**Consequences and limits.** The rooted pair (8) is now proved for its
stated range \(\alpha\ge4\): for \(\alpha\ge5\), apply the forest theorem
to \(T-w\); for \(\alpha=4\), use the same small-forest argument and lone
star exception as in the one-leaf case above. In the live
\(\alpha\in\{17,18,19\}\) window with \(b_2=0\),

\[
\frac{b_4}{e_4}\le\frac15<\frac7{27}\le c_\alpha,
\]

which closes this subcase of the joint endpoint bound (1), and therefore its
depth-three inequality using the previously established Pascal reserve.
The \(b_2>0\) subcase remains open. This is not a proof or refutation of
Erdős Problem #993, and no novelty claim is made here.

**Small-base replay.** A direct independent-set enumeration checked all
three trees with \(\alpha=4\) and \(b_1=b_2=0\) (all such trees have at
most eight vertices). All nine eligible rooted pendant-pair extensions
satisfied both bounds. Exactly one used the exceptional star base, giving
the displayed margins \((9,1)\); none had residual \(K_{1,4}\), though
the written proof covers that possibility without relying on its absence.
The exceptional star base was also replayed with one through five attached
leaves at both root types (centre and leaf). Every eligible case passed;
the three-leaf leaf-root case gave margins \((25,5)\). These read-only
checks reused `independent_masks` and `profile` from
`scripts/audit_pendant_p2_boundary_20260904.py` and completed in under one
second. The proof is complete as a written argument, but remains without
independent review or Lean verification. Stop this pass at that checkpoint.
