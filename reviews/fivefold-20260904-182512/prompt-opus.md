Review the mathematical correctness of the finite-forest fivefold theorem in Section 12 of the attached artifact, including its supporting arguments in Sections 6 and 9--11. Audience: a mathematician deciding whether the displayed proof establishes its stated theorem. Independently check the argument; report what is wrong with it and what changes are needed. Do not use other reviews or the authors' conversation. Do not infer validity from numerical scans. An exact counterexample is welcome if a claim is false.

Keep this a focused proof audit, not literature research or a manuscript-style review. Return at most 700 words: a one-sentence statement of the claim, specific findings ranked by severity with exact hypotheses/steps, and a verdict VALID AS WRITTEN / VALID WITH EXPLICIT REPAIR / GAP / FALSE / UNDETERMINED. Explain the decisive reasoning. If citing outside results, distinguish those independently checked from those merely recalled. No novelty judgment is requested.

Frozen artifact follows verbatim:

# Joint-density analytic sprint (2026-09-04)

## Status

This note follows the adverse multicell lift experiment in
`arxiv_astra_transfer_2026-09-04.md`. Sections 1--11 retain the chronological
proof-search record. Section 12 now assembles a complete written proof that
every forest with \(\alpha\ge5\) and \(b_2=0\) satisfies
\(5b_3\le e_3\) and \(5b_4\le e_4\). The proof has been checked here but
has not received independent review or Lean verification. It closes the
\(b_2=0\) subcase of the live joint-density target; the full conditional
joint bridge, the depth-three rung, and Erdős Problem #993 remain open.

## 1. Normalized target

Write

\[
C=2e_3b_3-e_2b_4-e_4b_2,
\qquad
D=(b_3^2-b_2b_4)+C.
\]

The proved Pascal reserve reduces the live branch \(D<0\) to

\[
e_2b_4+e_4b_2+b_2b_4
\le c_\alpha e_2e_4,
\qquad
c_\alpha=\frac{5\alpha+17}{27(\alpha-3)}. \tag{1}
\]

Equivalently, with \(x_d=b_d/e_d\),

\[
x_2+x_4+x_2x_4\le c_\alpha. \tag{2}
\]

On the live \(\alpha\in\{17,18,19\}\) window, \(c_\alpha\) decreases to
\(c_{19}=7/27\).

## 2. Exact \(b_2=0\) residual lemma

Let \(S\) be an independent set which is not contained in a maximum
independent set.

**Lemma 1.** If \(b_2=0\), every blocked defect-three set is maximal.

**Proof.** If a blocked set \(S\) of size \(\alpha-3\) were not maximal, then
\(S\cup\{v\}\) would be an independent set of size \(\alpha-2\). The hypothesis
\(b_2=0\) puts that set inside a maximum independent set, which would also
contain \(S\), a contradiction. ∎

**Lemma 2.** If \(b_2=0\) and \(S\) is blocked at defect four, then

\[
\alpha(T-N[S])\le1.
\]

Consequently, because the residual is a forest, it is one of
\(\varnothing\), \(K_1\), or \(K_2\).

**Proof.** Two independent residual vertices could be added to \(S\), giving
an independent defect-two set. That set would extend to a maximum independent
set, again contradicting that \(S\) is blocked. A forest with independence
number at most one has at most two vertices, and two vertices must be
adjacent. ∎

Thus, under \(b_2=0\), \(b_3\) counts maximal defect-three sets and every set
counted by \(b_4\) is at most one vertex short of maximality. This is a
tree-specific structural restriction; it is substantially stronger than the
bare nonpure-complex description of the blocked faces.

## 3. A proof-sized sublemma

In the \(b_2=0\) branch, (1) becomes \(b_4/e_4\le c_\alpha\). The single
stronger statement

\[
5b_4\le e_4 \tag{3}
\]

would close that branch throughout the live window, since
\(1/5<7/27\le c_\alpha\). The companion inequality

\[
5b_3\le e_3 \tag{4}
\]

is useful for induction under pendant matched-pair deletion.

The exact audit in `../scripts/audit_joint_density_subcases_20260904.py`
enumerates the 149,239 nonisomorphic trees through \(n\le18\) having
\(\alpha\ge4\) and matching deficiency at most five. Its relevant counts are:

- 586 trees have \(\alpha\ge5\) and \(b_2=0\);
- 234 of these have \(b_3>0\), and 422 have \(b_4>0\);
- there are zero failures of (3) or (4);
- equality in (3) occurs exactly once in this corpus, at the five-leaf star
  `EiP?`, where \((e_2,e_3,e_4)=(10,10,5)\) and
  \((b_2,b_3,b_4)=(0,0,1)\);
- 387 of the 419 trees with \(D<0\) have \(b_2=0\); the remaining 32 have
  \(b_2>0\);
- there is no \(D<0\) case with \(\alpha\ge5\) and deficiency zero or one.

The run took 116.45 seconds on one core. Exact certificate:
`../results/joint_density_subcases_n18_20260904.json`.

These are bounded facts, not proofs of (3), (4), or the low-deficiency
nonadversity statement.

## 4. The obvious local exchange injection fails

For a blocked defect-four set \(S\) and a maximum set \(M\), put
\(A=S\setminus M\) and \(B=M\setminus S\). Since
\(|B|=|A|+4\) and \(A\ne\varnothing\), replacing \(A\) by an equally large
subset of \(B\) produces at least five extendable defect-four sets. This
suggests a fivefold injection proving (3).

The suggestion is false in this local form, even if every maximum \(M\) is
allowed. The first exact obstruction is the 9-vertex tree
<code>Hs&#96;?GGC</code>, with

\[
(e_2,e_3,e_4)=(21,35,35),
\qquad
(b_2,b_3,b_4)=(0,2,6).
\]

Each of its six blocked defect-four sets has five candidates under the stated
repair rule, but the resulting bipartite demand problem has maximum flow 26,
short of the 30 distinct repair slots required. The global inequality still
holds: \(5b_4=30<35=e_4\). Therefore any proof of (3) must use repairs outside
that immediate exchange shadow, allow a weighted/global charge, or use a
recurrence.

Replay:

```bash
python3 scripts/verify_b2_zero_exchange_obstruction_20260904.py
```

Certificate: `../results/b2_zero_exchange_obstruction_20260904.json`.

## 5. What sequential Cohen--Macaulay theory does and does not give

There is an exact dictionary to the nonpure face-number literature. If
\(\Delta=I(T)\), then the pure \((\alpha-1)\)-skeleton is the subcomplex
generated by the maximum independent sets. Hence its faces of size
\(\alpha-d\) are counted by \(e_d\), the faces of \(\Delta\) of that size are
counted by \(s_d\), and their difference is \(b_d\). The
\(\widetilde h\)-triangle records the \(h\)-vectors of all pure skeleta, so it
organizes exactly the layers from which \(e_d\) and \(b_d\) are recovered.

This does not supply (3). Woodroofe proves that the independence complex of a
chordal graph is vertex decomposable, hence shellable and sequentially
Cohen--Macaulay. Adiprasito--Björner--Goodarzi numerically characterize the
possible \(h\)-triangles of sequentially Cohen--Macaulay complexes, while
emphasizing that pairwise Macaulay comparisons do not suffice in the nonpure
case. More decisively, the desired inequality is false even inside the
chordal graph class.

For \(m\ge1\), let \(G_m\) be the split graph with independent part
\(A=\{a_1,\ldots,a_6\}\), clique part \(X=\{x_1,\ldots,x_m\}\), and cross
edges \(x_i a_j\) exactly when \(j\ne1\). Its independence-complex facets are
\(A\) and the \(m\) edges \(\{a_1,x_i\}\). Thus \(\alpha=6\), \(b_2=0\),
\(e_4=\binom62=15\), and \(b_4=m\). Equation (3) fails for every \(m\ge4\),
although \(G_m\) is split, hence chordal, and its independence complex has all
the properties above.

Sources:

- Russ Woodroofe, [*Vertex decomposable graphs and obstructions to
  shellability*](https://arxiv.org/abs/0810.0311), especially Corollary 7.
- Karim Adiprasito, Anders Björner, and Afshin Goodarzi,
  [*Face numbers of sequentially Cohen--Macaulay complexes and Betti numbers
  of componentwise linear ideals*](https://arxiv.org/abs/1502.01183),
  especially the definitions of the pure skeleta and \(h\)-triangle and
  Proposition 2.

The transfer verdict is therefore precise: the \(h\)-triangle is useful
bookkeeping, but neither shellability nor sequential Cohen--Macaulayness is
the missing inequality. A successful proof must retain the original graph's
tree geometry.

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

## 7. Resulting frontier

The most concrete next proof task is no longer the full three-variable
density bound. It is the following rooted pair:

> For a rooted tree \((T,w)\) with \(b_1=b_2=h_2=0\), prove
> \(5h_3\le e_2+q_3\) and \(5h_4\le q_4\), or refute either statement exactly.

The residual lemma from Section 2 should be applied inside \(T-w\); the sharp
\(h_4\) example should be treated as the equality model. A proof would give a
pendant-pair induction for (3)--(4), while a counterexample would kill that
induction without affecting (1).

The \(b_2>0\) branch remains separate. It accounts for only 32 of 419 adverse
small-corpus trees but appears in actual-rung lifts, so it cannot be
discarded. There the correct target remains the joint density (2); neither
blocked-profile log-concavity nor the fivefold \(b_4\) bound alone is
sufficient.

No further blind scaling is recommended. The next useful work is a bounded
analytic attack on (8), with a small rooted census serving only as a kill
test for proposed injections or recurrences.

## 8. Reproduction and hashes

Commands:

```bash
python3 scripts/audit_joint_density_subcases_20260904.py --max-n 18 \
  --output results/joint_density_subcases_n18_20260904.json
python3 scripts/verify_b2_zero_exchange_obstruction_20260904.py
python3 scripts/audit_pendant_p2_boundary_20260904.py --max-n 14 \
  --output results/pendant_p2_boundary_n14_20260904.json
```

SHA-256:

- `scripts/audit_joint_density_subcases_20260904.py`:
  `915f206cb2130d3333aa29f281f6dd520b8ecd19dea371daf2e58419905220f2`
- `results/joint_density_subcases_n18_20260904.json`:
  `702ec9e2e1180cbd020ca9475038cf9163d6d3d0b26831e9b9e7d2c0f89aef52`
- `scripts/verify_b2_zero_exchange_obstruction_20260904.py`:
  `89782bc5d2b00907b5fc8e13ee7e9c553f1779d1554706ebe52646afdbe5d01f`
- `results/b2_zero_exchange_obstruction_20260904.json`:
  `faf0dcbec56b6a69136382c18215f1e2937ae766cdef919375b7cd43d384f06e`
- `scripts/audit_pendant_p2_boundary_20260904.py`:
  `87ec4cb58f609ce703c4777339862c20ef0f9a0edac7de1596698c2d2ec9a87f`
- `results/pendant_p2_boundary_n14_20260904.json`:
  `ebdbc10a55e35ea981e819676edc91a722ec707ff5d573a72ca655e8f10e9a15`

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
