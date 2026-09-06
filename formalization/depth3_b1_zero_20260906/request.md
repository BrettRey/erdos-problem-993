# Frozen Lean request: the b1-zero depth-three window theorem

Prove `b1_zero_window : DepthThree.B1ZeroWindowTarget` using the unchanged
specification RequestProject/DepthThreeSpec.lean. The exact universal target is:
every finite TREE G with 33 <= |V| <= 38, 17 <= alpha(G) <= 19,
2alpha(G) <= |V|+5, and b1=0 satisfies s2*s4 < s3^2.
The counts e_d and b_d are independent sets of size alpha-d that respectively
are and are not contained in a MAXIMUM independent set. Count sets once, not
once per completion. Negative-size counts are zero via the supplied guard.
This is an auxiliary result, NOT the full depth-three window or Erdős #993.

The complete computer-assisted written proof is WRITTEN_PROOF.md Sections 1--7.
Supplied Python/JSON files document the finite calculations, not trusted proof
axioms. A simpler sound proof is permitted.

## Reuse the completed low-density project: preferred shorter route

This is a continuation of project be8de245-5246-41b8-aa5e-ca01f73d8008, whose
previous task returned both exact graph targets. The main agent downloaded and
fresh-built all source locally (8,049 jobs), checked the unchanged specification,
read all six new modules, and independently printed the axioms of the targets
and key bridges: only propext, Classical.choice, Quot.sound.
The API status COMPLETE_WITH_ERRORS disagreed with its COMPLETE report, but
the actual local source replay passed. Keep all six completed modules.

In particular, `DepthThree.low_density_of_forest : DepthThree.LowDensityTarget`
is now proved, and the actual graph/matching-code count bridge is available.
Use this theorem to shorten the b1-zero proof: prove only the graph density
`3*(alpha-3)*b4 <= (alpha-7)*e4` in the b1-zero window, then apply it.

- r=0: all blocked sets vanish.
- 1<=r<=3: the new kernel-checked `all_small_forbidden_density_rays` proves
  every cone ray lies below `(alpha-7)/(3*(alpha-3))`. Thus the structural
  reduction, the p_j deletion inequality, and the cone bridge suffice.
- r=4: alpha=19, and the proved encoded bound is
  `52513/217404 < 1/4`, exactly below the low-density theorem's threshold.
  Prove the graph domination/encoding/evaluator bridges, then use low density.

This avoids proving b2 density bounds, the joint-endpoint scalar assembly, and
the separate special b3>=11b2 bound. Those are included as optional alternate
routes, not extra required targets. The frozen B1ZeroWindowTarget is unchanged.

## New locally kernel-checked supporting modules

The packet now includes five supporting modules, all checked locally with only
`propext`, `Classical.choice`, and `Quot.sound`:

- DepthThreeAlgebra: low-density, small-forbidden, and four-forbidden scalar assembly.
- ParameterCertificate: the five-ray monotone cone lemma, the exact joint
  rational bound, and the low-density bound on every ray in all twelve r<=3
  parameter cases (`all_small_forbidden_density_rays`).
- RootedForestCertificate: a complete recursive enumeration of PLANE rooted
  forests, a coverage theorem for every encoded forest of size10, a top-five
  reversed-coefficient recurrence, and the universal inequality
  `217404*numerator f <= 52513*denominator f` for every such encoded forest.
  The exhaustive Boolean proof uses `decide +kernel`, not native_decide.
- FivefoldBridge: proves `e_eq` and `b_eq` identifying DepthThree.e/b with the
  supplied finite-set eD/bD helpers, including maximum-set identification.
  This is not the matching-poset/projected-code bijection.
- PolynomialCertificate: proves coefficientwise attachment concentration via
  a subtraction-free semiring identity, and proves that the top-five evaluator
  agrees with the full reversed polynomial recurrence. The graph interpretation
  of that recurrence is still a separate obligation.

These are not the graph theorem. In particular, for RootedForestCertificate,
you must prove every rooted graph forest has a corresponding Forest encoding
and that its profile recurrence has the claimed graph-count meaning.
The numerator is 58p0+21p1+3p2+30q0+9q1+q2, and the denominator is
126p0+84p1+36p2+9p3+p4. They are the b4/e4 coefficients from equation(11)
for P=corona polynomial and Q=component-roots-deleted polynomial at base size10.
The encoded completeness theorem is already proved, so do not redo a canonical
1,842-type classification unless it simplifies the remaining graph bridge.

Four previously verified forest-count helpers are also supplied under
RequestProject/Fivefold (namespace FivefoldForest): Defs, Basic, Union, Paths.
They provide induced-set counts, positivity, component products, and acyclic
path lemmas. Only their import prefix was changed. The assembled input source
passed a local `lake build` (8,049 jobs). Reuse these helpers where useful;
identify their internal counts with DepthThree.e/b via the supplied FivefoldBridge
proofs. The matching-code count identification is still a separate obligation.

## Required proof obligations

1. For finite forests with b1=0, prove the structural reduction: allowed induced
graph H is well-covered, forced C are isolated in H, H-C has a pendant perfect
matching, |C|=|U|+delta, every forbidden vertex has >=3 forced neighbours, and
if |U|>0 then |U|<=delta-1. The matching-bag splicing proof is in Section 1.

2. Establish E=I(H) and the b4 bound of Section 3. For K=H-C of size2m,
its top coefficients obey (m-j)p_j <=2(j+1)p_(j+1) by matching-bag deletion
counts. Normalize by beta_j=choose(m,j)/2^j and use the nondecreasing cone.
The finite parameter calculation has 13 cases; all rays in the 12 with r<=3
satisfy the needed low-density threshold by the supplied theorem. Include
denominator positivity and r=0. The older joint bound is an alternate route.

3. For the last parameters n=33,alpha=19,delta=5,r=4, the forced/forbidden core
is a tree with 9 forced and4 forbidden vertices of forced degree3; K has ten
pendant pairs. Every K component attaches exactly once to U. Replace the core
by the central forced star shape using |N_C(S)|>=2|S|+1; then concentrate all
attachments at one forbidden vertex. The exact coefficientwise identity
new-old = ax(A_i-Q_i)(A_j-Q_j) >=0 proves concentration, with a=(1+x)^2.
K is a corona of a forest on ten base vertices; move each component root from
leaf to its mate if needed by the explicit size-preserving injection.

4. Prove the remaining coefficient inequality for EVERY rooted forest on ten
base vertices (one root per component). B*=aP((a+x)^3-a^3)+xQ(a+x)^3,
E=(1+x)^9P; b4/e4 <=52513/217404 <7/27. The Python generator lists all1842
unordered types and crosschecks with NetworkX rooted trees on11 vertices.
Lean must prove completeness AND the graph-to-representation/domination bridge,
not merely decide that all imported JSON rows pass. A fresh sound analytic bound
below7/27 is acceptable. Alternatively enumerate all ORDERED rooted forests of
size10 (equivalently plane rooted trees of size11, 16796 types): this larger
list may make a recursive completeness proof simpler than canonical isomorphism
classification. If using finite enumeration, use kernel-checkable `decide` or
another proof-producing method, not native_decide or an external oracle.

5. OPTIONAL alternate assembly, not needed with the preferred low-density route:
on the ORIGINAL tree prove b3>=11b2: b2=Σq_i0 and b3>=Σ(6q_i0+q_i1),
with q_i1>=5q_i0 by deletion counting in the ten original matching bags.
Also 16e3>=4e4. Do not transfer a lower bound from a dominating representative,
and do not assume coefficientwise domination preserves log-concavity.
Use the Pascal reserve to obtain full margin >=
(7/27-v)e2e4+(9/2-v)e4b2+b3^2>0, v=52513/217404.

6. Already proved in the previous task: GraphBridge and ExtendableCount supply
maximum matching existence, card Bag=indepNum, maximum-set identification,
and the extendable-set↔projected-word cardinality-preserving bijection.
`TreeMatching.e_eq_erasure` transports the reserve to the exact graph counts.
Reuse those declarations; do not replace them by assumptions or redo the bridge.
All source is pinned to Lean/Mathlib4.28.0.

## Graded deliverable and integrity

COMPLETE requires the exact B1ZeroWindowTarget graph theorem with all bridges.
Useful graded partials: structural theorem; graph-count/code bridge; r<=3
closure; universal rooted-forest bound with proved completeness; r=4 closure;
final assembly. Name every completed declaration and every remaining obligation.
Refutation or an exact missing inference is useful; an honest PARTIAL with real
proved steps is preferable to a conditional theorem presented as COMPLETE.

No sorry/admit, new axioms, opaque unsupported hypotheses, implemented_by,
native_decide, external oracle, or fixed-profile arithmetic substituting for the
universal graph theorem. Any internal representation needs a proved bridge.
Include final modules in the default Lake build, a statement/declaration map,
completion grading, and #print axioms results. Local verification will rebuild,
inspect exact target statements and check axioms. Leave supplied dependency
sources intact when practical; add new modules.

Source note SHA-256: 4d77890b084b1a7aec056217b84d2606f17c78b4ec5772aff2401426d300a8bb.
Repository HEAD: c1abc1ac7b4ebfedfc28d482ec7440e2da5a3bec (new source uncommitted).
