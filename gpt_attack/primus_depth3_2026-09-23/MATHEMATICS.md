# Exact mathematical brief

## 1. Objects and the primary target

All graphs here are finite, simple, and undirected. A tree is connected and
acyclic; a forest need not be connected. An independent set contains no adjacent
vertices. A **maximum** independent set has the greatest possible size alpha;
do not confuse maximum with inclusion-maximal.

Write

\[
I(T;x)=\sum_{k=0}^{\alpha}i_k(T)x^k,
\qquad s_d=i_{\alpha-d}(T).
\]

The target, denoted W, is

\[
\boxed{s_2s_4\le s_3^2}
\tag{W}
\]

for every tree satisfying

\[
33\le n\le38,\qquad 17\le\alpha\le19,\qquad 2\alpha\le n+5.
\]

For trees, delta=2 alpha-n is the matching deficiency, hence nonnegative. The
displayed conditions leave twelve (alpha,delta) cells: alpha=17 with delta=0,1;
alpha=18 with delta=0,1,2,3; alpha=19 with delta=0,1,2,3,4,5.

This is one local log-concavity inequality near the top of a polynomial in a
bounded order window. It neither establishes all coefficient positions nor
arbitrary tree orders. A failure of W is not automatically a failure of
unimodality: a positive sequence can be unimodal without being log-concave.

## 2. Intrinsic counts; no multiplicities of completions

An independent set S is **extendable** if S is a subset of at least one maximum
independent set of T. Otherwise it is **blocked**. Let e_d and b_d count these
two types at size alpha-d. Set the counts to zero when d>alpha. Then

\[
s_d=e_d+b_d,\qquad b_0=0.
\]

Each set is counted once, not once for every maximum completion. The definitions
in `lean/RequestProject/DepthThreeSpec.lean` are the formal source of truth.
These counts do not depend on a chosen matching or root.

Define the signed correction over the integers:

\[
D=(b_3^2-b_2b_4)+(2e_3b_3-e_2b_4-e_4b_2).
\]

Expanding s=e+b gives the exact identity

\[
s_3^2-s_2s_4=(e_3^2-e_2e_4)+D.\tag{1}
\]

In Lean, cast to integers before subtracting; natural-number subtraction would
erase the negative correction that this task is about.

## 3. Available results

The supplied source and local project records establish the following. The
historical verification date is 6 September 2026; packaging checks are listed
separately in `CHECKS.md`.

### A. Extendable Pascal reserve

For the forest counts above with alpha>=4,

\[
27(\alpha-3)e_3^2\ge32(\alpha-2)e_2e_4,
\]

and therefore

\[
e_3^2-e_2e_4\ge c_\alpha e_2e_4,
\qquad c_\alpha=\frac{5\alpha+17}{27(\alpha-3)}>0.\tag{2}
\]

Also e_d>0 for 0<=d<=alpha. The representation/count bridge is proved; this is
not an assumption that a graph happens to have a suitable code.

Relevant source: `RequestProject/PascalBridge.lean`,
`RequestProject/TreeCodeBridge.lean`, `RequestProject/ExtendableCount.lean`, and
the use of `MatchingBag.TreeMatching.erasure_depth_three_reserve` in
`RequestProject/LowDensity.lean`.

### B. Blocked shadow and extendable incidence

For every finite forest with alpha>=4,

\[
6b_3\ge(\alpha-4)b_2,\qquad
(\alpha-3)e_3\ge4e_4.\tag{3}
\]

The graph-level blocked-shadow declaration is
`DepthThree.blocked_shadow_of_forest : DepthThree.BlockedShadowTarget`.
The extendable incidence declaration is
`MatchingBag.TreeMatching.extendable_incidence`.

### C. Low-density closure

For every finite forest with alpha in {17,18,19},

\[
3(\alpha-3)b_4\le(\alpha-7)e_4
\quad\Longrightarrow\quad s_2s_4<s_3^2.\tag{4}
\]

The b_4/e_4 thresholds at alpha=17,18,19 are respectively 5/21, 11/45,
and 1/4. No b_1=0 condition is needed.

Declaration: `DepthThree.low_density_of_forest : DepthThree.LowDensityTarget`.

### D. The b_1=0 window

Every tree in W's window with b_1=0 satisfies the strict inequality
s_2s_4<s_3^2.

Declaration: `DepthThree.b1_zero_window : DepthThree.B1ZeroWindowTarget`.
The proof establishes the actual graph theorem, including a structural
reduction and a finite enumeration coverage argument. It does not merely
assume a list of coefficient rows is complete.

The formal proof uses a sufficient one-quarter density bound in its exceptional
case. Do not attribute the sharper historical written bound 52513/217404 or
the special estimate b_3>=11b_2 to this Lean theorem; they are not needed here.

## 4. What is left

If D>=0, equations (1) and (2) settle W. Results C and D settle the low-density
and b_1=0 branches. Thus it suffices to handle trees in the same window with

\[
\boxed{b_1>0,\quad
3(\alpha-3)b_4>(\alpha-7)e_4,\quad D<0.}\tag{R}
\]

Do not assume R is empty. A historical bounded search that failed to find an
example in R is not an emptiness proof. Proving R impossible would be a valid
solution, but requires a structural or coverage argument.

One sufficient target is

\[
(5\alpha+17)e_2e_4+27(\alpha-3)D\ge0.\tag{S}
\]

It implies W by (1)-(2), but is potentially stronger than W. An exact failure
of S with positive full margin refutes only this proposed route. A still
stronger endpoint estimate discards the favorable terms 2e_3b_3+b_3^2:

\[
(5\alpha+17)e_2e_4\ge
27(\alpha-3)(e_2b_4+e_4b_2+b_2b_4).\tag{E}
\]

E is not known universally and historical certificates contain failures outside
the adverse regime. Do not silently replace S by E or require E everywhere.

`lean/PrimusSpec.lean` defines W, the residual target under R, and S as open
propositions. It contains no proof of any of them. W's non-strict conclusion
is intentional; the established subcases happen to give strict positivity.

## 5. Computational representation supplied for checking

Choose a maximum matching of a forest. Each matched edge forms a two-vertex
bag; each unmatched vertex forms a singleton bag. There are alpha bags. A
maximum independent set selects exactly one vertex per bag, subject to the
original adjacency constraints. Contracting the matched edges gives a forest
of bags.

An extendable partial independent set corresponds to a partial bag assignment
that admits a full assignment. The supplied `extendable_defects` routine uses
dynamic programming over the bag forest, retaining the **set of feasible root
choices**, not a count of completions. Empty bags contribute one unit of
defect. Grouping by feasible-choice sets avoids counting the same partial set
several times. The tests compare this routine against direct independent-set
enumeration on small examples.

The independence polynomial is computed separately by the ordinary rooted-tree
recurrence: excluding a root multiplies the total child polynomials; including
it contributes x times the product of the child-excluded polynomials. The
blocked counts are s_d-e_d. All decisive replay calculations use arbitrary-
precision integers. These numerical tools are evidence and counterexample
checkers, not a substitute for a universal proof.

## 6. What would count as an advance

A complete W proof, an exact W counterexample, or a new proved lemma that
genuinely constrains R would count. A useful partial result could be an
explicit exchange/counting bound with proved multiplicities and constants,
or a smaller exhaustive model with a proved coverage map.

For any new representation, prove that it counts these intrinsic objects
without multiplicity and preserves the hypotheses and target. For any finite
certificate, separate exact arithmetic from completeness of its enumeration.
For any external theorem, state and discharge its actual hypotheses.
