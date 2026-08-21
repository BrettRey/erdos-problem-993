# The maximum layer is a forest-poset code

Date: 2026-08-20. Status updated 2026-08-21: the abstract poset/code layer and
the complete tree/maximum-matching-to-forest-poset bridge are independently
formalized. The blocked correction remains open. Companion exact checks:
`verify_erasure_shadow_poset_20260820.py`. Formalization:
`formalization/matching_bag_poset_aristotle/`.

The Lean development proves the projection/induced-ideal identity (3), the `B(P)`
identity (5), constant-coordinate multiplication, the antichain formulas (8)–(9),
coefficient formula (10), and the matching-cover pigeonhole lemma. It now also derives
König's theorem from Mathlib Hall, constructs the matching bags and forest comparison
poset for every finite forest, proves the solution/order-ideal and tree-code equivalences,
and closes equations (2), (6), (7), and (7a). It does not formalize the blocked
correction of Section 4.

## 1. From a maximum matching to a poset

Let $T$ be a tree with bipartition $L\sqcup R$, and fix a maximum
matching

\[
M=\{L_iR_i:1\le i\le q\}.
\]

There are \(\delta=n-2q\) unmatched vertices and
\(\alpha=q+\delta\). A minimum vertex cover has size \(q\). Since its
vertices must already meet the \(q\) disjoint edges of \(M\), it contains
exactly one endpoint of every matching edge and no unmatched vertex.
Write

\[
x_i=1\quad\Longleftrightarrow\quad L_i\text{ belongs to the cover}
\]

(and hence $x_i=0$ means that $R_i$ belongs to it). A nonmatching
edge $L_iR_j$ imposes

\[
x_i=1\text{ or }x_j=0,
\qquad\text{equivalently}\qquad x_j\le x_i. \tag{1}
\]

An edge from an unmatched $L$-vertex to $R_j$ forces $x_j=0$;
an edge from $L_i$ to an unmatched $R$-vertex forces $x_i=1$.

Contracting the matching edges gives a tree. Consequently the
undirected graph underlying the inequalities (1) is a forest, so the
inequalities define a poset. Propagate the unary values and remove all
forced variables. The remaining poset $P$ still has forest cover
graph. Indeed, if a directed comparison path between two free variables
passed through a forced variable $w$, then $w=0$ would force every
variable below it to 0, while $w=1$ would force every variable above it
to 1. At least one endpoint would not be free.

Thus the maximum independent sets of $T$, after complementing to
minimum covers and harmlessly flipping coordinates, form

\[
\Omega(T)=\Omega(P)\times\{\text{one constant word on }c\text{
coordinates}\}, \tag{2}
\]

where $P$ is a forest poset and $c\ge\delta$ counts the singleton bags
and forced matching bags. This is substantially narrower than an
arbitrary binary code.

## 2. Projection counts are induced order-ideal counts

For a binary code \(\Omega\subseteq\{0,1\}^m\), put

\[
p_k(\Omega)=\sum_{|K|=k}|\pi_K(\Omega)|.
\]

For a code on $m$ coordinates, the erasure-shadow count in the
matching-bag note is $e_d=p_{m-d}$.

For the order-ideal code \(\Omega(P)\), a partial assignment on
$K\subseteq P$ extends precisely when it is order-preserving on the
induced poset $P[K]$. Every order ideal of an induced subposet extends:
take its down-closure in $P$. Therefore

\[
p_k(P)=\sum_{\substack{K\subseteq P\\|K|=k}}|J(P[K])|. \tag{3}
\]

There is a useful graph translation. Let $B(P)$ have two vertices
$i^1,i^0$ for each $i\in P$, with an edge

\[
i^1j^0\quad\Longleftrightarrow\quad i\le_P j. \tag{4}
\]

The reflexive edges $i^1i^0$ prevent assigning both values to one
coordinate; every other edge forbids an order violation. Independent
sets of $B(P)$ are consequently the extendable partial assignments, and

\[
\boxed{\quad \sum_k p_k(P)t^k=I(B(P);t).\quad} \tag{5}
\]

Adding a constant coordinate multiplies this polynomial by $1+t$.
Combining (2) and (5), the extendable erasure polynomial for an arbitrary
tree is

\[
(1+t)^c I(B(P);t), \tag{6}
\]

read in reverse degree. The graph $B(P)$ is the standard bipartite
graph attached to a poset; its independence complex is pure and
Cohen--Macaulay. Here its poset has the additional strong restriction
that the cover graph is a forest.

## 3. The erasure-shadow lemma is now proved universally

The clean sufficient statement was:

> **Forest-poset shadow lemma.** If $P$ has forest cover graph, then
> the coefficient sequence of $I(B(P);t)$ is log-concave in the needed
> upper layers.

This is now a theorem in considerably greater generality. The elementary
Pascal-smoothing argument in
`pascal_smoothing_shadow_lemma_2026-08-20.md` applies to every
downward-closed family. It proves every upper-tail inequality through
defect depth eight
at every rank, and the entire transform through rank 33. Put $r=|P|$ and
$M=r+c=\alpha(T)$. For the present depth-3 frontier only the single
inequality

\[
E_3^2\ge E_2E_4,\qquad
E_d=[t^{M-d}](1+t)^cI(B(P);t), \tag{7}
\]

is needed, and only for $M\le19$. In fact the argument gives the reserve

\[
E_3^2\ge \frac{32(M-2)}{27(M-3)}E_2E_4. \tag{7a}
\]

There is also an exact $h$-vector formulation, but the relevant
$h$-vector is **not** a $P$-Eulerian descent polynomial. Let

\[
A_P(u)=\sum_j a_j u^j
\]

be the antichain polynomial of $P$. Group a partial isotone assignment
by its set $O$ of assigned 1s. Its assigned 0s may be any subset of
$P\setminus\uparrow O$, so

\[
I(B(P);t)=\sum_{O\subseteq P}t^{|O|}(1+t)^{r-|\uparrow O|}. \tag{8}
\]

Put $u=t/(1+t)$ and group $O$ by its minimal antichain $A$. For a
fixed $A$, write $O=A\sqcup S$ with
$S\subseteq\uparrow A\setminus A$. The binomial sum over $S$ is
exactly $u^{|A|}$. Hence

\[
\boxed{\quad I(B(P);t)=(1+t)^r
 A_P\!\left(\frac{t}{1+t}\right).\quad} \tag{9}
\]

Equivalently,

\[
E_d=\sum_j a_j\binom{M-j}{d}. \tag{10}
\]

Formula (9) is also the standard whiskering/corona identity: it is the
independence polynomial of the graph obtained from the comparability
graph of $P$ by appending one leaf at every vertex. A claim of full
log-concavity at arbitrary rank would lie inside the long-studied
log-concavity problem for fully whiskered (very well-covered) graphs.
The bounded-rank and upper-tail statement used here is much weaker and
follows directly from normalized face-density monotonicity; it does not
settle that open problem.

## 4. Relation to the alternating-port obstruction

The bipartition gauge removes the port parity from the earlier matching-
bag picture. An alternating-port path is exactly a directed comparison
path in $P$. A blocked partial assignment is one that respects the
cover relations whose internal bags remain occupied but violates a
transitive comparison whose internal bags were erased. Thus

\[
s_d=e_d+b_d
\]

compares the forest constraint graph with its transitive closure $B(P)$.
At depth at most four, $b_d$ is the union of order inversions whose
unique comparison paths have at most four erased internal bags. This is
the same finite obstruction list as before, now with all port labels
gauged away.

The maximum-matching condition has an equally concise interpretation:
there is no alternating path between two singleton bags, because such a
path would be $M$-augmenting. In the poset gauge this is exactly the
consistency of the unary boundary values.

## 5. A necessary caution: ultra-log-concavity is false

The tempting stronger conjecture that every binary erasure-shadow profile
is ultra-log-concave is false. For the seven-bit code

\[
\Omega=\{24,47,65,67,86,93,97,99\},
\]

where integers encode binary words, the retained-coordinate profile is

\[
(p_0,\ldots,p_7)=(1,14,80,187,220,145,52,8).
\]

At $k=6$, normalized log-concavity would require

\[
52^2\binom75\binom77
\ge145\cdot8\binom76^2,
\]

but the two sides are $56{,}784$ and $56{,}840$. Ordinary
log-concavity still holds. Exact exhaustive checks find ordinary and
normalized log-concavity for every nonempty binary code through length
four, so the seven-bit example is a useful kill of the stronger route,
not evidence against the needed ordinary inequality.

## 6. Next attack

The extendable term is finished. The immediate targets are:

1. express the explicitly classified inversion-path correction $b_d$ in
   the poset gauge;
2. compare it with the quantitative reserve (7a);
3. use $M\le19$ and $\delta\le5$ to exclude a negative corrected margin.

This route exposes the algebraic object behind the matching-bag model and
separates a proved universal smoothing lemma from the genuinely
tree-specific blocked-path correction.
