# The b1-zero forbidden-core reduction (5 September 2026)

## Result and limits

There is now a written proof, with an exhaustive exact finite certificate,
of the depth-three inequality
\(s_3^2>s_2s_4\) for every tree in the current window
\(33\le n\le38\), \(\alpha\in\{17,18,19\}\),
\(\delta=2\alpha-n\le5\), for which every independent
\((\alpha-1)\)-set extends to a maximum independent set (\(b_1=0\)).

Sections 1--4 close the cases with at most three forbidden vertices.
Section 7 closes the last case, \((n,\alpha,\delta,r)=(33,19,5,4)\),
by coefficientwise domination and a complete 1,842-case calculation.
Thus the whole \(b_1=0\) part of this window is covered.
The separate \(b_1>0\) branch remains open; this does not close the full
depth-three window or Erdős #993.

Status: written proof plus exact computational checks, not independently
reviewed or Lean-verified. The Pascal reserve used in the final assembly is
an existing formalized dependency. No novelty claim or manuscript change.

## 1. Definitions and structural lemma

For a finite forest \(F\), let \(\alpha=\alpha(F)\), \(n=|V(F)|\), and
\(\delta=2\alpha-n\). An independent set is *extendable* if it is contained
in a maximum independent set; otherwise it is *blocked*. At defect \(d\),
write \(e_d,b_d,s_d=e_d+b_d\) for their counts at size \(\alpha-d\).
Let \(A\) be the union of the maximum independent sets, \(C\) their
intersection, \(U=V(F)\setminus A\), and \(r=|U|\).
Thus \(A,C,U\) are the allowed, forced, and forbidden vertices, respectively.

Assume \(b_1=0\). Then:

1. \(|C|=r+\delta\).
2. The induced forest \(H=F[A]\) is well-covered: every independent set
   extends to one of size \(\alpha\). Its forced vertices are isolated,
   and its other vertices admit a perfect matching of pendant edges.
3. Every forbidden vertex has at least three forced neighbours.
4. If \(r>0\), then \(r\le\delta-1\).

Here a pendant edge has an endpoint of degree one in \(H\).

### Proof

Fix a maximum matching \(M\). Since a forest is bipartite, König's equality
gives \(\alpha=n-|M|\). Partition the vertices into the matched pairs and
the unmatched singletons. There are \(\alpha\) bags, and every maximum
independent set chooses exactly one vertex from every bag.
The \(\delta\) singleton vertices are therefore forced. In a matched bag,
either both endpoints are allowed, or one endpoint is forbidden and its
mate is forced. This proves \(|C|=r+\delta\). A forced vertex has no
allowed neighbour, since any maximum set containing that neighbour also
has to contain the forced vertex. Thus \(C\) is isolated in \(H\).

Consider a matching edge \(xy\) with both endpoints allowed. Suppose it is
not pendant in \(H\). Choose allowed neighbours \(x'\ne y\) of \(x\) and
\(y'\ne x\) of \(y\). They belong to different components of
\(F-\{x,y\}\), since otherwise there is a cycle. Every such component is
a union of whole remaining matching bags. The restriction of any maximum
independent set to it therefore contains exactly its number of bags.

In the component containing \(x'\), use the restriction of a maximum set
containing \(x'\). Do the same with \(y'\), and choose arbitrary maximum-set
restrictions on all other components. Their union is independent and has
size \(\alpha-1\), but it blocks both endpoints of the missing bag \(xy\).
It cannot extend to a maximum set, contradicting \(b_1=0\).

Consequently every allowed matching pair has a leaf endpoint in \(H\).
To extend any independent set of \(H\), add all missing forced isolates and
a leaf endpoint from every unoccupied pair. These additions cannot
conflict. This produces a maximum independent set of \(F\), proving (2).

Fix \(v\in U\). Acyclicity permits at most one neighbour of \(v\) in each
component of \(H\). A nonforced neighbour can be avoided by a maximum set
of its component; a forced neighbour is an isolated vertex. Writing
\(d_C(v)=|N(v)\cap C|\), we can therefore choose an independent set of
\(H-N(v)\) of size \(\alpha-d_C(v)\). Adding \(v\) gives size
\(\alpha+1-d_C(v)\). Values \(d_C(v)=0,1,2\) respectively contradict
maximal size, the fact that \(v\) is forbidden, and \(b_1=0\).
Hence \(d_C(v)\ge3\).

If \(r>0\), the forest induced by \(U\cup C\) has at least \(3r\) edges
and \(r+|C|=2r+\delta\) vertices. Thus
\(3r\le2r+\delta-1\), which proves (4). ∎

In particular, \(b_1=0\) and \(\delta\le1\) force \(U=\varnothing\) and
all blocked counts to vanish. In the current window, \(r\le4\).

## 2. What the reduction adds to the existing unary formula

Every independent set avoiding \(U\) now extends, so the whole extendable
polynomial is exactly \(I(H;x)\). There are no remaining pair obstructions.
Order \(U=(v_1,\ldots,v_r)\). The existing first-forbidden-vertex formula
from [the 4 September note](arxiv_astra_transfer_2026-09-04.md) becomes

\[
B_F(x)=x\sum_{i=1}^r
I\bigl(F-N[v_i]-\{v_j:j<i\};x\bigr). \tag{1}
\]

The formula itself is not new here: the new information is that it accounts
for *all* blocked sets and contains at most four summands in this window.
Each residual forest has independence number at most \(\alpha-3\), since
otherwise adding \(v_i\) would give a blocked set of size at least
\(\alpha-1\). If \(r>0\) in the window, \(\delta\ge2\), hence \(n\le36\);
each residual has at most 32 vertices.

Individual summand log-concavity does not prove the required aggregate
inequality. The previously refuted cross-summand route remains closed.
Instead, we use the coefficient bounds below.

## 3. A coefficientwise density bound

Put \(c=|C|=r+\delta\). By Section 1,

\[
I(H;x)=(1+x)^c P(x),\qquad P(x)=I(K;x),
\]

where \(K=H-C\) has \(2m\) vertices, a perfect matching, and
\(m=\alpha-c\). Write \(p_j=[x^{m-j}]P(x)\), taking out-of-range
coefficients to be zero.

For a nonempty \(S\subseteq U\), let \(t=|S|\). The bipartite forest
with vertex set \(S\cup N_C(S)\) has at least \(3t\) edges. Consequently

\[
|N_C(S)|\ge2t+1. \tag{2}
\]

For sets whose intersection with \(U\) is exactly \(S\), ignore internal
edges of \(S\) and all restrictions that \(S\) imposes on \(K\). Retain
only the forced vertices excluded by (2). Their generating polynomial is
coefficientwise bounded by
\(x^t(1+x)^{c-2t-1}P(x)\).
The exponent is nonnegative because \(c\ge2r+1\). This polynomial has
degree \(\alpha-t-1\), so at defect \(d\) only \(t\le d-1\) contributes.
Summing over \(S\) gives

\[
b_2\le r p_0, \tag{3}
\]

\[
\begin{aligned}
b_4\le{}&r\left[p_2+(c-3)p_1+\binom{c-3}{2}p_0\right]\\
&+\binom r2\left[p_1+(c-5)p_0\right]+\binom r3p_0.
\end{aligned} \tag{4}
\]

Terms for \(t>r\) are absent; their binomial prefactors are zero.

We also have, for \(0\le j<m\),

\[
(m-j)p_j\le2(j+1)p_{j+1}. \tag{5}
\]

Indeed, count one-vertex deletions from independent \((m-j)\)-sets of
\(K\). An independent \((m-j-1)\)-set occupies that many matching bags,
leaving \(j+1\) empty pairs. There are at most \(2(j+1)\) possible
one-vertex extensions. No log-concavity assumption is used.

Let \(\beta_j=2^{-j}\binom mj\). In the nonempty-forbidden cases of the
current window, \(m\ge10\), and (5) says
\(t_j=p_j/\beta_j\) is nonnegative and nondecreasing for \(0\le j\le4\).
Equation (3) and
\(e_2=p_2+cp_1+\binom c2p_0\) imply

\[
\frac{b_2}{e_2}\le
u:=\frac{r}{\beta_2+c\beta_1+\binom c2}. \tag{6}
\]

Write the right-hand side of (4) as \(\sum_{j=0}^4 A_jp_j\), where

\[
\begin{aligned}
A_0&=r\binom{c-3}{2}+\binom r2(c-5)+\binom r3,\\
A_1&=r(c-3)+\binom r2,\qquad A_2=r,\qquad A_3=A_4=0.
\end{aligned}
\]

Since \(e_4=\sum_{j=0}^4\binom c{4-j}p_j\), we obtain

\[
\frac{b_4}{e_4}\le
v:=\max_{0\le k\le4}
\frac{\sum_{j=k}^4 A_j\beta_j}
     {\sum_{j=k}^4\binom c{4-j}\beta_j}. \tag{7}
\]

To justify (7), express the nondecreasing sequence \(t_j\) as a
nonnegative combination of step sequences \(1_{j\ge k}\), with weights
\(t_0,t_1-t_0,\ldots,t_4-t_3\). On each step sequence, numerator is at
most \(v\) times denominator. Summing preserves that inequality. All
five denominators are positive.

## 4. Exact finite certificate and depth-three consequence

The existing Pascal reserve is

\[
e_3^2-e_2e_4\ge c_\alpha e_2e_4,\qquad
c_\alpha=\frac{5\alpha+17}{27(\alpha-3)}\ge\frac7{27}.
\]

Its source and formalization are recorded under `extendable_pascal_reserve`
in [the proof frontier](../proof_graph/erdos993_frontier.json).
Writing \(x_d=b_d/e_d\), the sufficient joint bound is
\(x_2+x_4+x_2x_4<c_\alpha\). Equations (6)--(7) bound its left side by
\(u+v+uv\).

There are exactly 13 admissible integer triples \((\alpha,\delta,r)\)
with \(r>0\), \(r\le\delta-1\), and the stated window constraints.
The exact rational calculation closes all 12 with \(r\le3\).
The largest upper bound among them occurs at \((19,4,3)\):

\[
u=\frac2{53},\qquad v=\frac{2632}{13039},\qquad
u+v+uv=\frac{170838}{691067}<\frac7{27}.
\]

The last comparison has positive cross-multiplied slack
\(7\cdot691067-27\cdot170838=224843\).
All 13 parameter records, including their five cone-ray ratios and
nonnegative ray slacks, are in the
[exact certificate](../results/b1_zero_forbidden_core_20260905.json).
The case \(r=0\) has \(b_d=0\) and follows directly from the reserve.

Finally,

\[
\begin{aligned}
s_3^2-s_2s_4
={}&(e_3^2-e_2e_4)+2e_3b_3+b_3^2\\
&-e_2b_4-e_4b_2-b_2b_4>0.
\end{aligned}
\]

This proof does not require the combined correction \(D\) to be negative;
the joint bound holds throughout this particular structural subcase.
It must not be promoted to an unconditional bound on other trees.

## 5. The four-forbidden family: the initial residual case

The sole parameter triple not closed by Section 4 is
\((\alpha,\delta,r)=(19,5,4)\).
Here \(c=9\), \(m=10\), and the coarse bound gives

\[
u=\frac{16}{369},\qquad v=\frac{2264}{8793},\qquad
u+v+uv=\frac{1012328}{3244617}>\frac7{27}.
\]

This is failure of a relaxation, not an actual tree counterexample.
The structural edge count is tight: the graph on \(U\cup C\) is a
13-vertex bipartite tree, every forbidden vertex has exactly three forced
neighbours, and there are no edges within \(U\). The remaining 20 vertices
form a well-covered forest \(K\) with ten matching pairs.
For an original *tree*, every component of \(K\) attaches to this core by
exactly one edge, ending at a forbidden vertex. Connectivity supplies an
attachment; two attachments would create a cycle, and a forced vertex
cannot have an allowed neighbour.

The follow-up in Section 7 retains the attachment restrictions in an exact
rooted-component calculation. The coarse joint bound still fails on some
actual members of this family, so the final assembly also restores the
positive defect-three term that Section 4 discarded.

## 6. Checks and rejected shortcut

Run:

```bash
python3 scripts/audit_b1_zero_forbidden_core_20260905.py --max-n 14
python3 test_all.py
```

The first command completed in 13.624 seconds and recorded:

- 5,446 nonisomorphic trees of orders 2--14, including 594 with \(b_1=0\);
- 105 disjoint-union checks from pairs of tree components of orders 1--6,
  including 36 with \(b_1=0\);
- 109 distinct archived graph6 profiles replayed, including 82 with
  \(b_1=0\);
- zero structural or coefficient-bound failures;
- four explicit trees at \(\alpha=19\) attaining
  \(r=\delta-1\), for \(r=1,2,3,4\);
- the 13 exact rational parameter certificates described above.

Small-graph profiles were checked against enumeration of every independent
set. The new matching-bag dynamic program counts a partial set once by its
unique feasible endpoint domain, not once per maximum completion. It also
reproduced every archived profile. The repository test suite passed all
55 tests. These finite checks support the implementation; the universal
structural assertion rests on the written proof.

An initial shortcut, \(b_1>0\Rightarrow D\ge0\), failed on five archived
trees. The first has \((n,\alpha,\delta)=(33,19,5)\),
\((b_1,b_2,b_3,b_4)=(8,117,994,7141)\), and
\(D=-56227048\). Its exact graph6 and profiles are retained in the JSON.
This kills the shortcut, not the depth-three target or Erdős #993.

## 7. Closing the four-forbidden case

The assumptions throughout this section are
\((n,\alpha,\delta,r)=(33,19,5,4)\) and \(b_1=0\).
In particular, \(|C|=9\) and the flexible forest \(K\) has ten matching
pairs. We prove

\[
\frac{b_4}{e_4}\le\frac{52513}{217404}
<\frac7{27},\qquad b_3\ge11b_2. \tag{8}
\]

These are enough even though the stronger joint endpoint bound can fail.

### 7.1 A dominating core and concentration of attachments

Keep the flexible components and their attachment roots fixed. For each
\(S\subseteq U\), count the independent sets with forbidden intersection
exactly \(S\). For nonempty \(S\), the forced contribution is
\((1+x)^{9-|N_C(S)|}\), which is coefficientwise at most
\((1+x)^{8-2|S|}\) by (2). Equality for every nonempty \(S\) is attained by
the core with one forced vertex adjacent to all four forbidden vertices,
and two private forced leaves at each forbidden vertex. Replacing the
core by this shape therefore increases every blocked coefficient while
leaving the extendable polynomial \((1+x)^9I(K;x)\) unchanged.

Group the flexible components by their attaching forbidden vertex. Let
\(A_i\) be the independence polynomial of group \(i\), and \(Q_i\) that of
the group with its attachment roots deleted. An empty group has
\(A_i=Q_i=1\). Put \(a=(1+x)^2\) and \(P=\prod_i A_i=I(K;x)\).
For the dominating core the blocked polynomial is exactly

\[
B(x)=\prod_{i=1}^4(aA_i+xQ_i)-a^4P. \tag{9}
\]

Indeed, the central forced vertex is absent from every blocked set;
each forbidden vertex is either absent, leaving its two private leaves
free, or present, excluding those leaves and its flexible attachment roots.
Subtract the term with all four forbidden vertices absent.

Move all the components of group \(j\) to group \(i\). For the two
affected factors, new minus old is

\[
\begin{aligned}
&(aA_iA_j+xQ_iQ_j)(a+x)
 -(aA_i+xQ_i)(aA_j+xQ_j)\\
&\hspace{2em}=ax(A_i-Q_i)(A_j-Q_j)\succeq0. \tag{10}
\end{aligned}
\]

Here \(\succeq0\) means coefficientwise nonnegative. Each difference
\(A_i-Q_i\) counts independent sets containing at least one attachment
root. The other two factors are nonnegative, and \(a^4P\) is unchanged.
Repeated merging proves that concentrating all flexible components at one
forbidden vertex increases every blocked coefficient.

### 7.2 Reduction to rooted forests on ten vertices

Choose a leaf endpoint from each pendant matching edge of \(K\), and call
its mate a base vertex. All edges between matching pairs join base vertices.
Thus \(K\) is obtained from a forest \(G\) on ten vertices by adding one
private leaf at every vertex: the corona of \(G\).

If an attachment root is a private leaf \(\ell\) with mate \(v\), replacing
it by \(v\) increases the root-deleted independence polynomial:
\(I(K-v;x)\succeq I(K-\ell;x)\). To see this, take an independent set
avoiding \(\ell\). If it contains \(v\), replace \(v\) by \(\ell\);
otherwise leave it unchanged. This is an injective, cardinality-preserving
map into the independent sets avoiding \(v\). Apply this separately to
each component.

It is therefore enough to consider all forests \(G\) on ten vertices,
with one distinguished base vertex in each component, and attach every
component at its distinguished vertex to the same forbidden vertex.
There are exactly 1,842 rooted-forest types. Completeness can be checked
in two independent ways: enumerate unordered multisets of rooted trees,
or add a new root joined to the component roots and enumerate rooted trees
on eleven vertices. The verifier implements both and compares their
canonical-code sets, not merely their counts.

Write \(W\) for these distinguished vertices, \(P=I(K;x)\), and
\(Q=I(K-W;x)\). Equation (9) reduces to

\[
B^*(x)=aP\bigl((a+x)^3-a^3\bigr)+xQ(a+x)^3,
\qquad E(x)=(1+x)^9P. \tag{11}
\]

The calculation checks (8)'s first bound on all 1,842 types. Its maximum
is attained when \(G\) is the ten-vertex star rooted at its centre. The
corresponding actual 33-vertex tree has

\[
(e_2,e_3,e_4)=(50841,217392,652212),\qquad
(b_2,b_3,b_4)=(2051,26672,157539).
\]

Thus the maximum ratio is \(157539/652212=52513/217404\).
Every original tree is coefficientwise dominated at the blocked levels by
one of these representatives with the same extendable polynomial, so the
bound transfers. This use of domination is only for \(b_4/e_4\); we do
not assume that coefficientwise domination preserves log-concavity.

### 7.3 The positive defect-three term pays for b2

For the original tree, let \(K_i\) be \(K\) with all flexible neighbours of
the \(i\)-th forbidden vertex deleted, and put
\(q_{i,j}=[x^{10-j}]I(K_i;x)\). Each such vertex has at most one neighbour
in each component of \(K\). These neighbours can all be avoided by
maximum sets of their components, so \(\alpha(K_i)=10\).

Any independent 9-set of \(K_i\) occupies nine of the ten original
matching bags. It admits at most two one-vertex extensions. Double-counting
deletions from independent 10-sets gives \(q_{i,1}\ge5q_{i,0}\).

A defect-two blocked set contains exactly one forbidden vertex: two would
exclude at least five forced vertices by (2), giving defect at least three.
For a fixed forbidden vertex, its three forced neighbours are absent and
the other six forced vertices are free. Consequently

\[
b_2=\sum_iq_{i,0},\qquad
b_3\ge\sum_i(6q_{i,0}+q_{i,1})\ge11b_2. \tag{12}
\]

This is a bound on the original tree, not a lower bound borrowed from its
dominating representative.

Every extendable independent 15-set has at least four extendable
one-vertex extensions, from any maximum set containing it. Every
extendable 16-set has sixteen one-vertex deletions. Hence
\(16e_3\ge4e_4\), or \(e_3\ge e_4/4\).
Combining the Pascal reserve, (8), and (12), with
\(v=52513/217404\), gives

\[
\begin{aligned}
s_3^2-s_2s_4
&\ge \frac7{27}e_2e_4-e_2b_4
       +b_2(22e_3-e_4-b_4)+b_3^2\\
&\ge \left(\frac7{27}-v\right)e_2e_4
       +\left(\frac92-v\right)e_4b_2+b_3^2>0.
\end{aligned} \tag{13}
\]

The first coefficient is exactly \(3851/217404>0\). This closes the
four-forbidden case and, with Sections 1--4, the whole \(b_1=0\) part of
the current window.

### 7.4 Exact verification and scope

Run:

```bash
python3 scripts/audit_four_forbidden_concentration_20260905.py
```

The 3.253-second run checked all 1,842 representatives in three ways:
subset enumeration on the ten base vertices, independence-polynomial DP
on the actual 33-vertex tree, and the feasible-domain matching-bag count
of extendable sets. All agree. It also checked coefficientwise domination
and (12) on 2,692 small decorated instances and 1,000 deterministic
33-vertex decorations, with zero failures. The small instances use all
four core types, every component attachment owner, and both root ports
through three flexible pairs; they are not a census at order 33.

The stronger joint endpoint bound fails on 17 concentrated representatives.
Those are not failures of (8) or (13). In particular, the largest-ratio
representative above has positive combined correction
\(D=2637714832\), despite failing the discarded stronger bound.

The [finite certificate](../results/four_forbidden_concentration_20260905.json)
stores a canonical code and exact coefficient row for every representative,
plus the graph6 certificate of the extremizer. This is a computer-assisted
written proof, not an independent review, a Lean formalization, a novelty
claim, or a solution of Erdős #993. All 55 repository tests still pass.

## 8. The remaining b1-positive branch: two sign shortcuts are false

The archived 33-vertex witness in Section 6 has \(b_1=8\), all from pair
obstructions, with no unary defect-one obstruction. Moving its leaf 3
from vertex 2 to vertex 10 preserves \((n,\alpha,\delta)=(33,19,5)\) and
produces a tree whose blocked profile is entirely unary:

\[
(b_1,b_2,b_3,b_4)=(8,132,2739,22070),\qquad
D=-178212783.
\]

Its extendable profile at defects 2--4 is
\((138240,562176,1601856)\), its full depth-three margin is
\(94423068753>0\), and its independence polynomial is unimodal.
This kills the further shortcut that a unary defect-one obstruction
forces \(D\ge0\); it was found on the ninth admissible leaf reattachment.
Thus the remaining \(b_1>0\) argument must handle both unary and pair
obstructions, without inferring the sign of the combined correction
from their presence.

Both original and mutated graph6 certificates, complete polynomials,
and the exact unary/pair split are retained in
[the sign-obstruction certificate](../results/b1_positive_sign_obstructions_20260905.json).
Reproduce them with:

```bash
python3 scripts/replay_b1_positive_sign_obstructions_20260905.py
```

## 9. A simpler sufficient bound for the remaining branch

The existing forest-poset representation supplies one more useful lemma,
without assuming \(b_1=0\):

\[
6b_3\ge(\alpha-4)b_2. \tag{14}
\]

For completeness, every blocked independent set contains a blocked subset
of size at most two. In the maximum-assignment code, a nonextendable
partial assignment either selects a value forbidden on a forced
coordinate, or violates one poset comparison between two specified free
coordinates. If neither happens, the specified free values extend to an
order ideal. This is the projection characterization in Section 2 of
[the forest-poset note](matching_bag_poset_reduction_2026-08-20.md).

Take a blocked independent \((\alpha-2)\)-set and such a one- or two-vertex
certificate. Deleting any vertex outside the certificate preserves
blockedness, giving at least \(\alpha-4\) blocked one-vertex deletions.
Conversely, an independent \((\alpha-3)\)-set occupies all but three
matching bags, and admits at most six one-vertex extensions. Counting
these incidences proves (14). This does not assume that *every* deletion
of a blocked set is blocked; deleting part of its certificate can restore
extendability.

The extendable incidence count used in Section 7 also gives, generally,
\((\alpha-3)e_3\ge4e_4\). Put

\[
h_\alpha=\frac{\alpha-7}{3(\alpha-3)},\qquad v=\frac{b_4}{e_4}.
\]

Keeping the positive \(2e_3b_3\) term, (14) and the Pascal reserve give

\[
s_3^2-s_2s_4
\ge(c_\alpha-v)e_2e_4+(h_\alpha-v)e_4b_2+b_3^2. \tag{15}
\]

For \(\alpha\in\{17,18,19\}\),
\(c_\alpha-h_\alpha=4(20-\alpha)/(27(\alpha-3))>0\).
Thus \(v\le h_\alpha\) proves the desired strict inequality, regardless
of \(b_1\). The thresholds are \(5/21,11/45,1/4\), respectively.

Combining this with the completed \(b_1=0\) argument, the only remaining
regime that needs work in the current window has all three properties:

\[
b_1>0,\qquad D<0,\qquad \frac{b_4}{e_4}>h_\alpha. \tag{16}
\]

Both incidence bounds passed on all 5,438 nonisomorphic trees with
\(\alpha\ge4\) through order 14, as well as the eligible disconnected
and archived checks in the updated forbidden-core verifier. These tests
support the implementation; (14)--(15) follow from the written argument.

A bounded adaptive search then targeted (16), starting from the five
archived actual-window \(b_1>0,D<0\) seeds and the new unary witness.
Of 5,000 mutations, 3,369 stayed in the window, with 3,085 distinct labeled
graph6 profiles and 254 visits to \(b_1>0,D<0\). There were no high-density
adverse cases, negative full depth-three margins, or non-unimodal trees.
The largest adverse density found was

\[
\frac{b_4}{e_4}=\frac{16283}{706415}\approx0.02305,
\]

at \((n,\alpha,\delta)=(33,19,5)\), with \(D=-394650\). This is only
\(65132/706415\), about 9.22%, of the proved \(1/4\) threshold, but the
search is adaptive and bounded, not a census or proof that (16) is empty.
Its full depth-three margin is \(15532731494>0\). The run took 3.815
seconds; a fresh graph6 replay reproduced the best record exactly.

```bash
python3 scripts/probe_high_b4_adverse_20260905.py --steps 5000 --seed 20260905
```

Record: [high-density adverse probe](../results/high_b4_adverse_probe_20260905.json).
The mathematical frontier is now (16), not a request to scale the same
sampler or an inference that all \(b_1>0\) corrections are nonnegative.
