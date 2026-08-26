# Attack packet: the matroid-representation number of a graph's independence system

**For:** an independent model (Ox Alpha via `ocx`, or any strong reasoner).
**Date frozen:** 2026-08-26.
**Disclosure:** safe to send externally. Everything here is either public
(Erdős Problem #993, arXiv:2608.23262) or a new elementary lemma the author is
happy to have checked in the open. No manuscript text, no referee material.

---

## 0. What I want from you, in one sentence

Attack the lemma in §2 until it breaks or you are convinced it holds, then tell
me **what this invariant is already called in the literature** and what exact
strings I should search to find it, without inventing a single citation.

**Refutation is a success mode.** If the lemma is false, or true only under a
definition I did not state, saying so precisely is the most valuable outcome
available. An honest "T1 REFUTED, here is the counterexample" beats a confident
survey.

---

## 1. Background, just enough to make the question well posed

Erdős Problem #993 asks whether the independent-set sequence of every tree is
unimodal. Write $i_k(G)$ for the number of independent sets of size $k$, and
call the sequence **log-concave (LC)** if $i_k^2 \ge i_{k-1} i_{k+1}$ for all
$k$, and **ultra log-concave (ULC)** if $i_k / \binom{n}{k}$ is log-concave.

A recent preprint (Schweitzer, arXiv:2608.23262, Aug 2026) proves:

- **(S1)** For any *two* matroids on a common $n$-element ground set, the numbers
  $I_k$ of common independent sets of size $k$ form a ULC sequence.
- **(S2)** This fails at *three*: he gives three partition matroids whose common
  independent sets have size sequence $(1, m+3, 3, 1)$, which is not LC for
  $m \ge 7$.

I want to know exactly where trees sit relative to that 2-vs-3 boundary. That
requires knowing how many matroids a graph's independence system actually needs.

---

## 2. The lemma under attack

Let $G$ be a **finite simple** graph with vertex set $V$, and let
$$\mathrm{IS}(G) = \{\, S \subseteq V : S \text{ contains no edge of } G \,\}.$$
For a matroid $M$ on ground set $V$, write $\mathrm{IN}(M)$ for its independent
sets. Define the **matroid-representation number**
$$\mathrm{mi}(G) \;=\; \min\Big\{\, k \;:\; \exists\ \text{matroids } M_1,\dots,M_k
\text{ on ground set } V \text{ with } \mathrm{IS}(G) = \bigcap_{i=1}^{k} \mathrm{IN}(M_i) \,\Big\}.$$
(Convention: the empty intersection is $2^V$, so $\mathrm{mi}(G) = 0$ iff $G$ has
no edges.)

Call a spanning subgraph of $G$ an **equivalence subgraph** if it is a disjoint
union of cliques, i.e. it is induced by a partition of $V$ every class of which
is a clique of $G$. Define $\mathrm{eqcov}(G)$ to be the least number of
equivalence subgraphs of $G$ whose union covers $E(G)$.

> **LEMMA.**
> **(a)** $\mathrm{mi}(G) = \mathrm{eqcov}(G)$ for every finite simple graph $G$.
> **(b)** If $G$ is triangle-free with at least one edge, $\mathrm{mi}(G) = \chi'(G)$, the chromatic index.
> **(c)** If $T$ is a tree with at least one edge, $\mathrm{mi}(T) = \Delta(T)$, the maximum degree.

### My proof, stated so you can find the hole

**Lower bound, $\mathrm{mi}(G) \ge \mathrm{eqcov}(G)$.** Suppose
$\mathrm{IS}(G) = \bigcap_i \mathrm{IN}(M_i)$.

1. Every singleton $\{v\}$ is in $\mathrm{IS}(G)$, hence in every
   $\mathrm{IN}(M_i)$, so no $M_i$ has a loop.
2. Fix $uv \in E(G)$. Then $\{u,v\} \notin \mathrm{IS}(G)$, so $\{u,v\} \notin
   \mathrm{IN}(M_i)$ for some $i$. Being a dependent 2-set in a loopless matroid,
   $\{u,v\}$ is a circuit, i.e. $u$ and $v$ are **parallel** in $M_i$.
3. Parallelism is an equivalence relation on a loopless matroid, so the parallel
   classes of $M_i$ partition $V$.
4. Conversely, if $u,v$ are parallel in some $M_i$ then $\{u,v\} \notin
   \mathrm{IN}(M_i)$, so $\{u,v\} \notin \mathrm{IS}(G)$, so $uv \in E(G)$.
   Hence every parallel class of $M_i$ is a **clique of $G$**.
5. By 2 and 4, the $k$ parallel-class partitions are equivalence subgraphs of $G$
   covering $E(G)$. So $k \ge \mathrm{eqcov}(G)$.

Note step 5 uses nothing about circuits of size $\ge 3$: the bound holds for
arbitrary matroids, not just partition matroids.

**Upper bound, $\mathrm{mi}(G) \le \mathrm{eqcov}(G)$.** Given equivalence
subgraphs $P_1,\dots,P_k$ covering $E(G)$, let $M_i$ be the partition matroid of
partition $P_i$ with capacity 1 on every class. A set $S$ lies in every
$\mathrm{IN}(M_i)$ iff no two of its vertices share a class, iff $S$ contains no
covered edge, iff (by covering) $S$ contains no edge at all, iff
$S \in \mathrm{IS}(G)$.

**(b)** Triangle-free means every clique has $\le 2$ vertices, so an equivalence
subgraph is a matching plus isolated vertices, and covering $E(G)$ by $k$
matchings is a proper $k$-edge-colouring.

**(c)** Trees are bipartite, so $\chi'(T) = \Delta(T)$ by König's edge-colouring
theorem. Direct lower-bound argument, avoiding König: at a vertex $v$ of degree
$\Delta$, no two neighbours of $v$ are adjacent (no triangles), so in any single
partition $v$ shares its class with at most one neighbour, hence each equivalence
subgraph covers at most one edge at $v$; so $\Delta$ of them are needed.

### What I already verified computationally
- The crude version ($n-1$ partition matroids, one per edge, block $\{u,v\}$)
  reproduces $\mathrm{IS}(T)$ on 32 trees with $n \le 9$.
- The sharp version: on **all 200 trees with $2 \le n \le 10$**, a proper
  $\Delta$-edge-colouring exists, each colour class is a matching, and the
  $\Delta$ induced partition matroids have common independent sets exactly
  $\mathrm{IS}(T)$.

Neither check touches the *lower* bound (that fewer matroids cannot work), which
is where I most want adversarial pressure.

---

## 3. Targets, graded

**T1 (highest value) — break the lemma or certify it.**
Attack specifically:
- **Ground set.** My definition forces all $M_i$ onto ground set $V$ itself. If
  auxiliary elements were allowed (represent $\mathrm{IS}(G)$ as the restriction,
  or a projection/contraction, of an intersection of $k$ matroids on a larger
  ground set), **can $k$ drop below $\mathrm{eqcov}(G)$?** Give a graph and an
  explicit construction, or an argument that it cannot. This is the gap I am
  least sure about, and it decides whether "$\mathrm{mi}$" is even the right
  invariant to name.
- **Higher circuits.** Step 5 claims circuits of size $\ge 3$ cannot help. Is
  there a matroid whose 3-element circuits kill an edge in a way my parallel-class
  argument misses? (I believe not, since a *2-set* being dependent is what an
  edge demands, but say so precisely or break it.)
- **Degenerate cases.** $k = 0$; $G$ edgeless; isolated vertices; disconnected
  $G$; $\Delta = 1$; $n = 1$; $G$ with a loop (excluded by "simple", but confirm
  nothing silently breaks).
- **Finiteness.** Does anything require finiteness beyond the obvious?

Verdict format: `T1: HOLDS-AS-STATED` / `HOLDS-WITH-REPAIR (state repair)` /
`REFUTED (give the counterexample explicitly, with the graph and the matroids)`.

**T2 (high value) — name it, do not cite it.**
This smells classical. I suspect $\mathrm{eqcov}$ is the **equivalence covering
number** studied in the 1980s, and that "every independence system is an
intersection of finitely many matroids, and the least number is an invariant" is
older still, possibly in the greedy-approximation literature on $k$-systems
(Korte–Hausmann, Jenkyns, Fisher–Nemhauser–Wolsey).

**Hard rule for T2: give me no bibliographic data you cannot stand behind.**
Do NOT produce author/year/title/venue tuples from memory; that is precisely the
failure mode I am guarding against. Instead give me, in order of usefulness:
1. **Concept names** the thing is likely to travel under, each with a confidence
   label (`CERTAIN` / `LIKELY` / `GUESS`).
2. **Exact search strings** I should run against Google Scholar, arXiv, zbMATH,
   MathSciNet — verbatim, ready to paste, ideally 8–15 of them, aimed at
   different vocabularies (matroid theory, graph covering, combinatorial
   optimization, approximation algorithms).
3. Any theorem *statements* you are confident exist, phrased as statements, with
   the label `STATEMENT-I-BELIEVE-EXISTS-BUT-CANNOT-CITE`.

If you are confident about a specific named theorem, name the theorem, not a
reference string.

**T3 (medium value) — is the tree corollary the sharpest reading?**
Given (S1)/(S2), my reading is: Schweitzer's theorem covers exactly the trees with
$\Delta \le 2$ (paths), and subcubic trees ($\Delta = 3$) are the first rung
past it, which is also the rung where his counterexample lives. Empirically,
**every one of the 1,230 known LC-failing trees ($n \le 32$) has $\Delta \ge 4$,
and subcubic trees are LC exhaustively through $n = 33$** (2,076,217,086 trees,
zero failures).

Questions: (i) Is there a sharper or different invariant that reads Schweitzer's
boundary onto trees better than $\Delta$? (ii) Does anything in the two-matroid
case (S1) give a *non-trivial* consequence for trees, given that only paths
qualify? (iii) Schweitzer's counterexample is a complete split graph, which is
triangle-rich. Is there a **triangle-free** independence system at the
three-matroid rung that fails LC? If yes, subcubicity alone cannot be what saves
trees, and acyclicity has to carry the argument. This last one is a concrete
finite search (subcubic bipartite graphs) and a genuine fork in the road.

---

## 4. Rules of engagement

- **No fabricated citations.** See T2. A wrong DOI or a plausible-sounding
  author/year is worse than nothing here and will cost more to detect than the
  answer is worth.
- **Exact arithmetic** for any numeric claim. If you compute an independence
  polynomial, use integers and show the coefficient list.
- **Explicit counterexamples**, not existence claims. "Some graph with property
  X would break this" is not a refutation; a named graph with its matroids is.
- **Mark every uncertainty inline.** I would much rather have a hedged true
  statement than a crisp false one.
- **Do not restate my proof back to me.** Assume I know it. Spend the tokens on
  where it fails, or on T2/T3.

## 5. Self-grading bar (required section in your output)

End with a short block, no more than 10 lines:

```
T1: <HOLDS-AS-STATED | HOLDS-WITH-REPAIR | REFUTED>  — one line why
T2: <how many concept names CERTAIN, how many search strings, 0 citations invented?>
T3: <ANSWERED | PARTIAL | NOT-REACHED> — one line
Biggest thing I am unsure about in my own answer: <one line>
What I would attack next if I had another hour: <one line>
```

An honest PARTIAL with one real proved step beats a grandiose complete-looking
answer. If you spent your effort on T1 and never reached T3, say so plainly.
