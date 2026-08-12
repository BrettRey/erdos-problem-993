# Atlas truncation: killed at the bottom level

Date: 2026-08-12. Companion script:
`verify_atlas_bottom_level_20260812.py` (exact arithmetic; run it, exit
0 expected). Produced with generative-AI assistance; the construction
below is instantiated from the full texts of arXiv:2203.01533 and
arXiv:2110.10740 in the research vault, not from memory.

## The question

The prefix-certification survey (2026-08-12) left one framework as
"open research": truncate the Chan-Pak combinatorial atlas's local
conditions by coefficient index, hoping to certify log-concavity of
tree independence sequences below the Levit-Mandrescu threshold. The
stated reopening trigger was a cheap computation: build the atlas's
local matrices for small tree independence complexes and see WHERE the
local conditions fail. If the failures had localized to head-level
structures, the lane would have come alive. This note reports that
computation. The lane is dead, with a clean mechanism.

## The instantiation

The atlas for a matroid and certifying index k places its entire
direct verification burden on sink vertices `(alpha, 0, 1)`, whose
matrix is `A(alpha, 1)`: rows and columns indexed by ground elements
plus a border element `null`, with `A[x][y] = 1` when `alpha+x+y` stays
in the language, `A[x][null] = 1` for every extendable `x`, and
`A[null][null] = 1` (their eq. 4.2 with m = 1). The required property
(Hyp) is equivalent to having at most one strictly positive eigenvalue
(their Lemma 3.5). Their Lemma 4.6 proves every sink hyperbolic by
introducing the conflict relation `x ~ y iff alpha+x+y` leaves the
language and observing that its TRANSITIVITY follows from the matroid
exchange property; the support-restricted matrix is then a bordered
complete multipartite matrix, which has exactly one positive
eigenvalue.

The same construction reads verbatim for a graph independence complex.
But at the empty word the conflict relation is graph ADJACENCY, which
is transitive precisely when the graph is a disjoint union of cliques.
No tree on three or more vertices is.

## The computation (exact, replayable)

Harness controls: the free matroid and a two-parallel-class matroid
(bordered complete multipartite `K_{3,3}`) give exactly one strictly
positive eigenvalue, as Lemma 4.6 requires. (A first version of the
count was wrong in an instructive way: sympy's `count_roots(0, oo)`
includes the endpoint, and these bordered matrices are rank-deficient,
so zero eigenvalues inflated the count. The verifier strips zero roots
exactly before counting.)

Tree results for `A(empty, 1)`, strictly positive eigenvalue counts:

- path `P6`: **3** (fails (Hyp))
- claw `K_{1,3}`: **2** (fails)
- spider (2,2,2): **2** (fails)
- exhaustive, all trees on 3 to 9 vertices: **every one fails**, and
  (Hyp)-holding coincides exactly with the conflict graph being a
  cluster graph across the sweep.

## Why this kills truncation specifically

1. **The failing vertex is in every atlas.** The empty-word sink
   `(empty, 0, 1)` is reachable from the top vertex of the atlas for
   EVERY certifying index k (follow null-edges down). Truncating by
   coefficient index removes no part of the failure.
2. **The failure anticorrelates with the pathology.** The path - real-
   rooted, ultra-log-concave, the best-behaved tree there is - fails
   with MORE positive eigenvalues (3) than the claw or the spider (2).
   The local-condition failure is not tracking the head-window
   geometry that organizes everything else in this problem; it is
   tracking distance-2 structure in the tree, which every tree has.
3. **The mechanism is the exchange property, localized.** Bottom-level
   transitivity of the conflict relation IS the exchange property in
   the atlas's own proof. This turns the survey's structural
   observation ("the engine is exchange, not index-uniformity") into a
   verified computation: the atlas certificate for trees does not fail
   at some deep level a truncation could avoid; it fails at the root,
   for every tree, for the exact reason matroids differ from graphs.

## Standing conclusion

Atlas truncation is closed as a lane for the tree problem. Any
atlas-like approach to graph independence complexes would need a
replacement bottom-level certificate for non-cluster conflict graphs,
which is a new local-global theory, not a modification. The reopening
trigger named in the survey has fired negative. Proof effort stays on
the ordered queue: LAW-V/LAW-W (packet frozen at
`gpt_attack/law_v_packet_2026-08-12/`), the U head strip, and B's body
log-concavity.
