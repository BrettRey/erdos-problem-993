# Aristotle formalization queue

Date: 2026-08-21.

This queue sends only frozen, definition-complete results. Every returned project must
be replayed locally with `lake build`, searched for proof escape hatches, checked with
`#print axioms`, and compared declaration-by-declaration with its mathematical packet.

## Active

None.

## Completed and locally replayed

0. **Matroid-representation lower bound.** If matroids `M₁..M_k` on `V(G)` have `IS(G)` as
   their common independent sets, then `deg v ≤ k` at every vertex whose neighbourhood is
   independent (`MatroidRep.le_card_of_represents`), with the triangle-free and tree
   corollaries (`le_card_of_represents_cliqueFree`, `le_card_of_represents_isTree`) and the
   ground-set-enlargement reduction (`represents_comap`,
   `le_card_of_represents_of_embedding`, `le_card_of_represents_minor` for arbitrary minors
   `(P i / C i) ↾ R i`). Project `546ab2b8-189c-4a4b-8316-3bbe6bfac967`, task
   `d8843604-e4ab-4ac2-ad3c-3f7b706446cd`. **Replayed locally 2026-08-26**: `lake build`
   clean, 8,027 jobs; escape-hatch grep empty; `#print axioms` reports only `propext`,
   `Classical.choice`, `Quot.sound` on all 13 declarations (`indep_singleton_of_represents`
   uses none); statements read and compared against the packet. Aristotle improved on the
   packet by dropping the ground-set hypothesis, having noticed the singleton argument
   already forces every vertex into every `(M i).E`. Source:
   `formalization/matroid_representation_lower_bound_aristotle/`.

1. **Full `C₂` assembly and abstract-tree wrapper.** The direct LLYZ coefficient bridge,
   arbitrary connector-state partition, odd/even parity assembly, recognition of the
   explicit models, and their transport to arbitrary finite trees are proved. In
   particular, every finite tree with at most two vertices of degree at least three has a
   log-concave independence polynomial. Project
   `76e84c3b-9ce0-4fc6-9074-45ea38dbb576`, final task
   `0c8f0365-f9a8-4a04-a709-81e14affcfc1`; local build: 8,078 jobs. A proof-escape scan
   is empty, and the principal axiom audits report only `propext`, `Classical.choice`, and
   `Quot.sound`.
2. **Clan normalization continuation: arbitrary even arms.** The explicit finite-family
   component model proves `W(image)=2^e(z+z⁻¹)`, the full `p=2` block identity
   `A(1,2^e;z)`, and `2^e≥1`. Project
   `76e84c3b-9ce0-4fc6-9074-45ea38dbb576`, task
   `f2b3e260-0547-4694-b2a5-b103ab16d201`; local build: 8,038 jobs.
3. **Universal Pascal smoothing.** Formalized
   `notes/pascal_smoothing_shadow_lemma_2026-08-20.md`, including the downward-closed
   incidence inequality, the `8/9` Pascal-transform bound, and its three stated
   consequences. Project `6519b82f-e80b-4c41-b633-990b72c7e3b9`, task
   `e6301287-ecb0-4d75-8b65-b6711945038f`.
4. **Order-ideal projection and graph-transform identities.** Formalize the frozen
   poset layer in `notes/matching_bag_poset_reduction_2026-08-20.md`, especially
   equations (3), (5), (8), (9), and (10), independently of the unresolved blocked
   correction. Project `2a2d691b-91bb-4204-bb00-ad47a51861ba`.
5. **Poset/Pascal bridge.** The transported antichain family is proved downward closed
   and profile-preserving; equations (7)/(7a), depth-eight log-concavity, and the full
   `M≤33` result now follow inside the poset project. Project
   `2a2d691b-91bb-4204-bb00-ad47a51861ba`, task
   `3dfa9612-a379-45a7-9af2-40b697b9aba4`.
6. **Even-connector Laurent block.** The uniform block is centrally unimodal for every
   `r,s≥1` and rational `c,d≥1`; Aristotle replaced the source's awkward range split
   with a uniform Catalan bound. Project `f9cd8c48-bd2e-4aa8-bc67-3a9316657b61`, task
   `dc1051b8-55b2-473f-9f3d-62a5f5e55e00`; local build: 8,029 jobs.
7. **Exhaustive adjacent-two-hub global partition.** The arbitrary-arm degree-`N` map
   space is partitioned by explicit idempotent representatives into blocks of size 4/2/1;
   disjointness, exhaustion, unique membership, total-order preservation, collision
   resolution, exact derived weights, and final degreewise normalized two-row
   nonnegativity are proved. Project `76e84c3b-9ce0-4fc6-9074-45ea38dbb576`, task
   `b78ed6c4-7595-43fd-b8d5-c2f236926a25`; local build: 8,051 jobs. Principal axiom
   audits report only `propext`, `Classical.choice`, and `Quot.sound`.
8. **Tree maximum-matching to forest-poset bridge.** König is derived from Mathlib Hall;
   matching contraction, forest comparison poset, solution/ideal equivalence, equations
   (2) and (6), and all tree-level Pascal inequalities are proved, in fact for arbitrary
   finite forests. Project `2a2d691b-91bb-4204-bb00-ad47a51861ba`, task
   `fcbe1949-5ec5-4c84-8e73-a085b2352297`. All new bridge modules replay directly and
   principal axiom audits report only the three standard Lean axioms. The unrelated old
   `ExhaustiveSmallCodes.lean` `native_decide` target was not awaited.
9. **LAW-V with one trivial hub, internal Lean layers.** Path coefficients, LR convolution,
   `R ≤lr Q`, the shifted cross term, and conditional arbitrary-arm assembly are proved;
   paths with at most two arms and all stars are unconditional. The arbitrary-spider
   log-concavity theorem remains an explicit hypothesis in Lean, although the mathematical
   note closes it by the published LLYZ theorem. Project
   `059ba8f5-39e6-43f0-9944-8d1268d90740`, task
   `9a6c6951-03eb-497b-bdb9-ed9f9013a1b7`; local build: 8,032 jobs.

## Not ready

- the unresolved blocked-correction terms `b₂,b₃,b₄`;
- broad LAW-W or head-strip conjectures;
- Conjecture A or the full Erdős 993 conjecture;
- any quarantined broad STP2 closure statement.
