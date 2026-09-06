# Independent semantic audit of the B1-zero Lean checkpoint

Checkpoint: `formalization/depth3_b1_zero_20260906/checkpoint-timeout-1/`.

Verdict: **PASS for source semantics of the frozen `DepthThree.B1ZeroWindowTarget`, conditional on the separately conducted fresh build and axiom audit succeeding.** I found no missing graph/count, structural, representation, or finite-coverage assumption in the final source proof. The proof takes a valid alternative route through the low-density theorem. It does not establish every sharper intermediate claim in the September 5 note, and the completion map contains two stale scope statements described below.

This is an independent source review, not a second compilation result. At the coordinating agent's request, I ran neither a Lean build nor either large finite kernel checker. I did not treat the completion report or its displayed axiom output as verification. I preserved the checkpoint, the original reviews, and the manifest.

## Inspection scope

I read the frozen specification, the complete target module, the declaration map, all 21 modules under `RequestProject/B1Zero/`, and the relevant complete common modules: `GraphBridge.lean`, `ExtendableCount.lean`, `BlockedCertificate.lean`, `BlockedShadow.lean`, `LowDensity.lean`, `DepthThreeAlgebra.lean`, and `RootedForestCertificate.lean`. I traced the final proof through their actual definitions and proof bodies, including the set injections and recursive representation.

A byte comparison against `formalization/depth3_b1_zero_20260906/input/RequestProject/` found all 31 existing input Lean modules present and byte-identical, including `DepthThreeSpec.lean`. Thus the checkpoint has not weakened the frozen specification or changed the supplied common dependencies to obtain its result. The older matching/poset and Pascal development was inspected in the initial audience-reader audit; this follow-up checks the new use of that development and the new bridges, rather than independently re-proving Mathlib or the unchanged imported library.

All file locations below are relative to the checkpoint's `RequestProject/` directory unless otherwise stated.

## Exact target and actual graph counts

`B1Zero/Target.lean:83` declares exactly:

```lean
theorem b1_zero_window : DepthThree.B1ZeroWindowTarget := by
```

The definition at `DepthThreeSpec.lean:40` universally quantifies over a finite vertex type and a simple graph, and requires precisely: `G.IsTree`, `33 ≤ |V| ≤ 38`, `17 ≤ G.indepNum ≤ 19`, `2 * G.indepNum ≤ |V| + 5`, and `b G 1 = 0`. The conclusion is `s G 2 * s G 4 < s G 3 ^ 2`. There is no supplied decomposition, matching, count profile, representability predicate, or certificate-coverage premise in that target.

The specification's `e` and `b` filter finite **vertex subsets** by the same size and independence conditions, and respectively by containment or noncontainment in a Mathlib maximum independent set. They count each subset once. Their sum `s` therefore has the intended independent-set-count meaning. The defect guard `d ≤ G.indepNum` is unchanged. The deficiency condition is expressed without potentially misleading natural-number subtraction.

The previously missing common bridges are now substantive proofs:

- `GraphBridge.lean:206`, `indepNum_eq_card_Bag`, proves both independence-number inequalities. It does not merely rename `|V|−|M|` as the independence number.
- `GraphBridge.lean:221`, `isMaximumIndepSet_iff_mem_maxIndepSets`, identifies the matching library's maximum-size collection with Mathlib maximum independent sets. `exists_treeMatching_of_isAcyclic` constructs the matching structure from an arbitrary finite forest after finite maximization over matchings.
- `ExtendableCount.lean:110`, `card_extSets_fiber`, gives a genuine finite-set bijection. It records occupied bag support and the selected endpoint word; injectivity recovers the vertex subset, and surjectivity takes the selected vertices of a maximum completion. The counting object is a set of projected words, not a list of completions.
- `ExtendableCount.lean:230`, `e_eq_erasure`, has the necessary `d ≤ card D.Bag` hypothesis. The specification-to-`extSets` and specification-to-`blkSets` identities carry the same guard. The erroneous unguarded equality beyond rank is not used.

The conversion of `b G 1 = 0` to `D.B1Zero` in the final module is also explicit: a hypothetical blocked independent `(α−1)`-set would belong to a finite collection whose cardinality is zero. No stronger `b₂ = 0` assumption appears.

## Structural and small-forbidden branch

The definitions in `B1Zero/Allowed.lean` mean some/every/no maximum independent set, respectively; `FlexV` is allowed minus forced. Nonemptiness of the maximum-set collection prevents vacuous forcedness. The bag classification is proved, and `BagCount.lean` derives `|V| + |ForcedV| = 2|Bag| + |ForbiddenV|` by summing identities for actual bag intersections.

`Pendant.lean` proves the crucial classification from `D.B1Zero`. Its spliced independent set is constructed from restrictions of two maximum sets, and bag-closure supplies the exact `(α−1)` cardinality. The endpoints' external allowed neighbours cannot be joined avoiding the matched edge, by acyclicity. `WellCovered.lean` then chooses one representative per bag, preserving an input independent allowed set and using pendant choices for unoccupied flexible bags. Thus the equivalence between extendability and avoiding forbidden vertices is a conclusion, not an assumed graph property.

`ForcedDegree.lean` constructs representatives avoiding a forbidden vertex's neighbours, and derives a forbidden-containing `(α−1)`-set if the forced-neighbour count were below three. `ForestCount.lean` proves the edge bound used to obtain `|N_C(S)| ≥ 2|S|+1`. Its nonempty-set and disjoint-side hypotheses are supplied at the application. The subsequent `r≤4` argument in `Target.lean` follows from these inequalities and the stated deficiency bound.

For `r≤3`, `Decompose.lean` counts forced and flexible parts by a finite-set bijection. `EBounds.lean` partitions blocked sets by their forbidden intersection and injects each fibre into an allowed residual family, giving an upper bound in the correct direction. `FlexCone.lean` double-counts one-vertex deletion incidences; a size-`l` flexible independent set leaves `m−l` bags available and therefore at most `2(m−l)` possible new vertices.

`Numeric.lean:19`, `density_numeric`, is an explicit universally quantified arithmetic implication about this size-indexed coefficient sequence. The proof splits the bounded integer parameters and uses the proved deletion inequalities and vanishing above `m`. `DensitySmall.lean` supplies these hypotheses from the graph. This is a new direct proof of the low-density condition; it need not reproduce the note's original joint `b₂/e₂`, `b₄/e₄` bound. Its use of natural subtraction does not discard a needed contribution: the fibre injection handles actual feasible cardinalities, and extension of the outer sum only adds nonnegative terms.

## Four-forbidden branch: what the certificate now means

This branch proves `4 b₄ ≤ e₄`, sufficient at `α=19`. It does not attempt to formalize the note's sharper ratio by the same concentration argument.

`B1Zero/FourStructure.lean:86` proves the actual extendable count identity

\[
e_4=126p_0+84p_1+36p_2+9p_3+p_4.
\]

Here `p_j = cntRev G FlexV 10 j`; `cntRev` is the cardinality of independent subsets `S` of the specified vertex set satisfying `S.card + j = n` (`Convolution.lean:32`). This exact equation also gives zero automatically beyond the reference rank; it avoids an off-by-one or truncated-subtraction reinterpretation of coefficients.

`FourStructure.lean:196` proves

\[
b_4\le73p_0+24p_1+3p_2+15q_0+6q_1+q_2,
\]

with `q_j` counting flexible independent sets avoiding **all** flexible vertices adjacent to forbidden vertices. The key inequality at `FourStructure.lean:147` is a correct union-bound rearrangement:

\[
\sum_{u\in U}q_{u,j}\le3p_j+q_j.
\]

Indeed, any independent set that does not avoid the full attachment set fails to avoid the neighbours of at least one forbidden vertex, so it contributes at most three times to the left side; a set avoiding all attachments contributes four times. The proof realizes this via finite-set differences and a covering union. Single-forbidden fibres contribute `15q_{u,0}+6q_{u,1}+q_{u,2}`; the six two-forbidden fibres contribute at most `4p₀+p₁` each, and the four three-forbidden fibres at most `p₀` each. Empty and four-forbidden fibres vanish. These terms give exactly `(73,24,3;15,6,1)`.

The representation and domination steps are also supplied:

- `FlexCorona.lean` constructs `CoronaData` from the proved pendant matching and proves that its full corona equals `FlexV`. Private leaves are only required to be private inside that corona, so attachments to forbidden vertices do not contradict the structure's meaning.
- `exists_attach_reach` uses the original graph's connectedness to show that every flexible component meets the attachment set. `exists_root_transversal` at `FlexCorona.lean:255` chooses one base vertex per base component, with that vertex or its pendant in the attachment set. Exact uniqueness of graph attachments is not required by this alternative argument.
- `CoronaMove.lean:25`, `cntRev_move_roots`, gives an explicit cardinality-preserving injection from sets avoiding the original attachment set to sets avoiding the chosen base roots. Only selected base vertices whose pendants were originally excluded are replaced; the source therefore cannot already contain those pendants. The proof includes an explicit recovery formula, and the inequality direction enlarges `q` as needed for its nonnegative numerator coefficients.
- `CoronaCount.lean` proves the profile recursions from disjoint finite-set partitions, splitting on root included, pendant included, or neither. At reference rank `|X|`, the “neither” case shifts the defect by one. The five-coordinate truncation is sufficient and correct because the product at depth at most four only involves depths at most four.
- `CoronaRepr.lean:28`, `exists_forest_repr`, is a universal theorem for arbitrary finite graph coronas with an acyclic base and a root transversal. Strong induction on the base set constructs an actual inductive `Forest` term with the correct size and **equal actual graph profiles**. Peeling a component and its root produces the child and remainder transversals required by the recursive evaluator.

Consequently `FlexCorona.lean:284`, `flex_certificate`, really transports the computed inequality to graph counts. It applies the representation theorem, applies root-moving domination for the `q` coordinates, and preserves equality for the `p` coordinates. It does not assume the graph is already one of the listed codes.

## Finite coverage and computational boundary

The finite object is a **plane rooted forest**, with constructors `nil` and `cons children rest`. It deliberately includes different orderings of an unordered rooted forest. `RootedForestCertificate.enumeration_complete` at line 47 proves by induction that every term of size `n` occurs in the fuelled enumeration for `n`. The graph-to-term theorem above supplies the term of size ten. There is no dependence on an unproved NetworkX census or on the JSON rows from the original written proof.

The relevant new predicate is `simpleWithinBound` in `B1Zero/SimpleCertificate.lean:23`, exactly `4 * simpleNumerator f ≤ denominator f`. `all_simple_codes_checked` at line 28 uses `decide +kernel`, and `rooted_forest_simple_bound` at line 33 combines that Boolean result with the proved enumeration membership. The list's 16,796 entries replace the 1,842 unordered representatives for this proof. Canonical uniqueness up to graph isomorphism is unnecessary for a universal inequality; extra ordered representations cannot create an omission.

I inspected the evaluator and its graph bridge rather than trusting the comment that its coordinates are graph coefficients. The equation `profile f = (C.vec X ∅, C.vec X Y)` proves the needed meaning. The older sharper `rooted_forest_bound` remains in the imported input module but is not the bound used by this final argument.

I did not execute the 16,796-case reduction. Successful kernel elaboration and the final axiom output remain the coordinating agent's separate verification responsibility. The source contains the correct proposition to check and a complete semantic chain from that proposition to the graph target.

## Final assembly and scope limitations

`FourStructure.lean:345` combines the count equalities, union bound, and transported certificate to prove `4 * card(blkSets 15) ≤ card(extSets 15)`. `Target.lean` derives `α=19`, `|ForcedV|=9`, and ten flexible bags when `r=4`, converts this to the general density condition, and joins it to the `r≤3` branch.

`LowDensity.lean:126` applies the actual-graph erasure bridge to the existing Pascal reserve, proves positive `e₂,e₄`, supplies the actual blocked-shadow and extendable-incidence inequalities, and invokes the scalar algebra lemma. All casts of `α−k` to real numbers have explicit lower-bound justifications. The final strict comparison is transported back to the natural graph counts. The general blocked-shadow proof uses a proved blocked vertex subset of cardinality at most two; it does not assume that every deletion of a blocked set stays blocked.

No semantic deficiency in the target proof was found. The following reporting limits should be retained:

1. **Documentation scope error:** `DECLARATION_MAP.md:111` says nothing in the development assumes connectedness, and line 113 says no vertex-count, deficiency, or `b₁` restriction is used. Those bullets must be labelled as statements about the general blocked-shadow/low-density layer. Read as claims about the whole completed checkpoint, they are false: the B1-zero target and its proof expressly use those hypotheses. This does not weaken or invalidate the actual declaration.
2. **Different intermediate result:** this checkpoint does not, through its final dependency chain, establish the note's sharp graph inequality `217404 b₄ ≤ 52513 e₄` or the special original-tree bound `11b₂ ≤ b₃`. It proves the weaker quarter-density bound and uses the general blocked-shadow theorem instead. The declaration map's later explanation of this route is accurate; completion of the frozen target must not be reported as formalization of every displayed intermediate equation in the note.
3. **Bounded theorem:** the verified mathematical meaning remains the stated B1-zero tree window. It neither closes the B1-positive window nor proves Erdős #993.

**Source-semantic conclusion:** the graph/count and representation/coverage gaps identified in the initial review have been discharged by explicit source proofs in this checkpoint. Subject to the separate clean build and axiom audit, `DepthThree.b1_zero_window` establishes the exact frozen theorem, not a conditional surrogate with missing bridges hidden in its hypotheses.
