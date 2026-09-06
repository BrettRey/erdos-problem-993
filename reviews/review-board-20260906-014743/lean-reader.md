# Lean audience-reader report

## Initial reader record

**Ten-second impression (title and result):** This looks like a reduction from a graph-theoretic condition, `b₁ = 0`, to a small exact certificate proving strict depth-three log-concavity in a bounded parameter window. The result box clearly limits the claim to this branch and explicitly denies Lean verification. I expect substantial work connecting graph counts to the certificate and to the existing reserve theorem.

**Sixty-second orientation (before detailed source and Lean inspection):**

| Question | Initial status | What I can recover immediately |
| --- | --- | --- |
| Exact quantifiers and count definitions | PARTLY | Trees in the stated `(n, α, δ)` window; extendable and blocked independent sets at size `α − d`. An explicit packaged pair of Lean target statements is absent. |
| Scope of each theorem | PARTLY | A general forest structural lemma and a bounded tree inequality are distinguishable; the request's “two new graph-count theorems” is not explicitly named in the opening. |
| Dependencies already proved | PARTLY | The Pascal reserve is identified by a frontier label, but its exact declaration, hypotheses, and count interface are not provided here. |
| Computation and completeness obligations | PARTLY | Thirteen parameter triples and a final 1,842-case calculation are announced, but the bridge from arbitrary forests to the finite objects is not visible in the orientation read. |
| Viable Lean dependency order | PARTLY | Structure → coefficient bounds → exact certificates → reserve → strict inequality is recognizable, but concrete reusable graph infrastructure is unresolved. |

## Outcome after the full read

The two desired graph-count conclusions are reconstructible, and the existing Lean projects contain useful, checked infrastructure. They do **not** already supply the graph-count reserve interface needed to assemble either conclusion. The most important missing common step is a cardinality-preserving correspondence between actual extendable independent sets and projected matching-bag words. The large additional obligation for the `b₁ = 0` theorem is the universal structural reduction and the completeness/correctness bridge to the 1,842 finite representatives.

I read the entire frozen principal note, the frozen forest-poset note, both frozen verifier programs, and parsed both complete JSON certificates. I inspected the declarations and relevant proof bodies in both named Lean projects. I did not read other reviews. The computations and Lean checks actually run are distinguished below from work that remains.

| Question | Final status | Resolution |
| --- | --- | --- |
| Exact quantifiers and count definitions | ANSWERED, with an explicit convention | Use the existing graph counts, with defect counts zero when `d > α`; take the window theorem over finite trees. State the auxiliary general incidence lemma with `α ≥ 4`. |
| Scope of each theorem | PARTLY in the source; reconstructible | The opening names only the `b₁ = 0` window result. Section 9 supplies a second strict log-concavity criterion and its blocked-shadow lemma, without packaging them as named targets. |
| Dependencies already proved | PARTLY | The code/poset reserve and the exact public graph-count definitions exist. Their connecting graph-count theorem does not occur in the inspected projects. |
| Computation and completeness obligations | PARTLY | Both exact certificates replay. The formal graph-to-representative reduction, enumeration coverage, and evaluator correctness remain to be proved. |
| Viable Lean dependency order | ANSWERED | Finish the shared graph/count bridge, derive the incidence and reserve inequalities, close the low-density theorem, then do the structural and finite work for the `b₁ = 0` theorem. |

## Exact theorem reconstruction

The phrase “two new graph-count theorems” is not a label used by the principal note. My initial orientation could therefore identify the main result and the Section 9 incidence lemma as the pair. For a complete verification deliverable, the Section 9 **low-density strict log-concavity conclusion** should be the second public theorem, with the incidence lemma a mandatory auxiliary theorem. This is a statement-packaging ambiguity, not evidence against either argument.

Use `α = G.indepNum`, `n = Fintype.card V`, `e_d = FivefoldForest.extendableCount G d`, `b_d = FivefoldForest.blockedCount G d`, and `s_d = FivefoldForest.iD G Finset.univ d`. The existing bridges identify `s_d = e_d + b_d`. All of these count actual vertex subsets, once each.

1. **Window theorem:** For every finite simple graph `T` that is acyclic and connected, if `33 ≤ n ≤ 38`, `α ∈ {17,18,19}`, `2α − n ≤ 5`, and `b₁ = 0`, then `s₂ s₄ < s₃²`. This is the opening result and Sections 1–7. In a natural-number statement, `2α ≤ n + 5` avoids silently interpreting the deficiency through truncated subtraction. The forest lower bound `n ≤ 2α` is also available. The connectedness hypothesis is used in the four-forbidden attachment reduction and must remain in this public theorem unless a separate forest argument is supplied.
2. **Low-density theorem:** For every finite forest `F`, if `α ∈ {17,18,19}` and `3(α−3)b₄ ≤ (α−7)e₄`, then `s₂ s₄ < s₃²`. This is Section 9, equations (14)–(15). It imposes no `b₁` or correction-sign condition. The written argument does not use the bounds on `n` or `δ`, or connectedness; retaining the narrower current-window scope would also be faithful. The cross-multiplied hypothesis is exactly `b₄/e₄ ≤ h_α` because `e₄ > 0` in this rank range.

Mandatory auxiliary statements are `6b₃ ≥ (α−4)b₂` for finite forests with `α ≥ 4`, and `4e₄ ≤ (α−3)e₃` for finite graphs with `α ≥ 4`. For the four-forbidden case, the exact graph-count targets are `217404 b₄ ≤ 52513 e₄` and `11b₂ ≤ b₃`, under the Section 7 tree hypotheses. They are intermediate results on the original tree, not merely inequalities on its dominating representative.

For both strict conclusions, use the existing positivity lemma to establish `e₂e₄ > 0`. In Section 9, even equality at the density threshold is allowed: `c_α−h_α > 0` supplies strictness. An implementation should not replace the threshold's `≤` by `<`.

## Exact declaration and path map

The following path prefixes are exact repository-relative paths:

- `FF/` = `formalization/fivefold_forest_20260904/output-final_aristotle/RequestProject/`
- `MB/` = `formalization/matching_bag_poset_aristotle/RequestProject/`
- `SRC/` = `reviews/review-board-20260906-014743/source/`

The table names existing declarations unless it explicitly says **missing**.

| Required interface | Declaration and exact file location | What it actually supplies |
| --- | --- | --- |
| Actual graph counts | `FivefoldForest.Extendable`, `extendableCount`, `blockedCount`; `FF/Main.lean:21`, `:25`, `:34` | Maximum-independent-set containment and defect counts, guarded by `d ≤ G.indepNum`. |
| Counts on vertex subsets | `FivefoldForest.indFam`, `alpha`, `Ext`, `icnt`, `ecnt`, `bcnt`, `eD`, `bD`; `FF/Defs.lean:29` onward | A useful fixed-ambient-graph interface for deletion and induced subgraphs. |
| Public/internal graph bridge | `FivefoldForest.alpha_univ`, `isMaximumIndepSet_iff`, `extendable_iff`, `extendableCount_eq`, `blockedCount_eq`; `FF/Main.lean:52`, `:64`, `:74`, `:84`, `:96` | These connect the public graph definitions to the internal counting framework. |
| Full count and partition | `FivefoldForest.iD`, `eD_add_bD`, `icnt_eq_iD`; `FF/Union.lean:30`, `:33`, `:40` | The required `s_d` and `s_d = e_d+b_d`, including the out-of-range convention. |
| Positivity and hereditary extendability | `FivefoldForest.eD_pos`, `Ext.mono`; `FF/Basic.lean:128`, `:57` | Nonzero extendable layers below `α`; subsets of an extendable set extend. |
| Component products | `FivefoldForest.Split`, `alpha_split`, `ext_split_iff`, `icnt_split`, `ecnt_split`, `iD_split`, `eD_split` (all in namespace `FivefoldForest`); `FF/Union.lean:24`, `:117`, `:182`, `:175`, `:208`, `:275`, `:283` | Reusable componentwise independence-number, extension, and convolution machinery. The requisite graph decomposition is still an input to these lemmas. |
| Matching setup | `MatchingBag.TreeMatching`; `MB/TreeMatching.lean:38`; `MatchingBag.exists_treeMatching`; `:59` | The structure contains acyclicity, a proper two-colouring, and a maximum matching. The existence theorem still requires an explicit maximum-matching witness and its maximal-cardinality proof. Choosing one from a finite set is additional wrapper work. |
| Matching/cover facts | `MatchingBag.konig`; `MB/KonigHall.lean:137`; `MatchingBag.TreeMatching.exists_minCover`, `maxIndepSets_eq_image_compl`, `minCover_structure`, `indep_not_both`; `MB/TreeMatching.lean:327`, `:332`, `:361`, `:600` | Cover existence, complements, bag occupancy for extrema, and at most one endpoint per matching edge for independent sets. |
| Rank of the code | `MatchingBag.TreeMatching.card_Bag`; `MB/TreeMatching.lean:271` | Its actual conclusion is `card D.Bag = card V − D.M.card`. It does **not** state equality with `D.G.indepNum`. |
| Maximum-set definition | `MatchingBag.TreeMatching.maxIndepSets`; `MB/TreeMatching.lean:317` | Defined by independence and size `card V − D.M.card`; an explicit equivalence with Mathlib `IsMaximumIndepSet` is **missing** in this project. |
| Free/forced coordinate representation | `MatchingBag.TreeMatching.Forced`, `Free`; `MB/BagPoset.lean:41`, `:47`; `treeCode_eq_codeRelabel`; `MB/TreeCodeBridge.lean:169` | A proved equivalence of the full code with a relabelled/flipped ideal code plus constants. `Forced` is a property of a matching coordinate, not yet the vertex set `C` of the new note. |
| Partial ideal extension | `MatchingBag.codeProj_idealCode`; `MB/PosetCode.lean:62` | Exactly characterizes extendable partial assignments on free coordinates by the induced order-ideal condition. This is the right reuse point for one- or two-vertex blocked certificates, after the graph/partial-word translation is proved. |
| Projection counts and complement gauge | `MatchingBag.codeP`; `MB/Codes.lean:32`; `MatchingBag.TreeMatching.codeP_maxIndepCode`; `MB/TreeCodeBridge.lean:211` | Counts distinct projected words; cover and maximum-independent-set codes have the same profile. It does not count arbitrary partial vertex subsets directly. |
| Code/poset profile bridge | `MatchingBag.TreeMatching.erasure`, `erasure_eq_erasureProfile`; `MB/TreeCodeBridge.lean:253`, `:256` | Defines `erasure d = codeP D.treeCode (card D.Bag − d)` and relates it to the poset profile when `d ≤ card D.Bag`. |
| Quantitative reserve | `MatchingBag.TreeMatching.erasure_depth_three_reserve`; `MB/TreeCodeBridge.lean:306` | Proves `32(M−2)(erasure 2·erasure 4) ≤ 27(M−3)(erasure 3)²`, with `M = card D.Bag` and `4 ≤ M`. This is the precise existing dependency. |
| General graph-count reserve | **Missing:** rank/maximum-set bridge and `extendableCount D.G d = D.erasure d` for `d ≤ D.G.indepNum` | Needed before substituting actual `e₂,e₃,e₄` into the reserve. No declaration in either inspected project supplies this correspondence. |
| Incidence inequalities needed here | **Missing:** blocked shadow and the lower extension-count bound for extendable sets | `FivefoldForest.ecnt_incidence` at `FF/Basic.lean:137` has conclusion `(k+1)ecnt(k+1) ≤ |S|ecnt(k)`. Its direction and constant do not give `4e₄ ≤ (α−3)e₃`; its finite-incidence proof pattern is reusable. |
| Existing fivefold theorem | `FivefoldForest.fivefold_of_forest`; `FF/Main.lean:115` | Requires `α ≥ 5` and `blockedCount G 2 = 0`. The new `b₁ = 0` hypothesis does not imply `b₂ = 0`; this theorem cannot close the general new branch. |
| Forbidden roles, allowed pendant matching, corona/rooted representation, domination, exhaustive representative coverage | **Missing** from the inspected Lean projects | These are new graph-theoretic and executable-combinatorics obligations. The local `FivefoldForest.Pendant` structure in `FF/Pendant.lean:39` describes a supplied pendant configuration; it is not a proof of the new universal pendant-matching classification. |

The code count's out-of-range behaviour deserves an explicit warning in the interface: `D.erasure d` uses natural subtraction without a guard, so beyond the number of bags it evaluates at retained size zero. The public graph counts instead return zero beyond `α`. The proposed equality must carry `d ≤ α` or use a separately guarded code profile. This causes no difficulty at depths 2–4 in either target, but an unqualified equality for every `d` is false.

## Friction log

| Location | Reader friction and consequence | Smallest useful handoff improvement |
| --- | --- | --- |
| Opening and Section 9 | The headline presents one theorem; the later density conclusion changes what a request for “two theorems” means. | Give the two public statements together, then name the general blocked-shadow lemma separately. |
| Section 1, lines 46–72 | “Every component is a union of whole remaining bags” and componentwise gluing do substantial work. A formalizer must prove that deleting the two endpoints of one matched edge leaves whole bags in each component, and that the glued set has exactly `α−1` vertices. | Break out the bag/component restriction and gluing lemma; do not supply the desired pendant matching as an assumption to the final theorem. |
| Section 1, lines 74–85 | Simultaneously avoiding nonforced neighbours and the forest edge bound depend on componentwise maximum-set facts. | State a lemma that a forbidden vertex meets each allowed component at most once, then prove the avoidance construction and `|N_C(v)| ≥ 3`. |
| Section 4, lines 202–210 | “Existing formalized dependency” can be read as an immediately importable graph-count reserve. The existing declaration instead concerns code erasures. | Name `erasure_depth_three_reserve` and list the graph-count translation as outstanding. |
| Section 3, lines 165–198 | The five-ray argument is a usable mathematical proof, but the ratios are not directly executable with the noncomputable graph counts. | Prove the monotone-sequence cone lemma over rational or real numbers, keep graph counts natural, and explicitly establish positivity of the denominators. Handle absent `t > r` terms before casting expressions such as `c−5`. |
| Sections 5 and 7.2 | The handoff changes from arbitrary trees to a tight 13-vertex core, coronas of forests, and rooted component types. Each change needs an isomorphism or count-preserving/dominating map. | Package a structural witness carrying the base forest, pendant pairing, exactly one attachment per component, and the relevant count equalities. Prove that every input tree supplies it. |
| Section 7.1, equations (9)–(10) | Polynomial identities and nonnegativity do not by themselves identify the polynomials with graph counts. | Prove the forced-intersection partition and group product formula first; then use the displayed algebra for concentration. |
| Section 7.2, enumeration paragraph | Agreement between two generated code sets is compelling executable crosschecking, but is not a theorem that every rooted forest has one of those codes. | Prove encoding/decoding and coverage of the recursive finite enumeration. The universal inequality needs coverage and invariant evaluation; a full uniqueness-up-to-isomorphism theorem is unnecessary unless the exact type count is itself a formal target. |
| Section 7.3 | The positive `b₃` lower bound must hold on the original tree. A coefficientwise upper representative cannot establish it. | Keep the original-tree incidence argument as its own theorem. The source correctly makes this distinction. |
| Section 9, lines 524–539 | The ideal-code projection theorem concerns partial words, while the claimed certificate is a subset of actual vertices. The translation must preserve containment and cardinality, including singleton and forced matched bags. | Prove that every blocked independent vertex set contains a blocked subset of cardinality at most two, then perform the deletion-incidence count on vertex sets. |
| Frozen scripts and links | The snapshot is not a standalone executable packet: both scripts infer `REPO` from `__file__`, import unfrozen repository helpers, and one loads three archived result files. Several relative links still point to the original note layout. | Include a dependency manifest and a replay command for the frozen layout, or freeze the helper modules and repair the paths in a later packet. No source repair was made in this review. |
| Lean project entry points | Both projects use the package/module root `RequestProject` and mathlib `v4.28.0`. The matching project's `Main.lean` only imports Mathlib and sets options; checking it alone checks none of the matching development. | Import `TreeCodeBridge` explicitly and choose a single integration package or unambiguous module layout. Preserve the public graph-count definitions. |

## What was actually checked

**Exact Python replay.** I loaded the frozen verifier modules without invoking their writing `main()` functions, set their repository root to the real checkout, and supplied the frozen first verifier under the import name used by the second. The imported polynomial, independent-mask, and matching-bag helper modules and archived inputs came from the current checkout; those dependencies are not all in the frozen source directory.

I ran `run(14)` from the first frozen verifier. It reproduced the 5,446 small-tree checks, 594 small `b₁ = 0` trees, 5,438 eligible general-shadow checks, 105 disconnected checks, 109 archive replays, and all 13 parameter records. After normal JSON serialization (which turns integer dictionary keys into strings) and excluding only `elapsed_seconds`, the entire result equals the frozen JSON.

I ran `run(1000, 3)` from the second frozen verifier. The complete result, excluding only `elapsed_seconds`, equals its frozen JSON: all 1,842 certificate rows, both canonical-code sets, all 2,692 small decoration checks, 1,000 deterministic order-33 decorations, the extremizer, the sharp ratio `52513/217404`, and the 17 failures of the stronger discarded joint bound. The programs' exact assertions passed. These are Python checks, not Lean proofs of the universal reductions. I did not rerun `test_all.py` or the separate Section 8/9 sign/probe programs, and do not promote the note's stated test results for them into new verification claims.

**Lean replay.** `lake env lean RequestProject/Main.lean` succeeded in the fivefold project's build directory. Importing its compiled `Main` and printing axioms for `FivefoldForest.fivefold_of_forest` and `FivefoldForest.extendableCount_eq` reported only `propext`, `Classical.choice`, and `Quot.sound`. This re-elaborated `Main` against the available compiled dependencies; it was not a fresh recompilation of every fivefold module.

For the matching reserve I freshly compiled, in dependency order, `Codes`, `MatchingCover`, `KonigHall`, `ForestLemmas`, `TreeMatching`, `PosetCode`, `BagPoset`, `CodeInvariance`, `Antichain`, `PascalSmoothing`, `PascalBridge`, and `TreeCodeBridge` into an isolated temporary build cache against the checkout's mathlib. All twelve succeeded. Printing axioms for `MatchingBag.TreeMatching.erasure_depth_three_reserve`, `MatchingBag.TreeMatching.treeCode_eq_codeRelabel`, and `MatchingBag.codeProj_idealCode` reported only the same three standard axioms.

A scan of all Lean sources in the two requested project directories found no `sorry` or `admit`. `MB/ExhaustiveSmallCodes.lean:61` uses `native_decide`, with its additional trust dependency documented in that file. It is outside the reserve's import closure; the printed reserve axiom set contains no `Lean.ofReduceBool`. A later finite-certificate formalization must state its own computational trust boundary rather than inherit a general claim that everything in the project has the same one.

## A viable dependency order

1. **Unify the graph-count interface.** Reuse the fivefold definitions and partition/positivity lemmas. Choose a maximum matching; prove `D.G.indepNum = card D.Bag` and identify `D.maxIndepSets` with Mathlib maximum independent sets.
2. **Prove the partial-set/code bijection.** Map an independent vertex set to its occupied bags and selected endpoints; prove occupied-bag cardinality equals vertex-set cardinality. Prove extension to a maximum set is equivalent to membership in the corresponding projected code. Count projected words once, independently of the number of completions. This gives the guarded actual-graph erasure equality.
3. **Export the actual-graph Pascal reserve.** Apply the checked code reserve through that equality. Keep its cross-multiplied natural-number form until subtraction or division is necessary.
4. **Close the low-density theorem first.** Transfer the ideal projection characterization to a blocked vertex certificate of size at most two. Prove the at-most-six extension bound using the three empty bags, prove the blocked deletion incidence and the extendable lower extension incidence, and combine them with the reserve and positivity. This target does not need the 1,842-case computation or the `b₁ = 0` classification.
5. **Formalize the `b₁ = 0` structural theorem.** Define allowed, forced, and forbidden vertices using maximum sets. Establish the pendant matching, extension of all allowed independent sets, forced-neighbour bound, and `r ≤ δ−1`. Derive the blocked coefficient bounds and coefficient-ratio cone lemma. Prove the finite parameter split and the twelve small cases by exact arithmetic.
6. **Formalize the four-forbidden reduction.** Prove the tight core, exactly one attachment per flexible component, corona representation, core domination, concentration, and the leaf-root-to-base-root injection. Separately derive `11b₂ ≤ b₃` for the original tree.
7. **Certify the finite coverage and evaluation.** Implement finite rooted-forest codes of total base order ten with a coverage proof. Prove the corona subset formula used by `SRC/audit_four_forbidden_concentration_20260905.py:106` computes the corresponding `P` and `Q`, then establish the concentrated graph-count formula. Check `217404b₄ ≤ 52513e₄` on every encoded representative. Checking arithmetic on the stored rows alone leaves both the evaluator and coverage bridges open.
8. **Assemble the window theorem on actual graph counts.** Split `r=0`, `1≤r≤3`, and `r=4`; use the appropriate inequalities and `e₂e₄>0`. The public statement should quantify over the graph and stated window hypotheses, with no assumed structural decomposition, representative-coverage property, or unproved count identity.

## Largest remaining formalization obligations

The largest shared obligation is the actual-vertex-set to partial-word bijection and its rank/maximum-set interfaces. Without it, the existing reserve cannot yet be cited as a Lean theorem about these graph counts, and the Section 9 blocked-certificate argument has not reached actual vertex subsets.

The largest obligation specific to the window theorem is the universal graph-to-corona-and-rooted-forest reduction followed by a proved complete enumeration and correct evaluator. The exact 1,842-row replay supplies useful finite data; it does not discharge those universal obligations. The original-tree lower bound `b₃ ≥ 11b₂` must also remain a separate proof. Closing either public target while assuming any of these bridges would verify a conditional surrogate, not either requested graph-count theorem.
