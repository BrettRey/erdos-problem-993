import RequestProject.TreeCodeBridge
import RequestProject.DepthThreeSpec

/-!
# Bridge between the matching–bag library and the frozen graph specification

This module connects `MatchingBag.TreeMatching` (a finite forest together with a maximum
matching, oriented by a proper 2-colouring) with the frozen specification in
`RequestProject/DepthThreeSpec.lean`.

Main contents:

* `MatchingBag.exists_maximum_matching`: a finite graph has a maximum matching.
* `MatchingBag.exists_treeMatching_of_isAcyclic`: every finite forest underlies a
  `TreeMatching`.
* `MatchingBag.TreeMatching.bagSet` / `bagOf`: the vertex set of a bag, and the bag of a
  vertex; bags partition the vertex set and have at most two elements.
* `MatchingBag.TreeMatching.indepNum_eq_card_Bag`: `α(G) = #bags`.
* `MatchingBag.TreeMatching.isMaximumIndepSet_iff_mem_maxIndepSets`: `D.maxIndepSets` is
  exactly the set of maximum independent sets of `D.G` in Mathlib's sense.
-/

open Finset SimpleGraph

namespace MatchingBag

attribute [local instance 1] Classical.propDecidable

variable {V : Type*} [Fintype V] [DecidableEq V]

/-- Every finite graph on a finite vertex type has a maximum matching. -/
theorem exists_maximum_matching (G : SimpleGraph V) :
    ∃ M₀ : Finset (V × V), IsMatchingSet G M₀ ∧
      ∀ S : Finset (V × V), IsMatchingSet G S → S.card ≤ M₀.card := by
  classical
  set T := (Finset.univ : Finset (Finset (V × V))).filter (fun S => IsMatchingSet G S) with hT
  have hne : T.Nonempty := ⟨∅, by simp [hT, IsMatchingSet]⟩
  obtain ⟨M₀, hM₀, hmax⟩ := T.exists_max_image Finset.card hne
  rw [hT, Finset.mem_filter] at hM₀
  refine ⟨M₀, hM₀.2, fun S hS => hmax S ?_⟩
  rw [hT, Finset.mem_filter]
  exact ⟨Finset.mem_univ _, hS⟩

/-- Every finite forest underlies a `TreeMatching`. -/
theorem exists_treeMatching_of_isAcyclic (G : SimpleGraph V) (hG : G.IsAcyclic) :
    ∃ D : TreeMatching V, D.G = G := by
  obtain ⟨M₀, hM₀, hmax⟩ := exists_maximum_matching G
  obtain ⟨D, hD, -⟩ := exists_treeMatching G hG M₀ hM₀ hmax
  exact ⟨D, hD⟩

omit [Fintype V] [DecidableEq V] in
/-- Independence of a finset, unfolded. -/
lemma isIndepSet_coe_iff {G : SimpleGraph V} {S : Finset V} :
    G.IsIndepSet (S : Set V) ↔ ∀ u ∈ S, ∀ v ∈ S, ¬ G.Adj u v := by
  constructor
  · intro h u hu v hv huv
    exact h (by exact_mod_cast hu) (by exact_mod_cast hv) huv.ne huv
  · intro h u hu v hv _ huv
    exact h u (by exact_mod_cast hu) v (by exact_mod_cast hv) huv

namespace TreeMatching

variable (D : TreeMatching V)

/-! ### Bags as vertex sets -/

/-- The set of vertices in a bag: the two endpoints of a matched edge, or a single
unmatched vertex. -/
noncomputable def bagSet : D.Bag → Finset V :=
  Sum.elim (fun i => {D.Lv i, D.Rv i}) (fun v => {(v : V)})

variable {D}

@[simp] lemma mem_bagSet_inl {i : D.Idx} {v : V} :
    v ∈ D.bagSet (Sum.inl i) ↔ v = D.Lv i ∨ v = D.Rv i := by
  simp [bagSet]

@[simp] lemma mem_bagSet_inr {w : D.Unm} {v : V} :
    v ∈ D.bagSet (Sum.inr w) ↔ v = (w : V) := by
  simp [bagSet]

lemma card_bagSet_le (b : D.Bag) : (D.bagSet b).card ≤ 2 := by
  cases b with
  | inl i => exact (Finset.card_insert_le _ _).trans (by simp)
  | inr w => simp [bagSet]

lemma notMem_bagSet_of_ne {b b' : D.Bag} (h : b ≠ b') {v : V} (hv : v ∈ D.bagSet b) :
    v ∉ D.bagSet b' := by
  cases b with
  | inl i =>
      cases b' with
      | inl j =>
          have hij : i ≠ j := fun hh => h (by rw [hh])
          obtain ⟨d1, d2, d3, d4⟩ := bags_disjoint hij
          rw [mem_bagSet_inl] at hv ⊢
          rcases hv with rfl | rfl <;> rintro (hh | hh)
          exacts [d1 hh, d2 hh, d3 hh, d4 hh]
      | inr w =>
          rw [mem_bagSet_inl] at hv
          rw [mem_bagSet_inr]
          rintro rfl
          have := mem_unmatchedVerts_iff.1 w.2 i
          rcases hv with hh | hh
          exacts [this.1 hh, this.2 hh]
  | inr w =>
      rw [mem_bagSet_inr] at hv
      subst hv
      cases b' with
      | inl j =>
          rw [mem_bagSet_inl]
          have := mem_unmatchedVerts_iff.1 w.2 j
          rintro (hh | hh)
          exacts [this.1 hh, this.2 hh]
      | inr w' =>
          rw [mem_bagSet_inr]
          intro hh
          exact h (by rw [Subtype.ext hh])

variable (D)

/-- The bag containing a given vertex. -/
noncomputable def bagOf (v : V) : D.Bag :=
  if h : v ∈ D.unmatchedVerts then Sum.inr ⟨v, h⟩
  else Sum.inl (Classical.choose (mem_matchedVerts_iff.1
    (by simpa [unmatchedVerts, Finset.mem_sdiff] using h)))

variable {D}

lemma mem_bagSet_bagOf (v : V) : v ∈ D.bagSet (D.bagOf v) := by
  unfold bagOf
  split
  · next h => simp
  · next h =>
      have hs := Classical.choose_spec (mem_matchedVerts_iff.1
        (by simpa [unmatchedVerts, Finset.mem_sdiff] using h))
      rw [mem_bagSet_inl]
      exact hs

lemma bagOf_eq_iff {v : V} {b : D.Bag} : D.bagOf v = b ↔ v ∈ D.bagSet b := by
  constructor
  · rintro rfl; exact mem_bagSet_bagOf v
  · intro hv
    by_contra hne
    exact notMem_bagSet_of_ne hne (mem_bagSet_bagOf v) hv

lemma bagOf_eq_of_mem {v : V} {b : D.Bag} (hv : v ∈ D.bagSet b) : D.bagOf v = b :=
  bagOf_eq_iff.2 hv

/-- Two distinct vertices of an independent set lie in distinct bags. -/
lemma bagOf_injOn_indep {S : Finset V} (hS : ∀ u ∈ S, ∀ v ∈ S, ¬ D.G.Adj u v) :
    Set.InjOn D.bagOf (S : Set V) := by
  intro u hu v hv huv
  have hu' : u ∈ S := by exact_mod_cast hu
  have hv' : v ∈ S := by exact_mod_cast hv
  have hu2 : u ∈ D.bagSet (D.bagOf v) := by rw [← huv]; exact mem_bagSet_bagOf u
  have hv2 : v ∈ D.bagSet (D.bagOf v) := mem_bagSet_bagOf v
  cases hb : D.bagOf v with
  | inl i =>
      rw [hb, mem_bagSet_inl] at hu2 hv2
      rcases hu2 with rfl | rfl <;> rcases hv2 with h2 | h2
      · rw [h2]
      · exact absurd (h2 ▸ adj_Lv_Rv i) (hS _ hu' _ hv')
      · exact absurd (h2 ▸ (adj_Lv_Rv i).symm) (hS _ hu' _ hv')
      · rw [h2]
  | inr w =>
      rw [hb, mem_bagSet_inr] at hu2 hv2
      rw [hu2, hv2]

lemma card_image_bagOf {S : Finset V} (hS : ∀ u ∈ S, ∀ v ∈ S, ¬ D.G.Adj u v) :
    (S.image D.bagOf).card = S.card :=
  Finset.card_image_of_injOn (bagOf_injOn_indep hS)

/-! ### The independence number is the number of bags -/

lemma maxIndepSets_indep {I : Finset V} (hI : I ∈ D.maxIndepSets) :
    ∀ u ∈ I, ∀ v ∈ I, ¬ D.G.Adj u v := by
  rw [maxIndepSets, Finset.mem_filter] at hI; exact hI.2.1

lemma maxIndepSets_card {I : Finset V} (hI : I ∈ D.maxIndepSets) :
    I.card = Fintype.card V - D.M.card := by
  rw [maxIndepSets, Finset.mem_filter] at hI; exact hI.2.2

lemma maxIndepSets_nonempty (D : TreeMatching V) : D.maxIndepSets.Nonempty := by
  obtain ⟨C, hC⟩ := exists_minCover D
  refine ⟨Cᶜ, ?_⟩
  rw [maxIndepSets_eq_image_compl]
  exact Finset.mem_image_of_mem _ hC

/-- Any independent set has at most `#bags` elements. -/
lemma card_le_card_Bag {S : Finset V} (hS : ∀ u ∈ S, ∀ v ∈ S, ¬ D.G.Adj u v) :
    S.card ≤ Fintype.card D.Bag := by
  have hcov : IsVertexCover D.G Sᶜ := by
    intro u v huv
    by_contra hc
    push_neg at hc
    simp only [Finset.mem_compl, not_not] at hc
    exact hS u hc.1 v hc.2 huv
  have h := card_matching_le_cover D.isMatching hcov
  rw [Finset.card_compl] at h
  have h2 : S.card ≤ Fintype.card V := by
    simpa [Finset.card_univ] using Finset.card_le_card (Finset.subset_univ S)
  rw [card_Bag]
  omega

/-- **The independence number equals the number of bags.** -/
theorem indepNum_eq_card_Bag (D : TreeMatching V) : D.G.indepNum = Fintype.card D.Bag := by
  obtain ⟨I, hI⟩ := maxIndepSets_nonempty D
  have hcard : I.card = Fintype.card D.Bag := by
    rw [maxIndepSets_card hI, card_Bag]
  have hle : Fintype.card D.Bag ≤ D.G.indepNum := by
    rw [← hcard]
    exact SimpleGraph.IsIndepSet.card_le_indepNum
      (isIndepSet_coe_iff.2 (maxIndepSets_indep hI))
  obtain ⟨J, hJ⟩ := SimpleGraph.maximumIndepSet_exists (G := D.G)
  have hJcard := SimpleGraph.maximumIndepSet_card_eq_indepNum J hJ
  have hJle : J.card ≤ Fintype.card D.Bag :=
    card_le_card_Bag (isIndepSet_coe_iff.1 hJ.isIndepSet)
  omega

/-- `D.maxIndepSets` is exactly the collection of maximum independent sets of `D.G`. -/
theorem isMaximumIndepSet_iff_mem_maxIndepSets {I : Finset V} :
    D.G.IsMaximumIndepSet I ↔ I ∈ D.maxIndepSets := by
  constructor
  · intro hI
    have hcard := SimpleGraph.maximumIndepSet_card_eq_indepNum I hI
    rw [maxIndepSets, Finset.mem_filter]
    refine ⟨Finset.mem_univ _, isIndepSet_coe_iff.1 hI.isIndepSet, ?_⟩
    rw [hcard, indepNum_eq_card_Bag, card_Bag]
  · intro hI
    refine ⟨isIndepSet_coe_iff.2 (maxIndepSets_indep hI), fun t ht => ?_⟩
    have h1 : t.card ≤ Fintype.card D.Bag := card_le_card_Bag (isIndepSet_coe_iff.1 ht)
    have h2 : I.card = Fintype.card D.Bag := by rw [maxIndepSets_card hI, card_Bag]
    omega

/-- Extendability, expressed with `D.maxIndepSets`. -/
theorem extendable_iff {S : Finset V} :
    DepthThree.Extendable D.G S ↔ ∃ I ∈ D.maxIndepSets, S ⊆ I := by
  constructor
  · rintro ⟨I, hI, hSI⟩
    exact ⟨I, isMaximumIndepSet_iff_mem_maxIndepSets.1 hI, hSI⟩
  · rintro ⟨I, hI, hSI⟩
    exact ⟨I, isMaximumIndepSet_iff_mem_maxIndepSets.2 hI, hSI⟩

/-! ### The chosen vertex of a bag -/

variable (D)

/-- The vertex of a bag selected by a Boolean value: for a matched bag, the `L`-endpoint
if the bit is `true` and the `R`-endpoint otherwise; for a singleton bag, its vertex. -/
noncomputable def pick : D.Bag → Bool → V :=
  Sum.elim (fun i t => if t then D.Lv i else D.Rv i) (fun v _ => (v : V))

variable {D}

lemma pick_mem_bagSet (b : D.Bag) (t : Bool) : D.pick b t ∈ D.bagSet b := by
  cases b with
  | inl i => cases t <;> simp [pick]
  | inr w => simp [pick]

@[simp] lemma bagOf_pick (b : D.Bag) (t : Bool) : D.bagOf (D.pick b t) = b :=
  bagOf_eq_of_mem (pick_mem_bagSet b t)

/-- A maximum independent set meets every bag: the vertex it contains is the one selected by
its code word. -/
lemma pick_coverWord_mem {I : Finset V} (hI : I ∈ D.maxIndepSets) (b : D.Bag) :
    D.pick b (D.coverWord I b) ∈ I := by
  rw [maxIndepSets_eq_image_compl] at hI
  obtain ⟨C, hC, rfl⟩ := Finset.mem_image.1 hI
  obtain ⟨hone, hcov, hunm⟩ := minCover_structure hC
  cases b with
  | inl i =>
      by_cases h : D.Lv i ∈ Cᶜ
      · have hw : D.coverWord Cᶜ (Sum.inl i) = true := by
          simp only [coverWord, Sum.elim_inl, decide_eq_true_eq]; exact h
        rw [hw]
        simpa [pick] using h
      · have hw : D.coverWord Cᶜ (Sum.inl i) = false := by
          simp only [coverWord, Sum.elim_inl, decide_eq_false_iff_not]; exact h
        rw [hw]
        have hLC : D.Lv i ∈ C := by simpa using h
        have hRC : D.Rv i ∉ C := fun hh => hone i ⟨hLC, hh⟩
        simpa [pick] using hRC
  | inr w =>
      simpa [pick] using hunm _ w.2

/-- The unique element of a maximum independent set inside a given bag. -/
lemma eq_pick_coverWord {I : Finset V} (hI : I ∈ D.maxIndepSets) {u : V} (hu : u ∈ I) :
    u = D.pick (D.bagOf u) (D.coverWord I (D.bagOf u)) := by
  have h1 : D.pick (D.bagOf u) (D.coverWord I (D.bagOf u)) ∈ I := pick_coverWord_mem hI _
  have h2 : D.bagOf (D.pick (D.bagOf u) (D.coverWord I (D.bagOf u))) = D.bagOf u := bagOf_pick _ _
  exact (bagOf_injOn_indep (maxIndepSets_indep hI) (by exact_mod_cast hu)
    (by exact_mod_cast h1) h2.symm)

end TreeMatching

end MatchingBag
