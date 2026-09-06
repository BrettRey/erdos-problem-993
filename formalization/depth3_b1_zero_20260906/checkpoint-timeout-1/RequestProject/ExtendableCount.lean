import RequestProject.GraphBridge

/-!
# The count of extendable independent sets is the erasure profile of the tree code

This module proves the cardinality-preserving bijection between extendable independent
vertex sets of a given size and projected partial words of the code of maximum independent
sets.  The consequence is `MatchingBag.TreeMatching.e_eq_erasure`:

`DepthThree.e D.G d = D.erasure d`  for `d ≤ #bags = α(D.G)`,

which is the interface obligation making `TreeCodeBridge.erasure_depth_three_reserve` a
statement about the graph counts of `DepthThreeSpec`.
-/

open Finset SimpleGraph

namespace MatchingBag

attribute [local instance 1] Classical.propDecidable

variable {V : Type*} [Fintype V] [DecidableEq V]

namespace TreeMatching

variable {D : TreeMatching V}

/-! ### The word of a partial set -/

/-- On the bags met by `S`, the code word of `S` agrees with the code word of any maximum
independent set containing `S`. -/
lemma coverWord_agree {S I : Finset V} (hI : I ∈ D.maxIndepSets) (hSI : S ⊆ I)
    {b : D.Bag} (hb : b ∈ S.image D.bagOf) :
    D.coverWord S b = D.coverWord I b := by
  obtain ⟨u, hu, hbu⟩ := Finset.mem_image.1 hb
  have hub : u ∈ D.bagSet b := bagOf_eq_iff.1 hbu
  have huI : u ∈ I := hSI hu
  cases b with
  | inl i =>
      rw [mem_bagSet_inl] at hub
      simp only [coverWord, Sum.elim_inl]
      rcases hub with rfl | rfl
      · simp [hu, huI]
      · have h1 : D.Lv i ∉ S := fun hh =>
          maxIndepSets_indep hI _ (hSI hh) _ huI (adj_Lv_Rv i)
        have h2 : D.Lv i ∉ I := fun hh =>
          maxIndepSets_indep hI _ hh _ huI (adj_Lv_Rv i)
        simp [h1, h2]
  | inr w =>
      rw [mem_bagSet_inr] at hub
      subst hub
      simp [coverWord, hu, huI]

variable (D)

/-- The vertex set chosen inside the bags of `K` by a word `f`. -/
noncomputable def pickSet (K : Finset D.Bag) (f : D.Bag → Bool) : Finset V :=
  K.image (fun b => D.pick b (f b))

variable {D}

@[simp] lemma image_bagOf_pickSet (K : Finset D.Bag) (f : D.Bag → Bool) :
    (D.pickSet K f).image D.bagOf = K := by
  rw [pickSet, Finset.image_image]
  have hid : (D.bagOf ∘ fun b => D.pick b (f b)) = id := by
    funext b; simp
  rw [hid, Finset.image_id]

lemma card_pickSet (K : Finset D.Bag) (f : D.Bag → Bool) :
    (D.pickSet K f).card = K.card := by
  refine Finset.card_image_of_injOn ?_
  intro b _ b' _ hbb
  have := congrArg D.bagOf hbb
  rwa [bagOf_pick, bagOf_pick] at this

lemma pickSet_subset_of_maxIndep {I : Finset V} (hI : I ∈ D.maxIndepSets) (K : Finset D.Bag) :
    D.pickSet K (D.coverWord I) ⊆ I := by
  intro v hv
  obtain ⟨b, -, rfl⟩ := Finset.mem_image.1 hv
  exact pick_coverWord_mem hI b

/-! ### The extendable sets of a given size -/

variable (D)

/-- The extendable independent sets of a given size. -/
noncomputable def extSets (k : ℕ) : Finset (Finset V) :=
  Finset.univ.filter (fun S => S.card = k ∧ (∀ u ∈ S, ∀ v ∈ S, ¬ D.G.Adj u v) ∧
    ∃ I ∈ D.maxIndepSets, S ⊆ I)

/-- The blocked independent sets of a given size. -/
noncomputable def blkSets (k : ℕ) : Finset (Finset V) :=
  Finset.univ.filter (fun S => S.card = k ∧ (∀ u ∈ S, ∀ v ∈ S, ¬ D.G.Adj u v) ∧
    ¬ ∃ I ∈ D.maxIndepSets, S ⊆ I)

variable {D}

@[simp] lemma mem_extSets {k : ℕ} {S : Finset V} :
    S ∈ D.extSets k ↔ S.card = k ∧ (∀ u ∈ S, ∀ v ∈ S, ¬ D.G.Adj u v) ∧
      ∃ I ∈ D.maxIndepSets, S ⊆ I := by
  simp [extSets]

@[simp] lemma mem_blkSets {k : ℕ} {S : Finset V} :
    S ∈ D.blkSets k ↔ S.card = k ∧ (∀ u ∈ S, ∀ v ∈ S, ¬ D.G.Adj u v) ∧
      ¬ ∃ I ∈ D.maxIndepSets, S ⊆ I := by
  simp [blkSets]

/-- **The fibrewise bijection**: the extendable independent sets whose bag support is `K`
are in bijection with the projections of the code of maximum independent sets onto `K`. -/
theorem card_extSets_fiber (K : Finset D.Bag) :
    ((D.extSets K.card).filter (fun S => S.image D.bagOf = K)).card
      = (codeProj K (D.maxIndepSets.image D.coverWord)).card := by
  refine Finset.card_bij (fun S _ => restrictTo K (D.coverWord S)) ?_ ?_ ?_
  · -- well defined
    intro S hS
    rw [Finset.mem_filter, mem_extSets] at hS
    obtain ⟨⟨-, -, I, hI, hSI⟩, hKS⟩ := hS
    refine Finset.mem_image.2 ⟨D.coverWord I, Finset.mem_image_of_mem _ hI, ?_⟩
    funext b
    by_cases hb : b ∈ K
    · have hb' : b ∈ S.image D.bagOf := by rw [hKS]; exact hb
      simp [restrictTo, hb, coverWord_agree hI hSI hb']
    · simp [restrictTo, hb]
  · -- injective
    intro S hS S' hS' hww
    rw [Finset.mem_filter, mem_extSets] at hS hS'
    obtain ⟨⟨-, hind, -⟩, hKS⟩ := hS
    obtain ⟨⟨-, hind', -⟩, hKS'⟩ := hS'
    have key : ∀ (A B : Finset V), A.image D.bagOf = K → B.image D.bagOf = K →
        (∀ u ∈ A, ∀ v ∈ A, ¬ D.G.Adj u v) →
        restrictTo K (D.coverWord A) = restrictTo K (D.coverWord B) → A ⊆ B := by
      intro A B hA hB hindA hw u hu
      have hbu : D.bagOf u ∈ K := by rw [← hA]; exact Finset.mem_image_of_mem _ hu
      have hwb := congrFun hw (D.bagOf u)
      simp only [restrictTo, if_pos hbu, Option.some.injEq] at hwb
      have hub : u ∈ D.bagSet (D.bagOf u) := mem_bagSet_bagOf u
      obtain ⟨u', hu', hbu'⟩ : ∃ u' ∈ B, D.bagOf u' = D.bagOf u := by
        have : D.bagOf u ∈ B.image D.bagOf := by rw [hB]; exact hbu
        obtain ⟨u', hu', h⟩ := Finset.mem_image.1 this
        exact ⟨u', hu', h⟩
      have hu'b : u' ∈ D.bagSet (D.bagOf u) := by
        rw [← hbu']; exact mem_bagSet_bagOf u'
      cases hbb : D.bagOf u with
      | inl i =>
          rw [hbb] at hwb hub hu'b
          simp only [coverWord, Sum.elim_inl, decide_eq_decide] at hwb
          rw [mem_bagSet_inl] at hub hu'b
          rcases hub with rfl | rfl
          · exact hwb.1 hu
          · have hLA : D.Lv i ∉ A := fun hh => hindA _ hh _ hu (adj_Lv_Rv i)
            have hLB : D.Lv i ∉ B := fun hh => hLA (hwb.2 hh)
            rcases hu'b with rfl | rfl
            · exact absurd hu' hLB
            · exact hu'
      | inr w =>
          rw [hbb] at hub hu'b
          rw [mem_bagSet_inr] at hub hu'b
          rw [hub, ← hu'b]
          exact hu'
    exact Finset.Subset.antisymm (key S S' hKS hKS' hind hww)
      (key S' S hKS' hKS hind' hww.symm)
  · -- surjective
    intro w hw
    obtain ⟨f, hf, rfl⟩ := Finset.mem_image.1 hw
    obtain ⟨I, hI, rfl⟩ := Finset.mem_image.1 hf
    refine ⟨D.pickSet K (D.coverWord I), ?_, ?_⟩
    · rw [Finset.mem_filter, mem_extSets]
      have hsub := pickSet_subset_of_maxIndep hI K
      exact ⟨⟨card_pickSet K _, fun u hu v hv =>
        maxIndepSets_indep hI u (hsub hu) v (hsub hv), I, hI, hsub⟩,
        image_bagOf_pickSet K _⟩
    · funext b
      by_cases hb : b ∈ K
      · have hb' : b ∈ (D.pickSet K (D.coverWord I)).image D.bagOf := by
          rw [image_bagOf_pickSet]; exact hb
        simp [restrictTo, hb,
          coverWord_agree hI (pickSet_subset_of_maxIndep hI K) hb']
      · simp [restrictTo, hb]

/-- **The number of extendable independent sets of size `k` is `p_k` of the tree code.** -/
theorem card_extSets (D : TreeMatching V) (k : ℕ) :
    (D.extSets k).card = codeP D.treeCode k := by
  rw [← codeP_maxIndepCode, codeP]
  have hmaps : Set.MapsTo (fun S : Finset V => S.image D.bagOf)
      (D.extSets k : Set (Finset V)) ((Finset.univ : Finset D.Bag).powersetCard k) := by
    intro S hS
    rw [Finset.mem_coe, mem_extSets] at hS
    rw [Finset.mem_coe, Finset.mem_powersetCard]
    exact ⟨Finset.subset_univ _, by rw [card_image_bagOf hS.2.1, hS.1]⟩
  rw [Finset.card_eq_sum_card_fiberwise hmaps]
  refine Finset.sum_congr rfl fun K hK => ?_
  rw [Finset.mem_powersetCard] at hK
  rw [← hK.2]
  exact card_extSets_fiber K

/-! ### Identification with the specification counts -/

/-- The extendable count of the specification is the number of extendable independent sets
of size `α - d`. -/
theorem e_eq_card_extSets (D : TreeMatching V) {d : ℕ} (hd : d ≤ Fintype.card D.Bag) :
    DepthThree.e D.G d = (D.extSets (Fintype.card D.Bag - d)).card := by
  have hα : D.G.indepNum = Fintype.card D.Bag := indepNum_eq_card_Bag D
  rw [DepthThree.e, if_pos (by rw [hα]; exact hd)]
  congr 1
  ext S
  simp only [Finset.mem_filter, Finset.mem_univ, true_and, mem_extSets, hα]
  constructor
  · rintro ⟨h1, h2, h3⟩
    exact ⟨h1, isIndepSet_coe_iff.1 h2, extendable_iff.1 h3⟩
  · rintro ⟨h1, h2, h3⟩
    exact ⟨h1, isIndepSet_coe_iff.2 h2, extendable_iff.2 h3⟩

/-- The blocked count of the specification is the number of blocked independent sets of
size `α - d`. -/
theorem b_eq_card_blkSets (D : TreeMatching V) {d : ℕ} (hd : d ≤ Fintype.card D.Bag) :
    DepthThree.b D.G d = (D.blkSets (Fintype.card D.Bag - d)).card := by
  have hα : D.G.indepNum = Fintype.card D.Bag := indepNum_eq_card_Bag D
  rw [DepthThree.b, if_pos (by rw [hα]; exact hd)]
  congr 1
  ext S
  simp only [Finset.mem_filter, Finset.mem_univ, true_and, mem_blkSets, hα]
  constructor
  · rintro ⟨h1, h2, h3⟩
    exact ⟨h1, isIndepSet_coe_iff.1 h2, fun hc => h3 (extendable_iff.2 hc)⟩
  · rintro ⟨h1, h2, h3⟩
    exact ⟨h1, isIndepSet_coe_iff.2 h2, fun hc => h3 (extendable_iff.1 hc)⟩

/-- **The interface obligation**: the specification's extendable count is the erasure
profile of the tree code. -/
theorem e_eq_erasure (D : TreeMatching V) {d : ℕ} (hd : d ≤ Fintype.card D.Bag) :
    DepthThree.e D.G d = D.erasure d := by
  rw [e_eq_card_extSets D hd, card_extSets, erasure]

end TreeMatching

end MatchingBag
