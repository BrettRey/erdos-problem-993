import RequestProject.B1Zero.Pendant

/-!
# `b₁ = 0` makes the allowed subgraph well covered

Assuming `b₁ = 0`, every independent set of allowed vertices extends to a maximum
independent set.  Consequently an independent set is extendable **iff** it avoids the
forbidden vertices.
-/

open Finset SimpleGraph

namespace MatchingBag

attribute [local instance 1] Classical.propDecidable

variable {V : Type*} [Fintype V] [DecidableEq V]

namespace TreeMatching

variable {D : TreeMatching V}

lemma mem_maxIndepSets_of_card {I : Finset V} (hind : ∀ u ∈ I, ∀ v ∈ I, ¬ D.G.Adj u v)
    (hcard : I.card = Fintype.card D.Bag) : I ∈ D.maxIndepSets := by
  rw [maxIndepSets, Finset.mem_filter]
  exact ⟨Finset.mem_univ _, hind, by rw [hcard, card_Bag]⟩

/-- **Well-coveredness.**  Assuming `b₁ = 0`, every independent set of allowed vertices is
contained in a maximum independent set. -/
theorem exists_maxIndepSet_superset (hb1 : D.B1Zero) {S : Finset V}
    (hSA : S ⊆ D.AllowedV) (hSind : ∀ u ∈ S, ∀ v ∈ S, ¬ D.G.Adj u v) :
    ∃ I ∈ D.maxIndepSets, S ⊆ I := by
  classical
  -- for each bag, choose a representative
  have hchoice : ∀ b : D.Bag, ∃ v, v ∈ D.bagSet b ∧ v ∈ D.AllowedV ∧
      (∀ u ∈ S, u ∈ D.bagSet b → u = v) ∧
      (v ∉ S → ∀ w, D.G.Adj v w → w ∈ D.AllowedV → w ∈ D.bagSet b) := by
    intro b
    by_cases hS : ∃ u ∈ S, u ∈ D.bagSet b
    · obtain ⟨u, huS, hub⟩ := hS
      refine ⟨u, hub, hSA huS, ?_, ?_⟩
      · intro u' hu'S hu'b
        by_contra hne
        exact hSind u' hu'S u huS (bagSet_adj hu'b hub hne)
      · intro h; exact absurd huS h
    · push_neg at hS
      by_cases hF : ∃ f ∈ D.bagSet b, f ∈ D.ForcedV
      · obtain ⟨f, hfb, hfF⟩ := hF
        refine ⟨f, hfb, Forced_subset_Allowed hfF, ?_, ?_⟩
        · intro u huS hub; exact absurd hub (hS u huS)
        · intro _hnS w hadj hwA
          exact absurd (forbidden_of_adj_forced hfF hadj)
            (Finset.disjoint_left.1 Allowed_disjoint_Forbidden hwA)
      · push_neg at hF
        obtain ⟨-, hflex⟩ := bagSet_flex_of_no_forced hF
        obtain ⟨x, hxb, hxleaf⟩ := exists_leaf_in_flex_bag hb1 hflex
        exact ⟨x, hxb, (mem_Flex.1 (hflex x hxb)).1,
          fun u huS hub => absurd hub (hS u huS), fun _ => hxleaf⟩
  choose sel hselb hselA hselS hselLeaf using hchoice
  refine ⟨Finset.univ.image sel, ?_, ?_⟩
  · refine mem_maxIndepSets_of_card ?_ ?_
    · intro u hu v hv hadj
      obtain ⟨b1, -, rfl⟩ := Finset.mem_image.1 hu
      obtain ⟨b2, -, rfl⟩ := Finset.mem_image.1 hv
      have hb1b : D.bagOf (sel b1) = b1 := bagOf_eq_of_mem (hselb b1)
      have hb2b : D.bagOf (sel b2) = b2 := bagOf_eq_of_mem (hselb b2)
      have hne : b1 ≠ b2 := by
        rintro rfl
        exact hadj.ne rfl
      by_cases h2 : sel b2 ∈ S
      · by_cases h1 : sel b1 ∈ S
        · exact hSind _ h1 _ h2 hadj
        · have h3 := hselLeaf b1 h1 (sel b2) hadj (hselA b2)
          exact hne ((bagOf_eq_of_mem h3).symm.trans hb2b)
      · have h3 := hselLeaf b2 h2 (sel b1) hadj.symm (hselA b1)
        exact hne (hb1b.symm.trans (bagOf_eq_of_mem h3))
    · rw [Finset.card_image_of_injective _ ?_, Finset.card_univ]
      intro b1 b2 h
      rw [← bagOf_eq_of_mem (hselb b1), ← bagOf_eq_of_mem (hselb b2), h]
  · intro u huS
    refine Finset.mem_image.2 ⟨D.bagOf u, Finset.mem_univ _, ?_⟩
    exact (hselS (D.bagOf u) u huS (mem_bagSet_bagOf u)).symm

/-- Assuming `b₁ = 0`, an independent set is extendable exactly when it avoids the forbidden
vertices. -/
theorem extendable_iff_subset_allowed (hb1 : D.B1Zero) {S : Finset V}
    (hSind : ∀ u ∈ S, ∀ v ∈ S, ¬ D.G.Adj u v) :
    (∃ I ∈ D.maxIndepSets, S ⊆ I) ↔ S ⊆ D.AllowedV := by
  constructor
  · rintro ⟨I, hI, hSI⟩ v hv
    exact mem_Allowed.2 ⟨I, hI, hSI hv⟩
  · intro hSA
    exact exists_maxIndepSet_superset hb1 hSA hSind

end TreeMatching

end MatchingBag
