import RequestProject.B1Zero.CoronaRepr

/-!
# Moving the avoided vertices to the base

If `Y` is a set of base vertices each of which is itself avoided, or whose pendant is
avoided, then avoiding `Y` instead of the original set can only increase the number of
independent sets.  The injection replaces each avoided base vertex of an independent set by
its pendant.
-/

open Finset SimpleGraph

namespace MatchingBag

namespace B1Zero

namespace CoronaData

attribute [local instance 1] Classical.propDecidable

variable {V : Type*} [Fintype V] [DecidableEq V] {G : SimpleGraph V} {C : CoronaData G}

/-- **Root moving.** -/
theorem cntRev_move_roots {X R Y : Finset V} (hX : X ⊆ C.base) (hY : Y ⊆ X)
    (hmove : ∀ y ∈ Y, y ∈ R ∨ C.pend y ∈ R) (n j : ℕ) :
    cntRev G (C.cor X \ R) n j ≤ cntRev G (C.cor X \ Y) n j := by
  classical
  set Y₂ := Y \ R with hY₂def
  set ψ : Finset V → Finset V := fun S => (S \ Y₂) ∪ (S ∩ Y₂).image C.pend with hψ
  have hY₂Y : Y₂ ⊆ Y := Finset.sdiff_subset
  have hY₂base : Y₂ ⊆ C.base := hY₂Y.trans (hY.trans hX)
  have hpendR : ∀ y ∈ Y₂, C.pend y ∈ R := by
    intro y hy
    rw [hY₂def, Finset.mem_sdiff] at hy
    rcases hmove y hy.1 with h | h
    · exact absurd h hy.2
    · exact h
  -- basic facts about a member of the source
  have hsrc : ∀ S ∈ revSets G (C.cor X \ R) n j,
      S ⊆ C.cor X ∧ IndepOn G S ∧ S.card + j = n ∧ (∀ x ∈ S, x ∉ R) := by
    intro S hS
    rw [mem_revSets] at hS
    exact ⟨fun x hx => (Finset.mem_sdiff.1 (hS.1 hx)).1, hS.2.1, hS.2.2,
      fun x hx => (Finset.mem_sdiff.1 (hS.1 hx)).2⟩
  have hpendnotS : ∀ S ∈ revSets G (C.cor X \ R) n j, ∀ y ∈ Y₂, C.pend y ∉ S := by
    intro S hS y hy hmem
    exact (hsrc S hS).2.2.2 _ hmem (hpendR y hy)
  have hmemY : ∀ S ∈ revSets G (C.cor X \ R) n j, ∀ y ∈ Y, y ∈ S → y ∈ Y₂ := by
    intro S hS y hyY hyS
    rw [hY₂def, Finset.mem_sdiff]
    exact ⟨hyY, (hsrc S hS).2.2.2 y hyS⟩
  have hcardψ : ∀ S ∈ revSets G (C.cor X \ R) n j, (ψ S).card = S.card := by
    intro S hS
    have hdisj : Disjoint (S \ Y₂) ((S ∩ Y₂).image C.pend) := by
      rw [Finset.disjoint_right]
      intro z hz hz'
      obtain ⟨y, hy, rfl⟩ := Finset.mem_image.1 hz
      exact hpendnotS S hS y (Finset.mem_inter.1 hy).2 (Finset.mem_sdiff.1 hz').1
    have himg : ((S ∩ Y₂).image C.pend).card = (S ∩ Y₂).card := by
      refine Finset.card_image_of_injOn ?_
      intro a ha b hb hab
      rw [Finset.mem_coe, Finset.mem_inter] at ha hb
      exact C.pend_inj a (hY₂base ha.2) b (hY₂base hb.2) hab
    have hsum : (S ∩ Y₂).card + (S \ Y₂).card = S.card :=
      Finset.card_inter_add_card_sdiff S Y₂
    rw [hψ]
    simp only
    rw [Finset.card_union_of_disjoint hdisj, himg]
    omega
  refine Finset.card_le_card_of_injOn ψ ?_ ?_
  · intro S hS
    obtain ⟨hsub, hind, hcard, hnotR⟩ := hsrc S hS
    rw [Finset.mem_coe] at hS
    rw [Finset.mem_coe, mem_revSets]
    refine ⟨?_, ?_, ?_⟩
    · intro x hx
      rw [hψ] at hx
      simp only [Finset.mem_union, Finset.mem_image, Finset.mem_sdiff, Finset.mem_inter] at hx
      rw [Finset.mem_sdiff]
      rcases hx with ⟨hxS, hxY₂⟩ | ⟨y, ⟨hyS, hyY₂⟩, rfl⟩
      · exact ⟨hsub hxS, fun hxY => hxY₂ (hmemY S hS x hxY hxS)⟩
      · exact ⟨pend_mem_cor (hY (hY₂Y hyY₂)), fun hmem => C.pend_notMem y (hY₂base hyY₂)
          ((hY.trans hX) hmem)⟩
    · intro a ha b hb
      rw [hψ] at ha hb
      simp only [Finset.mem_union, Finset.mem_image, Finset.mem_sdiff, Finset.mem_inter] at ha hb
      rcases ha with ⟨haS, haY₂⟩ | ⟨y, ⟨hyS, hyY₂⟩, rfl⟩
      · rcases hb with ⟨hbS, hbY₂⟩ | ⟨y', ⟨hy'S, hy'Y₂⟩, rfl⟩
        · exact hind a haS b hbS
        · intro hadj
          have := C.pend_nbr y' (hY₂base hy'Y₂) a (cor_mono (hX) (hsub haS)) hadj.symm
          exact haY₂ (this ▸ hy'Y₂)
      · rcases hb with ⟨hbS, hbY₂⟩ | ⟨y', ⟨hy'S, hy'Y₂⟩, rfl⟩
        · intro hadj
          have := C.pend_nbr y (hY₂base hyY₂) b (cor_mono (hX) (hsub hbS)) hadj
          exact hbY₂ (this ▸ hyY₂)
        · intro hadj
          have := C.pend_nbr y (hY₂base hyY₂) (C.pend y')
            (cor_mono hX (pend_mem_cor (hY (hY₂Y hy'Y₂)))) hadj
          exact C.pend_notMem y' (hY₂base hy'Y₂) (this ▸ hY₂base hyY₂)
    · rw [hcardψ S hS]
      exact hcard
  · -- injectivity, by an explicit recovery formula
    have recover : ∀ S ∈ revSets G (C.cor X \ R) n j,
        S = (ψ S \ Y₂.image C.pend) ∪ Y₂.filter (fun y => C.pend y ∈ ψ S) := by
      intro S hS
      obtain ⟨hsub, hind, hcard, hnotR⟩ := hsrc S hS
      ext x
      simp only [hψ, Finset.mem_union, Finset.mem_sdiff, Finset.mem_image, Finset.mem_inter,
        Finset.mem_filter]
      constructor
      · intro hx
        by_cases hxY₂ : x ∈ Y₂
        · exact Or.inr ⟨hxY₂, Or.inr ⟨x, ⟨hx, hxY₂⟩, rfl⟩⟩
        · refine Or.inl ⟨Or.inl ⟨hx, hxY₂⟩, ?_⟩
          rintro ⟨y, hy, rfl⟩
          exact hnotR _ hx (hpendR y hy)
      · rintro (⟨hx, hnx⟩ | ⟨hxY₂, hpx⟩)
        · rcases hx with ⟨hxS, -⟩ | ⟨y, ⟨hyS, hyY₂⟩, rfl⟩
          · exact hxS
          · exact absurd ⟨y, hyY₂, rfl⟩ hnx
        · rcases hpx with ⟨hpS, -⟩ | ⟨y, ⟨hyS, hyY₂⟩, hpy⟩
          · exact absurd hpS (hpendnotS S hS x hxY₂)
          · exact (C.pend_inj y (hY₂base hyY₂) x (hY₂base hxY₂) hpy) ▸ hyS
    intro S hS S' hS' heq
    rw [Finset.mem_coe] at hS hS'
    rw [recover S hS, recover S' hS', heq]

end CoronaData

end B1Zero

end MatchingBag
