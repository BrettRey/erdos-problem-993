import RequestProject.DepthThreeSpec
import RequestProject.Fivefold.Basic

/-! Identification of the supplied finite-set helpers with the frozen graph
counts. This does not yet identify either count with the matching-poset code.
The proofs are adapted from the independently replayed fivefold project's Main.
-/
namespace DepthThree.FivefoldBridge

open Finset FivefoldForest
open scoped Classical

variable {V : Type*} [Fintype V] (G : SimpleGraph V)

lemma alpha_univ : alpha G Finset.univ = G.indepNum := by
  refine le_antisymm ?_ ?_
  · obtain ⟨M, hM, hcard⟩ := exists_max_indFam G Finset.univ
    rw [← hcard]
    exact SimpleGraph.IsIndepSet.card_le_indepNum (mem_indFam.1 hM).2
  · obtain ⟨s, hs⟩ := G.exists_isNIndepSet_indepNum
    have hmem : s ∈ indFam G Finset.univ :=
      mem_indFam.2 ⟨Finset.subset_univ _, hs.isIndepSet⟩
    have := card_le_alpha hmem
    rw [hs.card_eq] at this
    exact this

lemma isMaximumIndepSet_iff {M : Finset V} :
    G.IsMaximumIndepSet M ↔ G.IsIndepSet (M : Set V) ∧ M.card = G.indepNum := by
  constructor
  · intro h
    exact ⟨h.isIndepSet, G.maximumIndepSet_card_eq_indepNum M h⟩
  · rintro ⟨hind, hcard⟩
    refine ⟨hind, fun t ht => ?_⟩
    rw [hcard]
    exact SimpleGraph.IsIndepSet.card_le_indepNum ht

lemma extendable_iff {A : Finset V} :
    DepthThree.Extendable G A ↔ Ext G Finset.univ A := by
  constructor
  · rintro ⟨M, hM, hAM⟩
    rw [isMaximumIndepSet_iff] at hM
    exact ⟨M, mem_indFam.2 ⟨Finset.subset_univ _, hM.1⟩,
      by rw [alpha_univ]; exact hM.2, hAM⟩
  · rintro ⟨M, hM, hcard, hAM⟩
    refine ⟨M, ?_, hAM⟩
    rw [isMaximumIndepSet_iff]
    exact ⟨(mem_indFam.1 hM).2, by rw [← alpha_univ]; exact hcard⟩

lemma e_eq (d : ℕ) : DepthThree.e G d = eD G Finset.univ d := by
  unfold DepthThree.e eD
  rw [alpha_univ]
  split
  · rw [ecnt]
    congr 1
    ext A
    simp only [Finset.mem_filter, mem_indFam, Finset.mem_univ, true_and,
      Finset.subset_univ, extendable_iff]
    tauto
  · rfl

lemma b_eq (d : ℕ) : DepthThree.b G d = bD G Finset.univ d := by
  unfold DepthThree.b bD
  rw [alpha_univ]
  split
  · rw [bcnt]
    congr 1
    ext A
    simp only [Finset.mem_filter, mem_indFam, Finset.mem_univ, true_and,
      Finset.subset_univ, extendable_iff]
    tauto
  · rfl

#print axioms e_eq
#print axioms b_eq

end DepthThree.FivefoldBridge
