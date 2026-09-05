import RequestProject.Induction

/-!
# Fivefold bounds for blocked sets in finite forests

This file states the requested theorem exactly in the form of the problem statement,
using the definitions supplied there (`Extendable`, `extendableCount`, `blockedCount`,
`FivefoldConclusion`), and derives it from the induction carried out in
`RequestProject.Induction`.
-/

open Finset

namespace FivefoldForest

open scoped Classical

variable {V : Type*} [Fintype V]

/-- Containment in a maximum independent set, not merely a maximal one. -/
def Extendable (G : SimpleGraph V) (S : Finset V) : Prop :=
  ∃ M : Finset V, G.IsMaximumIndepSet M ∧ S ⊆ M

/-- The number of independent sets of size `α - d` contained in a maximum independent set. -/
noncomputable def extendableCount (G : SimpleGraph V) (d : ℕ) : ℕ := by
  classical
  exact if d ≤ G.indepNum then
    (Finset.univ.filter fun S : Finset V =>
      S.card = G.indepNum - d ∧ G.IsIndepSet (S : Set V) ∧ Extendable G S).card
  else 0

/-- The number of independent sets of size `α - d` *not* contained in a maximum
independent set. -/
noncomputable def blockedCount (G : SimpleGraph V) (d : ℕ) : ℕ := by
  classical
  exact if d ≤ G.indepNum then
    (Finset.univ.filter fun S : Finset V =>
      S.card = G.indepNum - d ∧ G.IsIndepSet (S : Set V) ∧ ¬ Extendable G S).card
  else 0

/-- Specification only; this definition is not a proof of the assertion. -/
def FivefoldConclusion (G : SimpleGraph V) : Prop :=
  5 * blockedCount G 3 ≤ extendableCount G 3 ∧
  5 * blockedCount G 4 ≤ extendableCount G 4

/-! ### Bridge to the internal definitions -/

section Bridge

variable (G : SimpleGraph V)

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

lemma extendable_iff {A : Finset V} : Extendable G A ↔ Ext G Finset.univ A := by
  constructor
  · rintro ⟨M, hM, hAM⟩
    rw [isMaximumIndepSet_iff] at hM
    exact ⟨M, mem_indFam.2 ⟨Finset.subset_univ _, hM.1⟩, by rw [alpha_univ]; exact hM.2, hAM⟩
  · rintro ⟨M, hM, hcard, hAM⟩
    refine ⟨M, ?_, hAM⟩
    rw [isMaximumIndepSet_iff]
    exact ⟨(mem_indFam.1 hM).2, by rw [← alpha_univ]; exact hcard⟩

lemma extendableCount_eq (d : ℕ) : extendableCount G d = eD G Finset.univ d := by
  unfold extendableCount eD
  rw [alpha_univ]
  split
  · rw [ecnt]
    congr 1
    ext A
    simp only [Finset.mem_filter, mem_indFam, Finset.mem_univ, true_and,
      Finset.subset_univ, extendable_iff]
    tauto
  · rfl

lemma blockedCount_eq (d : ℕ) : blockedCount G d = bD G Finset.univ d := by
  unfold blockedCount bD
  rw [alpha_univ]
  split
  · rw [bcnt]
    congr 1
    ext A
    simp only [Finset.mem_filter, mem_indFam, Finset.mem_univ, true_and,
      Finset.subset_univ, extendable_iff]
    tauto
  · rfl

end Bridge

/-- **Fivefold bounds for blocked sets in finite forests.**
Let `G` be a finite forest with independence number at least five.  If every independent
set of size `α - 2` is contained in a maximum independent set, then the number of blocked
independent sets at defects three and four is at most a fifth of the number of extendable
ones. -/
theorem fivefold_of_forest
    {V : Type*} [Fintype V] (G : SimpleGraph V)
    (hforest : G.IsAcyclic)
    (halpha : 5 ≤ G.indepNum)
    (hb2 : blockedCount G 2 = 0) :
    FivefoldConclusion G := by
  have hF := forestLike_of_isAcyclic hforest
  rw [blockedCount_eq] at hb2
  have h5 : 5 ≤ alpha G Finset.univ := by rw [alpha_univ]; exact halpha
  obtain ⟨h3, h4⟩ := fivefold hF h5 hb2
  constructor
  · rw [blockedCount_eq, extendableCount_eq]; exact h3
  · rw [blockedCount_eq, extendableCount_eq]; exact h4

end FivefoldForest
