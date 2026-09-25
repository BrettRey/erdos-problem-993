import Mathlib.Combinatorics.SimpleGraph.Acyclic
import Mathlib.Combinatorics.SimpleGraph.Finite
import Mathlib.Data.Set.Card
import Mathlib.Tactic

/-!
# Edge counts in a finite forest

* `SimpleGraph.IsAcyclic.ncard_edgeSet_add_one_le`: a finite acyclic graph on a nonempty
  vertex type has at most `|V| - 1` edges.  (Extend the forest to a spanning tree of the
  complete graph and use `SimpleGraph.IsTree.card_edgeFinset`.)
* `SimpleGraph.IsAcyclic.card_adjPairs_le`: for two disjoint finsets `A`, `B` of vertices of
  a finite forest, the number of edges between `A` and `B` is at most `|A| + |B| - 1`.
-/

open Finset

namespace SimpleGraph

/-- A finite forest on a nonempty vertex type has at most `|V| - 1` edges. -/
theorem IsAcyclic.ncard_edgeSet_add_one_le {V : Type*} [Fintype V] [Nonempty V]
    {G : SimpleGraph V} (hG : G.IsAcyclic) :
    G.edgeSet.ncard + 1 ≤ Fintype.card V := by
  classical
  obtain ⟨F, hGF, hF⟩ := SimpleGraph.exists_maximal_isAcyclic_of_le_isAcyclic (G := ⊤) le_top hG
  have hreach := SimpleGraph.reachable_eq_of_maximal_isAcyclic F hF
  have hpre : F.Preconnected := by
    intro u v
    have h : (⊤ : SimpleGraph V).Reachable u v := by
      rcases eq_or_ne u v with rfl | h
      · exact Reachable.refl _
      · exact (SimpleGraph.top_adj (V := V) u v |>.2 h).reachable
    rw [hreach]; exact h
  have hFtree : F.IsTree := ⟨⟨hpre⟩, hF.prop.2⟩
  have hcard := hFtree.card_edgeFinset
  rw [SimpleGraph.edgeFinset_card] at hcard
  have hle : G.edgeSet.ncard ≤ F.edgeSet.ncard :=
    Set.ncard_le_ncard (SimpleGraph.edgeSet_mono hGF) (Set.toFinite _)
  have h2 : F.edgeSet.ncard = Fintype.card F.edgeSet := by
    rw [Set.ncard_eq_toFinset_card', Set.toFinset_card]
  omega

/-- **Forest incidence bound.**  In a finite forest, the number of edges between two disjoint
sets `A`, `B` of vertices is at most `|A| + |B| - 1`. -/
theorem IsAcyclic.card_adjPairs_le {V : Type*} [Fintype V] [DecidableEq V] {G : SimpleGraph V}
    [DecidableRel G.Adj] (hG : G.IsAcyclic) (A B : Finset V) (hAB : Disjoint A B)
    (hne : (A ∪ B).Nonempty) :
    ((A ×ˢ B).filter (fun p : V × V => G.Adj p.1 p.2)).card + 1 ≤ A.card + B.card := by
  classical
  set W := A ∪ B with hW
  set s : Set V := (W : Set V) with hs
  obtain ⟨w0, hw0⟩ := hne
  have hw0s : w0 ∈ s := by simpa [hs] using hw0
  haveI : Nonempty ↥s := ⟨⟨w0, hw0s⟩⟩
  have hb : (G.induce s).edgeSet.ncard + 1 ≤ W.card := by
    have h := (hG.induce s).ncard_edgeSet_add_one_le
    have hcards : Fintype.card ↥s = W.card := by simp [hs]
    omega
  set ι : V → ↥s := fun v => if h : v ∈ s then ⟨v, h⟩ else ⟨w0, hw0s⟩ with hι
  have hιval : ∀ v ∈ s, (ι v : V) = v := by intro v hv; simp [hι, hv]
  set P := (A ×ˢ B).filter (fun p : V × V => G.Adj p.1 p.2) with hP
  have hmem : ∀ p ∈ P, p.1 ∈ s ∧ p.2 ∈ s ∧ G.Adj p.1 p.2 ∧ p.1 ∈ A ∧ p.2 ∈ B := by
    intro p hp
    rw [hP, Finset.mem_filter, Finset.mem_product] at hp
    exact ⟨by simp [hs, hW, hp.1.1], by simp [hs, hW, hp.1.2], hp.2, hp.1.1, hp.1.2⟩
  have hkey : P.card ≤ (G.induce s).edgeSet.ncard := by
    rw [Set.ncard_eq_toFinset_card']
    refine Finset.card_le_card_of_injOn (fun p => s(ι p.1, ι p.2)) ?_ ?_
    · intro p hp
      obtain ⟨hp1, hp2, hadj, -, -⟩ := hmem p hp
      have hadj2 : (G.induce s).Adj (ι p.1) (ι p.2) := by
        show G.Adj _ _
        simp only [Function.Embedding.coe_subtype]
        rw [hιval _ hp1, hιval _ hp2]
        exact hadj
      simpa using hadj2
    · intro p hp q hq hpq
      rw [Finset.mem_coe] at hp hq
      obtain ⟨hp1, hp2, -, hpA, hpB⟩ := hmem p hp
      obtain ⟨hq1, hq2, -, hqA, hqB⟩ := hmem q hq
      simp only at hpq
      rw [Sym2.eq_iff] at hpq
      rcases hpq with ⟨h1, h2⟩ | ⟨h1, h2⟩
      · have e1 : p.1 = q.1 := by
          have := congrArg (Subtype.val) h1; rwa [hιval _ hp1, hιval _ hq1] at this
        have e2 : p.2 = q.2 := by
          have := congrArg (Subtype.val) h2; rwa [hιval _ hp2, hιval _ hq2] at this
        exact Prod.ext e1 e2
      · exfalso
        have e1 : p.1 = q.2 := by
          have := congrArg (Subtype.val) h1; rwa [hιval _ hp1, hιval _ hq2] at this
        exact Finset.disjoint_left.1 hAB hpA (e1 ▸ hqB)
  have hWcard : W.card = A.card + B.card := Finset.card_union_of_disjoint hAB
  omega

end SimpleGraph
