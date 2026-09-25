import RequestProject.B1Zero.CoronaCount

/-!
# Representing a graph corona by a plane rooted forest

`exists_forest_repr`: if the base of a corona inside a forest is equipped with a transversal
`Y` (one vertex in each connected component), then there is a term `f` of the inductive type
`RootedForestCertificate.Forest` with `size f = |X|` whose `profile` is exactly the pair of
codimension profiles of the corona over `X`, the second one avoiding `Y`.

This is the bridge that makes the finite enumeration certificate applicable to an actual
graph.
-/

open Finset SimpleGraph

namespace MatchingBag

namespace B1Zero

open DepthThree.RootedForestCertificate

attribute [local instance 1] Classical.propDecidable

variable {V : Type*} [Fintype V] [DecidableEq V] {G : SimpleGraph V}

/-- **Graph-to-forest representation.** -/
theorem exists_forest_repr (hac : G.IsAcyclic) (C : CoronaData G) (X : Finset V) :
    X ⊆ C.base → ∀ Y : Finset V, IsTransversal G X Y →
      ∃ f : Forest, size f = X.card ∧ profile f = (C.vec X ∅, C.vec X Y) := by
  induction X using Finset.strongInduction with
  | _ X ih =>
    intro hXbase Y hY
    rcases X.eq_empty_or_nonempty with rfl | hXne
    · refine ⟨Forest.nil, by simp [size], ?_⟩
      rw [hY.eq_empty rfl]
      simp only [profile, CoronaData.vec_empty]
    · obtain ⟨v, hv⟩ := hY.nonempty hXne
      have hvX : v ∈ X := hY.1 hv
      have hvb : v ∈ C.base := hXbase hvX
      set T := compIn G X v with hT
      set Ch := T.erase v with hCh
      set Rest := X \ T with hRest
      set N := Ch.filter (fun c => G.Adj v c) with hN
      have hTsub : T ⊆ X := compIn_subset
      have hvT : v ∈ T := self_mem_compIn hvX
      have hChsub : Ch ⊆ X := (Finset.erase_subset _ _).trans hTsub
      have hRestsub : Rest ⊆ X := Finset.sdiff_subset
      have hvCh : v ∉ Ch := fun h => (Finset.mem_erase.1 h).1 rfl
      have hTeq : T = insert v Ch := (Finset.insert_erase hvT).symm
      have hChbase : Ch ⊆ C.base := hChsub.trans hXbase
      have hTbase : T ⊆ C.base := hTsub.trans hXbase
      have hRbase : Rest ⊆ C.base := hRestsub.trans hXbase
      have hChss : Ch ⊂ X := ⟨hChsub, fun h => hvCh (h hvX)⟩
      have hRestss : Rest ⊂ X := ⟨hRestsub, fun h => (Finset.mem_sdiff.1 (h hvX)).2 hvT⟩
      obtain ⟨fC, hfCs, hfC⟩ :=
        ih Ch hChss hChbase N (transversal_children hac hY hv)
      obtain ⟨fR, hfRs, hfR⟩ :=
        ih Rest hRestss hRbase (Y.erase v) (transversal_erase hY hv)
      -- the vertex counts
      have hTcard : T.card = Ch.card + 1 := by
        rw [hTeq, Finset.card_insert_of_notMem hvCh]
      have hXcard : Rest.card + T.card = X.card := Finset.card_sdiff_add_card_eq_card hTsub
      -- the splitting of `X`
      have hXeq : X = T ∪ Rest := (Finset.union_sdiff_of_subset hTsub).symm
      have hdisj : Disjoint T Rest := Finset.disjoint_sdiff
      have hcross : ∀ u ∈ T, ∀ w ∈ Rest, ¬ G.Adj u w := by
        intro u hu w hw hadj
        exact (Finset.mem_sdiff.1 hw).2 (compIn_closed hu (Finset.mem_sdiff.1 hw).1 hadj)
      -- `v` is the only element of `Y` inside the corona over `T`
      have hYT : ∀ y ∈ Y, y ∈ C.cor T → y = v := by
        intro y hy hyc
        have hyX : y ∈ X := hY.1 hy
        have hyT : y ∈ T := by
          rcases CoronaData.mem_cor.1 hyc with h | ⟨b, hb, rfl⟩
          · exact h
          · exact absurd (hXbase hyX) (C.pend_notMem b (hTbase hb))
        exact (hY.2.2 v hv y hy (CoronaData.mem_cor.1 hyc |>.elim
          (fun _ => (mem_compIn.1 hyT).2) (fun _ => (mem_compIn.1 hyT).2))).symm ▸ rfl
      have hvcorRest : v ∉ C.cor Rest :=
        CoronaData.notMem_cor_of_base hRbase hvb (fun h => (Finset.mem_sdiff.1 h).2 hvT)
      have hcongrT : C.vec T Y = C.vec T {v} := by
        refine CoronaData.vec_congr ?_
        ext x
        simp only [Finset.mem_sdiff, Finset.mem_singleton]
        constructor
        · rintro ⟨hx, hxY⟩
          exact ⟨hx, fun hxv => hxY (hxv ▸ hv)⟩
        · rintro ⟨hx, hxv⟩
          exact ⟨hx, fun hxY => hxv (hYT x hxY hx)⟩
      have hcongrR : C.vec Rest Y = C.vec Rest (Y.erase v) := by
        refine CoronaData.vec_congr ?_
        ext x
        simp only [Finset.mem_sdiff, Finset.mem_erase]
        constructor
        · rintro ⟨hx, hxY⟩
          exact ⟨hx, fun h => hxY h.2⟩
        · rintro ⟨hx, hxY⟩
          refine ⟨hx, fun h => hxY ⟨?_, h⟩⟩
          rintro rfl
          exact hvcorRest hx
      have hvecX : C.vec X ∅ = mul (C.vec T ∅) (C.vec Rest ∅) := by
        rw [hXeq]
        exact CoronaData.vec_union hTbase hRbase hdisj hcross ∅
      have hvecXY : C.vec X Y = mul (C.vec T Y) (C.vec Rest Y) := by
        rw [hXeq]
        exact CoronaData.vec_union hTbase hRbase hdisj hcross Y
      have hpeel : C.vec T ∅
          = add (add (C.vec Ch ∅) (shift (C.vec Ch ∅))) (C.vec Ch N) := by
        rw [hTeq]
        exact CoronaData.vec_peel hChbase hvb hvCh
      have hpeelr : C.vec T {v} = add (C.vec Ch ∅) (shift (C.vec Ch ∅)) := by
        rw [hTeq]
        exact CoronaData.vec_peel_root hChbase hvb hvCh
      refine ⟨Forest.cons fC fR, ?_, ?_⟩
      · simp only [size, hfCs, hfRs]
        omega
      · simp only [profile, hfC, hfR]
        rw [hvecX, hvecXY, hpeel, hcongrT, hpeelr, hcongrR]

end B1Zero

end MatchingBag
