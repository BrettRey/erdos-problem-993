import RequestProject.GraphBridge
import RequestProject.BagPoset

/-!
# The unary-or-pair certificate for blocked independent sets

In the maximum-assignment code of a forest, a nonextendable partial assignment either
selects a value forbidden on a forced coordinate (a *unary* certificate), or violates one
poset comparison between two specified free coordinates (a *pair* certificate).  This is the
projection characterisation `MatchingBag.codeProj_idealCode`.

The graph consequence proved here is
`MatchingBag.TreeMatching.exists_small_blocked_subset`: every blocked independent vertex
set contains a blocked subset of cardinality at most two.
-/

open Finset SimpleGraph

namespace MatchingBag

attribute [local instance 1] Classical.propDecidable

/-- A partial assignment on `K` violating no order comparison inside `K` extends to the
indicator of an order ideal.  This is the surjectivity half of `codeProj_idealCode`. -/
lemma exists_idealIndicator_extending {P : Type*} [Fintype P] [DecidableEq P] [PartialOrder P]
    [DecidableRel ((· ≤ ·) : P → P → Prop)] (K : Finset P) (g : P → Bool)
    (hg : ∀ i ∈ K, ∀ j ∈ K, i ≤ j → g j = true → g i = true) :
    ∃ f : P → Bool, IsIdealIndicator f ∧ ∀ i ∈ K, f i = g i := by
  classical
  have hmem : (fun i : {x // x ∈ K} => g i.1) ∈ inducedIdeals K := by
    rw [inducedIdeals, Finset.mem_filter]
    exact ⟨Finset.mem_univ _, fun i j hij hj => hg i.1 i.2 j.1 j.2 hij hj⟩
  have h2 : extendPartial K (fun i : {x // x ∈ K} => g i.1) ∈ codeProj K (idealCode P) := by
    rw [codeProj_idealCode]
    exact Finset.mem_image_of_mem _ hmem
  rw [codeProj, Finset.mem_image] at h2
  obtain ⟨f, hf, hfe⟩ := h2
  refine ⟨f, by simpa [idealCode] using hf, fun i hi => ?_⟩
  have hc := congrFun hfe i
  simp only [restrictTo, extendPartial, dif_pos hi, if_pos hi, Option.some.injEq] at hc
  exact hc

variable {V : Type*} [Fintype V] [DecidableEq V]

namespace TreeMatching

variable {D : TreeMatching V}

/-! ### Extendability as a constraint on the assignment -/

/-- A set is contained in the independent set complementary to `coverOf x` exactly when the
assignment `x` takes the value prescribed by the set on every bag the set meets. -/
lemma subset_compl_coverOf_iff {S : Finset V} (x : D.Idx → Bool) :
    S ⊆ (D.coverOf x)ᶜ ↔
      ∀ i : D.Idx, (D.Lv i ∈ S → x i = false) ∧ (D.Rv i ∈ S → x i = true) := by
  constructor
  · intro h i
    constructor
    · intro hL
      have := h hL
      rw [Finset.mem_compl] at this
      simpa using fun hc => this (Lv_mem_coverOf.2 hc)
    · intro hR
      have := h hR
      rw [Finset.mem_compl] at this
      by_contra hc
      simp only [Bool.not_eq_true] at hc
      exact this (Rv_mem_coverOf.2 hc)
  · intro h v hv
    rw [Finset.mem_compl]
    rcases vertex_cases D v with ⟨i, rfl⟩ | ⟨i, rfl⟩ | hu
    · intro hc
      rw [Lv_mem_coverOf] at hc
      rw [(h i).1 hv] at hc
      exact Bool.noConfusion hc
    · intro hc
      rw [Rv_mem_coverOf] at hc
      rw [(h i).2 hv] at hc
      exact Bool.noConfusion hc
    · exact unmatched_not_mem_coverOf hu

/-- Extendability of an arbitrary vertex set, read on the constraint system. -/
theorem exists_maxIndep_superset_iff {S : Finset V} :
    (∃ I ∈ D.maxIndepSets, S ⊆ I) ↔
      ∃ x ∈ D.Sol, ∀ i : D.Idx, (D.Lv i ∈ S → x i = false) ∧ (D.Rv i ∈ S → x i = true) := by
  constructor
  · rintro ⟨I, hI, hSI⟩
    rw [maxIndepSets_eq_image_compl] at hI
    obtain ⟨C, hC, rfl⟩ := Finset.mem_image.1 hI
    refine ⟨D.solOf C, solOf_mem_Sol hC, ?_⟩
    rw [← subset_compl_coverOf_iff, coverOf_solOf hC]
    exact hSI
  · rintro ⟨x, hx, hprop⟩
    refine ⟨(D.coverOf x)ᶜ, ?_, (subset_compl_coverOf_iff x).2 hprop⟩
    rw [maxIndepSets_eq_image_compl]
    exact Finset.mem_image_of_mem _ (coverOf_mem_minCovers (mem_Sol.1 hx))

/-! ### The certificate -/

/-- **Unary-or-pair certificate.**  Every blocked independent set contains a blocked subset
of cardinality at most two: either a single vertex sitting on a forced coordinate with the
wrong value, or a pair of vertices violating one poset comparison. -/
theorem exists_small_blocked_subset {S : Finset V}
    (hind : ∀ u ∈ S, ∀ v ∈ S, ¬ D.G.Adj u v)
    (hblk : ¬ ∃ I ∈ D.maxIndepSets, S ⊆ I) :
    ∃ T ⊆ S, T.card ≤ 2 ∧ ¬ ∃ I ∈ D.maxIndepSets, T ⊆ I := by
  classical
  by_contra hcon
  push_neg at hcon
  -- every subset of `S` of size at most two is extendable
  have hall : ∀ T ⊆ S, T.card ≤ 2 → ∃ I ∈ D.maxIndepSets, T ⊆ I := by
    intro T hT hcard
    by_contra hc
    exact hc (hcon T hT hcard)
  -- singleton certificates give the value of the forced coordinates
  have hsingle : ∀ v ∈ S, ∃ x ∈ D.Sol,
      ∀ i : D.Idx, (D.Lv i = v → x i = false) ∧ (D.Rv i = v → x i = true) := by
    intro v hv
    obtain ⟨I, hI, hTI⟩ := hall {v} (Finset.singleton_subset_iff.2 hv) (by simp)
    obtain ⟨x, hx, hprop⟩ := exists_maxIndep_superset_iff.1 ⟨I, hI, hTI⟩
    refine ⟨x, hx, fun i => ⟨fun h => (hprop i).1 ?_, fun h => (hprop i).2 ?_⟩⟩
    · rw [h]; exact Finset.mem_singleton_self v
    · rw [h]; exact Finset.mem_singleton_self v
  -- the partial assignment determined by `S` on the free coordinates
  set g : D.Free → Bool := fun j => decide (D.Rv j.1 ∈ S) with hgdef
  set K : Finset D.Free :=
    Finset.univ.filter (fun j : D.Free => D.Lv j.1 ∈ S ∨ D.Rv j.1 ∈ S) with hKdef
  have hmemK : ∀ j : D.Free, j ∈ K ↔ (D.Lv j.1 ∈ S ∨ D.Rv j.1 ∈ S) := by
    intro j; simp [hKdef]
  -- no pair certificate: the assignment is an ideal on `K`
  have hideal : ∀ i ∈ K, ∀ j ∈ K, i ≤ j → g j = true → g i = true := by
    intro i hi j hj hij hgj
    by_contra hgi
    rw [hgdef] at hgi hgj
    simp only [decide_eq_true_eq] at hgj
    simp only [decide_eq_true_eq] at hgi
    have hLi : D.Lv i.1 ∈ S := ((hmemK i).1 hi).resolve_right hgi
    obtain ⟨I, hI, hTI⟩ := hall {D.Lv i.1, D.Rv j.1}
      (by
        intro z hz
        rcases Finset.mem_insert.1 hz with rfl | hz
        · exact hLi
        · rw [Finset.mem_singleton] at hz; rw [hz]; exact hgj)
      ((Finset.card_insert_le _ _).trans (by simp))
    obtain ⟨x, hx, hprop⟩ := exists_maxIndep_superset_iff.1 ⟨I, hI, hTI⟩
    have hxi : x i.1 = false := (hprop i.1).1 (by simp)
    have hxj : x j.1 = true := (hprop j.1).2 (by simp)
    have hy := isIdealIndicator_restrictFree hx i j hij
    rw [restrictFree, restrictFree] at hy
    have := hy hxj
    rw [hxi] at this
    exact Bool.noConfusion this
  obtain ⟨f, hf, hfK⟩ := exists_idealIndicator_extending K g hideal
  -- build a solution realising `S`
  refine hblk ?_
  refine exists_maxIndep_superset_iff.2 ⟨D.extendFree f, extendFree_mem_Sol hf, fun i => ?_⟩
  by_cases hFi : D.Forced i
  · rw [extendFree_of_forced hFi]
    constructor
    · intro hL
      obtain ⟨x, hx, hprop⟩ := hsingle _ hL
      rw [← forcedVal_spec hFi hx]
      exact (hprop i).1 rfl
    · intro hR
      obtain ⟨x, hx, hprop⟩ := hsingle _ hR
      rw [← forcedVal_spec hFi hx]
      exact (hprop i).2 rfl
  · rw [extendFree_of_not_forced hFi]
    have hjK : ∀ _ : (D.Lv i ∈ S ∨ D.Rv i ∈ S), (⟨i, hFi⟩ : D.Free) ∈ K := by
      intro h; exact (hmemK ⟨i, hFi⟩).2 h
    constructor
    · intro hL
      rw [hfK _ (hjK (Or.inl hL)), hgdef]
      have hR : D.Rv i ∉ S := fun hh => hind _ hL _ hh (adj_Lv_Rv i)
      simp [hR]
    · intro hR
      rw [hfK _ (hjK (Or.inr hR)), hgdef]
      simp [hR]

end TreeMatching

end MatchingBag
