import RequestProject.B1Zero.ForcedDegree

/-!
# Decomposing independent sets of allowed vertices

Forced vertices are isolated inside the allowed subgraph, so an independent set of allowed
vertices splits as an arbitrary set of forced vertices together with an independent set of
flexible vertices.  This gives the convolution formula

`#{independent S ⊆ F ∪ Flex, |S| = k} = ∑_{i+j=k} C(|F|, i) · p_j`,

where `p_j` counts the independent `j`-subsets of `Flex`.
-/

open Finset SimpleGraph

namespace MatchingBag

attribute [local instance 1] Classical.propDecidable

variable {V : Type*} [Fintype V] [DecidableEq V]

namespace TreeMatching

variable (D : TreeMatching V)

/-- Independent subsets of the flexible vertices of a given size. -/
noncomputable def flexInd (k : ℕ) : Finset (Finset V) :=
  (D.FlexV.powersetCard k).filter (fun S => ∀ u ∈ S, ∀ v ∈ S, ¬ D.G.Adj u v)

/-- The number of independent `k`-subsets of the flexible vertices. -/
noncomputable def pflex (k : ℕ) : ℕ := (D.flexInd k).card

/-- Independent subsets of `F ∪ Flex` of a given size, for `F` a set of forced vertices. -/
noncomputable def indSubsets (F : Finset V) (k : ℕ) : Finset (Finset V) :=
  ((F ∪ D.FlexV).powersetCard k).filter (fun S => ∀ u ∈ S, ∀ v ∈ S, ¬ D.G.Adj u v)

/-- The convolution `∑_{i+j=k} C(cc, i) · p_j`. -/
noncomputable def EE (cc k : ℕ) : ℕ :=
  ∑ ij ∈ Finset.antidiagonal k, Nat.choose cc ij.1 * D.pflex ij.2

variable {D}

lemma Flex_disjoint_Forced : Disjoint D.FlexV D.ForcedV := by
  rw [Finset.disjoint_left]
  intro v hv
  exact (mem_Flex.1 hv).2

lemma flexInd_subset {k : ℕ} {S : Finset V} (hS : S ∈ D.flexInd k) :
    S ⊆ D.FlexV ∧ S.card = k ∧ ∀ u ∈ S, ∀ v ∈ S, ¬ D.G.Adj u v := by
  rw [flexInd, Finset.mem_filter, Finset.mem_powersetCard] at hS
  exact ⟨hS.1.1, hS.1.2, hS.2⟩

lemma mem_flexInd {k : ℕ} {S : Finset V} :
    S ∈ D.flexInd k ↔ S ⊆ D.FlexV ∧ S.card = k ∧ ∀ u ∈ S, ∀ v ∈ S, ¬ D.G.Adj u v := by
  rw [flexInd, Finset.mem_filter, Finset.mem_powersetCard]
  tauto

lemma mem_indSubsets {F : Finset V} {k : ℕ} {S : Finset V} :
    S ∈ D.indSubsets F k ↔
      S ⊆ F ∪ D.FlexV ∧ S.card = k ∧ ∀ u ∈ S, ∀ v ∈ S, ¬ D.G.Adj u v := by
  rw [indSubsets, Finset.mem_filter, Finset.mem_powersetCard]
  tauto

/-- Independence of a set of allowed vertices only constrains its flexible part. -/
lemma indep_iff_flex_part {F : Finset V} (hF : F ⊆ D.ForcedV) {S : Finset V}
    (hS : S ⊆ F ∪ D.FlexV) :
    (∀ u ∈ S, ∀ v ∈ S, ¬ D.G.Adj u v) ↔
      ∀ u ∈ S ∩ D.FlexV, ∀ v ∈ S ∩ D.FlexV, ¬ D.G.Adj u v := by
  constructor
  · intro h u hu v hv
    exact h u (Finset.mem_inter.1 hu).1 v (Finset.mem_inter.1 hv).1
  · intro h u hu v hv hadj
    have hu' := hS hu
    have hv' := hS hv
    rw [Finset.mem_union] at hu' hv'
    rcases hu' with huF | huL
    · exact Finset.disjoint_left.1 Allowed_disjoint_Forbidden
        (rcases_allowed hv') (forbidden_of_adj_forced (hF huF) hadj)
    · rcases hv' with hvF | hvL
      · exact Finset.disjoint_left.1 Allowed_disjoint_Forbidden
          (mem_Flex.1 huL).1 (forbidden_of_adj_forced (hF hvF) hadj.symm)
      · exact h u (Finset.mem_inter.2 ⟨hu, huL⟩) v (Finset.mem_inter.2 ⟨hv, hvL⟩) hadj
where
  rcases_allowed {v : V} (h : v ∈ F ∨ v ∈ D.FlexV) : v ∈ D.AllowedV := by
    rcases h with h | h
    · exact Forced_subset_Allowed (hF h)
    · exact (mem_Flex.1 h).1

/-- **The convolution formula.** -/
theorem card_indSubsets {F : Finset V} (hF : F ⊆ D.ForcedV) (k : ℕ) :
    (D.indSubsets F k).card = D.EE F.card k := by
  classical
  have hdisj : Disjoint F D.FlexV :=
    (Flex_disjoint_Forced.mono_right hF).symm
  have hmaps : ∀ S ∈ D.indSubsets F k,
      ((S ∩ F).card, (S ∩ D.FlexV).card) ∈ Finset.antidiagonal k := by
    intro S hS
    rw [mem_indSubsets] at hS
    have hu : (S ∩ F) ∪ (S ∩ D.FlexV) = S := by
      rw [← Finset.inter_union_distrib_left]
      exact Finset.inter_eq_left.2 hS.1
    have hd : Disjoint (S ∩ F) (S ∩ D.FlexV) :=
      hdisj.mono Finset.inter_subset_right Finset.inter_subset_right
    have hc := Finset.card_union_of_disjoint hd
    rw [hu] at hc
    have hk := hS.2.1
    simp only [Finset.mem_antidiagonal]
    omega
  rw [Finset.card_eq_sum_card_fiberwise
    (f := fun S => ((S ∩ F).card, (S ∩ D.FlexV).card)) (t := Finset.antidiagonal k)
    (fun S hS => Finset.mem_coe.2 (hmaps S (Finset.mem_coe.1 hS))), EE]
  refine Finset.sum_congr rfl fun ij _ => ?_
  simp only [pflex, flexInd]
  rw [← Finset.card_powersetCard, ← Finset.card_product]
  refine Finset.card_bij (fun S _ => (S ∩ F, S ∩ D.FlexV)) ?_ ?_ ?_
  · intro S hS
    rw [Finset.mem_filter, mem_indSubsets] at hS
    obtain ⟨⟨hsub, hcard, hind⟩, hfib⟩ := hS
    rw [Finset.mem_product, Finset.mem_powersetCard, Finset.mem_filter, Finset.mem_powersetCard]
    refine ⟨⟨Finset.inter_subset_right, ?_⟩, ⟨Finset.inter_subset_right, ?_⟩, ?_⟩
    · exact congrArg Prod.fst hfib
    · exact congrArg Prod.snd hfib
    · intro u hu v hv
      exact hind u (Finset.mem_inter.1 hu).1 v (Finset.mem_inter.1 hv).1
  · intro S hS S' hS' heq
    rw [Finset.mem_filter, mem_indSubsets] at hS hS'
    have e1 : S ∩ F = S' ∩ F := congrArg Prod.fst heq
    have e2 : S ∩ D.FlexV = S' ∩ D.FlexV := congrArg Prod.snd heq
    have h1 : (S ∩ F) ∪ (S ∩ D.FlexV) = S := by
      rw [← Finset.inter_union_distrib_left]
      exact Finset.inter_eq_left.2 hS.1.1
    have h2 : (S' ∩ F) ∪ (S' ∩ D.FlexV) = S' := by
      rw [← Finset.inter_union_distrib_left]
      exact Finset.inter_eq_left.2 hS'.1.1
    rw [← h1, ← h2, e1, e2]
  · rintro ⟨A, B⟩ hAB
    rw [Finset.mem_product, Finset.mem_powersetCard, Finset.mem_filter,
      Finset.mem_powersetCard] at hAB
    obtain ⟨⟨hA, hAc⟩, ⟨hB, hBc⟩, hBind⟩ := hAB
    have hABd : Disjoint A B := hdisj.mono hA hB
    have hAF : A ∩ D.FlexV = ∅ := by
      rw [← Finset.disjoint_iff_inter_eq_empty]
      exact hdisj.mono_left hA
    have hBF : B ∩ F = ∅ := by
      rw [← Finset.disjoint_iff_inter_eq_empty]
      exact (hdisj.mono_right hB).symm
    refine ⟨A ∪ B, ?_, ?_⟩
    · rw [Finset.mem_filter, mem_indSubsets]
      have hsub : A ∪ B ⊆ F ∪ D.FlexV := Finset.union_subset_union hA hB
      have hAB1 : (A ∪ B) ∩ F = A := by
        rw [Finset.union_inter_distrib_right, Finset.inter_eq_left.2 hA, hBF,
          Finset.union_empty]
      have hAB2 : (A ∪ B) ∩ D.FlexV = B := by
        rw [Finset.union_inter_distrib_right, Finset.inter_eq_left.2 hB, hAF,
          Finset.empty_union]
      refine ⟨⟨hsub, ?_, ?_⟩, ?_⟩
      · rw [Finset.card_union_of_disjoint hABd, hAc, hBc]
        exact (Finset.mem_antidiagonal.1 (by assumption)).symm ▸ rfl
      · rw [indep_iff_flex_part hF hsub, hAB2]
        exact hBind
      · rw [hAB1, hAB2, hAc, hBc]
    · have hAB1 : (A ∪ B) ∩ F = A := by
        rw [Finset.union_inter_distrib_right, Finset.inter_eq_left.2 hA, hBF,
          Finset.union_empty]
      have hAB2 : (A ∪ B) ∩ D.FlexV = B := by
        rw [Finset.union_inter_distrib_right, Finset.inter_eq_left.2 hB, hAF,
          Finset.empty_union]
      exact Prod.ext hAB1 hAB2

/-- `EE` is monotone in the number of forced vertices. -/
lemma EE_mono {cc cc' k : ℕ} (h : cc ≤ cc') : D.EE cc k ≤ D.EE cc' k := by
  refine Finset.sum_le_sum fun ij _ => Nat.mul_le_mul_right _ ?_
  exact Nat.choose_le_choose _ h

end TreeMatching

end MatchingBag
