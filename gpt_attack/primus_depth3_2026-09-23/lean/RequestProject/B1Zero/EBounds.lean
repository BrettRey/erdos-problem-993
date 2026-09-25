import RequestProject.B1Zero.FlexCone

/-!
# The extendable count and the blocked bound under `b₁ = 0`

* `card_extSets_eq_EE`: the number of extendable independent `k`-sets is `EE c k`, where
  `c = |Forced|`.
* `card_blkSets_le`: the number of blocked independent `k`-sets is at most
  `∑_{t ≥ 1} C(r, t) · EE (c - (2t+1)) (k - t)`, where `r = |Forbidden|`.
-/

open Finset SimpleGraph

namespace MatchingBag

attribute [local instance 1] Classical.propDecidable

variable {V : Type*} [Fintype V] [DecidableEq V]

namespace TreeMatching

variable {D : TreeMatching V}

lemma Forced_union_Flex : D.ForcedV ∪ D.FlexV = D.AllowedV := by
  ext v
  simp only [Finset.mem_union, mem_Flex]
  constructor
  · rintro (h | h)
    · exact Forced_subset_Allowed h
    · exact h.1
  · intro h
    by_cases hf : v ∈ D.ForcedV
    · exact Or.inl hf
    · exact Or.inr ⟨h, hf⟩

lemma extSets_eq_indSubsets (hb1 : D.B1Zero) (k : ℕ) :
    D.extSets k = D.indSubsets D.ForcedV k := by
  ext S
  rw [mem_extSets, mem_indSubsets, Forced_union_Flex]
  constructor
  · rintro ⟨hcard, hind, hext⟩
    exact ⟨(extendable_iff_subset_allowed hb1 hind).1 hext, hcard, hind⟩
  · rintro ⟨hsub, hcard, hind⟩
    exact ⟨hcard, hind, (extendable_iff_subset_allowed hb1 hind).2 hsub⟩

/-- **The extendable count.** -/
theorem card_extSets_eq_EE (hb1 : D.B1Zero) (k : ℕ) :
    (D.extSets k).card = D.EE D.ForcedV.card k := by
  rw [extSets_eq_indSubsets hb1, card_indSubsets (subset_refl _)]

lemma mem_blkSets_iff (hb1 : D.B1Zero) {k : ℕ} {S : Finset V} :
    S ∈ D.blkSets k ↔ S.card = k ∧ (∀ u ∈ S, ∀ v ∈ S, ¬ D.G.Adj u v) ∧
      (S ∩ D.ForbiddenV).Nonempty := by
  rw [mem_blkSets]
  constructor
  · rintro ⟨hcard, hind, hnot⟩
    refine ⟨hcard, hind, ?_⟩
    by_contra hcon
    rw [Finset.not_nonempty_iff_eq_empty] at hcon
    refine hnot ((extendable_iff_subset_allowed hb1 hind).2 fun v hv => ?_)
    by_contra hvA
    have hvF : v ∈ D.ForbiddenV := by rw [Forbidden_eq_compl, Finset.mem_compl]; exact hvA
    have hmem : v ∈ S ∩ D.ForbiddenV := Finset.mem_inter.2 ⟨hv, hvF⟩
    rw [hcon] at hmem
    exact absurd hmem (Finset.notMem_empty v)
  · rintro ⟨hcard, hind, hne⟩
    refine ⟨hcard, hind, fun hext => ?_⟩
    obtain ⟨v, hv⟩ := hne
    rw [Finset.mem_inter] at hv
    exact Finset.disjoint_left.1 Allowed_disjoint_Forbidden
      ((extendable_iff_subset_allowed hb1 hind).1 hext hv.1) hv.2

/-- **The blocked bound.** -/
theorem card_blkSets_le (hb1 : D.B1Zero) (ha : 3 ≤ Fintype.card D.Bag) (k : ℕ) :
    (D.blkSets k).card ≤ ∑ t ∈ Finset.Icc 1 D.ForbiddenV.card,
      Nat.choose D.ForbiddenV.card t * D.EE (D.ForcedV.card - (2 * t + 1)) (k - t) := by
  classical
  set r := D.ForbiddenV.card with hr
  set c := D.ForcedV.card with hc
  set g : ℕ → ℕ := fun t => if t = 0 then 0 else D.EE (c - (2 * t + 1)) (k - t) with hg
  have hmaps : ∀ S ∈ D.blkSets k, S ∩ D.ForbiddenV ∈ D.ForbiddenV.powerset := fun S _ =>
    Finset.mem_powerset.2 Finset.inter_subset_right
  have hfiber : ∀ T ∈ D.ForbiddenV.powerset,
      ((D.blkSets k).filter (fun S => S ∩ D.ForbiddenV = T)).card ≤ g T.card := by
    intro T hT
    rw [Finset.mem_powerset] at hT
    rcases Finset.eq_empty_or_nonempty T with rfl | hTne
    · -- a blocked set meets the forbidden vertices, so this fibre is empty
      have : ((D.blkSets k).filter (fun S => S ∩ D.ForbiddenV = ∅)) = ∅ := by
        rw [Finset.eq_empty_iff_forall_notMem]
        intro S hS
        rw [Finset.mem_filter, mem_blkSets_iff hb1] at hS
        obtain ⟨⟨-, -, hne⟩, hemp⟩ := hS
        rw [hemp] at hne
        exact hne.ne_empty rfl
      rw [this]
      simp
    · have hgT : g T.card = D.EE (c - (2 * T.card + 1)) (k - T.card) := by
        rw [hg]
        simp only
        rw [if_neg (by simpa using hTne.card_pos.ne')]
      rw [hgT]
      set N := D.forcedNbrsSet T with hN
      have hNsub : N ⊆ D.ForcedV := Finset.filter_subset _ _
      have hNcard : 2 * T.card + 1 ≤ N.card := card_forcedNbrsSet_ge hb1 ha hT hTne
      have hFN : (D.ForcedV \ N).card = c - N.card := by
        rw [Finset.card_sdiff, Finset.inter_eq_left.2 hNsub]
      have hmono : D.EE (D.ForcedV \ N).card (k - T.card)
          ≤ D.EE (c - (2 * T.card + 1)) (k - T.card) := by
        refine EE_mono ?_
        rw [hFN]
        omega
      refine le_trans (le_trans ?_ (le_of_eq (card_indSubsets
        (D := D) (F := D.ForcedV \ N) (Finset.sdiff_subset) (k - T.card)))) hmono
      refine Finset.card_le_card_of_injOn (fun S => S \ T) ?_ ?_
      · intro S hS
        rw [Finset.mem_coe, Finset.mem_filter, mem_blkSets_iff hb1] at hS
        obtain ⟨⟨hcard, hind, -⟩, hST⟩ := hS
        have hTS : T ⊆ S := by
          rw [← hST]; exact Finset.inter_subset_left
        rw [Finset.mem_coe, mem_indSubsets]
        refine ⟨?_, ?_, ?_⟩
        · intro v hv
          rw [Finset.mem_sdiff] at hv
          have hvA : v ∈ D.AllowedV := by
            by_contra hvA
            have : v ∈ D.ForbiddenV := by
              rw [Forbidden_eq_compl, Finset.mem_compl]; exact hvA
            exact hv.2 (hST ▸ Finset.mem_inter.2 ⟨hv.1, this⟩)
          rw [← Forced_union_Flex, Finset.mem_union] at hvA
          rcases hvA with hvF | hvL
          · refine Finset.mem_union_left _ (Finset.mem_sdiff.2 ⟨hvF, ?_⟩)
            intro hvN
            obtain ⟨-, u, huT, hadj⟩ := mem_forcedNbrsSet.1 hvN
            exact hind u (hTS huT) v hv.1 hadj
          · exact Finset.mem_union_right _ hvL
        · rw [Finset.card_sdiff, Finset.inter_eq_left.2 hTS, hcard]
        · intro u hu v hv
          exact hind u (Finset.mem_sdiff.1 hu).1 v (Finset.mem_sdiff.1 hv).1
      · intro S hS S' hS' heq
        rw [Finset.mem_coe, Finset.mem_filter] at hS hS'
        have hTS : T ⊆ S := by rw [← hS.2]; exact Finset.inter_subset_left
        have hTS' : T ⊆ S' := by rw [← hS'.2]; exact Finset.inter_subset_left
        have h1 : T ∪ (S \ T) = S := Finset.union_sdiff_of_subset hTS
        have h2 : T ∪ (S' \ T) = S' := Finset.union_sdiff_of_subset hTS'
        rw [← h1, ← h2]
        simp only at heq
        rw [heq]
  have hstep : (D.blkSets k).card ≤ ∑ T ∈ D.ForbiddenV.powerset, g T.card := by
    rw [Finset.card_eq_sum_card_fiberwise
      (fun S hS => Finset.mem_coe.2 (hmaps S (Finset.mem_coe.1 hS)))]
    exact Finset.sum_le_sum hfiber
  have hgroup : ∑ T ∈ D.ForbiddenV.powerset, g T.card
      = ∑ t ∈ Finset.range (r + 1), Nat.choose r t * g t := by
    rw [Finset.sum_powerset]
    refine Finset.sum_congr rfl fun j _ => ?_
    have hcst : ∀ T ∈ Finset.powersetCard j D.ForbiddenV, g T.card = g j := by
      intro T hT; rw [(Finset.mem_powersetCard.1 hT).2]
    rw [Finset.sum_congr rfl hcst, Finset.sum_const, Finset.card_powersetCard, smul_eq_mul,
      ← hr]
  have hrange : Finset.range (r + 1) = insert 0 (Finset.Icc 1 r) := by
    ext x; simp only [Finset.mem_range, Finset.mem_insert, Finset.mem_Icc]; omega
  calc (D.blkSets k).card
      ≤ ∑ T ∈ D.ForbiddenV.powerset, g T.card := hstep
    _ = ∑ t ∈ Finset.range (r + 1), Nat.choose r t * g t := hgroup
    _ = ∑ t ∈ Finset.Icc 1 r, Nat.choose r t * g t := by
        rw [hrange, Finset.sum_insert (by simp)]
        simp [hg]
    _ = ∑ t ∈ Finset.Icc 1 r, Nat.choose r t * D.EE (c - (2 * t + 1)) (k - t) := by
        refine Finset.sum_congr rfl fun t ht => ?_
        rw [Finset.mem_Icc] at ht
        simp only [hg, if_neg (show t ≠ 0 by omega)]

end TreeMatching

end MatchingBag
