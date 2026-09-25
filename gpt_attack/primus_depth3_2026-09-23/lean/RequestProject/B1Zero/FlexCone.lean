import RequestProject.B1Zero.Decompose

/-!
# Flexible bags and the matching-bag deletion inequality

* `flexBags`: the bags containing no forced vertex; there are `α - |Forced|` of them and each
  consists of two flexible vertices.
* `pflex_eq_zero_of_gt`: `p_l = 0` for `l` larger than the number of flexible bags.
* `pflex_cone`: `(l+1) · p_{l+1} ≤ 2 (m - l) · p_l`, the deletion inequality (5) of the
  written proof.
-/

open Finset SimpleGraph

namespace MatchingBag

attribute [local instance 1] Classical.propDecidable

variable {V : Type*} [Fintype V] [DecidableEq V]

namespace TreeMatching

variable (D : TreeMatching V)

/-- The bags containing no forced vertex. -/
noncomputable def flexBags : Finset D.Bag :=
  Finset.univ.filter (fun b => D.ForcedV ∩ D.bagSet b = ∅)

variable {D}

lemma mem_flexBags {b : D.Bag} : b ∈ D.flexBags ↔ D.ForcedV ∩ D.bagSet b = ∅ := by
  simp [flexBags]

lemma flexBags_all_flex {b : D.Bag} (hb : b ∈ D.flexBags) :
    ∀ v ∈ D.bagSet b, v ∈ D.FlexV := by
  refine (bagSet_flex_of_no_forced ?_).2
  intro v hv hvF
  have := mem_flexBags.1 hb
  rw [Finset.eq_empty_iff_forall_notMem] at this
  exact this v (Finset.mem_inter.2 ⟨hvF, hv⟩)

lemma bagOf_mem_flexBags {v : V} (hv : v ∈ D.FlexV) : D.bagOf v ∈ D.flexBags := by
  rw [mem_flexBags, Finset.eq_empty_iff_forall_notMem]
  intro w hw
  rw [Finset.mem_inter] at hw
  exact ((bag_of_flex hv).2 w hw.2 |> mem_Flex.1).2 hw.1

lemma card_inter_flex_bagSet (b : D.Bag) :
    (D.FlexV ∩ D.bagSet b).card = if b ∈ D.flexBags then 2 else 0 := by
  obtain ⟨-, h2⟩ := bag_local_counts (D := D) b
  by_cases hb : b ∈ D.flexBags
  · rw [if_pos hb]
    have : (D.ForcedV ∩ D.bagSet b).card = 0 := by rw [mem_flexBags.1 hb]; simp
    omega
  · rw [if_neg hb]
    have hne : (D.ForcedV ∩ D.bagSet b).Nonempty := by
      rw [mem_flexBags] at hb
      exact Finset.nonempty_of_ne_empty hb
    have := Finset.card_pos.2 hne
    omega

/-- There are `α - |Forced|` flexible bags. -/
theorem card_flexBags (D : TreeMatching V) :
    D.flexBags.card + D.ForcedV.card = Fintype.card D.Bag := by
  classical
  have h1 : D.FlexV.card = 2 * D.flexBags.card := by
    rw [card_eq_sum_inter_bagSet (D := D) D.FlexV]
    rw [Finset.sum_congr rfl (fun b _ => card_inter_flex_bagSet (D := D) b)]
    rw [Finset.sum_ite_mem, Finset.univ_inter, Finset.sum_const, smul_eq_mul, mul_comm]
  have h2 := card_flex D
  omega

/-- An independent set of flexible vertices meets each flexible bag at most once. -/
lemma card_le_card_flexBags {S : Finset V} (hS : S ⊆ D.FlexV)
    (hind : ∀ u ∈ S, ∀ v ∈ S, ¬ D.G.Adj u v) : S.card ≤ D.flexBags.card := by
  classical
  have := Finset.card_image_of_injOn (bagOf_injOn_indep hind) (s := S)
  rw [← this]
  refine Finset.card_le_card ?_
  intro b hb
  obtain ⟨v, hv, rfl⟩ := Finset.mem_image.1 hb
  exact bagOf_mem_flexBags (hS hv)

lemma pflex_eq_zero_of_gt {l : ℕ} (hl : D.flexBags.card < l) : D.pflex l = 0 := by
  rw [pflex, Finset.card_eq_zero, Finset.eq_empty_iff_forall_notMem]
  intro S hS
  obtain ⟨hsub, hcard, hind⟩ := mem_flexInd.1 hS
  have := card_le_card_flexBags hsub hind
  omega

/-! ### The deletion inequality -/

/-- **Matching-bag deletion inequality.**  `(l+1) p_{l+1} ≤ 2 (m - l) p_l`, where `m` is the
number of flexible bags. -/
theorem pflex_cone (D : TreeMatching V) (l : ℕ) :
    (l + 1) * D.pflex (l + 1) ≤ 2 * (D.flexBags.card - l) * D.pflex l := by
  classical
  set m := D.flexBags.card with hm
  set E1 := D.flexInd (l + 1) with hE1
  set E0 := D.flexInd l with hE0
  set Inc : Finset (V × Finset V) := E1.biUnion (fun S => S.image (fun v => (v, S))) with hInc
  have hmemInc1 : ∀ p ∈ Inc, p.1 ∈ p.2 := by
    intro p hp
    rw [hInc, Finset.mem_biUnion] at hp
    obtain ⟨S, -, hp2⟩ := hp
    obtain ⟨v, hv, rfl⟩ := Finset.mem_image.1 hp2
    exact hv
  have hmemInc2 : ∀ p ∈ Inc, p.2 ∈ E1 := by
    intro p hp
    rw [hInc, Finset.mem_biUnion] at hp
    obtain ⟨S, hS, hp2⟩ := hp
    obtain ⟨v, -, rfl⟩ := Finset.mem_image.1 hp2
    exact hS
  have hdisj : ∀ S ∈ E1, ∀ S' ∈ E1, S ≠ S' →
      Disjoint (S.image (fun v => (v, S))) (S'.image (fun v => (v, S'))) := by
    intro S _ S' _ hne
    rw [Finset.disjoint_left]
    rintro p hp hp'
    obtain ⟨v, -, rfl⟩ := Finset.mem_image.1 hp
    obtain ⟨v', -, hv'⟩ := Finset.mem_image.1 hp'
    exact hne (congrArg Prod.snd hv').symm
  have hcount : Inc.card = (l + 1) * E1.card := by
    rw [hInc, Finset.card_biUnion hdisj]
    rw [Finset.sum_congr rfl (fun S hS => (Finset.card_image_of_injOn
      (fun x _ y _ hxy => congrArg Prod.fst hxy)).trans (mem_flexInd.1 hS).2.1)]
    rw [Finset.sum_const, smul_eq_mul, mul_comm]
  have hmaps : ∀ p ∈ Inc, p.2.erase p.1 ∈ E0 := by
    intro p hp
    have h2 := hmemInc2 p hp
    have h1 := hmemInc1 p hp
    obtain ⟨hsub, hcard, hind⟩ := mem_flexInd.1 h2
    refine mem_flexInd.2 ⟨(Finset.erase_subset _ _).trans hsub, ?_, ?_⟩
    · rw [Finset.card_erase_of_mem h1, hcard]
      omega
    · intro u hu v hv
      exact hind u (Finset.mem_of_mem_erase hu) v (Finset.mem_of_mem_erase hv)
  -- each fibre has at most `2 (m - l)` elements
  have hfib : ∀ T ∈ E0, (Inc.filter (fun p => p.2.erase p.1 = T)).card ≤ 2 * (m - l) := by
    intro T hT
    obtain ⟨hTsub, hTcard, hTind⟩ := mem_flexInd.1 hT
    -- the possible new vertices lie in flexible bags not met by `T`
    set W := D.FlexV.filter (fun v => D.bagOf v ∉ T.image D.bagOf) with hW
    have hTim : (T.image D.bagOf).card = l := by
      rw [Finset.card_image_of_injOn (bagOf_injOn_indep hTind), hTcard]
    have hsubb : T.image D.bagOf ⊆ D.flexBags := by
      intro b hb
      obtain ⟨v, hv, rfl⟩ := Finset.mem_image.1 hb
      exact bagOf_mem_flexBags (hTsub hv)
    have hWcard : W.card ≤ 2 * (m - l) := by
      have hsub : W ⊆ (D.flexBags \ T.image D.bagOf).biUnion (fun b => D.bagSet b) := by
        intro v hv
        rw [hW, Finset.mem_filter] at hv
        exact Finset.mem_biUnion.2 ⟨D.bagOf v,
          Finset.mem_sdiff.2 ⟨bagOf_mem_flexBags hv.1, hv.2⟩, mem_bagSet_bagOf v⟩
      refine le_trans (Finset.card_le_card hsub) (le_trans Finset.card_biUnion_le ?_)
      calc ∑ b ∈ D.flexBags \ T.image D.bagOf, (D.bagSet b).card
          ≤ ∑ _b ∈ D.flexBags \ T.image D.bagOf, 2 :=
            Finset.sum_le_sum fun b _ => card_bagSet_le b
        _ = 2 * (m - l) := by
            rw [Finset.sum_const, smul_eq_mul, mul_comm, Finset.card_sdiff,
              Finset.inter_eq_left.2 hsubb, hTim]
    refine le_trans (Finset.card_le_card_of_injOn (fun p => p.1) ?_ ?_) hWcard
    · intro p hp
      rw [Finset.mem_coe, Finset.mem_filter] at hp
      have h2 := hmemInc2 p hp.1
      have h1 := hmemInc1 p hp.1
      obtain ⟨hsub, hcard, hind⟩ := mem_flexInd.1 h2
      rw [hW, Finset.mem_coe, Finset.mem_filter]
      refine ⟨hsub h1, ?_⟩
      intro hcon
      obtain ⟨w, hw, hbw⟩ := Finset.mem_image.1 hcon
      have hwT : w ∈ p.2 := by
        rw [← hp.2] at hw
        exact Finset.mem_of_mem_erase hw
      have hwne : w ≠ p.1 := by
        rw [← hp.2] at hw
        exact Finset.ne_of_mem_erase hw
      exact hwne (bagOf_injOn_indep hind (by exact_mod_cast hwT) (by exact_mod_cast h1) hbw)
    · intro p hp q hq hpq
      rw [Finset.mem_coe, Finset.mem_filter] at hp hq
      have hpq' : p.1 = q.1 := hpq
      have h1 : p.2 = q.2 := by
        rw [← Finset.insert_erase (hmemInc1 p hp.1), ← Finset.insert_erase (hmemInc1 q hq.1),
          hp.2, hq.2, hpq']
      exact Prod.ext hpq' h1
  have hkey : Inc.card ≤ 2 * (m - l) * E0.card := by
    rw [Finset.card_eq_sum_card_fiberwise (t := E0) (fun p hp => hmaps p (Finset.mem_coe.1 hp))]
    calc ∑ T ∈ E0, (Inc.filter (fun p => p.2.erase p.1 = T)).card
        ≤ ∑ _T ∈ E0, 2 * (m - l) := Finset.sum_le_sum hfib
      _ = 2 * (m - l) * E0.card := by rw [Finset.sum_const, smul_eq_mul, mul_comm]
  rw [pflex, pflex, ← hE1, ← hE0]
  omega

end TreeMatching

end MatchingBag
