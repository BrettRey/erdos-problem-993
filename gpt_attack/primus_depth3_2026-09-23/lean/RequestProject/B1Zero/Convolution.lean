import RequestProject.B1Zero.EBounds

/-!
# Counting independent sets by codimension

For the four-forbidden-vertex case we need to count independent subsets of a vertex set `W`
by their *codimension* `j`, i.e. by the deficiency `|S| + j = n` from a reference size `n`.
This file introduces `cntRev` and proves the two facts used throughout:

* `cntRev_choose`: an independent vertex set of size `n` contributes binomial coefficients;
* `cntRev_union`: for a splitting of `W` into two parts with no edges between them, the
  codimension counts convolve.
-/

open Finset SimpleGraph

namespace MatchingBag

namespace B1Zero

attribute [local instance 1] Classical.propDecidable

variable {V : Type*} [Fintype V] [DecidableEq V]

/-- `S` spans no edge of `G`. -/
def IndepOn (G : SimpleGraph V) (S : Finset V) : Prop := ∀ u ∈ S, ∀ v ∈ S, ¬ G.Adj u v

lemma IndepOn.subset {G : SimpleGraph V} {S T : Finset V} (hS : IndepOn G S) (hTS : T ⊆ S) :
    IndepOn G T := fun u hu v hv => hS u (hTS hu) v (hTS hv)

/-- The independent subsets of `W` of codimension `j` below `n`. -/
noncomputable def revSets (G : SimpleGraph V) (W : Finset V) (n j : ℕ) : Finset (Finset V) :=
  W.powerset.filter (fun S => IndepOn G S ∧ S.card + j = n)

/-- The number of independent subsets `S ⊆ W` with `|S| + j = n`. -/
noncomputable def cntRev (G : SimpleGraph V) (W : Finset V) (n j : ℕ) : ℕ :=
  (revSets G W n j).card

lemma mem_revSets {G : SimpleGraph V} {W : Finset V} {n j : ℕ} {S : Finset V} :
    S ∈ revSets G W n j ↔ S ⊆ W ∧ IndepOn G S ∧ S.card + j = n := by
  rw [revSets, Finset.mem_filter, Finset.mem_powerset]

lemma revSets_mono {G : SimpleGraph V} {W W' : Finset V} (h : W ⊆ W') (n j : ℕ) :
    revSets G W n j ⊆ revSets G W' n j := by
  intro S hS
  rw [mem_revSets] at hS ⊢
  exact ⟨hS.1.trans h, hS.2⟩

lemma cntRev_mono {G : SimpleGraph V} {W W' : Finset V} (h : W ⊆ W') (n j : ℕ) :
    cntRev G W n j ≤ cntRev G W' n j :=
  Finset.card_le_card (revSets_mono h n j)

/-- On an independent vertex set the codimension counts are binomial coefficients. -/
theorem cntRev_choose {G : SimpleGraph V} {W : Finset V} (hW : IndepOn G W) (j : ℕ) :
    cntRev G W W.card j = W.card.choose j := by
  classical
  by_cases hj : j ≤ W.card
  · have hset : revSets G W W.card j = W.powersetCard (W.card - j) := by
      ext S
      rw [mem_revSets, Finset.mem_powersetCard]
      constructor
      · rintro ⟨hsub, -, hc⟩
        exact ⟨hsub, by omega⟩
      · rintro ⟨hsub, hc⟩
        exact ⟨hsub, hW.subset hsub, by
          have := Finset.card_le_card hsub
          omega⟩
    rw [cntRev, hset, Finset.card_powersetCard, Nat.choose_symm hj]
  · replace hj : W.card < j := by omega
    have hset : revSets G W W.card j = ∅ := by
      rw [Finset.eq_empty_iff_forall_notMem]
      intro S hS
      rw [mem_revSets] at hS
      omega
    rw [cntRev, hset, Finset.card_empty, Nat.choose_eq_zero_of_lt hj]

/-- **Convolution.**  If `W₁` and `W₂` are disjoint with no edges between them, and
independent subsets of `Wᵢ` have at most `nᵢ` elements, then the codimension counts of the
union are the convolution of those of the parts. -/
theorem cntRev_union {G : SimpleGraph V} {W₁ W₂ : Finset V} {n₁ n₂ : ℕ}
    (hdisj : Disjoint W₁ W₂)
    (hcross : ∀ u ∈ W₁, ∀ v ∈ W₂, ¬ G.Adj u v)
    (hb₁ : ∀ S ⊆ W₁, IndepOn G S → S.card ≤ n₁)
    (hb₂ : ∀ S ⊆ W₂, IndepOn G S → S.card ≤ n₂)
    (j : ℕ) :
    cntRev G (W₁ ∪ W₂) (n₁ + n₂) j
      = ∑ ij ∈ Finset.antidiagonal j, cntRev G W₁ n₁ ij.1 * cntRev G W₂ n₂ ij.2 := by
  classical
  have hsplit : ∀ S : Finset V, S ⊆ W₁ ∪ W₂ →
      (S ∩ W₁).card + (S ∩ W₂).card = S.card := by
    intro S hS
    have hu : (S ∩ W₁) ∪ (S ∩ W₂) = S := by
      rw [← Finset.inter_union_distrib_left]
      exact Finset.inter_eq_left.2 hS
    have hd : Disjoint (S ∩ W₁) (S ∩ W₂) :=
      hdisj.mono Finset.inter_subset_right Finset.inter_subset_right
    have := Finset.card_union_of_disjoint hd
    rw [hu] at this
    omega
  have hcard₁ : ∀ S ∈ revSets G (W₁ ∪ W₂) (n₁ + n₂) j, (S ∩ W₁).card ≤ n₁ := by
    intro S hS
    rw [mem_revSets] at hS
    exact hb₁ _ Finset.inter_subset_right (hS.2.1.subset Finset.inter_subset_left)
  have hcard₂ : ∀ S ∈ revSets G (W₁ ∪ W₂) (n₁ + n₂) j, (S ∩ W₂).card ≤ n₂ := by
    intro S hS
    rw [mem_revSets] at hS
    exact hb₂ _ Finset.inter_subset_right (hS.2.1.subset Finset.inter_subset_left)
  have hmaps : ∀ S ∈ revSets G (W₁ ∪ W₂) (n₁ + n₂) j,
      (n₁ - (S ∩ W₁).card, n₂ - (S ∩ W₂).card) ∈ Finset.antidiagonal j := by
    intro S hS
    have h1 := hcard₁ S hS
    have h2 := hcard₂ S hS
    rw [mem_revSets] at hS
    have h3 := hsplit S hS.1
    simp only [Finset.mem_antidiagonal]
    omega
  rw [cntRev, Finset.card_eq_sum_card_fiberwise
    (f := fun S => (n₁ - (S ∩ W₁).card, n₂ - (S ∩ W₂).card))
    (t := Finset.antidiagonal j)
    (fun S hS => Finset.mem_coe.2 (hmaps S (Finset.mem_coe.1 hS)))]
  refine Finset.sum_congr rfl fun ij hij => ?_
  obtain ⟨j1, j2⟩ := ij
  rw [Finset.mem_antidiagonal] at hij
  rw [cntRev, cntRev, ← Finset.card_product]
  refine Finset.card_bij (fun S _ => (S ∩ W₁, S ∩ W₂)) ?_ ?_ ?_
  · intro S hS
    rw [Finset.mem_filter] at hS
    obtain ⟨hSmem, hfib⟩ := hS
    have h1 := hcard₁ S hSmem
    have h2 := hcard₂ S hSmem
    have hj1 : n₁ - (S ∩ W₁).card = j1 := congrArg Prod.fst hfib
    have hj2 : n₂ - (S ∩ W₂).card = j2 := congrArg Prod.snd hfib
    rw [mem_revSets] at hSmem
    simp only [Finset.mem_product, mem_revSets]
    exact ⟨⟨Finset.inter_subset_right, hSmem.2.1.subset Finset.inter_subset_left, by omega⟩,
      ⟨Finset.inter_subset_right, hSmem.2.1.subset Finset.inter_subset_left, by omega⟩⟩
  · intro S hS S' hS' heq
    rw [Finset.mem_filter, mem_revSets] at hS hS'
    have e1 : S ∩ W₁ = S' ∩ W₁ := congrArg Prod.fst heq
    have e2 : S ∩ W₂ = S' ∩ W₂ := congrArg Prod.snd heq
    have h1 : (S ∩ W₁) ∪ (S ∩ W₂) = S := by
      rw [← Finset.inter_union_distrib_left]
      exact Finset.inter_eq_left.2 hS.1.1
    have h2 : (S' ∩ W₁) ∪ (S' ∩ W₂) = S' := by
      rw [← Finset.inter_union_distrib_left]
      exact Finset.inter_eq_left.2 hS'.1.1
    rw [← h1, ← h2, e1, e2]
  · rintro ⟨A, B⟩ hAB
    simp only [Finset.mem_product, mem_revSets] at hAB
    obtain ⟨⟨hA, hAind, hAc⟩, ⟨hB, hBind, hBc⟩⟩ := hAB
    have hABd : Disjoint A B := hdisj.mono hA hB
    have hAB2 : A ∩ W₂ = ∅ := by
      rw [← Finset.disjoint_iff_inter_eq_empty]
      exact hdisj.mono_left hA
    have hBA1 : B ∩ W₁ = ∅ := by
      rw [← Finset.disjoint_iff_inter_eq_empty]
      exact (hdisj.mono_right hB).symm
    have hind : IndepOn G (A ∪ B) := by
      intro u hu v hv
      rcases Finset.mem_union.1 hu with hu' | hu' <;> rcases Finset.mem_union.1 hv with hv' | hv'
      · exact hAind u hu' v hv'
      · exact hcross u (hA hu') v (hB hv')
      · exact fun h => hcross v (hA hv') u (hB hu') h.symm
      · exact hBind u hu' v hv'
    have e1 : (A ∪ B) ∩ W₁ = A := by
      rw [Finset.union_inter_distrib_right, Finset.inter_eq_left.2 hA, hBA1, Finset.union_empty]
    have e2 : (A ∪ B) ∩ W₂ = B := by
      rw [Finset.union_inter_distrib_right, Finset.inter_eq_left.2 hB, hAB2, Finset.empty_union]
    refine ⟨A ∪ B, ?_, ?_⟩
    · rw [Finset.mem_filter, mem_revSets]
      refine ⟨⟨Finset.union_subset_union hA hB, hind, ?_⟩, ?_⟩
      · rw [Finset.card_union_of_disjoint hABd]
        omega
      · show (n₁ - ((A ∪ B) ∩ W₁).card, n₂ - ((A ∪ B) ∩ W₂).card) = (j1, j2)
        rw [e1, e2]
        simp only [Prod.mk.injEq]
        omega
    · exact Prod.ext e1 e2

end B1Zero

end MatchingBag
