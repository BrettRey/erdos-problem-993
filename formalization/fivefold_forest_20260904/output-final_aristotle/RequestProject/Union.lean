import RequestProject.Basic

/-!
# Disjoint unions

If the vertex set `S` splits as `S₁ ∪ S₂` with no edges between the two parts,
then independent sets, extendable sets and blocked sets all factor as
convolutions.
-/

open Finset

namespace FivefoldForest

open scoped Classical

variable {V : Type*} {G : SimpleGraph V} {S S₁ S₂ A B : Finset V}

/-- There are no edges between `S₁` and `S₂`. -/
def NoEdges (G : SimpleGraph V) (S₁ S₂ : Finset V) : Prop :=
  ∀ x ∈ S₁, ∀ y ∈ S₂, ¬ G.Adj x y

/-- `S` is the disjoint union of the mutually non-adjacent sets `S₁` and `S₂`. -/
structure Split (G : SimpleGraph V) (S S₁ S₂ : Finset V) : Prop where
  eq_union : S = S₁ ∪ S₂
  disj : Disjoint S₁ S₂
  noEdges : NoEdges G S₁ S₂

/-- The number of independent sets of `S` at defect `d`. -/
noncomputable def iD (G : SimpleGraph V) (S : Finset V) (d : ℕ) : ℕ :=
  if d ≤ alpha G S then icnt G S (alpha G S - d) else 0

lemma eD_add_bD (G : SimpleGraph V) (S : Finset V) (d : ℕ) :
    eD G S d + bD G S d = iD G S d := by
  unfold eD bD iD
  split
  · exact ecnt_add_bcnt _ _ _
  · rfl

lemma icnt_eq_iD {d : ℕ} (h : d ≤ alpha G S) : icnt G S (alpha G S - d) = iD G S d :=
  (if_pos h).symm

lemma iD_zero (G : SimpleGraph V) (S : Finset V) : iD G S 0 = eD G S 0 := by
  have := eD_add_bD G S 0
  rw [bD_zero] at this
  omega

lemma iD_one (G : SimpleGraph V) (S : Finset V) : iD G S 1 = eD G S 1 + bD G S 1 :=
  (eD_add_bD G S 1).symm

namespace Split

lemma subset_left (h : Split G S S₁ S₂) : S₁ ⊆ S := by
  rw [h.eq_union]; exact Finset.subset_union_left

lemma subset_right (h : Split G S S₁ S₂) : S₂ ⊆ S := by
  rw [h.eq_union]; exact Finset.subset_union_right

lemma inter_left_mem (_h : Split G S S₁ S₂) (hA : A ∈ indFam G S) : A ∩ S₁ ∈ indFam G S₁ := by
  rw [mem_indFam] at hA ⊢
  exact ⟨Finset.inter_subset_right, indep_of_subset hA.2 Finset.inter_subset_left⟩

lemma inter_right_mem (_h : Split G S S₁ S₂) (hA : A ∈ indFam G S) : A ∩ S₂ ∈ indFam G S₂ := by
  rw [mem_indFam] at hA ⊢
  exact ⟨Finset.inter_subset_right, indep_of_subset hA.2 Finset.inter_subset_left⟩

lemma union_mem (h : Split G S S₁ S₂) (h1 : A ∈ indFam G S₁) (h2 : B ∈ indFam G S₂) :
    A ∪ B ∈ indFam G S := by
  rw [mem_indFam] at h1 h2 ⊢
  refine ⟨by rw [h.eq_union]; exact Finset.union_subset_union h1.1 h2.1, ?_⟩
  intro x hx y hy hxy
  simp only [Finset.coe_union, Set.mem_union, Finset.mem_coe] at hx hy
  rcases hx with hx | hx <;> rcases hy with hy | hy
  · exact h1.2 (by exact_mod_cast hx) (by exact_mod_cast hy) hxy
  · exact h.noEdges x (h1.1 hx) y (h2.1 hy)
  · exact fun hadj => h.noEdges y (h1.1 hy) x (h2.1 hx) hadj.symm
  · exact h2.2 (by exact_mod_cast hx) (by exact_mod_cast hy) hxy

lemma inter_union_self (h : Split G S S₁ S₂) (hA : A ⊆ S) : A ∩ S₁ ∪ A ∩ S₂ = A := by
  rw [← Finset.inter_union_distrib_left, ← h.eq_union]
  exact Finset.inter_eq_left.2 hA

lemma union_inter_left (h : Split G S S₁ S₂) (h1 : A ⊆ S₁) (h2 : B ⊆ S₂) : (A ∪ B) ∩ S₁ = A := by
  ext x
  simp only [Finset.mem_inter, Finset.mem_union]
  constructor
  · rintro ⟨hx | hx, hs⟩
    · exact hx
    · exact absurd (Finset.mem_inter.2 ⟨hs, h2 hx⟩) (by simp [Finset.disjoint_left.1 h.disj hs])
  · intro hx
    exact ⟨Or.inl hx, h1 hx⟩

lemma union_inter_right (h : Split G S S₁ S₂) (h1 : A ⊆ S₁) (h2 : B ⊆ S₂) : (A ∪ B) ∩ S₂ = B := by
  ext x
  simp only [Finset.mem_inter, Finset.mem_union]
  constructor
  · rintro ⟨hx | hx, hs⟩
    · exact absurd (Finset.disjoint_left.1 h.disj (h1 hx)) (by simp [hs])
    · exact hx
  · intro hx
    exact ⟨Or.inr hx, h2 hx⟩

lemma card_inter_add (h : Split G S S₁ S₂) (hA : A ⊆ S) :
    (A ∩ S₁).card + (A ∩ S₂).card = A.card := by
  rw [← Finset.card_union_of_disjoint, h.inter_union_self hA]
  refine Finset.disjoint_left.2 fun x hx1 hx2 => ?_
  rw [Finset.mem_inter] at hx1 hx2
  exact Finset.disjoint_left.1 h.disj hx1.2 hx2.2

lemma card_union (h : Split G S S₁ S₂) (h1 : A ⊆ S₁) (h2 : B ⊆ S₂) :
    (A ∪ B).card = A.card + B.card :=
  Finset.card_union_of_disjoint (Finset.disjoint_of_subset_left h1
    (Finset.disjoint_of_subset_right h2 h.disj))

end Split

lemma alpha_split (h : Split G S S₁ S₂) : alpha G S = alpha G S₁ + alpha G S₂ := by
  refine le_antisymm ?_ ?_
  · obtain ⟨M, hM, hcard⟩ := exists_max_indFam G S
    have h1 := card_le_alpha (h.inter_left_mem hM)
    have h2 := card_le_alpha (h.inter_right_mem hM)
    have h3 := h.card_inter_add (mem_indFam.1 hM).1
    omega
  · obtain ⟨M₁, hM₁, hc₁⟩ := exists_max_indFam G S₁
    obtain ⟨M₂, hM₂, hc₂⟩ := exists_max_indFam G S₂
    have := card_le_alpha (h.union_mem hM₁ hM₂)
    rwa [h.card_union (mem_indFam.1 hM₁).1 (mem_indFam.1 hM₂).1, hc₁, hc₂] at this

/-- The general convolution lemma: a property that factors through the two parts gives
a convolution identity for the corresponding layer counts. -/
lemma card_split_conv (h : Split G S S₁ S₂) (Q Q₁ Q₂ : Finset V → Prop)
    (hQ : ∀ A ∈ indFam G S, (Q A ↔ Q₁ (A ∩ S₁) ∧ Q₂ (A ∩ S₂)))
    (hQ' : ∀ A ∈ indFam G S₁, ∀ B ∈ indFam G S₂, Q₁ A → Q₂ B → Q (A ∪ B)) (k : ℕ) :
    ((indFam G S).filter (fun A => A.card = k ∧ Q A)).card
      = ∑ i ∈ range (k + 1), ((indFam G S₁).filter (fun A => A.card = i ∧ Q₁ A)).card
          * ((indFam G S₂).filter (fun A => A.card = k - i ∧ Q₂ A)).card := by
  classical
  rw [Finset.card_eq_sum_card_fiberwise
    (f := fun A => (A ∩ S₁).card) (t := range (k + 1)) (fun A hA => ?_)]
  · refine Finset.sum_congr rfl fun i hi => ?_
    rw [Finset.mem_range] at hi
    rw [← Finset.card_product]
    refine Finset.card_bij' (fun A _ => (A ∩ S₁, A ∩ S₂)) (fun q _ => q.1 ∪ q.2) ?_ ?_ ?_ ?_
    · intro A hA
      simp only [Finset.mem_filter, Finset.mem_product] at hA ⊢
      obtain ⟨⟨hA1, hcard, hQA⟩, hfib⟩ := hA
      have hsub := (mem_indFam.1 hA1).1
      have hadd := h.card_inter_add hsub
      exact ⟨⟨h.inter_left_mem hA1, hfib, ((hQ A hA1).1 hQA).1⟩,
        ⟨h.inter_right_mem hA1, by omega, ((hQ A hA1).1 hQA).2⟩⟩
    · rintro ⟨A, B⟩ hq
      simp only [Finset.mem_product, Finset.mem_filter] at hq
      obtain ⟨⟨hA, hcA, hQA⟩, ⟨hB, hcB, hQB⟩⟩ := hq
      have hAs := (mem_indFam.1 hA).1
      have hBs := (mem_indFam.1 hB).1
      simp only [Finset.mem_filter]
      refine ⟨⟨h.union_mem hA hB, ?_, hQ' A hA B hB hQA hQB⟩, ?_⟩
      · rw [h.card_union hAs hBs, hcA, hcB]; omega
      · rw [h.union_inter_left hAs hBs, hcA]
    · intro A hA
      rw [Finset.mem_filter, Finset.mem_filter] at hA
      exact h.inter_union_self (mem_indFam.1 hA.1.1).1
    · rintro ⟨A, B⟩ hq
      simp only [Finset.mem_product, Finset.mem_filter] at hq
      have hAs := (mem_indFam.1 hq.1.1).1
      have hBs := (mem_indFam.1 hq.2.1).1
      rw [Prod.mk.injEq]
      exact ⟨h.union_inter_left hAs hBs, h.union_inter_right hAs hBs⟩
  · simp only [Finset.mem_coe, Finset.mem_filter] at hA
    have hsub := (mem_indFam.1 hA.1).1
    have := h.card_inter_add hsub
    simp only [Finset.mem_coe, Finset.mem_range]
    omega

lemma icnt_split (h : Split G S S₁ S₂) (k : ℕ) :
    icnt G S k = ∑ i ∈ range (k + 1), icnt G S₁ i * icnt G S₂ (k - i) := by
  have := card_split_conv h (fun _ => True) (fun _ => True) (fun _ => True)
    (fun A _ => by tauto) (fun A _ B _ _ _ => trivial) k
  simpa [icnt] using this

/-- Extendability in a disjoint union factors through the two parts. -/
lemma ext_split_iff (h : Split G S S₁ S₂) (hA : A ∈ indFam G S) :
    Ext G S A ↔ Ext G S₁ (A ∩ S₁) ∧ Ext G S₂ (A ∩ S₂) := by
  constructor
  · rintro ⟨M, hM, hcard, hAM⟩
    have h1 := card_le_alpha (h.inter_left_mem hM)
    have h2 := card_le_alpha (h.inter_right_mem hM)
    have h3 := h.card_inter_add (mem_indFam.1 hM).1
    have h4 := alpha_split h
    constructor
    · exact ⟨M ∩ S₁, h.inter_left_mem hM, by omega,
        Finset.inter_subset_inter hAM (Finset.Subset.refl _)⟩
    · exact ⟨M ∩ S₂, h.inter_right_mem hM, by omega,
        Finset.inter_subset_inter hAM (Finset.Subset.refl _)⟩
  · rintro ⟨⟨M₁, hM₁, hc₁, h₁⟩, ⟨M₂, hM₂, hc₂, h₂⟩⟩
    refine ⟨M₁ ∪ M₂, h.union_mem hM₁ hM₂, ?_, ?_⟩
    · rw [h.card_union (mem_indFam.1 hM₁).1 (mem_indFam.1 hM₂).1, hc₁, hc₂, alpha_split h]
    · calc A = A ∩ S₁ ∪ A ∩ S₂ := (h.inter_union_self (mem_indFam.1 hA).1).symm
        _ ⊆ M₁ ∪ M₂ := Finset.union_subset_union h₁ h₂

/-- Extendable sets of the two parts combine to an extendable set of the union. -/
lemma ext_union_of_split (h : Split G S S₁ S₂) (hEA : Ext G S₁ A) (hEB : Ext G S₂ B) : Ext G S (A ∪ B) := by
  obtain ⟨M₁, hM₁, hc₁, h₁⟩ := hEA
  obtain ⟨M₂, hM₂, hc₂, h₂⟩ := hEB
  refine ⟨M₁ ∪ M₂, h.union_mem hM₁ hM₂, ?_, Finset.union_subset_union h₁ h₂⟩
  rw [h.card_union (mem_indFam.1 hM₁).1 (mem_indFam.1 hM₂).1, hc₁, hc₂, alpha_split h]

lemma ecnt_split (h : Split G S S₁ S₂) (k : ℕ) :
    ecnt G S k = ∑ i ∈ range (k + 1), ecnt G S₁ i * ecnt G S₂ (k - i) := by
  have := card_split_conv h (Ext G S) (Ext G S₁) (Ext G S₂)
    (fun A hA => ext_split_iff h hA)
    (fun _ _ _ _ hEA hEB => ext_union_of_split h hEA hEB) k
  simpa [ecnt] using this

/-- Reindexing a size-graded convolution into a defect-graded one. -/
private lemma conv_reindex (a₁ a₂ d : ℕ) (c₁ c₂ : ℕ → ℕ)
    (h₁ : ∀ k, a₁ < k → c₁ k = 0) (h₂ : ∀ k, a₂ < k → c₂ k = 0) (hd : d ≤ a₁ + a₂) :
    ∑ j ∈ range (a₁ + a₂ - d + 1), c₁ j * c₂ (a₁ + a₂ - d - j)
      = ∑ i ∈ range (d + 1), (if i ≤ a₁ then c₁ (a₁ - i) else 0)
          * (if d - i ≤ a₂ then c₂ (a₂ - (d - i)) else 0) := by
  classical
  set n := a₁ + a₂ - d with hn
  rw [← Finset.sum_subset (s₁ := (range (n + 1)).filter (fun j => j ≤ a₁ ∧ a₁ ≤ d + j))
      (Finset.filter_subset _ _), ← Finset.sum_subset
      (s₁ := (range (d + 1)).filter (fun i => i ≤ a₁ ∧ d - i ≤ a₂)) (Finset.filter_subset _ _)]
  · refine Finset.sum_nbij' (fun j => a₁ - j) (fun i => a₁ - i) ?_ ?_ ?_ ?_ ?_
    · intro j hj
      simp only [Finset.mem_filter, Finset.mem_range] at hj ⊢
      omega
    · intro i hi
      simp only [Finset.mem_filter, Finset.mem_range] at hi ⊢
      omega
    · intro j hj
      simp only [Finset.mem_filter, Finset.mem_range] at hj
      show a₁ - (a₁ - j) = j
      omega
    · intro i hi
      simp only [Finset.mem_filter, Finset.mem_range] at hi
      show a₁ - (a₁ - i) = i
      omega
    · intro j hj
      simp only [Finset.mem_filter, Finset.mem_range] at hj
      have e1 : a₁ - (a₁ - j) = j := by omega
      have e2 : d - (a₁ - j) ≤ a₂ := by omega
      have e3 : a₂ - (d - (a₁ - j)) = n - j := by omega
      show _ = (if a₁ - j ≤ a₁ then c₁ (a₁ - (a₁ - j)) else 0)
          * (if d - (a₁ - j) ≤ a₂ then c₂ (a₂ - (d - (a₁ - j))) else 0)
      rw [if_pos (by omega : a₁ - j ≤ a₁), if_pos e2, e1, e3]
  · intro i hi hni
    simp only [Finset.mem_filter, Finset.mem_range] at hi hni
    by_cases hle : i ≤ a₁
    · rw [if_neg (by omega : ¬ (d - i ≤ a₂)), mul_zero]
    · rw [if_neg hle, zero_mul]
  · intro j hj hnj
    simp only [Finset.mem_filter, Finset.mem_range] at hj hnj
    by_cases hle : j ≤ a₁
    · rw [h₂ _ (by omega), mul_zero]
    · rw [h₁ _ (by omega), zero_mul]

/-- The defect-graded form of a size-graded convolution. -/
private lemma conv_defect (a₁ a₂ d : ℕ) (c₁ c₂ : ℕ → ℕ)
    (h₁ : ∀ k, a₁ < k → c₁ k = 0) (h₂ : ∀ k, a₂ < k → c₂ k = 0) :
    (if d ≤ a₁ + a₂ then ∑ j ∈ range (a₁ + a₂ - d + 1), c₁ j * c₂ (a₁ + a₂ - d - j) else 0)
      = ∑ i ∈ range (d + 1), (if i ≤ a₁ then c₁ (a₁ - i) else 0)
          * (if d - i ≤ a₂ then c₂ (a₂ - (d - i)) else 0) := by
  rcases le_or_gt d (a₁ + a₂) with hd | hd
  · rw [if_pos hd]; exact conv_reindex a₁ a₂ d c₁ c₂ h₁ h₂ hd
  · rw [if_neg (by omega)]
    refine (Finset.sum_eq_zero fun i hi => ?_).symm
    rw [Finset.mem_range] at hi
    by_cases hle : i ≤ a₁
    · rw [if_neg (by omega : ¬ (d - i ≤ a₂)), mul_zero]
    · rw [if_neg hle, zero_mul]

lemma iD_split (h : Split G S S₁ S₂) (d : ℕ) :
    iD G S d = ∑ i ∈ range (d + 1), iD G S₁ i * iD G S₂ (d - i) := by
  have ha := alpha_split h
  simp only [iD, ha]
  rw [icnt_split h]
  exact conv_defect (alpha G S₁) (alpha G S₂) d (icnt G S₁) (icnt G S₂)
    (fun _ hk => icnt_eq_zero_of_gt hk) (fun _ hk => icnt_eq_zero_of_gt hk)

lemma eD_split (h : Split G S S₁ S₂) (d : ℕ) :
    eD G S d = ∑ i ∈ range (d + 1), eD G S₁ i * eD G S₂ (d - i) := by
  have ha := alpha_split h
  simp only [eD, ha]
  rw [ecnt_split h]
  exact conv_defect (alpha G S₁) (alpha G S₂) d (ecnt G S₁) (ecnt G S₂)
    (fun _ hk => ecnt_eq_zero_of_gt hk) (fun _ hk => ecnt_eq_zero_of_gt hk)

/-- Blocked layers of a part inject into blocked layers of a disjoint union. -/
lemma bD_le_of_split (h : Split G S S₁ S₂) (d : ℕ) : bD G S₁ d ≤ bD G S d := by
  classical
  rcases le_or_gt d (alpha G S₁) with hd | hd
  · obtain ⟨M₂, hM₂, hc₂⟩ := exists_max_indFam G S₂
    have hM₂s := (mem_indFam.1 hM₂).1
    have ha := alpha_split h
    rw [bD, bD, if_pos hd, if_pos (by omega)]
    unfold bcnt
    refine Finset.card_le_card_of_injOn (fun A => A ∪ M₂) ?_ ?_
    · intro A hA
      simp only [Finset.mem_coe, Finset.mem_filter] at hA ⊢
      obtain ⟨hA1, hcA, hnA⟩ := hA
      have hAs := (mem_indFam.1 hA1).1
      refine ⟨h.union_mem hA1 hM₂, ?_, ?_⟩
      · rw [h.card_union hAs hM₂s, hcA, hc₂]; omega
      · intro hext
        exact hnA (by
          have := ((ext_split_iff h (h.union_mem hA1 hM₂)).1 hext).1
          rwa [h.union_inter_left hAs hM₂s] at this)
    · intro A hA B hB hAB
      simp only [Finset.coe_filter, Set.mem_setOf_eq] at hA hB
      have hAs := (mem_indFam.1 hA.1).1
      have hBs := (mem_indFam.1 hB.1).1
      have := congrArg (fun T => T ∩ S₁) hAB
      simpa [h.union_inter_left hAs hM₂s, h.union_inter_left hBs hM₂s] using this
  · rw [bD, if_neg (by omega)]
    exact Nat.zero_le _

end FivefoldForest
