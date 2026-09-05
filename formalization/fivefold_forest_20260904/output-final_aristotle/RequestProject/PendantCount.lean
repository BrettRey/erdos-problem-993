import RequestProject.Pendant

/-!
# Layer identities for pendant configurations

The independent sets of `S` split into those avoiding the support vertex `p` — which are
exactly the independent sets of the disjoint union `pRest S p = pBase S p L ⊔ L` — and
those containing `p`, which are `{p} ∪ X` for `X` an independent set of `pCore S p w L`.
-/

open Finset

namespace FivefoldForest

open scoped Classical

variable {V : Type*} {G : SimpleGraph V} {S : Finset V} {p w : V} {L X : Finset V}

namespace Pendant

lemma p_notMem_base (S : Finset V) (p : V) (L : Finset V) : p ∉ pBase S p L := by
  simp [pBase]

lemma p_notMem_core (S : Finset V) (p w : V) (L : Finset V) : p ∉ pCore S p w L := fun h =>
  p_notMem_base S p L (Finset.mem_of_mem_erase h)

lemma core_subset (S : Finset V) (p w : V) (L : Finset V) : pCore S p w L ⊆ S :=
  (Finset.erase_subset _ _).trans (base_subset S p L)

lemma core_not_adj_p (hp : Pendant G S p w L) {x : V} (hx : x ∈ pCore S p w L) :
    ¬ G.Adj p x := by
  intro hadj
  simp only [pCore, pBase, Finset.mem_erase, Finset.mem_sdiff, Finset.mem_insert] at hx
  rcases hp.pNbrs x hx.2.1 hadj with rfl | hxL
  · exact hx.1 rfl
  · exact hx.2.2 (Or.inr hxL)

/-- Independent sets containing the support vertex are `{p}` together with an independent
set of the core. -/
lemma insert_p_mem_iff (hp : Pendant G S p w L) (hpX : p ∉ X) :
    insert p X ∈ indFam G S ↔ X ∈ indFam G (pCore S p w L) := by
  constructor
  · intro hA
    rw [mem_indFam] at hA ⊢
    have hXsub : X ⊆ S := fun x hx => hA.1 (Finset.mem_insert_of_mem hx)
    refine ⟨fun x hx => ?_, indep_of_subset hA.2 (Finset.subset_insert _ _)⟩
    have hxp : x ≠ p := fun h => hpX (h ▸ hx)
    have hpmem : p ∈ (↑(insert p X) : Set V) := by simp
    have hxmem : x ∈ (↑(insert p X) : Set V) := by simp [hx]
    have hnadj : ¬ G.Adj p x := hA.2 hpmem hxmem (Ne.symm hxp)
    simp only [pCore, pBase, Finset.mem_erase, Finset.mem_sdiff, Finset.mem_insert]
    refine ⟨?_, hXsub hx, ?_⟩
    · rintro rfl; exact hnadj hp.adj
    · rintro (rfl | hxL)
      · exact hxp rfl
      · exact hnadj (hp.leafAdj x hxL)
  · intro hX
    rw [mem_indFam] at hX ⊢
    constructor
    · intro x hx
      rcases Finset.mem_insert.1 hx with rfl | hx
      · exact hp.pmem
      · exact core_subset S p w L (hX.1 hx)
    · intro x hx y hy hxy hadj
      simp only [Finset.coe_insert, Set.mem_insert_iff, Finset.mem_coe] at hx hy
      rcases hx with rfl | hx <;> rcases hy with rfl | hy
      · exact hxy rfl
      · exact hp.core_not_adj_p (hX.1 hy) hadj
      · exact hp.core_not_adj_p (hX.1 hx) hadj.symm
      · exact hX.2 (by exact_mod_cast hx) (by exact_mod_cast hy) hxy hadj

lemma ext_of_mem_rest_iff (hp : Pendant G S p w L) {A : Finset V}
    (hA : A ∈ indFam G (pRest S p)) : Ext G S A ↔ Ext G (pRest S p) A :=
  hp.ext_rest_iff (mem_indFam.1 hA).1

/-- Splitting a layer count of `S` according to whether the support vertex is used. -/
lemma card_split_p (hp : Pendant G S p w L) (P : Finset V → Prop) (k : ℕ) :
    ((indFam G S).filter fun A => A.card = k ∧ P A).card
      = ((indFam G (pRest S p)).filter fun A => A.card = k ∧ P A).card
        + ((indFam G (pCore S p w L)).filter fun Y => Y.card + 1 = k ∧ P (insert p Y)).card := by
  classical
  have h1 : (((indFam G S).filter fun A => A.card = k ∧ P A).filter fun A => ¬ p ∈ A)
      = (indFam G (pRest S p)).filter fun A => A.card = k ∧ P A := by
    ext A
    simp only [Finset.mem_filter, mem_indFam, pRest, Finset.subset_erase]
    tauto
  have h2 : ((((indFam G S).filter fun A => A.card = k ∧ P A).filter fun A => p ∈ A)).card
      = ((indFam G (pCore S p w L)).filter fun Y => Y.card + 1 = k ∧ P (insert p Y)).card := by
    refine Finset.card_bij' (fun A _ => A.erase p) (fun Y _ => insert p Y) ?_ ?_ ?_ ?_
    · intro A hA
      simp only [Finset.mem_filter] at hA ⊢
      obtain ⟨⟨hA1, hcard, hPA⟩, hpA⟩ := hA
      have hins : insert p (A.erase p) = A := Finset.insert_erase hpA
      refine ⟨(hp.insert_p_mem_iff (Finset.notMem_erase p A)).1 (by rwa [hins]), ?_, ?_⟩
      · rw [Finset.card_erase_of_mem hpA]
        have : 1 ≤ A.card := Finset.card_pos.2 ⟨p, hpA⟩
        omega
      · rwa [hins]
    · intro Y hY
      simp only [Finset.mem_filter] at hY ⊢
      obtain ⟨hY1, hcard, hPY⟩ := hY
      have hpY : p ∉ Y := fun h => p_notMem_core S p w L ((mem_indFam.1 hY1).1 h)
      exact ⟨⟨(hp.insert_p_mem_iff hpY).2 hY1,
        by rw [Finset.card_insert_of_notMem hpY]; exact hcard, hPY⟩,
        Finset.mem_insert_self _ _⟩
    · intro A hA
      simp only [Finset.mem_filter] at hA
      exact Finset.insert_erase hA.2
    · intro Y hY
      simp only [Finset.mem_filter] at hY
      exact Finset.erase_insert (fun h => p_notMem_core S p w L ((mem_indFam.1 hY.1).1 h))
  rw [← Finset.card_filter_add_card_filter_not
    (s := (indFam G S).filter fun A => A.card = k ∧ P A) (p := fun A => p ∈ A), h2, h1,
    Nat.add_comm]

/-- The negated form of `card_split_p`, stated so that it can be rewritten with. -/
lemma card_split_p_not (hp : Pendant G S p w L) (P : Finset V → Prop) (k : ℕ) :
    ((indFam G S).filter fun A => A.card = k ∧ ¬ P A).card
      = ((indFam G (pRest S p)).filter fun A => A.card = k ∧ ¬ P A).card
        + ((indFam G (pCore S p w L)).filter
            fun Y => Y.card + 1 = k ∧ ¬ P (insert p Y)).card := by
  convert hp.card_split_p (fun A => ¬ P A) k using 3 <;> congr 1

/-- With at least two leaves, no maximum independent set contains the support vertex. -/
lemma no_p_in_max (hp : Pendant G S p w L) (hm : 2 ≤ L.card) {M : Finset V}
    (hM : M ∈ indFam G S) (hcard : M.card = alpha G S) : p ∉ M := by
  intro hpM
  obtain ⟨h1, h2⟩ := hp.lift_of_mem hM hpM
  have h3 := card_le_alpha h1
  rw [hp.alpha_rest] at h3
  omega

lemma not_ext_insert_p (hp : Pendant G S p w L) (hm : 2 ≤ L.card) (Y : Finset V) :
    ¬ Ext G S (insert p Y) := by
  rintro ⟨M, hM, hcard, hsub⟩
  exact hp.no_p_in_max hm hM hcard (hsub (Finset.mem_insert_self p Y))

section TwoLe

/-- With at least two leaves, the extendable layers ignore the support vertex. -/
lemma eD_two_le (hp : Pendant G S p w L) (hm : 2 ≤ L.card) (d : ℕ) :
    eD G S d = eD G (pRest S p) d := by
  have hcnt : ∀ k, ecnt G S k = ecnt G (pRest S p) k := by
    intro k
    rw [ecnt, hp.card_split_p (Ext G S) k, ecnt]
    have hz : ((indFam G (pCore S p w L)).filter
        fun Y => Y.card + 1 = k ∧ Ext G S (insert p Y)) = ∅ :=
      Finset.eq_empty_of_forall_notMem fun Y hY =>
        hp.not_ext_insert_p hm Y (Finset.mem_filter.1 hY).2.2
    rw [hz, Finset.card_empty, add_zero]
    exact congrArg Finset.card (Finset.filter_congr fun A hA => by
      rw [hp.ext_of_mem_rest_iff hA])
  simp only [eD, hp.alpha_rest, hcnt]

lemma bD_two_le (hp : Pendant G S p w L) (hm : 2 ≤ L.card) (d : ℕ) (hd : d + 1 ≤ alpha G S) :
    bD G S d = bD G (pRest S p) d + icnt G (pCore S p w L) (alpha G S - d - 1) := by
  have hcnt : ∀ k, bcnt G S k = bcnt G (pRest S p) k
      + ((indFam G (pCore S p w L)).filter fun Y => Y.card + 1 = k).card := by
    intro k
    rw [bcnt, hp.card_split_p_not (Ext G S) k, bcnt]
    congr 1
    · exact congrArg Finset.card (Finset.filter_congr fun A hA => by
        rw [hp.ext_of_mem_rest_iff hA])
    · exact congrArg Finset.card (Finset.filter_congr fun Y _ => by
        simp only [and_iff_left (hp.not_ext_insert_p hm Y)])
  simp only [bD, hp.alpha_rest]
  split_ifs with hs
  · rw [hcnt, icnt]
    congr 1
    exact congrArg Finset.card (Finset.filter_congr fun Y _ => by
      constructor
      · intro h1; omega
      · intro h1; omega)
  · omega

end TwoLe

section One

lemma alpha_base_add (hp : Pendant G S p w L) :
    alpha G S = alpha G (pBase S p L) + L.card := by
  rw [← hp.alpha_rest, alpha_split hp.split, hp.leaves_alpha]

/-- For a pendant `P₂` (one leaf), eligibility forces a maximum independent set of the base
to avoid `w`. -/
lemma alpha_core_one (hp : Pendant G S p w L) (hm : L.card = 1) (helig : Eligible G S) :
    alpha G (pCore S p w L) = alpha G (pBase S p L) := by
  by_contra hne
  have h1 : alpha G (pCore S p w L) ≤ alpha G (pBase S p L) := alpha_erase_le _ _ _
  have h2 : alpha G (pBase S p L) ≤ alpha G (pCore S p w L) + 1 := alpha_erase_succ_ge _ _ _
  have hS : alpha G S = alpha G (pBase S p L) + 1 := by rw [hp.alpha_base_add, hm]
  obtain ⟨Y, hY, hYc⟩ := exists_max_indFam G (pCore S p w L)
  have hpY : p ∉ Y := fun h => p_notMem_core S p w L ((mem_indFam.1 hY).1 h)
  have hmem : insert p Y ∈ indFam G S := (hp.insert_p_mem_iff hpY).2 hY
  have hcard : (insert p Y).card = alpha G S - 1 := by
    rw [Finset.card_insert_of_notMem hpY, hYc]; omega
  have hnotext : ¬ Ext G S (insert p Y) := by
    rintro ⟨M, hM, hMc, hsub⟩
    have hpM : p ∈ M := hsub (Finset.mem_insert_self _ _)
    have hEr : M.erase p ∈ indFam G (pCore S p w L) :=
      (hp.insert_p_mem_iff (Finset.notMem_erase _ _)).1 (by rwa [Finset.insert_erase hpM])
    have h4 := card_le_alpha hEr
    rw [Finset.card_erase_of_mem hpM] at h4
    omega
  have hpos : 1 ≤ bD G S 1 := by
    rw [bD, if_pos (by omega), bcnt]
    refine Finset.card_pos.2 ⟨insert p Y, ?_⟩
    simp only [Finset.mem_filter]
    exact ⟨hmem, hcard, hnotext⟩
  rw [helig.2.1] at hpos
  omega

/-- With one leaf and a core of full independence number, `{p} ∪ X` is extendable exactly
when `X` is extendable in the core. -/
lemma ext_insert_p_iff (hp : Pendant G S p w L)
    (halpha : alpha G (pCore S p w L) + 1 = alpha G S) {Y : Finset V}
    (hY : Y ∈ indFam G (pCore S p w L)) :
    Ext G S (insert p Y) ↔ Ext G (pCore S p w L) Y := by
  constructor
  · rintro ⟨M, hM, hMc, hsub⟩
    have hpM : p ∈ M := hsub (Finset.mem_insert_self _ _)
    have hEr : M.erase p ∈ indFam G (pCore S p w L) :=
      (hp.insert_p_mem_iff (Finset.notMem_erase _ _)).1 (by rwa [Finset.insert_erase hpM])
    refine ⟨M.erase p, hEr, ?_, ?_⟩
    · rw [Finset.card_erase_of_mem hpM, hMc]; omega
    · intro x hx
      exact Finset.mem_erase.2 ⟨fun h => p_notMem_core S p w L (h ▸ (mem_indFam.1 hY).1 hx),
        hsub (Finset.mem_insert_of_mem hx)⟩
  · rintro ⟨Z, hZ, hZc, hsub⟩
    have hpZ : p ∉ Z := fun h => p_notMem_core S p w L ((mem_indFam.1 hZ).1 h)
    refine ⟨insert p Z, (hp.insert_p_mem_iff hpZ).2 hZ, ?_, Finset.insert_subset_insert _ hsub⟩
    rw [Finset.card_insert_of_notMem hpZ, hZc]; omega

lemma eD_one (hp : Pendant G S p w L) (hm : L.card = 1)
    (hF : alpha G (pCore S p w L) = alpha G (pBase S p L)) (d : ℕ) :
    eD G S d = eD G (pRest S p) d + eD G (pCore S p w L) d := by
  have halpha : alpha G (pCore S p w L) + 1 = alpha G S := by
    rw [hF, hp.alpha_base_add, hm]
  have hcnt : ∀ k, ecnt G S k = ecnt G (pRest S p) k
      + ((indFam G (pCore S p w L)).filter
          fun Y => Y.card + 1 = k ∧ Ext G (pCore S p w L) Y).card := by
    intro k
    rw [ecnt, hp.card_split_p (Ext G S) k, ecnt]
    congr 1
    · exact congrArg Finset.card (Finset.filter_congr fun A hA => by
        rw [hp.ext_of_mem_rest_iff hA])
    · exact congrArg Finset.card (Finset.filter_congr fun Y hY => by
        rw [hp.ext_insert_p_iff halpha hY])
  simp only [eD, hp.alpha_rest]
  split_ifs with hs hc
  · rw [hcnt, ecnt]
    congr 1
    exact congrArg Finset.card (Finset.filter_congr fun Y _ => by
      constructor
      · rintro ⟨h1, h2⟩; exact ⟨by omega, h2⟩
      · rintro ⟨h1, h2⟩; exact ⟨by omega, h2⟩)
  · rw [hcnt]
    congr 1
    refine Finset.card_eq_zero.2 (Finset.eq_empty_of_forall_notMem fun Y hY => ?_)
    have := (Finset.mem_filter.1 hY).2.1
    omega
  · omega
  · omega

lemma bD_one (hp : Pendant G S p w L) (hm : L.card = 1)
    (hF : alpha G (pCore S p w L) = alpha G (pBase S p L)) (d : ℕ) :
    bD G S d = bD G (pRest S p) d + bD G (pCore S p w L) d := by
  have halpha : alpha G (pCore S p w L) + 1 = alpha G S := by
    rw [hF, hp.alpha_base_add, hm]
  have hcnt : ∀ k, bcnt G S k = bcnt G (pRest S p) k
      + ((indFam G (pCore S p w L)).filter
          fun Y => Y.card + 1 = k ∧ ¬ Ext G (pCore S p w L) Y).card := by
    intro k
    rw [bcnt, hp.card_split_p_not (Ext G S) k, bcnt]
    congr 1
    · exact congrArg Finset.card (Finset.filter_congr fun A hA => by
        rw [hp.ext_of_mem_rest_iff hA])
    · exact congrArg Finset.card (Finset.filter_congr fun Y hY => by
        rw [hp.ext_insert_p_iff halpha hY])
  simp only [bD, hp.alpha_rest]
  split_ifs with hs hc
  · rw [hcnt, bcnt]
    congr 1
    exact congrArg Finset.card (Finset.filter_congr fun Y _ => by
      constructor
      · rintro ⟨h1, h2⟩; exact ⟨by omega, h2⟩
      · rintro ⟨h1, h2⟩; exact ⟨by omega, h2⟩)
  · rw [hcnt]
    congr 1
    refine Finset.card_eq_zero.2 (Finset.eq_empty_of_forall_notMem fun Y hY => ?_)
    have := (Finset.mem_filter.1 hY).2.1
    omega
  · omega
  · omega

end One

end Pendant

end FivefoldForest
