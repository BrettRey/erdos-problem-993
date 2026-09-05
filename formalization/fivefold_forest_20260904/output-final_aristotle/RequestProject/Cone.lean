import RequestProject.Basic

/-!
# Cone identities

If every maximum independent set of `S` contains the vertex `w`, the extendable layers
of `S` are built from those of `S.erase w`.
-/

open Finset

namespace FivefoldForest

open scoped Classical

variable {V : Type*} {G : SimpleGraph V} {S : Finset V} {w : V}

/-- If there are no blocked sets at defect `d`, every independent set of that size extends. -/
lemma ext_of_bD_zero {d : ℕ} (hd : d ≤ alpha G S) (h : bD G S d = 0) {A : Finset V}
    (hA : A ∈ indFam G S) (hc : A.card = alpha G S - d) : Ext G S A := by
  by_contra hn
  rw [bD, if_pos hd, bcnt, Finset.card_eq_zero] at h
  have hmem : A ∈ ((indFam G S).filter fun B => B.card = alpha G S - d ∧ ¬ Ext G S B) :=
    Finset.mem_filter.2 ⟨hA, hc, hn⟩
  rw [h] at hmem
  exact absurd hmem (Finset.notMem_empty _)

/-- Deleting `w` from a maximum independent set of `S` gives a maximum independent set
of `S.erase w`. -/
lemma cone_max_erase (hall : ∀ M ∈ indFam G S, M.card = alpha G S → w ∈ M)
    (halpha : alpha G (S.erase w) + 1 = alpha G S) {M : Finset V} (hM : M ∈ indFam G S)
    (hc : M.card = alpha G S) :
    M.erase w ∈ indFam G (S.erase w) ∧ (M.erase w).card = alpha G (S.erase w) := by
  have hwM := hall M hM hc
  rw [mem_indFam] at hM
  refine ⟨mem_indFam.2 ⟨fun x hx => ?_, indep_of_subset hM.2 (Finset.erase_subset _ _)⟩, ?_⟩
  · exact Finset.mem_erase.2 ⟨(Finset.mem_erase.1 hx).1, hM.1 (Finset.mem_erase.1 hx).2⟩
  · rw [Finset.card_erase_of_mem hwM, hc]; omega

/-- Adding `w` to a maximum independent set of `S.erase w` gives a maximum independent set
of `S`. -/
lemma cone_max_insert (hall : ∀ M ∈ indFam G S, M.card = alpha G S → w ∈ M)
    (hb1 : bD G S 1 = 0) (halpha : alpha G (S.erase w) + 1 = alpha G S) {N : Finset V}
    (hN : N ∈ indFam G (S.erase w)) (hc : N.card = alpha G (S.erase w)) :
    insert w N ∈ indFam G S ∧ (insert w N).card = alpha G S := by
  have hNS : N ∈ indFam G S := indFam_mono (Finset.erase_subset _ _) hN
  have hwN : w ∉ N := fun h => (Finset.mem_erase.1 ((mem_indFam.1 hN).1 h)).1 rfl
  obtain ⟨M, hM, hMc, hNM⟩ := ext_of_bD_zero (by omega) hb1 hNS (by omega)
  have hwM := hall M hM hMc
  obtain ⟨h1, h2⟩ := cone_max_erase hall halpha hM hMc
  have hsub : N ⊆ M.erase w := fun x hx => Finset.mem_erase.2 ⟨fun h => hwN (h ▸ hx), hNM hx⟩
  have hNe : N = M.erase w := Finset.eq_of_subset_of_card_le hsub (by omega)
  rw [hNe, Finset.insert_erase hwM]
  exact ⟨hM, hMc⟩

/-- Extendability in a cone is controlled by extendability in the core. -/
lemma cone_ext_iff (hall : ∀ M ∈ indFam G S, M.card = alpha G S → w ∈ M)
    (hb1 : bD G S 1 = 0) (halpha : alpha G (S.erase w) + 1 = alpha G S) {A : Finset V} :
    Ext G S A ↔ Ext G (S.erase w) (A.erase w) := by
  constructor
  · rintro ⟨M, hM, hMc, hAM⟩
    obtain ⟨h1, h2⟩ := cone_max_erase hall halpha hM hMc
    exact ⟨M.erase w, h1, h2, Finset.erase_subset_erase w hAM⟩
  · rintro ⟨N, hN, hNc, hAN⟩
    obtain ⟨h1, h2⟩ := cone_max_insert hall hb1 halpha hN hNc
    refine ⟨insert w N, h1, h2, fun x hx => ?_⟩
    by_cases hxw : x = w
    · subst hxw; exact Finset.mem_insert_self _ _
    · exact Finset.mem_insert_of_mem (hAN (Finset.mem_erase.2 ⟨hxw, hx⟩))

/-- Counting extendable sets of a cone by whether they use the apex `w`. -/
lemma cone_ecnt (hall : ∀ M ∈ indFam G S, M.card = alpha G S → w ∈ M)
    (hb1 : bD G S 1 = 0) (halpha : alpha G (S.erase w) + 1 = alpha G S) (k : ℕ) :
    ecnt G S k = ecnt G (S.erase w) k
      + ((indFam G (S.erase w)).filter
          fun X => X.card + 1 = k ∧ Ext G (S.erase w) X).card := by
  classical
  have h1 : (((indFam G S).filter fun A => A.card = k ∧ Ext G S A).filter fun A => ¬ w ∈ A)
      = (indFam G (S.erase w)).filter fun A => A.card = k ∧ Ext G (S.erase w) A := by
    ext A
    simp only [Finset.mem_filter, mem_indFam, Finset.subset_erase]
    constructor
    · rintro ⟨⟨⟨hs, hind⟩, hcard, hext⟩, hwA⟩
      refine ⟨⟨⟨hs, hwA⟩, hind⟩, hcard, ?_⟩
      rw [cone_ext_iff hall hb1 halpha, Finset.erase_eq_of_notMem hwA] at hext
      exact hext
    · rintro ⟨⟨⟨hs, hwA⟩, hind⟩, hcard, hext⟩
      refine ⟨⟨⟨hs, hind⟩, hcard, ?_⟩, hwA⟩
      rw [cone_ext_iff hall hb1 halpha, Finset.erase_eq_of_notMem hwA]
      exact hext
  have h2 : ((((indFam G S).filter fun A => A.card = k ∧ Ext G S A).filter fun A => w ∈ A)).card
      = ((indFam G (S.erase w)).filter
          fun X => X.card + 1 = k ∧ Ext G (S.erase w) X).card := by
    refine Finset.card_bij' (fun A _ => A.erase w) (fun X _ => insert w X) ?_ ?_ ?_ ?_
    · intro A hA
      simp only [Finset.mem_filter] at hA ⊢
      obtain ⟨⟨hA1, hcard, hext⟩, hwA⟩ := hA
      refine ⟨?_, ?_, ?_⟩
      · rw [mem_indFam] at hA1 ⊢
        exact ⟨fun x hx => Finset.mem_erase.2
            ⟨(Finset.mem_erase.1 hx).1, hA1.1 (Finset.mem_erase.1 hx).2⟩,
          indep_of_subset hA1.2 (Finset.erase_subset _ _)⟩
      · rw [Finset.card_erase_of_mem hwA]
        have : 1 ≤ A.card := Finset.card_pos.2 ⟨w, hwA⟩
        omega
      · exact (cone_ext_iff hall hb1 halpha).1 hext
    · intro X hX
      simp only [Finset.mem_filter] at hX ⊢
      obtain ⟨hX1, hcard, hext⟩ := hX
      have hwX : w ∉ X := fun h => (Finset.mem_erase.1 ((mem_indFam.1 hX1).1 h)).1 rfl
      obtain ⟨N, hN, hNc, hXN⟩ := hext
      obtain ⟨hm1, hm2⟩ := cone_max_insert hall hb1 halpha hN hNc
      refine ⟨⟨indFam_subset hm1 (Finset.insert_subset_insert _ hXN), ?_, ?_⟩,
        Finset.mem_insert_self _ _⟩
      · rw [Finset.card_insert_of_notMem hwX]; exact hcard
      · exact ⟨insert w N, hm1, hm2, Finset.insert_subset_insert _ hXN⟩
    · intro A hA
      simp only [Finset.mem_filter] at hA
      exact Finset.insert_erase hA.2
    · intro X hX
      simp only [Finset.mem_filter] at hX
      exact Finset.erase_insert
        (fun h => (Finset.mem_erase.1 ((mem_indFam.1 hX.1).1 h)).1 rfl)
  rw [ecnt, ← Finset.card_filter_add_card_filter_not
    (s := (indFam G S).filter fun A => A.card = k ∧ Ext G S A) (p := fun A => w ∈ A),
    h2, h1, Nat.add_comm, ecnt]

/-- The cone identity at defect zero. -/
lemma cone_eD_zero
    (hall : ∀ M ∈ indFam G S, M.card = alpha G S → w ∈ M)
    (hb1 : bD G S 1 = 0) (halpha : alpha G (S.erase w) + 1 = alpha G S) :
    eD G S 0 = eD G (S.erase w) 0 := by
  have h := cone_ecnt hall hb1 halpha (alpha G S)
  rw [ecnt_eq_zero_of_gt (S := S.erase w) (by omega), zero_add] at h
  simp only [eD, Nat.sub_zero, if_pos (Nat.zero_le _)]
  rw [h, ecnt]
  exact congrArg Finset.card (Finset.filter_congr fun X _ => by
    constructor
    · rintro ⟨ha, hb⟩; exact ⟨by omega, hb⟩
    · rintro ⟨ha, hb⟩; exact ⟨by omega, hb⟩)

/-- The cone identity at positive defect. -/
lemma cone_eD_succ
    (hall : ∀ M ∈ indFam G S, M.card = alpha G S → w ∈ M)
    (hb1 : bD G S 1 = 0) (halpha : alpha G (S.erase w) + 1 = alpha G S) (d : ℕ) :
    eD G S (d + 1) = eD G (S.erase w) (d + 1) + eD G (S.erase w) d := by
  rcases le_or_gt (d + 1) (alpha G S) with hd | hd
  · have h := cone_ecnt hall hb1 halpha (alpha G S - (d + 1))
    rw [eD, if_pos hd, h]
    have e1 : eD G (S.erase w) d = ecnt G (S.erase w) (alpha G S - (d + 1)) := by
      rw [eD, if_pos (by omega)]
      congr 1
      omega
    have e2 : eD G (S.erase w) (d + 1)
        = ((indFam G (S.erase w)).filter
            fun X => X.card + 1 = alpha G S - (d + 1) ∧ Ext G (S.erase w) X).card := by
      rcases le_or_gt (d + 1) (alpha G (S.erase w)) with hd' | hd'
      · rw [eD, if_pos hd', ecnt]
        exact congrArg Finset.card (Finset.filter_congr fun X _ => by
          constructor
          · rintro ⟨ha, hb⟩; exact ⟨by omega, hb⟩
          · rintro ⟨ha, hb⟩; exact ⟨by omega, hb⟩)
      · rw [eD, if_neg (by omega)]
        refine (Finset.card_eq_zero.2 (Finset.eq_empty_of_forall_notMem fun X hX => ?_)).symm
        have := (Finset.mem_filter.1 hX).2.1
        omega
    rw [e1, e2, Nat.add_comm]
  · rw [eD_eq_zero_of_gt (by omega), eD_eq_zero_of_gt (by omega), eD_eq_zero_of_gt (by omega)]

/-- The core of a cone has no blocked sets at defect one. -/
lemma cone_bD_one
    (hall : ∀ M ∈ indFam G S, M.card = alpha G S → w ∈ M)
    (hb1 : bD G S 1 = 0) (hb2 : bD G S 2 = 0) (halpha : alpha G (S.erase w) + 1 = alpha G S) :
    bD G (S.erase w) 1 = 0 := by
  rw [bD]
  split
  · rename_i hd
    rw [bcnt]
    refine Finset.card_eq_zero.2 (Finset.eq_empty_of_forall_notMem fun X hX => ?_)
    obtain ⟨hX1, hXc, hXn⟩ := Finset.mem_filter.1 hX
    have hXS : X ∈ indFam G S := indFam_mono (Finset.erase_subset _ _) hX1
    have hwX : w ∉ X := fun h => (Finset.mem_erase.1 ((mem_indFam.1 hX1).1 h)).1 rfl
    have hnot : ¬ Ext G S X := by
      rw [cone_ext_iff hall hb1 halpha, Finset.erase_eq_of_notMem hwX]
      exact hXn
    rw [bD, if_pos (by omega), bcnt, Finset.card_eq_zero] at hb2
    have hmem : X ∈ ((indFam G S).filter fun A => A.card = alpha G S - 2 ∧ ¬ Ext G S A) :=
      Finset.mem_filter.2 ⟨hXS, by omega, hnot⟩
    rw [hb2] at hmem
    exact absurd hmem (Finset.notMem_empty _)
  · rfl

end FivefoldForest
