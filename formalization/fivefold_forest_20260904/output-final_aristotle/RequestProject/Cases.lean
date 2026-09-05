import RequestProject.PendantCount
import RequestProject.Cone

/-!
# The individual cases of the tree induction

Each lemma here takes the relevant structural configuration together with the
conclusion of the induction hypothesis for the strictly smaller pieces, and produces
the fivefold bounds for `S`.
-/

open Finset

namespace FivefoldForest

open scoped Classical

variable {V : Type*} {G : SimpleGraph V} {S : Finset V} {p w c : V} {L : Finset V}

/-! ### Stars -/

/-- `S` is a star with centre `c`: every other vertex is a leaf attached to `c`. -/
def IsStarAt (G : SimpleGraph V) (S : Finset V) (c : V) : Prop :=
  c ∈ S ∧ ∀ x ∈ S, x ≠ c → G.Adj c x ∧ ∀ y ∈ S, G.Adj x y → y = c

/-- The independent sets of a star: subsets of the leaf set, plus the centre alone. -/
lemma star_indFam (h : IsStarAt G S c) :
    indFam G S = insert {c} ((S.erase c).powerset) := by
  ext A
  simp only [mem_indFam, Finset.mem_insert, Finset.mem_powerset]
  constructor
  · rintro ⟨hsub, hind⟩
    by_cases hc : c ∈ A
    · left
      refine Finset.eq_singleton_iff_unique_mem.2 ⟨hc, fun x hx => ?_⟩
      by_contra hxc
      exact hind (by exact_mod_cast hc) (by exact_mod_cast hx) (Ne.symm hxc)
        (h.2 x (hsub hx) hxc).1
    · right
      intro x hx
      exact Finset.mem_erase.2 ⟨fun hxc => hc (hxc ▸ hx), hsub hx⟩
  · rintro (rfl | hsub)
    · refine ⟨by simpa using h.1, ?_⟩
      simp
    · refine ⟨hsub.trans (Finset.erase_subset _ _), ?_⟩
      intro x hx y hy _ hadj
      have hx' := Finset.mem_erase.1 (hsub (by exact_mod_cast hx))
      have hy' := Finset.mem_erase.1 (hsub (by exact_mod_cast hy))
      exact hy'.1 ((h.2 x hx'.2 hx'.1).2 y hy'.2 hadj)

lemma star_mem_indFam_of_subset (h : IsStarAt G S c) {A : Finset V} (hA : A ⊆ S.erase c) :
    A ∈ indFam G S := by
  rw [star_indFam h]; exact Finset.mem_insert_of_mem (Finset.mem_powerset.2 hA)

lemma star_erase_mem_indFam (h : IsStarAt G S c) : S.erase c ∈ indFam G S :=
  star_mem_indFam_of_subset h (Finset.Subset.refl _)

lemma star_alpha_small (h : IsStarAt G S c) (hm : (S.erase c).card ≤ 1) : alpha G S ≤ 1 := by
  refine Finset.sup_le fun A hA => ?_
  rw [star_indFam h, Finset.mem_insert, Finset.mem_powerset] at hA
  rcases hA with rfl | hA
  · simp
  · exact le_trans (Finset.card_le_card hA) hm

lemma star_alpha (h : IsStarAt G S c) (hm : 2 ≤ (S.erase c).card) :
    alpha G S = (S.erase c).card := by
  refine le_antisymm ?_ (card_le_alpha (star_erase_mem_indFam h))
  refine Finset.sup_le fun A hA => ?_
  rw [star_indFam h, Finset.mem_insert, Finset.mem_powerset] at hA
  rcases hA with rfl | hA
  · rw [Finset.card_singleton]; omega
  · exact Finset.card_le_card hA

/-- In a star with at least two leaves, the leaf set is the unique maximum independent set. -/
lemma star_ext_iff (h : IsStarAt G S c) (hm : 2 ≤ (S.erase c).card) {A : Finset V} :
    Ext G S A ↔ A ⊆ S.erase c := by
  constructor
  · rintro ⟨M, hM, hcard, hAM⟩
    rw [star_alpha h hm] at hcard
    rw [star_indFam h, Finset.mem_insert, Finset.mem_powerset] at hM
    rcases hM with rfl | hM
    · rw [Finset.card_singleton] at hcard; omega
    · exact hAM.trans hM
  · intro hA
    exact ⟨S.erase c, star_erase_mem_indFam h, (star_alpha h hm).symm, hA⟩

lemma star_ecnt (h : IsStarAt G S c) (hm : 2 ≤ (S.erase c).card) (k : ℕ) :
    ecnt G S k = (S.erase c).card.choose k := by
  have hset : ((indFam G S).filter fun A => A.card = k ∧ Ext G S A)
      = Finset.powersetCard k (S.erase c) := by
    ext A
    simp only [Finset.mem_filter, Finset.mem_powersetCard, star_ext_iff h hm]
    constructor
    · rintro ⟨-, hk, hs⟩; exact ⟨hs, hk⟩
    · rintro ⟨hs, hk⟩
      exact ⟨star_mem_indFam_of_subset h hs, hk, hs⟩
  rw [ecnt, hset, Finset.card_powersetCard]

lemma star_bcnt (h : IsStarAt G S c) (hm : 2 ≤ (S.erase c).card) (k : ℕ) :
    bcnt G S k = if k = 1 then 1 else 0 := by
  rcases eq_or_ne k 1 with rfl | hk
  · rw [if_pos rfl]
    have hset : ((indFam G S).filter fun A => A.card = 1 ∧ ¬ Ext G S A)
        = {({c} : Finset V)} := by
      ext A
      simp only [Finset.mem_filter, Finset.mem_singleton, star_ext_iff h hm]
      constructor
      · rintro ⟨hA, -, hnot⟩
        rw [star_indFam h, Finset.mem_insert, Finset.mem_powerset] at hA
        rcases hA with rfl | hA
        · rfl
        · exact absurd hA hnot
      · rintro rfl
        refine ⟨?_, Finset.card_singleton c, ?_⟩
        · rw [star_indFam h]; exact Finset.mem_insert_self _ _
        · intro hsub
          exact (Finset.mem_erase.1 (hsub (Finset.mem_singleton_self c))).1 rfl
    rw [bcnt, hset, Finset.card_singleton]
  · rw [if_neg hk, bcnt]
    refine Finset.card_eq_zero.2 (Finset.eq_empty_of_forall_notMem fun A hA => ?_)
    simp only [Finset.mem_filter, star_ext_iff h hm] at hA
    obtain ⟨hA1, hcard, hnot⟩ := hA
    rw [star_indFam h, Finset.mem_insert, Finset.mem_powerset] at hA1
    rcases hA1 with rfl | hA1
    · rw [Finset.card_singleton] at hcard; exact hk hcard.symm
    · exact hnot hA1

lemma star_eD (h : IsStarAt G S c) (hm : 2 ≤ (S.erase c).card) (d : ℕ) :
    eD G S d = (S.erase c).card.choose d := by
  rw [eD, star_alpha h hm]
  split
  · rename_i hd
    rw [star_ecnt h hm, Nat.choose_symm hd]
  · rename_i hd
    exact (Nat.choose_eq_zero_of_lt (by omega)).symm

lemma star_bD (h : IsStarAt G S c) (hm : 2 ≤ (S.erase c).card) (d : ℕ) :
    bD G S d = if d + 1 = (S.erase c).card then 1 else 0 := by
  rw [bD, star_alpha h hm]
  split
  · rename_i hd
    rw [star_bcnt h hm]
    by_cases h1 : d + 1 = (S.erase c).card
    · rw [if_pos (by omega), if_pos h1]
    · rw [if_neg (by omega), if_neg h1]
  · rename_i hd
    rw [if_neg (by omega)]

/-- Stars satisfy the fivefold bounds, except for `K₁,₄`, which has the exceptional profile. -/
lemma star_good_or_star (h : IsStarAt G S c) (helig : Eligible G S) :
    Good G S ∨ StarProf G S := by
  rcases le_or_gt (S.erase c).card 1 with hm | hm
  · have := star_alpha_small h hm
    exact Or.inl ⟨margin_nonneg_of_alpha_le (by omega), margin_nonneg_of_alpha_le (by omega)⟩
  · have hm2 : 2 ≤ (S.erase c).card := hm
    have hb := star_bD h hm2
    have he := star_eD h hm2
    have ha := star_alpha h hm2
    have hne2 : (S.erase c).card ≠ 2 := by
      intro hc
      have := helig.2.1
      rw [hb 1, if_pos (by omega)] at this
      exact absurd this one_ne_zero
    have hne3 : (S.erase c).card ≠ 3 := by
      intro hc
      have := helig.2.2
      rw [hb 2, if_pos (by omega)] at this
      exact absurd this one_ne_zero
    rcases eq_or_ne (S.erase c).card 4 with h4 | h4
    · right
      refine ⟨by rw [ha, h4], ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
      · rw [he 0, h4]; decide
      · rw [he 1, h4]; decide
      · rw [he 2, h4]; decide
      · rw [he 3, h4]; decide
      · rw [he 4, h4]; decide
      · rw [hb 3, if_pos (by omega)]
      · rw [hb 4, if_neg (by omega)]
    · left
      have hm5 : 5 ≤ (S.erase c).card := by omega
      have hb3 : bD G S 3 = 0 := by rw [hb 3, if_neg (by omega)]
      constructor
      · simp only [margin, hb3, Nat.cast_zero, mul_zero, sub_zero]
        exact Int.natCast_nonneg _
      · rcases eq_or_ne (S.erase c).card 5 with h5 | h5
        · have hb4 : bD G S 4 = 1 := by rw [hb 4, if_pos (by omega)]
          have he4 : eD G S 4 = 5 := by rw [he 4, h5]; decide
          simp only [margin, hb4, he4, Nat.cast_one, Nat.cast_ofNat]
          norm_num
        · have hb4 : bD G S 4 = 0 := by rw [hb 4, if_neg (by omega)]
          simp only [margin, hb4, Nat.cast_zero, mul_zero, sub_zero]
          exact Int.natCast_nonneg _

/-! ### The pendant `P₂` case (one leaf) -/

lemma pendant_one_good (hp : Pendant G S p w L) (hm : L.card = 1) (helig : Eligible G S)
    (hGs2 : 2 ≤ (pBase S p L).card)
    (IH : ∀ T : Finset V, T.card < S.card → Eligible G T → Good G T ∨ StarProf G T) :
    Good G S := by
  have hF : alpha G (pCore S p w L) = alpha G (pBase S p L) :=
    Pendant.alpha_core_one hp hm helig
  have hsplit := hp.split
  have hLne : L.Nonempty := hp.Lne
  have hGsne : (pBase S p L).Nonempty := ⟨w, hp.wmem_base⟩
  -- layer identities
  have hmar : ∀ d, margin G S d = margin G (pRest S p) d + margin G (pCore S p w L) d := by
    intro d
    simp only [margin, Pendant.eD_one hp hm hF d, Pendant.bD_one hp hm hF d]
    push_cast
    ring
  -- eligibility of the two pieces
  have hb1 : bD G (pRest S p) 1 = 0 ∧ bD G (pCore S p w L) 1 = 0 := by
    have := Pendant.bD_one hp hm hF 1
    have h0 := helig.2.1
    omega
  have hb2 : bD G (pRest S p) 2 = 0 ∧ bD G (pCore S p w L) 2 = 0 := by
    have := Pendant.bD_one hp hm hF 2
    have h0 := helig.2.2
    omega
  have hbase1 : bD G (pBase S p L) 1 = 0 := by
    have := bD_le_of_split hsplit 1
    omega
  have hbase2 : bD G (pBase S p L) 2 = 0 := by
    have := bD_le_of_split hsplit 2
    omega
  have heligGs : Eligible G (pBase S p L) := ⟨hGsne, hbase1, hbase2⟩
  have hFne : (pCore S p w L).Nonempty := by
    rw [pCore_eq, ← Finset.card_pos]
    have := Finset.pred_card_le_card_erase (s := pBase S p L) (a := w)
    omega
  have hbase := IH _ hp.base_card_lt heligGs
  have hcore := IH _ hp.core_card_lt ⟨hFne, hb1.2, hb2.2⟩
  -- the leaf is a single vertex
  have hL0 : eD G L 0 = 1 := by rw [hp.leaves_eD 0, hm]; decide
  have hL1 : eD G L 1 = 1 := by rw [hp.leaves_eD 1, hm]; decide
  have hL2 : eD G L 2 = 0 := by rw [hp.leaves_eD 2, hm]; decide
  have hL3 : eD G L 3 = 0 := by rw [hp.leaves_eD 3, hm]; decide
  have hL4 : eD G L 4 = 0 := by rw [hp.leaves_eD 4, hm]; decide
  have hmL3 : margin G L 3 = 0 := by simp [margin, hL3, hp.leaves_bD 3]
  have hmL4 : margin G L 4 = 0 := by simp [margin, hL4, hp.leaves_bD 4]
  have m3' := margin3_split hsplit heligGs hp.leaves_eligible
  have m4' := margin4_split hsplit heligGs hp.leaves_eligible
  rw [hmL3, hL0, hL2, hL1] at m3'
  rw [hmL3, hmL4, hL0, hL1, hL2] at m4'
  push_cast at m3' m4'
  -- assemble
  have e3 := hmar 3
  have e4 := hmar 4
  rw [m3'] at e3
  rw [m4'] at e4
  have hn2 : (0 : ℤ) ≤ (eD G (pBase S p L) 2 : ℤ) := Int.natCast_nonneg _
  have hcore2 : StarProf G (pCore S p w L) → (1 : ℤ) ≤ (eD G (pBase S p L) 2 : ℤ) := by
    intro s2
    have h4 : alpha G (pBase S p L) = 4 := hF ▸ s2.1
    have := eD_pos (G := G) (S := pBase S p L) (d := 2) (by rw [h4]; omega)
    exact_mod_cast this
  constructor
  · rcases hbase with ⟨g3, _⟩ | s1 <;> rcases hcore with ⟨f3, _⟩ | s2
    · linarith
    · rw [s2.margin3] at e3; have := hcore2 s2; linarith
    · rw [s1.margin3, s1.2.2.2.1] at e3; push_cast at e3; linarith
    · rw [s1.margin3, s1.2.2.2.1, s2.margin3] at e3; push_cast at e3; linarith
  · rcases hbase with ⟨g3, g4⟩ | s1 <;> rcases hcore with ⟨f3, f4⟩ | s2
    · linarith
    · rw [s2.margin4] at e4; linarith
    · rw [s1.margin3, s1.margin4] at e4; linarith
    · rw [s1.margin3, s1.margin4, s2.margin4] at e4; linarith

/-! ### Two or more leaves -/

lemma le_choose_two {m : ℕ} (h : 3 ≤ m) : m ≤ m.choose 2 := by
  rw [Nat.choose_two_right, Nat.le_div_iff_mul_le (by norm_num)]
  have h1 : 2 ≤ m - 1 := by omega
  exact Nat.mul_le_mul_left _ h1

/-- The margin of `S` at defect three, in terms of the base. -/
lemma pendant_many_margin3 (hp : Pendant G S p w L) (hm : 2 ≤ L.card)
    (heligGs : Eligible G (pBase S p L)) (hd : 4 ≤ alpha G S) :
    margin G S 3 = (eD G (pBase S p L) 0 : ℤ) * (L.card.choose 3 : ℤ)
      + margin G (pBase S p L) 3
      + (eD G (pBase S p L) 1 : ℤ) * (L.card.choose 2 : ℤ)
      + (eD G (pBase S p L) 2 : ℤ) * (L.card : ℤ)
      - 5 * (icnt G (pCore S p w L) (alpha G S - 3 - 1) : ℤ) := by
  have hsplit := hp.split
  have he := Pendant.eD_two_le hp hm 3
  have hb := Pendant.bD_two_le hp hm 3 (by omega)
  have m3' := margin3_split hsplit heligGs hp.leaves_eligible
  rw [hp.leaves_eD 0, hp.leaves_eD 1, hp.leaves_eD 2] at m3'
  have hmL3 : margin G L 3 = (L.card.choose 3 : ℤ) := by
    simp [margin, hp.leaves_eD 3, hp.leaves_bD 3]
  rw [hmL3] at m3'
  simp only [Nat.choose_zero_right, Nat.choose_one_right, Nat.cast_one, one_mul] at m3'
  simp only [margin, he, hb]
  push_cast
  rw [show margin G (pRest S p) 3 = (eD G (pRest S p) 3 : ℤ) - 5 * (bD G (pRest S p) 3 : ℤ) from rfl]
    at m3'
  rw [show margin G (pBase S p L) 3
      = (eD G (pBase S p L) 3 : ℤ) - 5 * (bD G (pBase S p L) 3 : ℤ) from rfl] at *
  linarith [m3']

/-- The margin of `S` at defect four, in terms of the base. -/
lemma pendant_many_margin4 (hp : Pendant G S p w L) (hm : 2 ≤ L.card)
    (heligGs : Eligible G (pBase S p L)) (hd : 5 ≤ alpha G S) :
    margin G S 4 = (eD G (pBase S p L) 0 : ℤ) * (L.card.choose 4 : ℤ)
      + margin G (pBase S p L) 4
      + (eD G (pBase S p L) 1 : ℤ) * (L.card.choose 3 : ℤ)
      + (L.card : ℤ) * margin G (pBase S p L) 3
      + (eD G (pBase S p L) 2 : ℤ) * (L.card.choose 2 : ℤ)
      - 5 * (icnt G (pCore S p w L) (alpha G S - 4 - 1) : ℤ) := by
  have hsplit := hp.split
  have he := Pendant.eD_two_le hp hm 4
  have hb := Pendant.bD_two_le hp hm 4 (by omega)
  have m4' := margin4_split hsplit heligGs hp.leaves_eligible
  rw [hp.leaves_eD 0, hp.leaves_eD 1, hp.leaves_eD 2] at m4'
  have hmL3 : margin G L 3 = (L.card.choose 3 : ℤ) := by
    simp [margin, hp.leaves_eD 3, hp.leaves_bD 3]
  have hmL4 : margin G L 4 = (L.card.choose 4 : ℤ) := by
    simp [margin, hp.leaves_eD 4, hp.leaves_bD 4]
  rw [hmL3, hmL4] at m4'
  simp only [Nat.choose_zero_right, Nat.choose_one_right, Nat.cast_one, one_mul] at m4'
  simp only [margin, he, hb]
  push_cast
  rw [show margin G (pRest S p) 4 = (eD G (pRest S p) 4 : ℤ) - 5 * (bD G (pRest S p) 4 : ℤ) from rfl]
    at m4'
  rw [show margin G (pBase S p L) 4
      = (eD G (pBase S p L) 4 : ℤ) - 5 * (bD G (pBase S p L) 4 : ℤ) from rfl] at *
  rw [show margin G (pBase S p L) 3
      = (eD G (pBase S p L) 3 : ℤ) - 5 * (bD G (pBase S p L) 3 : ℤ) from rfl] at *
  linarith [m4']

/-- Shared setup for a support vertex with at least two leaves. -/
lemma pendant_many_setup (hp : Pendant G S p w L) (hm : 2 ≤ L.card) (helig : Eligible G S) :
    Eligible G (pBase S p L) ∧ alpha G S = alpha G (pBase S p L) + L.card ∧ 3 ≤ L.card := by
  have hsplit := hp.split
  have hGsne : (pBase S p L).Nonempty := ⟨w, hp.wmem_base⟩
  have ha1 : 1 ≤ alpha G (pBase S p L) := one_le_alpha hGsne
  have halpha : alpha G S = alpha G (pBase S p L) + L.card := by
    rw [← hp.alpha_rest, alpha_split hsplit, hp.leaves_alpha]
  have hbS' : ∀ d, d + 1 ≤ alpha G S →
      bD G S d = bD G (pRest S p) d + icnt G (pCore S p w L) (alpha G S - d - 1) := fun d hd =>
    Pendant.bD_two_le hp hm d hd
  have hb1S' : bD G (pRest S p) 1 = 0 := by
    have h := hbS' 1 (by omega)
    have := helig.2.1
    omega
  have hb2S' : bD G (pRest S p) 2 = 0 := by
    have h := hbS' 2 (by omega)
    have := helig.2.2
    omega
  refine ⟨⟨hGsne, Nat.le_zero.mp (hb1S' ▸ bD_le_of_split hsplit 1),
    Nat.le_zero.mp (hb2S' ▸ bD_le_of_split hsplit 2)⟩, halpha, ?_⟩
  by_contra hcon
  have hm2 : L.card = 2 := by omega
  have h2 := hbS' 2 (by omega)
  rw [helig.2.2, hb2S'] at h2
  have hle : alpha G (pBase S p L) ≤ alpha G (pCore S p w L) + 1 := by
    rw [pCore]; exact alpha_erase_succ_ge G _ w
  have hpos : 1 ≤ icnt G (pCore S p w L) (alpha G S - 2 - 1) := icnt_pos (by omega)
  omega


/-- Four or more leaves. -/
lemma pendant_four_le_good (hp : Pendant G S p w L) (hm : 4 ≤ L.card) (helig : Eligible G S)
    (hforest : ∀ T : Finset V, T.card ≤ 2 * alpha G T)
    (IH : ∀ T : Finset V, T.card < S.card → Eligible G T → Good G T ∨ StarProf G T) :
    Good G S := by
  obtain ⟨heligGs, halpha, -⟩ := pendant_many_setup hp (by omega) helig
  have hbase := IH _ hp.base_card_lt heligGs
  have hGsne := heligGs.1
  have ha1 : 1 ≤ alpha G (pBase S p L) := one_le_alpha hGsne
  have hFsub : pCore S p w L ⊆ pBase S p L := Pendant.core_subset_base S p w L
  have hFa : alpha G (pCore S p w L) ≤ alpha G (pBase S p L) := alpha_mono hFsub
  -- values of the base profile
  have hn0 : (1 : ℤ) ≤ (eD G (pBase S p L) 0 : ℤ) := by
    exact_mod_cast eD_pos (G := G) (S := pBase S p L) (Nat.zero_le _)
  have hn1 : (0 : ℤ) ≤ (eD G (pBase S p L) 1 : ℤ) := Int.natCast_nonneg _
  have hn2 : (0 : ℤ) ≤ (eD G (pBase S p L) 2 : ℤ) := Int.natCast_nonneg _
  have hq0 : (eD G (pBase S p L) 0 : ℤ) ≤ 2 * (eD G (pBase S p L) 1 : ℤ) := by
    exact_mod_cast e0_le_two_e1 (G := G) (S := pBase S p L) hGsne (hforest _)
  have hq1 : (eD G (pBase S p L) 1 : ℤ) ≤ (eD G (pBase S p L) 0 : ℤ)
      + 6 * (eD G (pBase S p L) 2 : ℤ) := by
    exact_mod_cast e1_le_e0_add_six_e2 (G := G) (S := pBase S p L) hGsne (hforest _)
  have hmar3 := pendant_many_margin3 hp (by omega) heligGs (by omega)
  have hmar4 := pendant_many_margin4 hp (by omega) heligGs (by omega)
  -- the two correction terms
  have hNle : ∀ k, icnt G (pCore S p w L) k ≤ icnt G (pBase S p L) k := fun k =>
    icnt_mono hFsub k
  have hiD0 : icnt G (pBase S p L) (alpha G (pBase S p L)) = eD G (pBase S p L) 0 := by
    have := icnt_eq_iD (G := G) (S := pBase S p L) (d := 0) (Nat.zero_le _)
    rw [Nat.sub_zero] at this
    rw [this, iD_zero]
  have hiD1 : icnt G (pBase S p L) (alpha G (pBase S p L) - 1) = eD G (pBase S p L) 1 := by
    have := icnt_eq_iD (G := G) (S := pBase S p L) (d := 1) (by omega)
    rw [this, iD_one, heligGs.2.1, Nat.add_zero]
  rcases eq_or_lt_of_le hm with hm4 | hm5
  · -- exactly four leaves
    have hN3 : (icnt G (pCore S p w L) (alpha G S - 3 - 1) : ℤ) ≤ (eD G (pBase S p L) 0 : ℤ) := by
      have h : alpha G S - 3 - 1 = alpha G (pBase S p L) := by omega
      rw [h]
      exact_mod_cast (hNle _).trans (le_of_eq hiD0)
    have hN4 : (icnt G (pCore S p w L) (alpha G S - 4 - 1) : ℤ) ≤ (eD G (pBase S p L) 1 : ℤ) := by
      have h : alpha G S - 4 - 1 = alpha G (pBase S p L) - 1 := by omega
      rw [h]
      exact_mod_cast (hNle _).trans (le_of_eq hiD1)
    have hc1 : (L.card : ℤ) = 4 := by rw [← hm4]; norm_num
    have hc2 : (L.card.choose 2 : ℤ) = 6 := by rw [← hm4]; decide
    have hc3 : (L.card.choose 3 : ℤ) = 4 := by rw [← hm4]; decide
    have hc4 : (L.card.choose 4 : ℤ) = 1 := by rw [← hm4]; decide
    rw [hc1, hc2, hc3] at hmar3
    rw [hc1, hc2, hc3, hc4] at hmar4
    rcases hbase with ⟨g3, g4⟩ | s1
    · exact ⟨by rw [hmar3]; linarith, by rw [hmar4]; linarith⟩
    · simp only [s1.margin3] at hmar3 hmar4
      simp only [s1.margin4] at hmar4
      simp only [s1.2.1, s1.2.2.1, s1.2.2.2.1] at hmar3 hmar4 hN3 hN4
      push_cast at hmar3 hmar4 hN3 hN4
      exact ⟨by rw [hmar3]; linarith, by rw [hmar4]; linarith⟩
  · have hm5 : 5 ≤ L.card := hm5
    have hN3 : icnt G (pCore S p w L) (alpha G S - 3 - 1) = 0 := by
      refine icnt_eq_zero_of_gt ?_
      omega
    rcases eq_or_lt_of_le hm5 with hm5' | hm6
    · -- exactly five leaves
      have hN4 : (icnt G (pCore S p w L) (alpha G S - 4 - 1) : ℤ) ≤ (eD G (pBase S p L) 0 : ℤ) := by
        have h : alpha G S - 4 - 1 = alpha G (pBase S p L) := by omega
        rw [h]
        exact_mod_cast (hNle _).trans (le_of_eq hiD0)
      have hc1 : (L.card : ℤ) = 5 := by rw [← hm5']; norm_num
      have hc2 : (L.card.choose 2 : ℤ) = 10 := by rw [← hm5']; decide
      have hc3 : (L.card.choose 3 : ℤ) = 10 := by rw [← hm5']; decide
      have hc4 : (L.card.choose 4 : ℤ) = 5 := by rw [← hm5']; decide
      rw [hc1, hc2, hc3, hN3] at hmar3
      rw [hc1, hc2, hc3, hc4] at hmar4
      norm_num at hmar3
      rcases hbase with ⟨g3, g4⟩ | s1
      · exact ⟨by rw [hmar3]; linarith, by rw [hmar4]; linarith⟩
      · simp only [s1.margin3] at hmar3 hmar4
        simp only [s1.margin4] at hmar4
        simp only [s1.2.1, s1.2.2.1, s1.2.2.2.1] at hmar3 hmar4 hN4
        push_cast at hmar3 hmar4 hN4
        exact ⟨by rw [hmar3]; linarith, by rw [hmar4]; linarith⟩
    · -- six or more leaves
      have hm6 : 6 ≤ L.card := hm6
      have hN4 : icnt G (pCore S p w L) (alpha G S - 4 - 1) = 0 := by
        refine icnt_eq_zero_of_gt ?_
        omega
      rw [hN3] at hmar3
      rw [hN4] at hmar4
      norm_num at hmar3 hmar4
      have hc2 : (L.card : ℤ) ≤ (L.card.choose 2 : ℤ) := by
        exact_mod_cast le_choose_two (by omega)
      have hc3 : (0 : ℤ) ≤ (L.card.choose 3 : ℤ) := Int.natCast_nonneg _
      have hc4 : (0 : ℤ) ≤ (L.card.choose 4 : ℤ) := Int.natCast_nonneg _
      have hcm : (6 : ℤ) ≤ (L.card : ℤ) := by exact_mod_cast hm6
      rcases hbase with ⟨g3, g4⟩ | s1
      · have t1 : (0:ℤ) ≤ (eD G (pBase S p L) 0 : ℤ) * (L.card.choose 3 : ℤ) :=
          mul_nonneg (by linarith) hc3
        have t2 : (0:ℤ) ≤ (eD G (pBase S p L) 1 : ℤ) * (L.card.choose 2 : ℤ) :=
          mul_nonneg hn1 (by linarith)
        have t3 : (0:ℤ) ≤ (eD G (pBase S p L) 2 : ℤ) * (L.card : ℤ) :=
          mul_nonneg hn2 (by linarith)
        have t4 : (0:ℤ) ≤ (eD G (pBase S p L) 0 : ℤ) * (L.card.choose 4 : ℤ) :=
          mul_nonneg (by linarith) hc4
        have t5 : (0:ℤ) ≤ (eD G (pBase S p L) 1 : ℤ) * (L.card.choose 3 : ℤ) :=
          mul_nonneg hn1 hc3
        have t6 : (0:ℤ) ≤ (L.card : ℤ) * margin G (pBase S p L) 3 :=
          mul_nonneg (by linarith) g3
        have t7 : (0:ℤ) ≤ (eD G (pBase S p L) 2 : ℤ) * (L.card.choose 2 : ℤ) :=
          mul_nonneg hn2 (by linarith)
        exact ⟨by rw [hmar3]; linarith, by rw [hmar4]; linarith⟩
      · simp only [s1.margin3] at hmar3 hmar4
        simp only [s1.margin4] at hmar4
        simp only [s1.2.1, s1.2.2.1, s1.2.2.2.1] at hmar3 hmar4
        push_cast at hmar3 hmar4
        exact ⟨by rw [hmar3]; linarith, by rw [hmar4]; linarith⟩


/-- Exactly three leaves. -/
lemma pendant_three_good (hp : Pendant G S p w L) (hm : L.card = 3) (helig : Eligible G S)
    (hGs2 : 2 ≤ (pBase S p L).card)
    (hforest : ∀ T : Finset V, T.card ≤ 2 * alpha G T)
    (IH : ∀ T : Finset V, T.card < S.card → Eligible G T → Good G T ∨ StarProf G T) :
    Good G S := by
  obtain ⟨heligGs, halpha, -⟩ := pendant_many_setup hp (by omega) helig
  have hbase := IH _ hp.base_card_lt heligGs
  have hGsne := heligGs.1
  have ha1 : 1 ≤ alpha G (pBase S p L) := one_le_alpha hGsne
  have hFne : (pCore S p w L).Nonempty := by
    rw [pCore_eq, ← Finset.card_pos]
    have := Finset.pred_card_le_card_erase (s := pBase S p L) (a := w)
    omega
  have hle : alpha G (pBase S p L) ≤ alpha G (pCore S p w L) + 1 := by
    rw [pCore_eq]; exact alpha_erase_succ_ge G _ w
  have hFmono : alpha G (pCore S p w L) ≤ alpha G (pBase S p L) := by
    rw [pCore_eq]; exact alpha_erase_le G _ w
  -- the deleted vertex is essential
  have hzero : icnt G (pCore S p w L) (alpha G (pBase S p L)) = 0 := by
    have h2 := Pendant.bD_two_le hp (by omega) 2 (by omega)
    have hb2S' : bD G (pRest S p) 2 = 0 :=
      Nat.le_zero.mp (helig.2.2 ▸ (by omega : bD G (pRest S p) 2 ≤ bD G S 2))
    rw [helig.2.2] at h2
    have hidx : alpha G S - 2 - 1 = alpha G (pBase S p L) := by omega
    rw [hidx] at h2
    omega
  have haF : alpha G (pCore S p w L) + 1 = alpha G (pBase S p L) := by
    rcases le_or_gt (alpha G (pBase S p L)) (alpha G (pCore S p w L)) with h | h
    · exact absurd (icnt_pos (G := G) (S := pCore S p w L) h) (by omega)
    · omega
  -- every maximum independent set of the base contains `w`
  have hall : ∀ M ∈ indFam G (pBase S p L), M.card = alpha G (pBase S p L) → w ∈ M := by
    intro M hM hcard
    by_contra hw
    have hM' : M ∈ indFam G (pCore S p w L) := by
      rw [pCore_eq, mem_indFam] at *
      exact ⟨Finset.subset_erase.2 ⟨hM.1, hw⟩, hM.2⟩
    have := card_le_alpha hM'
    omega
  -- cone identities
  have hc0 := cone_eD_zero (G := G) hall heligGs.2.1 (by rw [← pCore_eq] at *; omega)
  have hc1 := cone_eD_succ (G := G) hall heligGs.2.1 (by rw [← pCore_eq] at *; omega) 0
  have hc2 := cone_eD_succ (G := G) hall heligGs.2.1 (by rw [← pCore_eq] at *; omega) 1
  have hcb := cone_bD_one (G := G) hall heligGs.2.1 heligGs.2.2
    (by rw [← pCore_eq] at *; omega)
  rw [← pCore_eq] at hc0 hc1 hc2 hcb
  norm_num at hc1 hc2
  -- the two correction terms
  have hN3 : icnt G (pCore S p w L) (alpha G S - 3 - 1) = eD G (pCore S p w L) 0 := by
    have hidx : alpha G S - 3 - 1 = alpha G (pCore S p w L) - 0 := by omega
    rw [hidx, icnt_eq_iD (Nat.zero_le _), iD_zero]
  -- nonnegativity facts about the core profile
  have hf0 : (1 : ℤ) ≤ (eD G (pCore S p w L) 0 : ℤ) := by
    exact_mod_cast eD_pos (G := G) (S := pCore S p w L) (Nat.zero_le _)
  have hf1 : (0 : ℤ) ≤ (eD G (pCore S p w L) 1 : ℤ) := Int.natCast_nonneg _
  have hf2 : (0 : ℤ) ≤ (eD G (pCore S p w L) 2 : ℤ) := Int.natCast_nonneg _
  have hq0 : (eD G (pCore S p w L) 0 : ℤ) ≤ 2 * (eD G (pCore S p w L) 1 : ℤ) := by
    exact_mod_cast e0_le_two_e1 (G := G) (S := pCore S p w L) hFne (hforest _)
  have hkey : (0 : ℤ) ≤ (eD G (pCore S p w L) 0 : ℤ) - (eD G (pCore S p w L) 1 : ℤ)
      + 3 * (eD G (pCore S p w L) 2 : ℤ) := by
    rcases le_or_gt 3 (alpha G (pCore S p w L)) with h3 | h3
    · have := e1_le_three_e2 (G := G) (S := pCore S p w L) h3 (hforest _)
      have : (eD G (pCore S p w L) 1 : ℤ) ≤ 3 * (eD G (pCore S p w L) 2 : ℤ) := by
        exact_mod_cast this
      linarith
    · have h1 : 1 ≤ alpha G (pCore S p w L) := one_le_alpha hFne
      interval_cases h : alpha G (pCore S p w L)
      · have e1 : eD G (pCore S p w L) 1 = 1 := by
          have := eD_alpha_eq_one G (pCore S p w L); rwa [h] at this
        have e2 : eD G (pCore S p w L) 2 = 0 := eD_eq_zero_of_gt (by omega)
        rw [e1, e2]; push_cast; linarith
      · have e2 : eD G (pCore S p w L) 2 = 1 := by
          have := eD_alpha_eq_one G (pCore S p w L); rwa [h] at this
        have e1 : eD G (pCore S p w L) 1 ≤ 4 * eD G (pCore S p w L) 2 :=
          e1_le_four_e2 (G := G) (S := pCore S p w L) (by omega) (hforest _)
        rw [e2] at e1 ⊢
        have : (eD G (pCore S p w L) 1 : ℤ) ≤ 4 := by exact_mod_cast e1
        push_cast; linarith
  -- margin at defect three
  have hmar3 := pendant_many_margin3 hp (by omega) heligGs (by omega)
  have hcc1 : (L.card.choose 2 : ℤ) = 3 := by rw [hm]; decide
  have hcc2 : (L.card.choose 3 : ℤ) = 1 := by rw [hm]; decide
  have hcc0 : (L.card : ℤ) = 3 := by rw [hm]; norm_num
  rw [hcc0, hcc1, hcc2, hN3, hc0, hc1, hc2] at hmar3
  push_cast at hmar3
  constructor
  · rcases hbase with ⟨g3, g4⟩ | s1
    · rw [hmar3]; linarith
    · rw [s1.margin3] at hmar3
      have q0 := s1.2.1
      have q1 := s1.2.2.1
      have q2 := s1.2.2.2.1
      rw [hc0] at q0
      rw [hc1] at q1
      rw [hc2] at q2
      have p0 : (eD G (pCore S p w L) 0 : ℤ) = 1 := by exact_mod_cast q0
      have p1 : (eD G (pCore S p w L) 1 : ℤ) + (eD G (pCore S p w L) 0 : ℤ) = 4 := by
        exact_mod_cast q1
      have p2 : (eD G (pCore S p w L) 2 : ℤ) + (eD G (pCore S p w L) 1 : ℤ) = 6 := by
        exact_mod_cast q2
      rw [hmar3]; linarith
  · rcases le_or_gt (alpha G S) 4 with hsmall | hbig
    · exact margin_nonneg_of_alpha_le hsmall
    · have hmar4 := pendant_many_margin4 hp (by omega) heligGs (by omega)
      have hN4 : icnt G (pCore S p w L) (alpha G S - 4 - 1) = eD G (pCore S p w L) 1 := by
        have hidx : alpha G S - 4 - 1 = alpha G (pCore S p w L) - 1 := by omega
        rw [hidx, icnt_eq_iD (by omega), iD_one, hcb, Nat.add_zero]
      have hcc3 : (L.card.choose 4 : ℤ) = 0 := by rw [hm]; decide
      rw [hcc0, hcc1, hcc2, hcc3, hN4, hc0, hc1, hc2] at hmar4
      push_cast at hmar4
      rcases hbase with ⟨g3, g4⟩ | s1
      · rw [hmar4]; linarith
      · rw [s1.margin3, s1.margin4] at hmar4
        have q0 := s1.2.1
        have q1 := s1.2.2.1
        have q2 := s1.2.2.2.1
        rw [hc0] at q0
        rw [hc1] at q1
        rw [hc2] at q2
        have p0 : (eD G (pCore S p w L) 0 : ℤ) = 1 := by exact_mod_cast q0
        have p1 : (eD G (pCore S p w L) 1 : ℤ) + (eD G (pCore S p w L) 0 : ℤ) = 4 := by
          exact_mod_cast q1
        have p2 : (eD G (pCore S p w L) 2 : ℤ) + (eD G (pCore S p w L) 1 : ℤ) = 6 := by
          exact_mod_cast q2
        rw [hmar4]; linarith

end FivefoldForest
