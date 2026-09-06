import RequestProject.B1Zero.FlexCorona

/-!
# The case of four forbidden vertices

In the `b₁ = 0` window the number `r` of forbidden vertices is at most four, and `r = 4`
forces `(|V|, α, |Forced|, #flexBags) = (33, 19, 9, 10)`.  This file proves the density bound
`3 (α-3) b₄ ≤ (α-7) e₄` in that case, i.e. `4 b₄ ≤ e₄`:

* `card_extSets_eq_convolution`: `e₄ = 126 p₀ + 84 p₁ + 36 p₂ + 9 p₃ + p₄`;
* `card_blkSets_le_union_bound`:
  `b₄ ≤ 73 p₀ + 24 p₁ + 3 p₂ + 15 q₀ + 6 q₁ + q₂`, by fibring the blocked sets over their
  forbidden part, using the expansion `|N_C(T)| ≥ 2|T|+1` and a union bound over the four
  forbidden vertices;
* `density_four`, obtained by combining these with `flex_certificate`.
-/

open Finset SimpleGraph

namespace MatchingBag

open B1Zero

attribute [local instance 1] Classical.propDecidable

variable {V : Type*} [Fintype V] [DecidableEq V]

namespace TreeMatching

variable {D : TreeMatching V}

/-- The flexible neighbours of a set of vertices. -/
noncomputable def flexNbrsSet (D : TreeMatching V) (T : Finset V) : Finset V :=
  D.FlexV.filter (fun z => ∃ u ∈ T, D.G.Adj u z)

lemma mem_flexNbrsSet {T : Finset V} {z : V} :
    z ∈ D.flexNbrsSet T ↔ z ∈ D.FlexV ∧ ∃ u ∈ T, D.G.Adj u z := by
  rw [flexNbrsSet, Finset.mem_filter]

lemma flexAttach_eq : D.flexAttach = D.flexNbrsSet D.ForbiddenV := rfl

/-! ### Elementary independence facts -/

lemma indepOn_ForcedV : IndepOn D.G D.ForcedV := by
  intro u hu v hv hadj
  exact Finset.disjoint_left.1 Allowed_disjoint_Forbidden (Forced_subset_Allowed hv)
    (forbidden_of_adj_forced hu hadj)

lemma no_cross_forced_flex : ∀ u ∈ D.ForcedV, ∀ w ∈ D.FlexV, ¬ D.G.Adj u w := by
  intro u hu w hw hadj
  exact Finset.disjoint_left.1 Allowed_disjoint_Forbidden (mem_Flex.1 hw).1
    (forbidden_of_adj_forced hu hadj)

lemma card_le_flexBags_of_subset {W : Finset V} (hW : W ⊆ D.FlexV) :
    ∀ S ⊆ W, IndepOn D.G S → S.card ≤ D.flexBags.card :=
  fun _ hS hind => card_le_card_flexBags (hS.trans hW) hind

/-- Splitting the independent sets of a forced part and a flexible part. -/
theorem cntRev_split {F W : Finset V} (hF : F ⊆ D.ForcedV) (hW : W ⊆ D.FlexV) (j : ℕ) :
    cntRev D.G (F ∪ W) (F.card + D.flexBags.card) j
      = ∑ ij ∈ Finset.antidiagonal j, (F.card).choose ij.1 * cntRev D.G W D.flexBags.card ij.2 := by
  have hdisj : Disjoint F W :=
    (Flex_disjoint_Forced.mono (hW) hF).symm
  have hcross : ∀ u ∈ F, ∀ w ∈ W, ¬ D.G.Adj u w :=
    fun u hu w hw => no_cross_forced_flex u (hF hu) w (hW hw)
  rw [cntRev_union hdisj hcross (fun S hS _ => Finset.card_le_card hS)
    (card_le_flexBags_of_subset hW) j]
  refine Finset.sum_congr rfl fun ij _ => ?_
  rw [cntRev_choose (indepOn_ForcedV.subset hF)]

/-! ### The extendable count -/

lemma card_extSets_eq_cntRev (hb1 : D.B1Zero) {k j n : ℕ} (h : k + j = n) :
    (D.extSets k).card = cntRev D.G (D.ForcedV ∪ D.FlexV) n j := by
  rw [extSets_eq_indSubsets hb1, cntRev]
  congr 1
  ext S
  rw [mem_indSubsets, mem_revSets]
  constructor
  · rintro ⟨hsub, hcard, hind⟩
    exact ⟨hsub, hind, by omega⟩
  · rintro ⟨hsub, hind, hcard⟩
    exact ⟨hsub, by omega, hind⟩

/-- **The extendable count in the four-forbidden case.** -/
theorem card_extSets_eq_convolution (hb1 : D.B1Zero) (hc : D.ForcedV.card = 9)
    (hm : D.flexBags.card = 10) :
    (D.extSets 15).card
      = 126 * cntRev D.G D.FlexV 10 0 + 84 * cntRev D.G D.FlexV 10 1
        + 36 * cntRev D.G D.FlexV 10 2 + 9 * cntRev D.G D.FlexV 10 3
        + cntRev D.G D.FlexV 10 4 := by
  have h := card_extSets_eq_cntRev (k := 15) (j := 4) (n := 19) hb1 rfl
  rw [h, show (19 : ℕ) = D.ForcedV.card + D.flexBags.card by omega,
    cntRev_split (Finset.Subset.refl _) (Finset.Subset.refl _) 4]
  simp only [Finset.Nat.sum_antidiagonal_eq_sum_range_succ_mk, Finset.sum_range_succ,
    Finset.sum_range_zero, hc, hm]
  norm_num [Nat.choose]
  ring

/-! ### The blocked count -/

/-- Every fibre of the blocked sets over their forbidden part injects into the independent
sets of the surviving allowed vertices. -/
theorem card_fibre_le (hb1 : D.B1Zero) {k n j : ℕ} {T W₁ W₂ : Finset V}
    (hW₁ : D.ForcedV \ D.forcedNbrsSet T ⊆ W₁)
    (hW₂ : D.FlexV \ D.flexNbrsSet T ⊆ W₂)
    (hnj : (k - T.card) + j = n) :
    ((D.blkSets k).filter (fun S => S ∩ D.ForbiddenV = T)).card
      ≤ cntRev D.G (W₁ ∪ W₂) n j := by
  rw [cntRev]
  refine Finset.card_le_card_of_injOn (fun S => S \ T) ?_ ?_
  · intro S hS
    rw [Finset.mem_coe, Finset.mem_filter, mem_blkSets_iff hb1] at hS
    obtain ⟨⟨hcard, hind, -⟩, hST⟩ := hS
    have hTS : T ⊆ S := by rw [← hST]; exact Finset.inter_subset_left
    rw [Finset.mem_coe, mem_revSets]
    refine ⟨?_, ?_, ?_⟩
    · intro x hx
      rw [Finset.mem_sdiff] at hx
      have hxA : x ∈ D.ForcedV ∪ D.FlexV := by
        rw [Forced_union_Flex]
        by_contra hxA
        have hxF : x ∈ D.ForbiddenV := by
          rw [Forbidden_eq_compl, Finset.mem_compl]; exact hxA
        exact hx.2 (hST ▸ Finset.mem_inter.2 ⟨hx.1, hxF⟩)
      rcases Finset.mem_union.1 hxA with hxF | hxL
      · refine Finset.mem_union_left _ (hW₁ (Finset.mem_sdiff.2 ⟨hxF, ?_⟩))
        intro hxN
        obtain ⟨-, u, huT, hadj⟩ := mem_forcedNbrsSet.1 hxN
        exact hind u (hTS huT) x hx.1 hadj
      · refine Finset.mem_union_right _ (hW₂ (Finset.mem_sdiff.2 ⟨hxL, ?_⟩))
        intro hxN
        obtain ⟨-, u, huT, hadj⟩ := mem_flexNbrsSet.1 hxN
        exact hind u (hTS huT) x hx.1 hadj
    · exact fun u hu v hv => hind u (Finset.mem_sdiff.1 hu).1 v (Finset.mem_sdiff.1 hv).1
    · rw [Finset.card_sdiff, Finset.inter_eq_left.2 hTS, hcard]
      exact hnj
  · intro S hS S' hS' heq
    rw [Finset.mem_coe, Finset.mem_filter] at hS hS'
    have hTS : T ⊆ S := by rw [← hS.2]; exact Finset.inter_subset_left
    have hTS' : T ⊆ S' := by rw [← hS'.2]; exact Finset.inter_subset_left
    rw [← Finset.union_sdiff_of_subset hTS, ← Finset.union_sdiff_of_subset hTS']
    simp only at heq
    rw [heq]

/-- **The union bound over the four forbidden vertices.** -/
theorem sum_flexNbrs_le (hr : D.ForbiddenV.card = 4) (j : ℕ) :
    ∑ u ∈ D.ForbiddenV, cntRev D.G (D.FlexV \ D.flexNbrsSet {u}) D.flexBags.card j
      ≤ 3 * cntRev D.G D.FlexV D.flexBags.card j
        + cntRev D.G (D.FlexV \ D.flexAttach) D.flexBags.card j := by
  classical
  simp only [cntRev]
  set m := D.flexBags.card with hm
  set 𝒮 := revSets D.G D.FlexV m j with h𝒮
  set 𝒜 := revSets D.G (D.FlexV \ D.flexAttach) m j with h𝒜
  have hsub : ∀ u : V, revSets D.G (D.FlexV \ D.flexNbrsSet {u}) m j ⊆ 𝒮 :=
    fun u => revSets_mono Finset.sdiff_subset m j
  have hsub' : 𝒜 ⊆ 𝒮 := revSets_mono Finset.sdiff_subset m j
  have hcover : 𝒮 \ 𝒜 ⊆ D.ForbiddenV.biUnion
      (fun u => 𝒮 \ revSets D.G (D.FlexV \ D.flexNbrsSet {u}) m j) := by
    intro S hS
    rw [Finset.mem_sdiff, mem_revSets] at hS
    obtain ⟨⟨hSsub, hind, hcard⟩, hSA⟩ := hS
    have : ¬ (S ⊆ D.FlexV \ D.flexAttach) := by
      intro h
      exact hSA (mem_revSets.2 ⟨h, hind, hcard⟩)
    obtain ⟨x, hxS, hx⟩ : ∃ x ∈ S, x ∉ D.FlexV \ D.flexAttach := by
      by_contra hcon
      push_neg at hcon
      exact this hcon
    have hxA : x ∈ D.flexAttach := by
      rw [Finset.mem_sdiff] at hx
      push_neg at hx
      exact hx (hSsub hxS)
    rw [flexAttach_eq, mem_flexNbrsSet] at hxA
    obtain ⟨hxF, u, hu, hadj⟩ := hxA
    refine Finset.mem_biUnion.2 ⟨u, hu, Finset.mem_sdiff.2 ⟨mem_revSets.2 ⟨hSsub, hind, hcard⟩, ?_⟩⟩
    intro hmem
    have := (mem_revSets.1 hmem).1 hxS
    rw [Finset.mem_sdiff] at this
    exact this.2 (mem_flexNbrsSet.2 ⟨hxF, u, Finset.mem_singleton_self u, hadj⟩)
  have hcard1 : (𝒮 \ 𝒜).card
      ≤ ∑ u ∈ D.ForbiddenV, (𝒮 \ revSets D.G (D.FlexV \ D.flexNbrsSet {u}) m j).card :=
    le_trans (Finset.card_le_card hcover) Finset.card_biUnion_le
  have hstep : ∑ u ∈ D.ForbiddenV, ((revSets D.G (D.FlexV \ D.flexNbrsSet {u}) m j).card
      + (𝒮 \ revSets D.G (D.FlexV \ D.flexNbrsSet {u}) m j).card) = 4 * 𝒮.card := by
    rw [Finset.sum_congr rfl (fun u _ => ?_), Finset.sum_const, hr, smul_eq_mul]
    rw [add_comm, Finset.card_sdiff_add_card_eq_card (hsub u)]
  rw [Finset.sum_add_distrib] at hstep
  have hle : 𝒜.card ≤ 𝒮.card := Finset.card_le_card hsub'
  have hsd : (𝒮 \ 𝒜).card = 𝒮.card - 𝒜.card := by
    rw [Finset.card_sdiff, Finset.inter_eq_left.2 hsub']
  omega

/-- **The blocked count in the four-forbidden case.** -/
theorem card_blkSets_le_union_bound (hb1 : D.B1Zero) (ha : 17 ≤ Fintype.card D.Bag)
    (hc : D.ForcedV.card = 9) (hr : D.ForbiddenV.card = 4) (hm : D.flexBags.card = 10) :
    (D.blkSets 15).card
      ≤ 73 * cntRev D.G D.FlexV 10 0 + 24 * cntRev D.G D.FlexV 10 1
        + 3 * cntRev D.G D.FlexV 10 2
        + 15 * cntRev D.G (D.FlexV \ D.flexAttach) 10 0
        + 6 * cntRev D.G (D.FlexV \ D.flexAttach) 10 1
        + cntRev D.G (D.FlexV \ D.flexAttach) 10 2 := by
  classical
  set f : Finset V → ℕ :=
    fun T => ((D.blkSets 15).filter (fun S => S ∩ D.ForbiddenV = T)).card with hf
  have hfib : (D.blkSets 15).card = ∑ T ∈ D.ForbiddenV.powerset, f T :=
    Finset.card_eq_sum_card_fiberwise
      (fun S _ => Finset.mem_powerset.2 Finset.inter_subset_right)
  -- the expansion bound `|N_C(T)| ≥ 2|T|+1`
  have hexp : ∀ T ⊆ D.ForbiddenV, T.Nonempty →
      (D.ForcedV \ D.forcedNbrsSet T).card ≤ 8 - 2 * T.card := by
    intro T hT hTne
    have h1 := card_forcedNbrsSet_ge hb1 (by omega) hT hTne
    have h2 : D.forcedNbrsSet T ⊆ D.ForcedV := Finset.filter_subset _ _
    have h3 : (D.ForcedV \ D.forcedNbrsSet T).card
        = D.ForcedV.card - (D.forcedNbrsSet T).card := by
      rw [Finset.card_sdiff, Finset.inter_eq_left.2 h2]
    omega
  -- a superset of the surviving forced vertices of the prescribed size
  have hbig : ∀ T ⊆ D.ForbiddenV, T.Nonempty → ∀ s, 8 - 2 * T.card ≤ s → s ≤ 9 →
      ∃ F, D.ForcedV \ D.forcedNbrsSet T ⊆ F ∧ F ⊆ D.ForcedV ∧ F.card = s := by
    intro T hT hTne s hs hs'
    refine Finset.exists_subsuperset_card_eq Finset.sdiff_subset ?_ (by omega)
    exact le_trans (hexp T hT hTne) hs
  -- the fibre bounds
  have hone : ∀ u ∈ D.ForbiddenV,
      f {u} ≤ 15 * cntRev D.G (D.FlexV \ D.flexNbrsSet {u}) 10 0
        + 6 * cntRev D.G (D.FlexV \ D.flexNbrsSet {u}) 10 1
        + cntRev D.G (D.FlexV \ D.flexNbrsSet {u}) 10 2 := by
    intro u hu
    simp only [hf]
    obtain ⟨F, hF1, hF2, hF3⟩ := hbig {u} (by simpa using hu) ⟨u, Finset.mem_singleton_self u⟩ 6
      (by simp) (by norm_num)
    have hle := card_fibre_le (k := 15) (n := 16) (j := 2) (T := {u}) (W₁ := F)
      (W₂ := D.FlexV \ D.flexNbrsSet {u}) hb1 hF1 (Finset.Subset.refl _) (by simp)
    rw [show (16 : ℕ) = F.card + D.flexBags.card by omega,
      cntRev_split hF2 Finset.sdiff_subset 2] at hle
    simp only [Finset.Nat.sum_antidiagonal_eq_sum_range_succ_mk, Finset.sum_range_succ,
      Finset.sum_range_zero, hF3, hm] at hle
    norm_num [Nat.choose] at hle
    omega
  have htwo : ∀ T ⊆ D.ForbiddenV, T.card = 2 →
      f T ≤ 4 * cntRev D.G D.FlexV 10 0 + cntRev D.G D.FlexV 10 1 := by
    intro T hT hT2
    simp only [hf]
    obtain ⟨F, hF1, hF2, hF3⟩ := hbig T hT (Finset.card_pos.1 (by omega)) 4
      (by omega) (by norm_num)
    have hle := card_fibre_le (k := 15) (n := 14) (j := 1) (T := T) (W₁ := F)
      (W₂ := D.FlexV) hb1 hF1 Finset.sdiff_subset (by omega)
    rw [show (14 : ℕ) = F.card + D.flexBags.card by omega,
      cntRev_split hF2 (Finset.Subset.refl _) 1] at hle
    simp only [Finset.Nat.sum_antidiagonal_eq_sum_range_succ_mk, Finset.sum_range_succ,
      Finset.sum_range_zero, hF3, hm] at hle
    norm_num [Nat.choose] at hle
    omega
  have hthree : ∀ T ⊆ D.ForbiddenV, T.card = 3 → f T ≤ cntRev D.G D.FlexV 10 0 := by
    intro T hT hT3
    simp only [hf]
    obtain ⟨F, hF1, hF2, hF3⟩ := hbig T hT (Finset.card_pos.1 (by omega)) 2
      (by omega) (by norm_num)
    have hle := card_fibre_le (k := 15) (n := 12) (j := 0) (T := T) (W₁ := F)
      (W₂ := D.FlexV) hb1 hF1 Finset.sdiff_subset (by omega)
    rw [show (12 : ℕ) = F.card + D.flexBags.card by omega,
      cntRev_split hF2 (Finset.Subset.refl _) 0] at hle
    simp only [Finset.Nat.sum_antidiagonal_eq_sum_range_succ_mk, Finset.sum_range_succ,
      Finset.sum_range_zero, hF3, hm] at hle
    norm_num [Nat.choose] at hle
    omega
  have hfour : ∀ T ⊆ D.ForbiddenV, T.card = 4 → f T = 0 := by
    intro T hT hT4
    simp only [hf]
    have hTne : T.Nonempty := Finset.card_pos.1 (by omega)
    have h1 := card_forcedNbrsSet_ge hb1 (by omega) hT hTne
    have h2 : D.forcedNbrsSet T ⊆ D.ForcedV := Finset.filter_subset _ _
    have hempty : D.ForcedV \ D.forcedNbrsSet T = ∅ := by
      rw [Finset.sdiff_eq_empty_iff_subset]
      exact Finset.eq_of_subset_of_card_le h2 (by omega) ▸ Finset.Subset.refl _
    have hle := card_fibre_le (k := 15) (n := 11) (j := 0) (T := T) (W₁ := ∅)
      (W₂ := D.FlexV) hb1 (by rw [hempty]) Finset.sdiff_subset (by omega)
    have hzero : cntRev D.G ((∅ : Finset V) ∪ D.FlexV) 11 0 = 0 := by
      rw [Finset.empty_union, cntRev, Finset.card_eq_zero, Finset.eq_empty_iff_forall_notMem]
      intro S hS
      rw [mem_revSets] at hS
      have := card_le_card_flexBags hS.1 hS.2.1
      omega
    omega
  have hzeroe : f ∅ = 0 := by
    simp only [hf]
    rw [Finset.card_eq_zero, Finset.eq_empty_iff_forall_notMem]
    intro S hS
    rw [Finset.mem_filter, mem_blkSets_iff hb1] at hS
    obtain ⟨⟨-, -, hne⟩, hemp⟩ := hS
    rw [hemp] at hne
    exact hne.ne_empty rfl
  -- assemble the sum over the powerset
  rw [hfib, Finset.sum_powerset, hr]
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, Nat.zero_add]
  have e0 : ∑ T ∈ Finset.powersetCard 0 D.ForbiddenV, f T = 0 := by
    rw [Finset.powersetCard_zero, Finset.sum_singleton, hzeroe]
  have e1 : ∑ T ∈ Finset.powersetCard 1 D.ForbiddenV, f T
      ≤ 15 * (3 * cntRev D.G D.FlexV 10 0 + cntRev D.G (D.FlexV \ D.flexAttach) 10 0)
        + 6 * (3 * cntRev D.G D.FlexV 10 1 + cntRev D.G (D.FlexV \ D.flexAttach) 10 1)
        + (3 * cntRev D.G D.FlexV 10 2 + cntRev D.G (D.FlexV \ D.flexAttach) 10 2) := by
    rw [Finset.powersetCard_one, Finset.sum_map]
    simp only [Function.Embedding.coeFn_mk]
    have hb := Finset.sum_le_sum hone
    rw [Finset.sum_add_distrib, Finset.sum_add_distrib, ← Finset.mul_sum,
      ← Finset.mul_sum] at hb
    have h0 := sum_flexNbrs_le hr 0
    have h1 := sum_flexNbrs_le hr 1
    have h2 := sum_flexNbrs_le hr 2
    rw [hm] at h0 h1 h2
    omega
  have e2 : ∑ T ∈ Finset.powersetCard 2 D.ForbiddenV, f T
      ≤ 6 * (4 * cntRev D.G D.FlexV 10 0 + cntRev D.G D.FlexV 10 1) := by
    have hb : ∑ T ∈ Finset.powersetCard 2 D.ForbiddenV, f T
        ≤ ∑ _T ∈ Finset.powersetCard 2 D.ForbiddenV,
            (4 * cntRev D.G D.FlexV 10 0 + cntRev D.G D.FlexV 10 1) := by
      refine Finset.sum_le_sum fun T hT => ?_
      rw [Finset.mem_powersetCard] at hT
      exact htwo T hT.1 hT.2
    rw [Finset.sum_const, Finset.card_powersetCard, hr, smul_eq_mul] at hb
    norm_num [Nat.choose] at hb
    omega
  have e3 : ∑ T ∈ Finset.powersetCard 3 D.ForbiddenV, f T ≤ 4 * cntRev D.G D.FlexV 10 0 := by
    have hb : ∑ T ∈ Finset.powersetCard 3 D.ForbiddenV, f T
        ≤ ∑ _T ∈ Finset.powersetCard 3 D.ForbiddenV, cntRev D.G D.FlexV 10 0 := by
      refine Finset.sum_le_sum fun T hT => ?_
      rw [Finset.mem_powersetCard] at hT
      exact hthree T hT.1 hT.2
    rw [Finset.sum_const, Finset.card_powersetCard, hr, smul_eq_mul] at hb
    norm_num [Nat.choose] at hb
    omega
  have e4 : ∑ T ∈ Finset.powersetCard 4 D.ForbiddenV, f T = 0 := by
    rw [Finset.sum_eq_zero]
    intro T hT
    rw [Finset.mem_powersetCard] at hT
    exact hfour T hT.1 hT.2
  omega

/-! ### The density bound -/

/-- **The density bound for exactly four forbidden vertices.** -/
theorem density_four_aux (hb1 : D.B1Zero) (hconn : D.G.Connected)
    (ha : Fintype.card D.Bag = 19) (hc : D.ForcedV.card = 9) (hr : D.ForbiddenV.card = 4)
    (hm : D.flexBags.card = 10) :
    4 * (D.blkSets 15).card ≤ (D.extSets 15).card := by
  have hU : D.ForbiddenV.Nonempty := Finset.card_pos.1 (by omega)
  have hb := card_blkSets_le_union_bound hb1 (by omega) hc hr hm
  have he := card_extSets_eq_convolution hb1 hc hm
  have hcert := flex_certificate hb1 hconn hU hm
  omega

end TreeMatching

end MatchingBag
