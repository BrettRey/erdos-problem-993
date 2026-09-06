import RequestProject.BlockedShadow
import RequestProject.DepthThreeAlgebra

/-!
# The extendable incidence bound, positivity, and the low-density depth-three theorem

* `MatchingBag.TreeMatching.extendable_incidence`: `4 e₄ ≤ (α - 3) e₃`, by double counting
  extendable one-vertex extensions.
* `MatchingBag.TreeMatching.extSets_nonempty`: `e_d > 0` for `d ≤ α`.
* `MatchingBag.TreeMatching.depthThreeStrict_of_lowDensity`: the numerical assembly, using
  the graph Pascal reserve `erasure_depth_three_reserve`, the blocked shadow bound, the
  extendable incidence bound, and the scalar implication `DepthThree.low_density_algebra`.
-/

open Finset SimpleGraph

namespace MatchingBag

attribute [local instance 1] Classical.propDecidable

variable {V : Type*} [Fintype V] [DecidableEq V]

namespace TreeMatching

variable {D : TreeMatching V}

lemma card_of_maxIndepSets {I : Finset V} (hI : I ∈ D.maxIndepSets) :
    I.card = Fintype.card D.Bag := by
  rw [maxIndepSets_card hI, card_Bag]

/-- Every size below `α` is realised by an extendable independent set. -/
theorem extSets_nonempty (D : TreeMatching V) {k : ℕ} (hk : k ≤ Fintype.card D.Bag) :
    (D.extSets k).Nonempty := by
  classical
  obtain ⟨I, hI⟩ := maxIndepSets_nonempty D
  have hIcard : I.card = Fintype.card D.Bag := card_of_maxIndepSets hI
  obtain ⟨T, hTI, hTcard⟩ := Finset.exists_subset_card_eq (s := I) (n := k) (by omega)
  exact ⟨T, mem_extSets.2 ⟨hTcard, fun u hu v hv =>
    maxIndepSets_indep hI u (hTI hu) v (hTI hv), I, hI, hTI⟩⟩

/-- **Extendable incidence bound**: `4 e₄ ≤ (α - 3) e₃`. -/
theorem extendable_incidence (D : TreeMatching V) (ha : 4 ≤ Fintype.card D.Bag) :
    4 * (D.extSets (Fintype.card D.Bag - 4)).card
      ≤ (Fintype.card D.Bag - 3) * (D.extSets (Fintype.card D.Bag - 3)).card := by
  classical
  set a := Fintype.card D.Bag with hadef
  set E3 := D.extSets (a - 3) with hE3
  set E4 := D.extSets (a - 4) with hE4
  set Inc : Finset (V × Finset V) := E3.biUnion (fun S => S.image (fun v => (v, S))) with hInc
  have hmemInc1 : ∀ p ∈ Inc, p.1 ∈ p.2 := by
    intro p hp
    rw [hInc, Finset.mem_biUnion] at hp
    obtain ⟨S, -, hp2⟩ := hp
    obtain ⟨v, hv, rfl⟩ := Finset.mem_image.1 hp2
    exact hv
  have hmemInc2 : ∀ p ∈ Inc, p.2 ∈ E3 := by
    intro p hp
    rw [hInc, Finset.mem_biUnion] at hp
    obtain ⟨S, hS, hp2⟩ := hp
    obtain ⟨v, -, rfl⟩ := Finset.mem_image.1 hp2
    exact hS
  have hdisj : ∀ S ∈ E3, ∀ S' ∈ E3, S ≠ S' →
      Disjoint (S.image (fun v => (v, S))) (S'.image (fun v => (v, S'))) := by
    intro S _ S' _ hne
    rw [Finset.disjoint_left]
    rintro p hp hp'
    obtain ⟨v, -, rfl⟩ := Finset.mem_image.1 hp
    obtain ⟨v', -, hv'⟩ := Finset.mem_image.1 hp'
    exact hne (congrArg Prod.snd hv').symm
  have hcards : ∀ S ∈ E3, S.card = a - 3 := by
    intro S hS
    rw [hE3, mem_extSets] at hS
    exact hS.1
  have hcount : Inc.card = (a - 3) * E3.card := by
    rw [hInc, Finset.card_biUnion hdisj]
    rw [Finset.sum_congr rfl (fun S hS => (Finset.card_image_of_injOn
      (fun x _ y _ hxy => congrArg Prod.fst hxy)).trans (hcards S hS))]
    rw [Finset.sum_const, smul_eq_mul, mul_comm]
  have hmaps : ∀ p ∈ Inc, p.2.erase p.1 ∈ E4 := by
    intro p hp
    have h2 := hmemInc2 p hp
    have h1 := hmemInc1 p hp
    rw [hE3, mem_extSets] at h2
    obtain ⟨hc, hind, I, hI, hSI⟩ := h2
    rw [hE4, mem_extSets]
    refine ⟨?_, ?_, I, hI, (Finset.erase_subset _ _).trans hSI⟩
    · rw [Finset.card_erase_of_mem h1, hc]; omega
    · intro u hu v hv
      exact hind u (Finset.mem_of_mem_erase hu) v (Finset.mem_of_mem_erase hv)
  have hfib : ∀ T ∈ E4, 4 ≤ (Inc.filter (fun p => p.2.erase p.1 = T)).card := by
    intro T hT
    rw [hE4, mem_extSets] at hT
    obtain ⟨hTcard, hTind, I, hI, hTI⟩ := hT
    have hIcard : I.card = a := card_of_maxIndepSets hI
    have hdiff : (I \ T).card = 4 := by
      rw [Finset.card_sdiff_of_subset hTI, hIcard, hTcard]; omega
    have hsub : (I \ T).image (fun v => (v, insert v T))
        ⊆ Inc.filter (fun p => p.2.erase p.1 = T) := by
      intro p hp
      obtain ⟨v, hv, rfl⟩ := Finset.mem_image.1 hp
      rw [Finset.mem_sdiff] at hv
      have hvI : v ∈ I := hv.1
      have hvT : v ∉ T := hv.2
      have hins : insert v T ⊆ I := Finset.insert_subset hvI hTI
      have hinsE3 : insert v T ∈ E3 := by
        rw [hE3, mem_extSets]
        refine ⟨?_, fun u hu w hw => maxIndepSets_indep hI u (hins hu) w (hins hw), I, hI, hins⟩
        rw [Finset.card_insert_of_notMem hvT, hTcard]; omega
      rw [Finset.mem_filter]
      refine ⟨?_, ?_⟩
      · rw [hInc, Finset.mem_biUnion]
        exact ⟨insert v T, hinsE3, Finset.mem_image_of_mem _ (Finset.mem_insert_self v T)⟩
      · exact Finset.erase_insert hvT
    calc (4 : ℕ) = ((I \ T).image (fun v => (v, insert v T))).card := by
          rw [Finset.card_image_of_injOn (fun x _ y _ hxy => congrArg Prod.fst hxy), hdiff]
      _ ≤ _ := Finset.card_le_card hsub
  have hkey : 4 * E4.card ≤ Inc.card := by
    rw [Finset.card_eq_sum_card_fiberwise (t := E4) (fun p hp => hmaps p hp)]
    calc 4 * E4.card = ∑ _T ∈ E4, 4 := by rw [Finset.sum_const, smul_eq_mul, mul_comm]
      _ ≤ _ := Finset.sum_le_sum hfib
  omega

/-! ### Numerical assembly -/

/-- **The low-density depth-three inequality for a forest with a maximum matching.** -/
theorem depthThreeStrict_of_lowDensity (D : TreeMatching V)
    (h17 : 17 ≤ Fintype.card D.Bag) (h19 : Fintype.card D.Bag ≤ 19)
    (hden : 3 * (Fintype.card D.Bag - 3) * DepthThree.b D.G 4
      ≤ (Fintype.card D.Bag - 7) * DepthThree.e D.G 4) :
    DepthThree.DepthThreeStrict D.G := by
  classical
  set a := Fintype.card D.Bag with hadef
  have hE2 : DepthThree.e D.G 2 = D.erasure 2 := e_eq_erasure D (by omega)
  have hE3 : DepthThree.e D.G 3 = D.erasure 3 := e_eq_erasure D (by omega)
  have hE4 : DepthThree.e D.G 4 = D.erasure 4 := e_eq_erasure D (by omega)
  -- the graph Pascal reserve
  have hres : 32 * (a - 2) * (DepthThree.e D.G 2 * DepthThree.e D.G 4)
      ≤ 27 * (a - 3) * DepthThree.e D.G 3 ^ 2 := by
    rw [hE2, hE3, hE4]
    exact erasure_depth_three_reserve D (by omega)
  -- the blocked shadow bound
  have hblk : (a - 4) * DepthThree.b D.G 2 ≤ 6 * DepthThree.b D.G 3 := by
    rw [b_eq_card_blkSets D (show 2 ≤ a by omega), b_eq_card_blkSets D (show 3 ≤ a by omega)]
    exact blocked_shadow D (by omega)
  -- the extendable incidence bound
  have hext : 4 * DepthThree.e D.G 4 ≤ (a - 3) * DepthThree.e D.G 3 := by
    rw [e_eq_card_extSets D (show 4 ≤ a by omega), e_eq_card_extSets D (show 3 ≤ a by omega)]
    exact extendable_incidence D (by omega)
  -- positivity
  have hpos2 : 0 < DepthThree.e D.G 2 := by
    rw [e_eq_card_extSets D (show 2 ≤ a by omega)]
    exact Finset.card_pos.2 (extSets_nonempty D (by omega))
  have hpos4 : 0 < DepthThree.e D.G 4 := by
    rw [e_eq_card_extSets D (show 4 ≤ a by omega)]
    exact Finset.card_pos.2 (extSets_nonempty D (by omega))
  -- transport to the reals
  have hc2 : ((a - 2 : ℕ) : ℝ) = (a : ℝ) - 2 := by
    rw [Nat.cast_sub (by omega)]; norm_num
  have hc3 : ((a - 3 : ℕ) : ℝ) = (a : ℝ) - 3 := by
    rw [Nat.cast_sub (by omega)]; norm_num
  have hc4 : ((a - 4 : ℕ) : ℝ) = (a : ℝ) - 4 := by
    rw [Nat.cast_sub (by omega)]; norm_num
  have hc7 : ((a - 7 : ℕ) : ℝ) = (a : ℝ) - 7 := by
    rw [Nat.cast_sub (by omega)]; norm_num
  have key := DepthThree.low_density_algebra a h17 h19
    ((DepthThree.e D.G 2 : ℕ) : ℝ) ((DepthThree.e D.G 3 : ℕ) : ℝ) ((DepthThree.e D.G 4 : ℕ) : ℝ)
    ((DepthThree.b D.G 2 : ℕ) : ℝ) ((DepthThree.b D.G 3 : ℕ) : ℝ) ((DepthThree.b D.G 4 : ℕ) : ℝ)
    (by exact_mod_cast hpos2) (by positivity) (by exact_mod_cast hpos4)
    (by positivity) (by positivity) (by positivity)
    (by
      have := hres
      have hcast : ((32 * (a - 2) * (DepthThree.e D.G 2 * DepthThree.e D.G 4) : ℕ) : ℝ)
          ≤ ((27 * (a - 3) * DepthThree.e D.G 3 ^ 2 : ℕ) : ℝ) := by exact_mod_cast this
      push_cast [hc2, hc3] at hcast
      linarith)
    (by
      have hcast : (((a - 4) * DepthThree.b D.G 2 : ℕ) : ℝ)
          ≤ ((6 * DepthThree.b D.G 3 : ℕ) : ℝ) := by exact_mod_cast hblk
      push_cast [hc4] at hcast
      linarith)
    (by
      have hcast : ((4 * DepthThree.e D.G 4 : ℕ) : ℝ) ≤ (((a - 3) * DepthThree.e D.G 3 : ℕ) : ℝ) :=
        by exact_mod_cast hext
      push_cast [hc3] at hcast
      linarith)
    (by
      have hcast : ((3 * (a - 3) * DepthThree.b D.G 4 : ℕ) : ℝ)
          ≤ (((a - 7) * DepthThree.e D.G 4 : ℕ) : ℝ) := by exact_mod_cast hden
      push_cast [hc3, hc7] at hcast
      linarith)
  rw [DepthThree.DepthThreeStrict, DepthThree.s, DepthThree.s, DepthThree.s]
  exact_mod_cast key

end TreeMatching

end MatchingBag
