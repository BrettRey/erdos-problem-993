import RequestProject.Forest

/-!
# The strong induction

Every nonempty eligible vertex set of a forest either satisfies both fivefold bounds
or has the exceptional profile of the four-leaf star.
-/

open Finset

namespace FivefoldForest

open scoped Classical

variable {V : Type*} {G : SimpleGraph V}

/-- **Main induction.**  For a forest, every nonempty vertex set with no blocked sets at
defects one and two either satisfies both fivefold bounds, or is a four-leaf star. -/
theorem good_or_starProf (hF : ForestLike G) :
    ∀ (n : ℕ) (S : Finset V), S.card ≤ n → Eligible G S → Good G S ∨ StarProf G S := by
  intro n
  induction n with
  | zero =>
    intro S hcard helig
    exact absurd (Finset.card_pos.2 helig.1) (by omega)
  | succ n ih =>
    intro S hcard helig
    have IH : ∀ T : Finset V, T.card < S.card → Eligible G T → Good G T ∨ StarProf G T := by
      intro T hT helT
      exact ih T (by omega) helT
    rcases hF.decomp S helig.1 with ⟨c, hstar⟩ | ⟨S₁, S₂, hsp, hne1, hne2⟩ | ⟨p, w, L, hp, hGs2⟩
    · exact star_good_or_star hstar helig
    · -- disjoint union of two smaller parts
      have hsub1 : S₁ ⊆ S := by rw [hsp.eq_union]; exact Finset.subset_union_left
      have hsub2 : S₂ ⊆ S := by rw [hsp.eq_union]; exact Finset.subset_union_right
      have hcard1 : S₁.card < S.card := by
        have : S₁.card + S₂.card = S.card := by
          rw [hsp.eq_union, Finset.card_union_of_disjoint hsp.disj]
        have := Finset.card_pos.2 hne2
        omega
      have hcard2 : S₂.card < S.card := by
        have : S₁.card + S₂.card = S.card := by
          rw [hsp.eq_union, Finset.card_union_of_disjoint hsp.disj]
        have := Finset.card_pos.2 hne1
        omega
      have he1 : Eligible G S₁ :=
        ⟨hne1, Nat.le_zero.mp (helig.2.1 ▸ bD_le_of_split hsp 1),
          Nat.le_zero.mp (helig.2.2 ▸ bD_le_of_split hsp 2)⟩
      have hsp' : Split G S S₂ S₁ :=
        ⟨by rw [hsp.eq_union, Finset.union_comm], hsp.disj.symm,
          fun x hx y hy hadj => hsp.noEdges y hy x hx hadj.symm⟩
      have he2 : Eligible G S₂ :=
        ⟨hne2, Nat.le_zero.mp (helig.2.1 ▸ bD_le_of_split hsp' 1),
          Nat.le_zero.mp (helig.2.2 ▸ bD_le_of_split hsp' 2)⟩
      exact Or.inl (split_good hsp he1 he2 (hF.card_le _) (hF.card_le _)
        (IH S₁ hcard1 he1) (IH S₂ hcard2 he2))
    · -- a pendant configuration
      refine Or.inl ?_
      rcases lt_trichotomy L.card 2 with h1 | h2 | h3
      · have : L.card = 1 := by
          have := Finset.card_pos.2 hp.Lne
          omega
        exact pendant_one_good hp this helig hGs2 IH
      · -- two leaves are impossible
        obtain ⟨-, -, h3⟩ := pendant_many_setup hp (by omega) helig
        omega
      · rcases eq_or_lt_of_le (show 3 ≤ L.card by omega) with h3' | h4
        · exact pendant_three_good hp (by omega) helig hGs2 hF.card_le IH
        · exact pendant_four_le_good hp (by omega) helig hF.card_le IH

/-- Blocked sets of size `α - 1` cannot exist once those of size `α - 2` all extend
(Lemma 3 of the source note). -/
theorem bD_one_eq_zero_of_bD_two (hF : ForestLike G) {S : Finset V} (h4 : 4 ≤ alpha G S)
    (h2 : bD G S 2 = 0) : bD G S 1 = 0 := by
  classical
  by_contra hne
  rw [bD, if_pos (by omega : (1 : ℕ) ≤ alpha G S), bcnt] at hne
  obtain ⟨S₀, hS₀mem⟩ := Finset.card_ne_zero.1 hne
  obtain ⟨hS₀, hS₀card, hS₀not⟩ := Finset.mem_filter.1 hS₀mem
  have hS₀sub : S₀ ⊆ S := (mem_indFam.1 hS₀).1
  have hS₀ind : G.IsIndepSet (S₀ : Set V) := (mem_indFam.1 hS₀).2
  -- a blocked set of size `α - 1` is a maximal independent set
  have hmaximal : ∀ z ∈ S, z ∉ S₀ → ∃ u ∈ S₀, G.Adj z u := by
    intro z hz hzS₀
    by_contra hcon
    push_neg at hcon
    refine hS₀not ⟨insert z S₀, mem_indFam.2 ⟨Finset.insert_subset hz hS₀sub, ?_⟩, ?_,
      Finset.subset_insert _ _⟩
    · intro x hx y hy hxy
      simp only [Finset.coe_insert, Set.mem_insert_iff, Finset.mem_coe] at hx hy
      rcases hx with rfl | hx <;> rcases hy with rfl | hy
      · exact absurd rfl hxy
      · exact hcon y hy
      · exact fun hadj => hcon x hx hadj.symm
      · exact hS₀ind (by exact_mod_cast hx) (by exact_mod_cast hy) hxy
    · rw [Finset.card_insert_of_notMem hzS₀, hS₀card]; omega
  -- every vertex of the blocked set supports two private outside vertices
  have key : ∀ v : V, ∃ D : Finset V, v ∈ S₀ →
      2 ≤ D.card ∧ D ⊆ S \ S₀ ∧ (∀ z ∈ D, G.Adj z v) ∧
        (∀ z ∈ D, ∀ u ∈ S₀, G.Adj z u → u = v) := by
    intro v
    by_cases hv : v ∈ S₀
    swap
    · exact ⟨∅, fun h => absurd h hv⟩
    obtain ⟨M, hM, hMc, hsubM⟩ := ext_of_bD_zero (d := 2) (by omega) h2
      (indFam_subset hS₀ (Finset.erase_subset _ _))
      (by rw [Finset.card_erase_of_mem hv, hS₀card]; omega)
    have hvM : v ∉ M := by
      intro hvM
      refine hS₀not ⟨M, hM, hMc, fun x hx => ?_⟩
      by_cases hxv : x = v
      · exact hxv ▸ hvM
      · exact hsubM (Finset.mem_erase.2 ⟨hxv, hx⟩)
    have hMsub : M ⊆ S := (mem_indFam.1 hM).1
    have hMind : G.IsIndepSet (M : Set V) := (mem_indFam.1 hM).2
    have hnb : ∀ z ∈ M \ S₀, ∀ u ∈ S₀, G.Adj z u → u = v := by
      intro z hz u hu hadj
      by_contra huv
      obtain ⟨hzM, hzS₀⟩ := Finset.mem_sdiff.1 hz
      have huM : u ∈ M := hsubM (Finset.mem_erase.2 ⟨huv, hu⟩)
      exact hMind (by exact_mod_cast hzM) (by exact_mod_cast huM)
        (fun hh => hzS₀ (by rw [hh]; exact hu)) hadj
    refine ⟨M \ S₀, fun _ => ⟨?_, ?_, ?_, hnb⟩⟩
    · have hinter : M ∩ S₀ ⊆ S₀.erase v := by
        intro x hx
        obtain ⟨hxM, hxS₀⟩ := Finset.mem_inter.1 hx
        exact Finset.mem_erase.2 ⟨fun hh => hvM (hh ▸ hxM), hxS₀⟩
      have hsum : (M ∩ S₀).card + (M \ S₀).card = M.card :=
        Finset.card_inter_add_card_sdiff _ _
      have hle : (M ∩ S₀).card ≤ (S₀.erase v).card := Finset.card_le_card hinter
      rw [Finset.card_erase_of_mem hv, hS₀card] at hle
      omega
    · intro z hz
      obtain ⟨hzM, hzS₀⟩ := Finset.mem_sdiff.1 hz
      exact Finset.mem_sdiff.2 ⟨hMsub hzM, hzS₀⟩
    · intro z hz
      obtain ⟨hzM, hzS₀⟩ := Finset.mem_sdiff.1 hz
      obtain ⟨u, hu, hadj⟩ := hmaximal z (hMsub hzM) hzS₀
      exact (hnb z hz u hu hadj) ▸ hadj
  choose D hD using key
  have hdisj : (↑S₀ : Set V).PairwiseDisjoint D := by
    intro v hv v' hv' hvv'
    simp only [Finset.mem_coe] at hv hv'
    refine Finset.disjoint_left.2 fun z hz hz' => ?_
    have hzv := (hD v hv).2.2.1 z hz
    exact hvv' ((hD v' hv').2.2.2 z hz' v hv hzv)
  have hbU : S₀.biUnion D ⊆ S \ S₀ := by
    intro z hz
    obtain ⟨v, hv, hzv⟩ := Finset.mem_biUnion.1 hz
    exact (hD v hv).2.1 hzv
  have hcard1 : 2 * S₀.card ≤ (S₀.biUnion D).card := by
    rw [Finset.card_biUnion hdisj]
    calc 2 * S₀.card = ∑ _v ∈ S₀, 2 := by rw [Finset.sum_const, smul_eq_mul, mul_comm]
      _ ≤ ∑ v ∈ S₀, (D v).card := Finset.sum_le_sum fun v hv => (hD v hv).1
  have hcard2 : (S₀.biUnion D).card ≤ (S \ S₀).card := Finset.card_le_card hbU
  have hcard3 : (S \ S₀).card + S₀.card = S.card := Finset.card_sdiff_add_card_eq_card hS₀sub
  have hcard4 := hF.card_le S
  have hcard5 : alpha G S ≤ S.card := alpha_le_card G S
  omega

/-- **The fivefold bounds for forests.** -/
theorem fivefold (hF : ForestLike G) {S : Finset V} (h5 : 5 ≤ alpha G S)
    (h2 : bD G S 2 = 0) : 5 * bD G S 3 ≤ eD G S 3 ∧ 5 * bD G S 4 ≤ eD G S 4 := by
  have h1 : bD G S 1 = 0 := bD_one_eq_zero_of_bD_two hF (by omega) h2
  have hne : S.Nonempty := by
    rw [← Finset.card_pos]
    have := alpha_le_card G S
    omega
  rcases good_or_starProf hF S.card S le_rfl ⟨hne, h1, h2⟩ with ⟨g3, g4⟩ | s
  · constructor
    · have : (5 : ℤ) * (bD G S 3 : ℤ) ≤ (eD G S 3 : ℤ) := by
        simp only [margin] at g3; linarith
      exact_mod_cast this
    · have : (5 : ℤ) * (bD G S 4 : ℤ) ≤ (eD G S 4 : ℤ) := by
        simp only [margin] at g4; linarith
      exact_mod_cast this
  · exact absurd s.1 (by omega)

end FivefoldForest
