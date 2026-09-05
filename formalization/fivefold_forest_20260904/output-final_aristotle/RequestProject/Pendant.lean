import RequestProject.Closure

/-!
# Pendant configurations

A *pendant configuration* in `S` consists of a vertex `p ∈ S` whose neighbours inside
`S` are a nonempty set `L` of leaves of `S` together with exactly one further vertex `w`.

Notation used throughout:

* `pBase S p L = S \ insert p L`  — the tree that remains after deleting `p` and its leaves;
* `pRest S p = S.erase p`         — `S` with `p` deleted (the disjoint union of `pBase` and `L`);
* `pCore S p w L = (S \ insert p L).erase w` — the base with the attachment vertex removed.
-/

open Finset

namespace FivefoldForest

open scoped Classical

variable {V : Type*} {G : SimpleGraph V} {S A : Finset V} {p w : V} {L : Finset V}

/-- The base graph: what remains of `S` after deleting `p` and its leaves. -/
noncomputable def pBase (S : Finset V) (p : V) (L : Finset V) : Finset V := S \ insert p L

/-- `S` with the support vertex `p` deleted. -/
noncomputable def pRest (S : Finset V) (p : V) : Finset V := S.erase p

/-- The base graph with the attachment vertex `w` deleted. -/
noncomputable def pCore (S : Finset V) (p w : V) (L : Finset V) : Finset V :=
  (pBase S p L).erase w

lemma pCore_eq (S : Finset V) (p w : V) (L : Finset V) :
    pCore S p w L = (pBase S p L).erase w := rfl

/-- A pendant configuration: `p` is adjacent to the nonempty set `L` of leaves of `S`
and to exactly one other vertex `w`. -/
structure Pendant (G : SimpleGraph V) (S : Finset V) (p w : V) (L : Finset V) : Prop where
  pmem : p ∈ S
  wmem : w ∈ S
  adj : G.Adj p w
  wnotL : w ∉ L
  Lsub : L ⊆ S
  pnotL : p ∉ L
  Lne : L.Nonempty
  leafAdj : ∀ x ∈ L, G.Adj p x
  leafOnly : ∀ x ∈ L, ∀ y ∈ S, G.Adj x y → y = p
  pNbrs : ∀ y ∈ S, G.Adj p y → y = w ∨ y ∈ L

namespace Pendant

lemma wmem_base (hp : Pendant G S p w L) : w ∈ pBase S p L := by
  simp only [pBase, Finset.mem_sdiff, Finset.mem_insert]
  exact ⟨hp.wmem, by rintro (rfl | h); · exact hp.adj.ne' rfl
                     · exact hp.wnotL h⟩

lemma base_subset (S : Finset V) (p : V) (L : Finset V) : pBase S p L ⊆ S := Finset.sdiff_subset

lemma core_subset_base (S : Finset V) (p w : V) (L : Finset V) :
    pCore S p w L ⊆ pBase S p L := Finset.erase_subset _ _

lemma base_subset_erase (S : Finset V) (p : V) (L : Finset V) :
    pBase S p L ⊆ S.erase p := by
  intro x hx
  simp only [pBase, Finset.mem_sdiff, Finset.mem_insert] at hx
  exact Finset.mem_erase.2 ⟨fun h => hx.2 (Or.inl h), hx.1⟩

lemma base_card_lt (hp : Pendant G S p w L) : (pBase S p L).card < S.card := by
  refine lt_of_le_of_lt (Finset.card_le_card (base_subset_erase S p L)) ?_
  exact Finset.card_erase_lt_of_mem hp.pmem

lemma core_card_lt (hp : Pendant G S p w L) : (pCore S p w L).card < S.card :=
  lt_of_le_of_lt (Finset.card_le_card (core_subset_base S p w L)) hp.base_card_lt

/-- `S` without `p` splits as the base together with the leaves. -/
lemma split (hp : Pendant G S p w L) : Split G (pRest S p) (pBase S p L) L := by
  refine ⟨?_, ?_, ?_⟩
  · ext x
    simp only [pRest, pBase, Finset.mem_erase, Finset.mem_union, Finset.mem_sdiff,
      Finset.mem_insert]
    constructor
    · rintro ⟨hxp, hxS⟩
      by_cases hxL : x ∈ L
      · exact Or.inr hxL
      · exact Or.inl ⟨hxS, by rintro (rfl | hL); exacts [hxp rfl, hxL hL]⟩
    · rintro (⟨hxS, hx⟩ | hxL)
      · exact ⟨fun hh => hx (Or.inl hh), hxS⟩
      · exact ⟨fun hh => hp.pnotL (hh ▸ hxL), hp.Lsub hxL⟩
  · refine Finset.disjoint_left.2 fun x hx hxL => ?_
    simp only [pBase, Finset.mem_sdiff, Finset.mem_insert] at hx
    exact hx.2 (Or.inr hxL)
  · intro x hx y hy hadj
    simp only [pBase, Finset.mem_sdiff, Finset.mem_insert] at hx
    exact hx.2 (Or.inl (hp.leafOnly y hy x hx.1 hadj.symm))

/-- The leaves are pairwise non-adjacent. -/
lemma leaves_indep (hp : Pendant G S p w L) : G.IsIndepSet (L : Set V) := by
  intro x hx y hy _ hadj
  have hyp : y = p := hp.leafOnly x (by exact_mod_cast hx) y (hp.Lsub (by exact_mod_cast hy)) hadj
  exact hp.pnotL (hyp ▸ (by exact_mod_cast hy : y ∈ L))

lemma leaves_indFam (hp : Pendant G S p w L) : indFam G L = L.powerset := by
  ext A
  simp only [mem_indFam, Finset.mem_powerset]
  exact ⟨fun h => h.1, fun h => ⟨h, indep_of_subset hp.leaves_indep h⟩⟩

lemma leaves_mem (hp : Pendant G S p w L) : L ∈ indFam G L := by
  rw [hp.leaves_indFam]; exact Finset.mem_powerset_self L

/-- The leaves form an edgeless set. -/
lemma leaves_alpha (hp : Pendant G S p w L) : alpha G L = L.card :=
  le_antisymm (alpha_le_card G L) (card_le_alpha hp.leaves_mem)

lemma leaves_ext (hp : Pendant G S p w L) {A : Finset V} (hA : A ⊆ L) : Ext G L A :=
  ⟨L, hp.leaves_mem, hp.leaves_alpha.symm, hA⟩

lemma leaves_ecnt (hp : Pendant G S p w L) (k : ℕ) : ecnt G L k = L.card.choose k := by
  have hset : ((indFam G L).filter fun A => A.card = k ∧ Ext G L A) = Finset.powersetCard k L := by
    ext A
    simp only [Finset.mem_filter, Finset.mem_powersetCard, hp.leaves_indFam, Finset.mem_powerset]
    exact ⟨fun h => ⟨h.1, h.2.1⟩, fun h => ⟨h.1, h.2, hp.leaves_ext h.1⟩⟩
  rw [ecnt, hset, Finset.card_powersetCard]

lemma leaves_eD (hp : Pendant G S p w L) (d : ℕ) : eD G L d = L.card.choose d := by
  rw [eD, hp.leaves_alpha]
  split
  · rename_i hd
    rw [hp.leaves_ecnt, Nat.choose_symm hd]
  · rename_i hd
    exact (Nat.choose_eq_zero_of_lt (by omega)).symm

lemma leaves_bD (hp : Pendant G S p w L) (d : ℕ) : bD G L d = 0 := by
  rw [bD]
  split
  · rw [bcnt]
    refine Finset.card_eq_zero.2 (Finset.eq_empty_of_forall_notMem fun A hA => ?_)
    simp only [Finset.mem_filter, hp.leaves_indFam, Finset.mem_powerset] at hA
    exact hA.2.2 (hp.leaves_ext hA.1)
  · rfl

lemma leaves_eligible (hp : Pendant G S p w L) : Eligible G L := by
  refine ⟨hp.Lne, hp.leaves_bD 1, hp.leaves_bD 2⟩

/-- Replacing the support vertex by all of its leaves. -/
lemma lift_of_mem (hp : Pendant G S p w L) {M : Finset V} (hM : M ∈ indFam G S) (hpM : p ∈ M) :
    M.erase p ∪ L ∈ indFam G (pRest S p) ∧ (M.erase p ∪ L).card + 1 = M.card + L.card := by
  rw [mem_indFam] at hM
  simp only [pRest]
  have hdisj : Disjoint (M.erase p) L := by
    refine Finset.disjoint_left.2 fun x hx hxL => ?_
    obtain ⟨hxp, hxM⟩ := Finset.mem_erase.1 hx
    exact hM.2 (by exact_mod_cast hpM) (by exact_mod_cast hxM) (Ne.symm hxp) (hp.leafAdj x hxL)
  constructor
  · rw [mem_indFam]
    refine ⟨?_, ?_⟩
    · intro x hx
      rcases Finset.mem_union.1 hx with hx | hx
      · obtain ⟨hxp, hxM⟩ := Finset.mem_erase.1 hx
        exact Finset.mem_erase.2 ⟨hxp, hM.1 hxM⟩
      · exact Finset.mem_erase.2 ⟨fun hh => hp.pnotL (hh ▸ hx), hp.Lsub hx⟩
    · intro x hx y hy hxy hadj
      simp only [Finset.coe_union, Set.mem_union, Finset.mem_coe] at hx hy
      rcases hx with hx | hx <;> rcases hy with hy | hy
      · exact hM.2 (by exact_mod_cast (Finset.mem_erase.1 hx).2)
          (by exact_mod_cast (Finset.mem_erase.1 hy).2) hxy hadj
      · exact (Finset.mem_erase.1 hx).1
          (hp.leafOnly y hy x (hM.1 (Finset.mem_erase.1 hx).2) hadj.symm)
      · exact (Finset.mem_erase.1 hy).1
          (hp.leafOnly x hx y (hM.1 (Finset.mem_erase.1 hy).2) hadj)
      · exact hp.pnotL (hp.leafOnly x hx y (hp.Lsub hy) hadj ▸ hy)
  · rw [Finset.card_union_of_disjoint hdisj, Finset.card_erase_of_mem hpM]
    have hMpos : 1 ≤ M.card := Finset.card_pos.2 ⟨p, hpM⟩
    omega

/-- Any independent set of `S` can be pushed off the support vertex `p` without losing size. -/
lemma lift (hp : Pendant G S p w L) {M : Finset V} (hM : M ∈ indFam G S) :
    ∃ M' ∈ indFam G (pRest S p), M.card ≤ M'.card ∧ M.erase p ⊆ M' := by
  by_cases hpM : p ∈ M
  · obtain ⟨h1, h2⟩ := hp.lift_of_mem hM hpM
    have hL := Finset.card_pos.2 hp.Lne
    exact ⟨_, h1, by omega, Finset.subset_union_left⟩
  · rw [mem_indFam] at hM
    exact ⟨M, mem_indFam.2 ⟨fun x hx => Finset.mem_erase.2 ⟨fun hh => hpM (hh ▸ hx), hM.1 hx⟩,
      hM.2⟩, le_rfl, Finset.erase_subset _ _⟩

/-- Deleting `p` does not change the independence number. -/
lemma alpha_rest (hp : Pendant G S p w L) : alpha G (pRest S p) = alpha G S := by
  refine le_antisymm (alpha_mono (Finset.erase_subset _ _)) ?_
  obtain ⟨M, hM, hcard⟩ := exists_max_indFam G S
  obtain ⟨M', hM', hle, -⟩ := hp.lift hM
  calc alpha G S = M.card := hcard.symm
    _ ≤ M'.card := hle
    _ ≤ alpha G (pRest S p) := card_le_alpha hM'

/-- Independent sets avoiding `p` have the same extendability in `S` and in `S \ {p}`. -/
lemma ext_rest_iff (hp : Pendant G S p w L) (hA : A ⊆ pRest S p) :
    Ext G S A ↔ Ext G (pRest S p) A := by
  constructor
  · rintro ⟨M, hM, hcard, hAM⟩
    obtain ⟨M', hM', hle, hsub⟩ := hp.lift hM
    refine ⟨M', hM', le_antisymm (card_le_alpha hM') ?_, ?_⟩
    · rw [hp.alpha_rest, ← hcard]; exact hle
    · intro x hx
      exact hsub (Finset.mem_erase.2 ⟨(Finset.mem_erase.1 (hA hx)).1, hAM hx⟩)
  · rintro ⟨M, hM, hcard, hAM⟩
    exact ⟨M, indFam_mono (Finset.erase_subset _ _) hM, by rw [hcard, hp.alpha_rest], hAM⟩

end Pendant

end FivefoldForest
