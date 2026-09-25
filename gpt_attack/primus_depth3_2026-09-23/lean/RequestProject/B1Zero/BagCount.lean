import RequestProject.B1Zero.Allowed

/-!
# Counting vertices by bag type

Every bag is of exactly one of the following three kinds:

* it contains one forced vertex and no other allowed vertex (the other vertex, if any, is
  forbidden);
* it contains two flexible vertices.

Summing the two resulting per-bag identities over all bags gives

* `|V| + |Forced| = 2 α + |Forbidden|`  (`card_verts_add_card_forced`), i.e.
  `|Forced| = |Forbidden| + (2α - |V|)`;
* `|Flex| = 2 (α - |Forced|)`           (`card_flex`).
-/

open Finset SimpleGraph

namespace MatchingBag

attribute [local instance 1] Classical.propDecidable

variable {V : Type*} [Fintype V] [DecidableEq V]

namespace TreeMatching

variable {D : TreeMatching V}

/-- Any finset of vertices is partitioned by the bags. -/
lemma card_eq_sum_inter_bagSet (S : Finset V) :
    S.card = ∑ b : D.Bag, (S ∩ D.bagSet b).card := by
  classical
  rw [Finset.card_eq_sum_card_fiberwise (f := D.bagOf) (t := (Finset.univ : Finset D.Bag))
    (fun x _ => Finset.mem_univ _)]
  refine Finset.sum_congr rfl fun b _ => ?_
  congr 1
  ext v
  simp only [Finset.mem_filter, Finset.mem_inter]
  exact and_congr_right fun _ => bagOf_eq_iff

lemma card_univ_eq_sum_bagSet :
    Fintype.card V = ∑ b : D.Bag, (D.bagSet b).card := by
  have := card_eq_sum_inter_bagSet (D := D) (Finset.univ : Finset V)
  simpa [Finset.card_univ, Finset.univ_inter] using this

/-- The two per-bag identities. -/
lemma bag_local_counts (b : D.Bag) :
    (D.bagSet b).card + (D.ForcedV ∩ D.bagSet b).card
        = 2 + (D.ForbiddenV ∩ D.bagSet b).card ∧
      (D.FlexV ∩ D.bagSet b).card + 2 * (D.ForcedV ∩ D.bagSet b).card = 2 := by
  classical
  by_cases hb : ∃ f ∈ D.bagSet b, f ∈ D.ForcedV
  · obtain ⟨f, hfb, hfF⟩ := hb
    have hFsingle : D.ForcedV ∩ D.bagSet b = {f} := by
      ext v
      simp only [Finset.mem_inter, Finset.mem_singleton]
      constructor
      · rintro ⟨hvF, hvb⟩
        by_contra hvf
        exact Finset.disjoint_left.1 Allowed_disjoint_Forbidden
          (Forced_subset_Allowed hvF) (bagSet_eq_of_forced hfb hfF hvb hvf)
      · rintro rfl; exact ⟨hfF, hfb⟩
    have hFlexEmpty : D.FlexV ∩ D.bagSet b = ∅ := by
      ext v
      simp only [Finset.mem_inter, Finset.notMem_empty, iff_false, not_and]
      intro hvFlex hvb
      have hvf : v ≠ f := by rintro rfl; exact (mem_Flex.1 hvFlex).2 hfF
      exact Finset.disjoint_left.1 Allowed_disjoint_Forbidden (mem_Flex.1 hvFlex).1
        (bagSet_eq_of_forced hfb hfF hvb hvf)
    have hFor : D.ForbiddenV ∩ D.bagSet b = (D.bagSet b).erase f := by
      ext v
      simp only [Finset.mem_inter, Finset.mem_erase]
      constructor
      · rintro ⟨hvB, hvb⟩
        refine ⟨?_, hvb⟩
        rintro rfl
        exact Finset.disjoint_left.1 Allowed_disjoint_Forbidden
          (Forced_subset_Allowed hfF) hvB
      · rintro ⟨hvf, hvb⟩
        exact ⟨bagSet_eq_of_forced hfb hfF hvb hvf, hvb⟩
    rw [hFsingle, hFlexEmpty, hFor, Finset.card_singleton, Finset.card_erase_of_mem hfb,
      Finset.card_empty]
    have := card_bagSet_pos (D := D) b
    omega
  · push_neg at hb
    obtain ⟨hcard, hflex⟩ := bagSet_flex_of_no_forced hb
    have hF : D.ForcedV ∩ D.bagSet b = ∅ := by
      ext v; simp only [Finset.mem_inter, Finset.notMem_empty, iff_false, not_and]
      intro hvF hvb; exact hb v hvb hvF
    have hB : D.ForbiddenV ∩ D.bagSet b = ∅ := by
      ext v; simp only [Finset.mem_inter, Finset.notMem_empty, iff_false, not_and']
      intro hvb hvB
      exact Finset.disjoint_left.1 Allowed_disjoint_Forbidden (mem_Flex.1 (hflex v hvb)).1 hvB
    have hL : D.FlexV ∩ D.bagSet b = D.bagSet b := by
      apply Finset.inter_eq_right.2
      intro v hv; exact hflex v hv
    rw [hF, hB, hL, hcard]
    simp

/-- `|V| + |Forced| = 2 α + |Forbidden|`. -/
theorem card_verts_add_card_forced (D : TreeMatching V) :
    Fintype.card V + D.ForcedV.card = 2 * Fintype.card D.Bag + D.ForbiddenV.card := by
  classical
  rw [card_univ_eq_sum_bagSet (D := D), card_eq_sum_inter_bagSet (D := D) D.ForcedV,
    card_eq_sum_inter_bagSet (D := D) D.ForbiddenV, ← Finset.sum_add_distrib]
  rw [show 2 * Fintype.card D.Bag = ∑ _b : D.Bag, 2 by
    simp [Finset.sum_const, Finset.card_univ, mul_comm], ← Finset.sum_add_distrib]
  exact Finset.sum_congr rfl fun b _ => (bag_local_counts b).1

/-- `|Flex| = 2 (α - |Forced|)`. -/
theorem card_flex (D : TreeMatching V) :
    D.FlexV.card + 2 * D.ForcedV.card = 2 * Fintype.card D.Bag := by
  classical
  rw [card_eq_sum_inter_bagSet (D := D) D.FlexV, card_eq_sum_inter_bagSet (D := D) D.ForcedV,
    Finset.mul_sum, ← Finset.sum_add_distrib]
  rw [show 2 * Fintype.card D.Bag = ∑ _b : D.Bag, 2 by
    simp [Finset.sum_const, Finset.card_univ, mul_comm]]
  exact Finset.sum_congr rfl fun b _ => (bag_local_counts b).2

lemma card_forced_le_bags (D : TreeMatching V) : D.ForcedV.card ≤ Fintype.card D.Bag := by
  have := card_flex D
  omega

lemma card_forbidden_le_forced (D : TreeMatching V) : D.ForbiddenV.card ≤ D.ForcedV.card := by
  have h1 := card_verts_add_card_forced D
  have h2 := card_univ_eq_sum_bagSet (D := D)
  have h3 : ∑ b : D.Bag, (D.bagSet b).card ≤ ∑ _b : D.Bag, 2 :=
    Finset.sum_le_sum fun b _ => card_bagSet_le b
  simp only [Finset.sum_const, Finset.card_univ, smul_eq_mul] at h3
  omega

end TreeMatching

end MatchingBag
