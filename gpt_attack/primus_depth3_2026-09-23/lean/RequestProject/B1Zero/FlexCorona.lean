import RequestProject.B1Zero.CoronaMove
import RequestProject.B1Zero.SimpleCertificate

/-!
# The flexible part of a `b₁ = 0` forest is a corona

Assuming `b₁ = 0`, every flexible bag contains a vertex with no allowed neighbour outside
its bag (`exists_leaf_in_flex_bag`).  Choosing such a vertex as the pendant of its bag turns
the flexible vertices into the corona over the remaining `|flexBags|` base vertices.

The main result is `flex_certificate`: with ten flexible bags, the codimension counts of the
flexible part satisfy the union-bound inequality certified in `SimpleCertificate`, where the
avoided set is the set `flexAttach` of flexible vertices adjacent to a forbidden vertex.
-/

open Finset SimpleGraph

namespace MatchingBag

open B1Zero DepthThree.RootedForestCertificate

attribute [local instance 1] Classical.propDecidable

variable {V : Type*} [Fintype V] [DecidableEq V]

namespace B1Zero

omit [Fintype V] [DecidableEq V] in
/-- A set closed under adjacency contains the endpoint of every walk starting in it. -/
lemma mem_of_walk_closed {G : SimpleGraph V} {K : Finset V}
    (hK : ∀ z ∈ K, ∀ w, G.Adj z w → w ∈ K) {a b : V} (p : G.Walk a b) : a ∈ K → b ∈ K := by
  induction p with
  | nil => exact fun ha => ha
  | cons h _ ih => exact fun ha => ih (hK _ ha _ h)

end B1Zero

namespace TreeMatching

variable {D : TreeMatching V}

/-- The chosen pendant vertex of a bag. -/
noncomputable def flexLeaf (D : TreeMatching V) (b : D.Bag) : V :=
  if h : ∃ x, x ∈ D.bagSet b ∧ ∀ w, D.G.Adj x w → w ∈ D.AllowedV → w ∈ D.bagSet b
  then h.choose else (bagSet_nonempty (D := D) b).choose

/-- The chosen base vertex of a bag. -/
noncomputable def flexBase (D : TreeMatching V) (b : D.Bag) : V :=
  if h : ((D.bagSet b).erase (D.flexLeaf b)).Nonempty then h.choose else D.flexLeaf b

lemma flexLeaf_spec (hb1 : D.B1Zero) {b : D.Bag} (hb : b ∈ D.flexBags) :
    D.flexLeaf b ∈ D.bagSet b ∧
      ∀ w, D.G.Adj (D.flexLeaf b) w → w ∈ D.AllowedV → w ∈ D.bagSet b := by
  have h : ∃ x, x ∈ D.bagSet b ∧ ∀ w, D.G.Adj x w → w ∈ D.AllowedV → w ∈ D.bagSet b := by
    obtain ⟨x, hx, hx'⟩ := exists_leaf_in_flex_bag hb1 (flexBags_all_flex hb)
    exact ⟨x, hx, hx'⟩
  rw [flexLeaf, dif_pos h]
  exact h.choose_spec

lemma card_bagSet_flexBags {b : D.Bag} (hb : b ∈ D.flexBags) : (D.bagSet b).card = 2 := by
  refine (bagSet_flex_of_no_forced (fun v hv hvF => ?_)).1
  have h := mem_flexBags.1 hb
  rw [Finset.eq_empty_iff_forall_notMem] at h
  exact h v (Finset.mem_inter.2 ⟨hvF, hv⟩)

lemma flexBase_spec (hb1 : D.B1Zero) {b : D.Bag} (hb : b ∈ D.flexBags) :
    D.flexBase b ∈ D.bagSet b ∧ D.flexBase b ≠ D.flexLeaf b := by
  have hleaf := (flexLeaf_spec hb1 hb).1
  have hcard := card_bagSet_flexBags hb
  have hne : ((D.bagSet b).erase (D.flexLeaf b)).Nonempty := by
    rw [← Finset.card_pos, Finset.card_erase_of_mem hleaf]
    omega
  rw [flexBase, dif_pos hne]
  have h := hne.choose_spec
  rw [Finset.mem_erase] at h
  exact ⟨h.2, h.1⟩

lemma bagOf_flexBase (hb1 : D.B1Zero) {b : D.Bag} (hb : b ∈ D.flexBags) :
    D.bagOf (D.flexBase b) = b := bagOf_eq_of_mem (flexBase_spec hb1 hb).1

lemma bagOf_flexLeaf (hb1 : D.B1Zero) {b : D.Bag} (hb : b ∈ D.flexBags) :
    D.bagOf (D.flexLeaf b) = b := bagOf_eq_of_mem (flexLeaf_spec hb1 hb).1

lemma bagSet_eq_pair (hb1 : D.B1Zero) {b : D.Bag} (hb : b ∈ D.flexBags) :
    D.bagSet b = {D.flexBase b, D.flexLeaf b} := by
  have hbase := (flexBase_spec hb1 hb).1
  have hne := (flexBase_spec hb1 hb).2
  have hleaf := (flexLeaf_spec hb1 hb).1
  have hsub : ({D.flexBase b, D.flexLeaf b} : Finset V) ⊆ D.bagSet b := by
    intro x hx
    rcases Finset.mem_insert.1 hx with rfl | hx
    · exact hbase
    · rw [Finset.mem_singleton] at hx
      exact hx ▸ hleaf
  refine (Finset.eq_of_subset_of_card_le hsub ?_).symm
  rw [Finset.card_insert_of_notMem (by simpa using hne), Finset.card_singleton,
    card_bagSet_flexBags hb]

lemma flexBase_mem_Flex (hb1 : D.B1Zero) {b : D.Bag} (hb : b ∈ D.flexBags) :
    D.flexBase b ∈ D.FlexV := flexBags_all_flex hb _ (flexBase_spec hb1 hb).1

lemma flexLeaf_mem_Flex (hb1 : D.B1Zero) {b : D.Bag} (hb : b ∈ D.flexBags) :
    D.flexLeaf b ∈ D.FlexV := flexBags_all_flex hb _ (flexLeaf_spec hb1 hb).1

/-- The base vertices of the flexible bags. -/
noncomputable def flexBaseSet (D : TreeMatching V) : Finset V := D.flexBags.image D.flexBase

lemma mem_flexBaseSet (hb1 : D.B1Zero) {x : V} :
    x ∈ D.flexBaseSet ↔ D.bagOf x ∈ D.flexBags ∧ x = D.flexBase (D.bagOf x) := by
  constructor
  · intro hx
    obtain ⟨b, hb, rfl⟩ := Finset.mem_image.1 hx
    rw [bagOf_flexBase hb1 hb]
    exact ⟨hb, rfl⟩
  · rintro ⟨hb, hx⟩
    exact Finset.mem_image.2 ⟨D.bagOf x, hb, hx.symm⟩

lemma card_flexBaseSet (hb1 : D.B1Zero) : D.flexBaseSet.card = D.flexBags.card := by
  refine Finset.card_image_of_injOn ?_
  intro b hb b' hb' heq
  rw [Finset.mem_coe] at hb hb'
  rw [← bagOf_flexBase hb1 hb, ← bagOf_flexBase hb1 hb', heq]

/-- The corona structure carried by the flexible vertices. -/
noncomputable def flexCorona (D : TreeMatching V) (hb1 : D.B1Zero) : CoronaData D.G where
  base := D.flexBaseSet
  pend := fun v => D.flexLeaf (D.bagOf v)
  pend_notMem := by
    intro b hb hmem
    rw [mem_flexBaseSet hb1] at hb hmem
    rw [bagOf_flexLeaf hb1 hb.1] at hmem
    exact (flexBase_spec hb1 hb.1).2 hmem.2.symm
  pend_inj := by
    intro b hb b' hb' heq
    rw [mem_flexBaseSet hb1] at hb hb'
    have h : D.bagOf b = D.bagOf b' := by
      rw [← bagOf_flexLeaf hb1 hb.1, ← bagOf_flexLeaf hb1 hb'.1, heq]
    rw [hb.2, hb'.2, h]
  adj_pend := by
    intro b hb
    have hb' := (mem_flexBaseSet hb1).1 hb
    refine bagSet_adj (b := D.bagOf b) (mem_bagSet_bagOf b) (flexLeaf_spec hb1 hb'.1).1 ?_
    intro heq
    exact (flexBase_spec hb1 hb'.1).2 (by rw [← hb'.2]; exact heq)
  pend_nbr := by
    intro b hb w hw hadj
    have hbb := (mem_flexBaseSet hb1).1 hb
    have hwA : w ∈ D.AllowedV := by
      rcases Finset.mem_union.1 hw with hw | hw
      · have hw' := (mem_flexBaseSet hb1).1 hw
        have hwF : w ∈ D.FlexV := by rw [hw'.2]; exact flexBase_mem_Flex hb1 hw'.1
        exact (mem_Flex.1 hwF).1
      · obtain ⟨b', hb', rfl⟩ := Finset.mem_image.1 hw
        exact (mem_Flex.1 (flexLeaf_mem_Flex hb1 ((mem_flexBaseSet hb1).1 hb').1)).1
    have hwb := (flexLeaf_spec hb1 hbb.1).2 w hadj hwA
    rw [bagSet_eq_pair hb1 hbb.1, Finset.mem_insert, Finset.mem_singleton] at hwb
    rcases hwb with rfl | rfl
    · exact hbb.2.symm
    · exact absurd hadj (D.G.irrefl)

lemma flexCorona_base (hb1 : D.B1Zero) : (D.flexCorona hb1).base = D.flexBaseSet := rfl

lemma flexCorona_pend (hb1 : D.B1Zero) (v : V) :
    (D.flexCorona hb1).pend v = D.flexLeaf (D.bagOf v) := rfl

/-- The corona over the base vertices is exactly the set of flexible vertices. -/
theorem cor_flexBaseSet (hb1 : D.B1Zero) :
    (D.flexCorona hb1).cor (D.flexCorona hb1).base = D.FlexV := by
  ext x
  rw [CoronaData.mem_cor]
  constructor
  · rintro (hx | ⟨b, hb, rfl⟩)
    · rw [flexCorona_base, mem_flexBaseSet hb1] at hx
      exact hx.2 ▸ flexBase_mem_Flex hb1 hx.1
    · rw [flexCorona_base, mem_flexBaseSet hb1] at hb
      rw [flexCorona_pend]
      exact flexLeaf_mem_Flex hb1 hb.1
  · intro hx
    have hb : D.bagOf x ∈ D.flexBags := bagOf_mem_flexBags hx
    have hmem : x ∈ D.bagSet (D.bagOf x) := mem_bagSet_bagOf x
    rw [bagSet_eq_pair hb1 hb, Finset.mem_insert, Finset.mem_singleton] at hmem
    rcases hmem with hx' | hx'
    · exact Or.inl (by rw [flexCorona_base, mem_flexBaseSet hb1]; exact ⟨hb, hx'⟩)
    · refine Or.inr ⟨D.flexBase (D.bagOf x), ?_, ?_⟩
      · rw [flexCorona_base, mem_flexBaseSet hb1, bagOf_flexBase hb1 hb]
        exact ⟨hb, rfl⟩
      · rw [flexCorona_pend, bagOf_flexBase hb1 hb]
        exact hx'.symm

/-- The flexible vertices adjacent to a forbidden vertex. -/
noncomputable def flexAttach (D : TreeMatching V) : Finset V :=
  D.FlexV.filter (fun z => ∃ u ∈ D.ForbiddenV, D.G.Adj u z)

/-- Every flexible vertex is joined inside the flexible part to an attachment vertex. -/
theorem exists_attach_reach (hconn : D.G.Connected) (hU : D.ForbiddenV.Nonempty)
    {x : V} (hx : x ∈ D.FlexV) :
    ∃ z ∈ D.flexAttach, ReachIn D.G D.FlexV x z := by
  by_contra hcon
  push_neg at hcon
  set K := D.FlexV.filter (fun z => ReachIn D.G D.FlexV x z) with hK
  have hxK : x ∈ K := by
    rw [hK, Finset.mem_filter]
    exact ⟨hx, ReachIn.refl⟩
  have hclosed : ∀ z ∈ K, ∀ w, D.G.Adj z w → w ∈ K := by
    intro z hz w hadj
    rw [hK, Finset.mem_filter] at hz
    have hzF : z ∈ D.FlexV := hz.1
    have hwnotForced : w ∉ D.ForcedV := by
      intro hwF
      exact Finset.disjoint_left.1 Allowed_disjoint_Forbidden (mem_Flex.1 hzF).1
        (forbidden_of_adj_forced hwF hadj.symm)
    have hwnotForbidden : w ∉ D.ForbiddenV := by
      intro hwFb
      exact hcon z (Finset.mem_filter.2 ⟨hzF, ⟨w, hwFb, hadj.symm⟩⟩) hz.2
    have hwA : w ∈ D.AllowedV := by
      by_contra hwA
      exact hwnotForbidden (by rw [Forbidden_eq_compl, Finset.mem_compl]; exact hwA)
    have hwF : w ∈ D.FlexV := mem_Flex.2 ⟨hwA, hwnotForced⟩
    rw [hK, Finset.mem_filter]
    exact ⟨hwF, hz.2.trans (ReachIn.single ⟨hzF, hwF, hadj⟩)⟩
  obtain ⟨u, hu⟩ := hU
  obtain ⟨p⟩ := hconn.preconnected x u
  have huK := B1Zero.mem_of_walk_closed hclosed p hxK
  rw [hK, Finset.mem_filter] at huK
  exact Finset.disjoint_left.1 Allowed_disjoint_Forbidden (mem_Flex.1 huK.1).1 hu

/-- Reachability in the corona projects to reachability in the base. -/
theorem reachIn_base_of_reachIn_cor (hb1 : D.B1Zero) {a c : V}
    (ha : a ∈ (D.flexCorona hb1).base) (h : ReachIn D.G D.FlexV a c) :
    ReachIn D.G (D.flexCorona hb1).base a ((D.flexCorona hb1).pinv c) := by
  set C := D.flexCorona hb1 with hC
  have hcor : C.cor C.base = D.FlexV := cor_flexBaseSet hb1
  induction h with
  | refl => rw [CoronaData.pinv_base ha]
  | @tail c d _ hstep ih =>
      have hcC : c ∈ C.cor C.base := hcor ▸ hstep.1
      have hdC : d ∈ C.cor C.base := hcor ▸ hstep.2.1
      have hpc : C.pinv c ∈ C.base := CoronaData.pinv_mem (Finset.Subset.refl _) hcC
      have hpd : C.pinv d ∈ C.base := CoronaData.pinv_mem (Finset.Subset.refl _) hdC
      rcases CoronaData.eq_or_eq_pend (Finset.Subset.refl _) hcC with hc | hc
      · rcases CoronaData.eq_or_eq_pend (Finset.Subset.refl _) hdC with hd | hd
        · refine ih.trans (ReachIn.single ⟨hpc, hpd, ?_⟩)
          rw [← hc, ← hd]
          exact hstep.2.2
        · have hcd : c = C.pinv d :=
            C.pend_nbr (C.pinv d) hpd c hcC (by rw [← hd]; exact hstep.2.2.symm)
          rw [← hcd, hc]
          exact ih
      · have hdc : d = C.pinv c :=
          C.pend_nbr (C.pinv c) hpc d hdC (by rw [← hc]; exact hstep.2.2)
        rw [hdc, CoronaData.pinv_base hpc]
        exact ih

/-- A transversal of the base whose vertices are attachment vertices, up to pendants. -/
theorem exists_root_transversal (hb1 : D.B1Zero) (hconn : D.G.Connected)
    (hU : D.ForbiddenV.Nonempty) :
    ∃ Y : Finset V, IsTransversal D.G (D.flexCorona hb1).base Y ∧
      ∀ y ∈ Y, y ∈ D.flexAttach ∨ (D.flexCorona hb1).pend y ∈ D.flexAttach := by
  set C := D.flexCorona hb1 with hC
  have hcor : C.cor C.base = D.FlexV := cor_flexBaseSet hb1
  set R := C.base.filter (fun y => y ∈ D.flexAttach ∨ C.pend y ∈ D.flexAttach) with hR
  have hmeet : ∀ x ∈ C.base, ∃ y ∈ R, ReachIn D.G C.base x y := by
    intro x hx
    have hxF : x ∈ D.FlexV := hcor ▸ CoronaData.subset_cor hx
    obtain ⟨z, hz, hreach⟩ := exists_attach_reach hconn hU hxF
    have hzC : z ∈ C.cor C.base := by
      rw [hcor]
      exact (Finset.mem_filter.1 hz).1
    refine ⟨C.pinv z, ?_, reachIn_base_of_reachIn_cor hb1 hx hreach⟩
    rw [hR, Finset.mem_filter]
    refine ⟨CoronaData.pinv_mem (Finset.Subset.refl _) hzC, ?_⟩
    rcases CoronaData.eq_or_eq_pend (Finset.Subset.refl _) hzC with h | h
    · exact Or.inl (h ▸ hz)
    · exact Or.inr (h ▸ hz)
  obtain ⟨Y, hYR, hY⟩ := exists_transversal_subset (Finset.filter_subset _ _) hmeet
  refine ⟨Y, hY, fun y hy => ?_⟩
  have hmem := hYR hy
  rw [Finset.mem_filter] at hmem
  exact hmem.2

/-- **The flexible certificate.**  With ten flexible bags, the union-bound numerator built
from the codimension counts of the flexible part is at most a quarter of the extendable
count denominator. -/
theorem flex_certificate (hb1 : D.B1Zero) (hconn : D.G.Connected) (hU : D.ForbiddenV.Nonempty)
    (hm : D.flexBags.card = 10) :
    4 * (73 * cntRev D.G D.FlexV 10 0 + 24 * cntRev D.G D.FlexV 10 1
        + 3 * cntRev D.G D.FlexV 10 2
        + 15 * cntRev D.G (D.FlexV \ D.flexAttach) 10 0
        + 6 * cntRev D.G (D.FlexV \ D.flexAttach) 10 1
        + cntRev D.G (D.FlexV \ D.flexAttach) 10 2)
      ≤ 126 * cntRev D.G D.FlexV 10 0 + 84 * cntRev D.G D.FlexV 10 1
        + 36 * cntRev D.G D.FlexV 10 2 + 9 * cntRev D.G D.FlexV 10 3
        + cntRev D.G D.FlexV 10 4 := by
  classical
  set C := D.flexCorona hb1 with hC
  have hcor : C.cor C.base = D.FlexV := cor_flexBaseSet hb1
  have hcard : C.base.card = 10 := by
    rw [hC, flexCorona_base, card_flexBaseSet hb1, hm]
  obtain ⟨Y, hY, hYmove⟩ := exists_root_transversal hb1 hconn hU
  obtain ⟨f, hsize, hprof⟩ :=
    exists_forest_repr D.acyclic C C.base (Finset.Subset.refl _) Y hY
  have hP : ∀ j, C.vecAt C.base ∅ j = cntRev D.G D.FlexV 10 j := by
    intro j
    rw [CoronaData.vecAt_empty_avoid, hcor, hcard]
  have hQ : ∀ j, cntRev D.G (D.FlexV \ D.flexAttach) 10 j ≤ C.vecAt C.base Y j := by
    intro j
    rw [CoronaData.vecAt, hcard, ← hcor]
    exact CoronaData.cntRev_move_roots (Finset.Subset.refl _) hY.1 hYmove 10 j
  have hcert := rooted_forest_simple_bound f (by omega)
  simp only [simpleNumerator, denominator, hprof, CoronaData.vec] at hcert
  simp only [hP] at hcert
  have h0 := hQ 0
  have h1 := hQ 1
  have h2 := hQ 2
  omega

end TreeMatching

end MatchingBag
