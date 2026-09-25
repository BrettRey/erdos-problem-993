import RequestProject.B1Zero.WellCovered
import RequestProject.B1Zero.ForestCount

/-!
# Every forbidden vertex has at least three forced neighbours

Assuming `b₁ = 0`, a forbidden vertex `v` has `d_C(v) ≥ 3` forced neighbours.  The proof
constructs an independent set of size `α - d_C(v)` avoiding `N(v)`; adjoining `v` gives an
independent set of size `α + 1 - d_C(v)`, which contradicts either maximality (if
`d_C(v) = 1`) or `b₁ = 0` (if `d_C(v) = 2`).

Together with the forest incidence bound this yields `|N_C(S)| ≥ 2|S| + 1` for every
nonempty set `S` of forbidden vertices.
-/

open Finset SimpleGraph

namespace MatchingBag

attribute [local instance 1] Classical.propDecidable

variable {V : Type*} [Fintype V] [DecidableEq V]

namespace TreeMatching

variable {D : TreeMatching V}

/-- In a forest, every walk between adjacent vertices uses the connecting edge. -/
lemma edge_mem_walk_edges {p q : V} (hpq : D.G.Adj p q) (w : D.G.Walk p q) :
    s(p, q) ∈ w.edges :=
  (SimpleGraph.isBridge_iff_adj_and_forall_walk_mem_edges.1
    ((SimpleGraph.isAcyclic_iff_forall_adj_isBridge.1 D.acyclic) hpq)).2 w

/-- A forest has no triangle. -/
lemma no_triangle {x y z : V} (hxy : D.G.Adj x y) (hyz : D.G.Adj y z) (hxz : D.G.Adj x z) :
    False := by
  have h := edge_mem_walk_edges hxz
    (SimpleGraph.Walk.cons hxy (SimpleGraph.Walk.cons hyz SimpleGraph.Walk.nil))
  simp only [SimpleGraph.Walk.edges_cons, SimpleGraph.Walk.edges_nil, List.mem_cons,
    List.not_mem_nil, or_false] at h
  rcases h with he | he
  · rw [Sym2.eq_iff] at he
    rcases he with ⟨-, h2⟩ | ⟨h1, -⟩
    · exact hyz.ne' h2
    · exact hxy.ne h1
  · rw [Sym2.eq_iff] at he
    rcases he with ⟨h1, -⟩ | ⟨-, h2⟩
    · exact hxy.ne h1
    · exact hyz.ne' h2

variable (D)

/-- The forced neighbours of a vertex. -/
noncomputable def forcedNbrs (v : V) : Finset V := D.ForcedV.filter (fun f => D.G.Adj v f)

/-- The forced neighbours of a set of vertices. -/
noncomputable def forcedNbrsSet (S : Finset V) : Finset V :=
  D.ForcedV.filter (fun f => ∃ u ∈ S, D.G.Adj u f)

variable {D}

@[simp] lemma mem_forcedNbrs {v f : V} :
    f ∈ D.forcedNbrs v ↔ f ∈ D.ForcedV ∧ D.G.Adj v f := by simp [forcedNbrs]

@[simp] lemma mem_forcedNbrsSet {S : Finset V} {f : V} :
    f ∈ D.forcedNbrsSet S ↔ f ∈ D.ForcedV ∧ ∃ u ∈ S, D.G.Adj u f := by simp [forcedNbrsSet]

lemma forcedNbrs_subset_forcedNbrsSet {S : Finset V} {u : V} (hu : u ∈ S) :
    D.forcedNbrs u ⊆ D.forcedNbrsSet S := by
  intro f hf
  rw [mem_forcedNbrs] at hf
  exact mem_forcedNbrsSet.2 ⟨hf.1, u, hu, hf.2⟩

/-- The mate of a forbidden vertex is a forced neighbour. -/
lemma forcedNbrs_nonempty_of_forbidden {v : V} (hv : v ∈ D.ForbiddenV) :
    (D.forcedNbrs v).Nonempty := by
  obtain ⟨f, hfb, hfF, hfv⟩ := bag_of_forbidden hv
  exact ⟨f, mem_forcedNbrs.2 ⟨hfF, bagSet_adj (mem_bagSet_bagOf v) hfb (Ne.symm hfv)⟩⟩

/-- **Three forced neighbours.**  Assuming `b₁ = 0`, every forbidden vertex has at least
three forced neighbours. -/
theorem three_le_card_forcedNbrs (hb1 : D.B1Zero) (ha : 3 ≤ Fintype.card D.Bag)
    {v : V} (hv : v ∈ D.ForbiddenV) : 3 ≤ (D.forcedNbrs v).card := by
  classical
  by_contra hlt
  push_neg at hlt
  set a := Fintype.card D.Bag with hadef
  set Nc := D.forcedNbrs v with hNc
  set NB := Nc.image D.bagOf with hNB
  have hNBcard : NB.card = Nc.card := by
    refine Finset.card_image_of_injOn ?_
    intro f hf g hg h
    refine bagOf_injOn_Forced ?_ ?_ h
    · exact_mod_cast (mem_forcedNbrs.1 hf).1
    · exact_mod_cast (mem_forcedNbrs.1 hg).1
  -- select a representative in each bag that avoids `N(v)`
  have hchoice : ∀ b : D.Bag, ∃ u, b ∉ NB →
      (u ∈ D.bagSet b ∧ u ∈ D.AllowedV ∧ ¬ D.G.Adj v u ∧
        ((∀ w, D.G.Adj u w → w ∈ D.AllowedV → w ∈ D.bagSet b) ∨
         (∃ x ∈ D.bagSet b, x ≠ u ∧ D.G.Adj v x ∧ D.G.Adj x u))) := by
    intro b
    by_cases hbNB : b ∈ NB
    · exact ⟨v, fun h => absurd hbNB h⟩
    by_cases hF : ∃ f ∈ D.bagSet b, f ∈ D.ForcedV
    · obtain ⟨f, hfb, hfF⟩ := hF
      refine ⟨f, fun _ => ⟨hfb, Forced_subset_Allowed hfF, ?_, Or.inl ?_⟩⟩
      · intro hadj
        exact hbNB (Finset.mem_image.2 ⟨f, mem_forcedNbrs.2 ⟨hfF, hadj⟩, bagOf_eq_of_mem hfb⟩)
      · intro w hadj hwA
        exact absurd (forbidden_of_adj_forced hfF hadj)
          (Finset.disjoint_left.1 Allowed_disjoint_Forbidden hwA)
    · push_neg at hF
      obtain ⟨hcard2, hflex⟩ := bagSet_flex_of_no_forced hF
      obtain ⟨x, hxb, hxleaf⟩ := exists_leaf_in_flex_bag hb1 hflex
      by_cases hvx : D.G.Adj v x
      · -- take the mate of the leaf
        obtain ⟨p, q, hpq, hbpq⟩ := Finset.card_eq_two.1 hcard2
        have hpb : p ∈ D.bagSet b := by rw [hbpq]; simp
        have hqb : q ∈ D.bagSet b := by rw [hbpq]; simp
        have hxpq : x = p ∨ x = q := by
          have := hxb; rw [hbpq] at this; simpa using this
        refine ⟨if x = p then q else p, fun _ => ?_⟩
        by_cases hxp : x = p
        · subst hxp
          rw [if_pos (rfl : x = x)]
          have hadjxq : D.G.Adj x q := bagSet_adj hpb hqb hpq
          refine ⟨hqb, (mem_Flex.1 (hflex q hqb)).1, ?_, Or.inr ⟨x, hxb, hpq, hvx, hadjxq⟩⟩
          intro hadj
          exact no_triangle hvx hadjxq hadj
        · simp only [if_neg hxp]
          have hxq : x = q := by tauto
          subst hxq
          have hadjxp : D.G.Adj x p := bagSet_adj hqb hpb (Ne.symm hpq)
          refine ⟨hpb, (mem_Flex.1 (hflex p hpb)).1, ?_, Or.inr ⟨x, hxb, Ne.symm hpq, hvx, hadjxp⟩⟩
          intro hadj
          exact no_triangle hvx hadjxp hadj
      · exact ⟨x, fun _ => ⟨hxb, (mem_Flex.1 (hflex x hxb)).1, hvx, Or.inl hxleaf⟩⟩
  choose sel hsel using hchoice
  set T := (Finset.univ \ NB).image sel with hT
  have hselb : ∀ b ∈ Finset.univ \ NB, sel b ∈ D.bagSet b := fun b hb =>
    (hsel b (Finset.mem_sdiff.1 hb).2).1
  have hinj : Set.InjOn sel ((Finset.univ \ NB : Finset D.Bag) : Set D.Bag) := by
    intro b1 hb1' b2 hb2' h
    rw [← bagOf_eq_of_mem (hselb b1 (by exact_mod_cast hb1')),
      ← bagOf_eq_of_mem (hselb b2 (by exact_mod_cast hb2')), h]
  have hTcard : T.card = a - Nc.card := by
    rw [hT, Finset.card_image_of_injOn hinj, Finset.card_sdiff, Finset.inter_univ,
      Finset.card_univ, hNBcard]
  have hvT : v ∉ T := by
    rw [hT]
    intro h
    obtain ⟨b, hb, hbv⟩ := Finset.mem_image.1 h
    have := (hsel b (Finset.mem_sdiff.1 hb).2).2.1
    rw [hbv] at this
    exact Finset.disjoint_left.1 Allowed_disjoint_Forbidden this hv
  -- `T ∪ {v}` is independent
  have hTind : ∀ u ∈ T, ∀ u' ∈ T, ¬ D.G.Adj u u' := by
    intro u hu u' hu' hadj
    obtain ⟨b1, hb1', rfl⟩ := Finset.mem_image.1 hu
    obtain ⟨b2, hb2', rfl⟩ := Finset.mem_image.1 hu'
    obtain ⟨h1b, h1A, -, h1d⟩ := hsel b1 (Finset.mem_sdiff.1 hb1').2
    obtain ⟨h2b, h2A, -, h2d⟩ := hsel b2 (Finset.mem_sdiff.1 hb2').2
    have hbb1 : D.bagOf (sel b1) = b1 := bagOf_eq_of_mem h1b
    have hbb2 : D.bagOf (sel b2) = b2 := bagOf_eq_of_mem h2b
    have hne : b1 ≠ b2 := by rintro rfl; exact hadj.ne rfl
    rcases h1d with h1 | ⟨x1, hx1b, hx1u, hvx1, hx1a⟩
    · exact hne ((bagOf_eq_of_mem (h1 _ hadj h2A)).symm.trans hbb2)
    rcases h2d with h2 | ⟨x2, hx2b, hx2u, hvx2, hx2a⟩
    · exact hne (hbb1.symm.trans (bagOf_eq_of_mem (h2 _ hadj.symm h1A)))
    -- the five-vertex configuration would close a cycle
    have hvu1 : v ≠ sel b1 := fun h =>
      Finset.disjoint_left.1 Allowed_disjoint_Forbidden (h ▸ h1A) hv
    have hvu2 : v ≠ sel b2 := fun h =>
      Finset.disjoint_left.1 Allowed_disjoint_Forbidden (h ▸ h2A) hv
    have hx1ne : x1 ≠ sel b2 := fun h => notMem_bagSet_of_ne hne hx1b (h ▸ h2b)
    have hx2ne : x2 ≠ sel b1 := fun h => notMem_bagSet_of_ne (Ne.symm hne) hx2b (h ▸ h1b)
    have hw := edge_mem_walk_edges hadj
      (SimpleGraph.Walk.cons hx1a.symm (SimpleGraph.Walk.cons hvx1.symm
        (SimpleGraph.Walk.cons hvx2 (SimpleGraph.Walk.cons hx2a SimpleGraph.Walk.nil))))
    simp only [SimpleGraph.Walk.edges_cons, SimpleGraph.Walk.edges_nil, List.mem_cons,
      List.not_mem_nil, or_false] at hw
    rcases hw with he | he | he | he
    · rw [Sym2.eq_iff] at he
      rcases he with ⟨-, h2⟩ | ⟨-, h2⟩
      · exact hx1ne h2.symm
      · exact hadj.ne h2.symm
    · rw [Sym2.eq_iff] at he
      rcases he with ⟨-, h2⟩ | ⟨h1, -⟩
      · exact hvu2 h2.symm
      · exact hvu1 h1.symm
    · rw [Sym2.eq_iff] at he
      rcases he with ⟨h1, -⟩ | ⟨h1, -⟩
      · exact hvu1 h1.symm
      · exact hx2ne h1.symm
    · rw [Sym2.eq_iff] at he
      rcases he with ⟨h1, -⟩ | ⟨h1, -⟩
      · exact hx2ne h1.symm
      · exact hadj.ne h1
  have hvadj : ∀ u ∈ T, ¬ D.G.Adj v u := by
    intro u hu
    obtain ⟨b, hb, rfl⟩ := Finset.mem_image.1 hu
    exact (hsel b (Finset.mem_sdiff.1 hb).2).2.2.1
  -- pick a subset of size `α - 1` containing `v`
  have hNcpos : 1 ≤ Nc.card := Finset.card_pos.2 (forcedNbrs_nonempty_of_forbidden hv)
  obtain ⟨T', hT'sub, hT'card⟩ :=
    Finset.exists_subset_card_eq (s := T) (n := a - 2) (by omega)
  set S := insert v T' with hS
  have hvT' : v ∉ T' := fun h => hvT (hT'sub h)
  have hScard : S.card = a - 1 := by
    rw [hS, Finset.card_insert_of_notMem hvT', hT'card]
    omega
  have hSind : ∀ u ∈ S, ∀ u' ∈ S, ¬ D.G.Adj u u' := by
    intro u hu u' hu'
    rw [hS, Finset.mem_insert] at hu hu'
    rcases hu with rfl | hu <;> rcases hu' with rfl | hu'
    · exact fun h => h.ne rfl
    · exact hvadj u' (hT'sub hu')
    · exact fun h => hvadj u (hT'sub hu) h.symm
    · exact hTind u (hT'sub hu) u' (hT'sub hu')
  obtain ⟨I, hI, hSI⟩ := hb1 S hScard hSind
  exact mem_Forbidden.1 hv I hI (hSI (Finset.mem_insert_self v T'))

/-- **Expansion.**  Assuming `b₁ = 0`, every nonempty set of forbidden vertices has at least
`2|S| + 1` forced neighbours. -/
theorem card_forcedNbrsSet_ge (hb1 : D.B1Zero) (ha : 3 ≤ Fintype.card D.Bag)
    {S : Finset V} (hS : S ⊆ D.ForbiddenV) (hne : S.Nonempty) :
    2 * S.card + 1 ≤ (D.forcedNbrsSet S).card := by
  classical
  set N := D.forcedNbrsSet S with hN
  have hdisj : Disjoint S N := by
    rw [Finset.disjoint_left]
    intro u hu huN
    exact Finset.disjoint_left.1 Allowed_disjoint_Forbidden
      (Forced_subset_Allowed (mem_forcedNbrsSet.1 huN).1) (hS hu)
  have hUne : (S ∪ N).Nonempty := hne.mono Finset.subset_union_left
  have hbound := D.acyclic.card_adjPairs_le S N hdisj hUne
  have hcount : ((S ×ˢ N).filter (fun p : V × V => D.G.Adj p.1 p.2)).card
      = ∑ u ∈ S, (N.filter (fun f => D.G.Adj u f)).card := by
    rw [Finset.card_filter, Finset.sum_product]
    exact Finset.sum_congr rfl fun u _ => (Finset.card_filter _ _).symm
  have hge : ∀ u ∈ S, 3 ≤ (N.filter (fun f => D.G.Adj u f)).card := by
    intro u hu
    have h1 : D.forcedNbrs u ⊆ N.filter (fun f => D.G.Adj u f) := by
      intro f hf
      rw [Finset.mem_filter]
      exact ⟨forcedNbrs_subset_forcedNbrsSet hu hf, (mem_forcedNbrs.1 hf).2⟩
    exact le_trans (three_le_card_forcedNbrs hb1 ha (hS hu)) (Finset.card_le_card h1)
  have hsum : 3 * S.card ≤ ∑ u ∈ S, (N.filter (fun f => D.G.Adj u f)).card := by
    calc 3 * S.card = ∑ _u ∈ S, 3 := by rw [Finset.sum_const, smul_eq_mul, mul_comm]
      _ ≤ _ := Finset.sum_le_sum hge
  omega

end TreeMatching

end MatchingBag
