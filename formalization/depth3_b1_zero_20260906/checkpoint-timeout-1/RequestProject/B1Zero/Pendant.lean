import RequestProject.B1Zero.BagCount

/-!
# The pendant lemma for `b₁ = 0`

If every independent set of size `α - 1` extends to a maximum independent set, then in every
bag consisting of two flexible vertices one of the two vertices has no allowed neighbour
outside its own bag ("it is a leaf of the allowed subgraph").

The proof is the matching-bag splicing argument: if both endpoints `x, y` of a flexible bag
had allowed neighbours `x', y'` outside the bag, then acyclicity separates `x'` from `y'` in
`G - {x, y}`, and gluing maximum independent sets across that separation produces an
independent `(α-1)`-set blocking both `x` and `y`.
-/

open Finset SimpleGraph

namespace MatchingBag

attribute [local instance 1] Classical.propDecidable

variable {V : Type*} [Fintype V] [DecidableEq V]

namespace TreeMatching

variable (D : TreeMatching V)

/-- `b₁ = 0`: every independent set of size `α - 1` extends to a maximum independent set. -/
def B1Zero : Prop :=
  ∀ S : Finset V, S.card = Fintype.card D.Bag - 1 →
    (∀ u ∈ S, ∀ v ∈ S, ¬ D.G.Adj u v) → ∃ I ∈ D.maxIndepSets, S ⊆ I

/-- The set of vertices reachable from `s` by a walk avoiding both `x` and `y`. -/
noncomputable def avoidReach (x y s : V) : Finset V :=
  Finset.univ.filter (fun v => ∃ w : D.G.Walk s v, x ∉ w.support ∧ y ∉ w.support)

variable {D}

lemma mem_avoidReach {x y s v : V} :
    v ∈ D.avoidReach x y s ↔ ∃ w : D.G.Walk s v, x ∉ w.support ∧ y ∉ w.support := by
  simp [avoidReach]

lemma self_mem_avoidReach {x y s : V} (hx : x ≠ s) (hy : y ≠ s) : s ∈ D.avoidReach x y s :=
  mem_avoidReach.2 ⟨SimpleGraph.Walk.nil, by simpa using hx, by simpa using hy⟩

lemma fst_notMem_avoidReach {x y s : V} : x ∉ D.avoidReach x y s := by
  intro h
  obtain ⟨w, hx, -⟩ := mem_avoidReach.1 h
  exact hx w.end_mem_support

lemma snd_notMem_avoidReach {x y s : V} : y ∉ D.avoidReach x y s := by
  intro h
  obtain ⟨w, -, hy⟩ := mem_avoidReach.1 h
  exact hy w.end_mem_support

lemma avoidReach_closure {x y s u v : V} (hu : u ∈ D.avoidReach x y s)
    (huv : D.G.Adj u v) (hvx : v ≠ x) (hvy : v ≠ y) : v ∈ D.avoidReach x y s := by
  obtain ⟨w, hx, hy⟩ := mem_avoidReach.1 hu
  refine mem_avoidReach.2 ⟨w.append (SimpleGraph.Walk.cons huv SimpleGraph.Walk.nil), ?_, ?_⟩
  · rw [SimpleGraph.Walk.support_append]
    simp only [List.mem_append, SimpleGraph.Walk.support_cons, SimpleGraph.Walk.support_nil,
      List.tail_cons, List.mem_singleton]
    rintro (h | h)
    · exact hx h
    · exact hvx h.symm
  · rw [SimpleGraph.Walk.support_append]
    simp only [List.mem_append, SimpleGraph.Walk.support_cons, SimpleGraph.Walk.support_nil,
      List.tail_cons, List.mem_singleton]
    rintro (h | h)
    · exact hy h
    · exact hvy h.symm

/-- **Separation from acyclicity.**  If `xy` is an edge, `x'` a neighbour of `x` other than
`y`, and `y'` a neighbour of `y` other than `x`, then `y'` is not reachable from `x'` while
avoiding `x` and `y`. -/
lemma notMem_avoidReach_of_acyclic {x y x' y' : V} (hxy : D.G.Adj x y)
    (hxx' : D.G.Adj x x') (hx'y : x' ≠ y) (hyy' : D.G.Adj y y') (hy'x : y' ≠ x) :
    y' ∉ D.avoidReach x y x' := by
  intro hmem
  obtain ⟨p, hpx, hpy⟩ := mem_avoidReach.1 hmem
  have hbridge : D.G.IsBridge s(x, y) :=
    (SimpleGraph.isAcyclic_iff_forall_adj_isBridge.1 D.acyclic) hxy
  obtain ⟨-, hall⟩ := SimpleGraph.isBridge_iff_adj_and_forall_walk_mem_edges.1 hbridge
  -- the walk  x → x' → … → y' → y  avoids the edge `xy`
  set q : D.G.Walk x y :=
    (SimpleGraph.Walk.cons hxx' p).append (SimpleGraph.Walk.cons hyy'.symm SimpleGraph.Walk.nil)
    with hq
  have hmem2 := hall q
  rw [hq, SimpleGraph.Walk.edges_append, SimpleGraph.Walk.edges_cons,
    SimpleGraph.Walk.edges_cons] at hmem2
  simp only [SimpleGraph.Walk.edges_nil, List.mem_append, List.mem_cons,
    List.not_mem_nil, or_false] at hmem2
  rcases hmem2 with (he | he) | he
  · rw [Sym2.eq_iff] at he
    rcases he with ⟨-, h2⟩ | ⟨-, h2⟩
    · exact hx'y h2.symm
    · exact hxy.ne' h2
  · exact hpx (p.fst_mem_support_of_mem_edges he)
  · rw [Sym2.eq_iff] at he
    rcases he with ⟨h1, -⟩ | ⟨h1, -⟩
    · exact hy'x h1.symm
    · exact hxy.ne h1

/-! ### Bag-closed sets -/

/-- If `Z` is a union of bags, every maximum independent set meets it in the same number of
vertices, namely the number of bags contained in `Z`. -/
lemma card_inter_of_bagClosed {Z : Finset V}
    (hZ : ∀ u v : V, u ∈ Z → D.bagOf u = D.bagOf v → v ∈ Z)
    {I : Finset V} (hI : I ∈ D.maxIndepSets) :
    (I ∩ Z).card = (Finset.univ.filter (fun b : D.Bag => D.bagSet b ⊆ Z)).card := by
  classical
  refine Finset.card_bij (fun u _ => D.bagOf u) ?_ ?_ ?_
  · intro u hu
    rw [Finset.mem_inter] at hu
    rw [Finset.mem_filter]
    refine ⟨Finset.mem_univ _, fun v hv => ?_⟩
    exact hZ u v hu.2 (bagOf_eq_of_mem hv).symm
  · intro u hu v hv huv
    rw [Finset.mem_inter] at hu hv
    exact bagOf_injOn_indep (maxIndepSets_indep hI) (by exact_mod_cast hu.1)
      (by exact_mod_cast hv.1) huv
  · intro b hb
    rw [Finset.mem_filter] at hb
    refine ⟨D.pick b (D.coverWord I b), ?_, ?_⟩
    · rw [Finset.mem_inter]
      exact ⟨pick_coverWord_mem hI b, hb.2 (pick_mem_bagSet b _)⟩
    · exact bagOf_pick b _


/-! ### The pendant lemma -/

/-- **Pendant lemma.**  Assuming `b₁ = 0`, in every bag of two flexible vertices one of the
two vertices has no allowed neighbour outside its bag. -/
theorem exists_leaf_in_flex_bag (hb1 : D.B1Zero) {b : D.Bag}
    (hb : ∀ v ∈ D.bagSet b, v ∈ D.FlexV) :
    ∃ x ∈ D.bagSet b, ∀ w, D.G.Adj x w → w ∈ D.AllowedV → w ∈ D.bagSet b := by
  classical
  by_contra hcon
  push_neg at hcon
  have hnf : ∀ v ∈ D.bagSet b, v ∉ D.ForcedV := fun v hv => (mem_Flex.1 (hb v hv)).2
  obtain ⟨hcard2, -⟩ := bagSet_flex_of_no_forced hnf
  obtain ⟨x, y, hxy, hbxy⟩ := Finset.card_eq_two.1 hcard2
  have hxb : x ∈ D.bagSet b := by rw [hbxy]; simp
  have hyb : y ∈ D.bagSet b := by rw [hbxy]; simp
  have hadjxy : D.G.Adj x y := bagSet_adj hxb hyb hxy
  obtain ⟨x', hxx', hx'A, hx'b⟩ := hcon x hxb
  obtain ⟨y', hyy', hy'A, hy'b⟩ := hcon y hyb
  have hx'x : x' ≠ x := by rintro rfl; exact hx'b hxb
  have hx'y : x' ≠ y := by rintro rfl; exact hx'b hyb
  have hy'x : y' ≠ x := by rintro rfl; exact hy'b hxb
  have hy'y : y' ≠ y := by rintro rfl; exact hy'b hyb
  set Z := D.avoidReach x y x' with hZdef
  have hxZ : x ∉ Z := fst_notMem_avoidReach
  have hyZ : y ∉ Z := snd_notMem_avoidReach
  have hx'Z : x' ∈ Z := self_mem_avoidReach (Ne.symm hx'x) (Ne.symm hx'y)
  have hy'Z : y' ∉ Z := notMem_avoidReach_of_acyclic hadjxy hxx' hx'y hyy' hy'x
  -- `Z` is a union of bags
  have hZclosed : ∀ u v : V, u ∈ Z → D.bagOf u = D.bagOf v → v ∈ Z := by
    intro u v hu huv
    by_cases hvu : v = u
    · rwa [hvu]
    have hvb : v ∈ D.bagSet (D.bagOf u) := by rw [huv]; exact mem_bagSet_bagOf v
    have hub : u ∈ D.bagSet (D.bagOf u) := mem_bagSet_bagOf u
    have hadj : D.G.Adj u v := bagSet_adj hub hvb (Ne.symm hvu)
    have hvx : v ≠ x := by
      rintro rfl
      have : D.bagOf u = b := by rw [huv]; exact bagOf_eq_of_mem hxb
      have : u ∈ D.bagSet b := by rw [← this]; exact mem_bagSet_bagOf u
      rw [hbxy] at this
      rcases Finset.mem_insert.1 this with rfl | h
      · exact hxZ hu
      · rw [Finset.mem_singleton] at h; subst h; exact hyZ hu
    have hvy : v ≠ y := by
      rintro rfl
      have : D.bagOf u = b := by rw [huv]; exact bagOf_eq_of_mem hyb
      have : u ∈ D.bagSet b := by rw [← this]; exact mem_bagSet_bagOf u
      rw [hbxy] at this
      rcases Finset.mem_insert.1 this with rfl | h
      · exact hxZ hu
      · rw [Finset.mem_singleton] at h; subst h; exact hyZ hu
    exact avoidReach_closure hu hadj hvx hvy
  obtain ⟨Ix, hIx, hx'Ix⟩ := mem_Allowed.1 hx'A
  obtain ⟨Iy, hIy, hy'Iy⟩ := mem_Allowed.1 hy'A
  set ZB := (Finset.univ.filter (fun b' : D.Bag => D.bagSet b' ⊆ Z)).card with hZB
  have hIxZ : (Ix ∩ Z).card = ZB := card_inter_of_bagClosed hZclosed hIx
  have hIyZ : (Iy ∩ Z).card = ZB := card_inter_of_bagClosed hZclosed hIy
  set a := Fintype.card D.Bag with ha
  have hIycard : Iy.card = a := card_of_maxIndepSets hIy
  -- the spliced set
  set J := (Iy \ Z) \ D.bagSet b with hJ
  have hsplit1 : (Iy ∩ Z).card + (Iy \ Z).card = Iy.card := Finset.card_inter_add_card_sdiff _ _
  have hinter : (Iy \ Z) ∩ D.bagSet b = Iy ∩ D.bagSet b := by
    ext v
    simp only [Finset.mem_inter, Finset.mem_sdiff]
    constructor
    · rintro ⟨⟨h1, -⟩, h2⟩; exact ⟨h1, h2⟩
    · rintro ⟨h1, h2⟩
      refine ⟨⟨h1, ?_⟩, h2⟩
      rw [hbxy] at h2
      rcases Finset.mem_insert.1 h2 with rfl | h
      · exact hxZ
      · rw [Finset.mem_singleton] at h; subst h; exact hyZ
  have hsplit2 : ((Iy \ Z) ∩ D.bagSet b).card + J.card = (Iy \ Z).card :=
    Finset.card_inter_add_card_sdiff _ _
  have hone : (Iy ∩ D.bagSet b).card = 1 := card_inter_bagSet hIy b
  have hJcard : J.card + ZB + 1 = a := by
    rw [hinter, hone] at hsplit2
    omega
  set W := (Ix ∩ Z) ∪ J with hW
  have hWdisj : Disjoint (Ix ∩ Z) J := by
    rw [Finset.disjoint_left]
    intro v hv hv'
    exact (Finset.mem_sdiff.1 (Finset.mem_sdiff.1 hv').1).2 (Finset.mem_inter.1 hv).2
  have hWcard : W.card = a - 1 := by
    rw [hW, Finset.card_union_of_disjoint hWdisj, hIxZ]
    omega
  have hWind : ∀ u ∈ W, ∀ v ∈ W, ¬ D.G.Adj u v := by
    have key : ∀ u ∈ Ix ∩ Z, ∀ v ∈ J, ¬ D.G.Adj u v := by
      intro u hu v hv hadj
      rw [Finset.mem_inter] at hu
      have hvZ : v ∉ Z := (Finset.mem_sdiff.1 (Finset.mem_sdiff.1 hv).1).2
      have hvb : v ∉ D.bagSet b := (Finset.mem_sdiff.1 hv).2
      have hvx : v ≠ x := by rintro rfl; exact hvb hxb
      have hvy : v ≠ y := by rintro rfl; exact hvb hyb
      exact hvZ (avoidReach_closure hu.2 hadj hvx hvy)
    intro u hu v hv
    rw [hW, Finset.mem_union] at hu hv
    rcases hu with hu | hu <;> rcases hv with hv | hv
    · exact maxIndepSets_indep hIx u (Finset.mem_inter.1 hu).1 v (Finset.mem_inter.1 hv).1
    · exact key u hu v hv
    · exact fun hadj => key v hv u hu hadj.symm
    · exact maxIndepSets_indep hIy u (Finset.mem_sdiff.1 (Finset.mem_sdiff.1 hu).1).1 v
        (Finset.mem_sdiff.1 (Finset.mem_sdiff.1 hv).1).1
  obtain ⟨I, hI, hWI⟩ := hb1 W hWcard hWind
  -- `x'` and `y'` both lie in `W`
  have hx'W : x' ∈ W := by
    rw [hW, Finset.mem_union]
    exact Or.inl (Finset.mem_inter.2 ⟨hx'Ix, hx'Z⟩)
  have hy'W : y' ∈ W := by
    rw [hW, Finset.mem_union]
    refine Or.inr (Finset.mem_sdiff.2 ⟨Finset.mem_sdiff.2 ⟨hy'Iy, hy'Z⟩, ?_⟩)
    rw [hbxy]
    simp [hy'x, hy'y]
  -- but `I` must contain `x` or `y`
  have hIb : (I ∩ D.bagSet b).card = 1 := card_inter_bagSet hI b
  obtain ⟨z, hz⟩ := Finset.card_eq_one.1 hIb
  have hzI : z ∈ I := (Finset.mem_inter.1 (hz ▸ Finset.mem_singleton_self z)).1
  have hzb : z ∈ D.bagSet b := (Finset.mem_inter.1 (hz ▸ Finset.mem_singleton_self z)).2
  rw [hbxy] at hzb
  rcases Finset.mem_insert.1 hzb with rfl | hz2
  · exact maxIndepSets_indep hI z hzI x' (hWI hx'W) hxx'
  · rw [Finset.mem_singleton] at hz2
    subst hz2
    exact maxIndepSets_indep hI z hzI y' (hWI hy'W) hyy'

end TreeMatching

end MatchingBag
