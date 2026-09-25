import RequestProject.LowDensity

/-!
# AllowedV, forced, forbidden and flexible vertices of a forest with a maximum matching

For a finite forest `G` with a maximum matching, organised into bags by
`MatchingBag.TreeMatching`, we split the vertices into

* `AllowedV`  — vertices lying in at least one maximum independent set;
* `ForcedV`   — vertices lying in **every** maximum independent set;
* `ForbiddenV`— vertices lying in **no** maximum independent set (the complement of `AllowedV`);
* `FlexV`     — allowed but not forced.

The main results of this file are the *bag classification*: every bag is either a single
forced vertex, or a forbidden vertex together with its forced mate, or a pair of flexible
vertices; and the resulting count
`|V| + |ForcedV| = 2 α + |ForbiddenV|`, i.e. `|ForcedV| = |ForbiddenV| + (2α - |V|)`.

None of this uses `b₁ = 0`.
-/

open Finset SimpleGraph

namespace MatchingBag

attribute [local instance 1] Classical.propDecidable

variable {V : Type*} [Fintype V] [DecidableEq V]

namespace TreeMatching

variable (D : TreeMatching V)

/-- Vertices lying in at least one maximum independent set. -/
noncomputable def AllowedV : Finset V :=
  Finset.univ.filter (fun v => ∃ I ∈ D.maxIndepSets, v ∈ I)

/-- Vertices lying in every maximum independent set. -/
noncomputable def ForcedV : Finset V :=
  Finset.univ.filter (fun v => ∀ I ∈ D.maxIndepSets, v ∈ I)

/-- Vertices lying in no maximum independent set. -/
noncomputable def ForbiddenV : Finset V :=
  Finset.univ.filter (fun v => ∀ I ∈ D.maxIndepSets, v ∉ I)

/-- AllowedV but not forced vertices. -/
noncomputable def FlexV : Finset V := D.AllowedV \ D.ForcedV

variable {D}

@[simp] lemma mem_Allowed {v : V} : v ∈ D.AllowedV ↔ ∃ I ∈ D.maxIndepSets, v ∈ I := by
  simp [AllowedV]

@[simp] lemma mem_Forced {v : V} : v ∈ D.ForcedV ↔ ∀ I ∈ D.maxIndepSets, v ∈ I := by
  simp [ForcedV]

@[simp] lemma mem_Forbidden {v : V} : v ∈ D.ForbiddenV ↔ ∀ I ∈ D.maxIndepSets, v ∉ I := by
  simp [ForbiddenV]

@[simp] lemma mem_Flex {v : V} : v ∈ D.FlexV ↔ v ∈ D.AllowedV ∧ v ∉ D.ForcedV := by
  simp [FlexV]

lemma Forbidden_eq_compl : D.ForbiddenV = D.AllowedVᶜ := by
  ext v; simp [not_exists]

lemma Forced_subset_Allowed : D.ForcedV ⊆ D.AllowedV := by
  intro v hv
  obtain ⟨I, hI⟩ := maxIndepSets_nonempty D
  exact mem_Allowed.2 ⟨I, hI, mem_Forced.1 hv I hI⟩

lemma Allowed_disjoint_Forbidden : Disjoint D.AllowedV D.ForbiddenV := by
  rw [Forbidden_eq_compl]; exact disjoint_compl_right

/-! ### Bags -/

/-- Two distinct vertices of the same bag are adjacent. -/
lemma bagSet_adj {b : D.Bag} {u v : V} (hu : u ∈ D.bagSet b) (hv : v ∈ D.bagSet b)
    (huv : u ≠ v) : D.G.Adj u v := by
  cases b with
  | inl i =>
      rw [mem_bagSet_inl] at hu hv
      rcases hu with rfl | rfl <;> rcases hv with rfl | rfl
      · exact absurd rfl huv
      · exact adj_Lv_Rv i
      · exact (adj_Lv_Rv i).symm
      · exact absurd rfl huv
  | inr w =>
      rw [mem_bagSet_inr] at hu hv
      exact absurd (hu.trans hv.symm) huv

lemma bagSet_nonempty (b : D.Bag) : (D.bagSet b).Nonempty := by
  cases b with
  | inl i => exact ⟨D.Lv i, by simp⟩
  | inr w => exact ⟨(w : V), by simp⟩

lemma card_bagSet_pos (b : D.Bag) : 0 < (D.bagSet b).card :=
  Finset.card_pos.2 (bagSet_nonempty b)

/-- Every maximum independent set meets every bag in exactly one vertex. -/
lemma card_inter_bagSet {I : Finset V} (hI : I ∈ D.maxIndepSets) (b : D.Bag) :
    (I ∩ D.bagSet b).card = 1 := by
  refine Finset.card_eq_one.2 ⟨D.pick b (D.coverWord I b), ?_⟩
  apply Finset.eq_singleton_iff_unique_mem.2
  refine ⟨Finset.mem_inter.2 ⟨pick_coverWord_mem hI b, pick_mem_bagSet b _⟩, ?_⟩
  intro u hu
  rw [Finset.mem_inter] at hu
  have h1 : D.bagOf u = b := bagOf_eq_of_mem hu.2
  have := eq_pick_coverWord hI hu.1
  rw [h1] at this
  exact this

/-- Every bag contains an allowed vertex. -/
lemma exists_allowed_mem_bagSet (b : D.Bag) : ∃ v ∈ D.bagSet b, v ∈ D.AllowedV := by
  obtain ⟨I, hI⟩ := maxIndepSets_nonempty D
  exact ⟨D.pick b (D.coverWord I b), pick_mem_bagSet b _,
    mem_Allowed.2 ⟨I, hI, pick_coverWord_mem hI b⟩⟩

/-- A forced vertex has no allowed neighbour. -/
lemma forbidden_of_adj_forced {u v : V} (hu : u ∈ D.ForcedV) (huv : D.G.Adj u v) :
    v ∈ D.ForbiddenV := by
  refine mem_Forbidden.2 fun I hI hvI => ?_
  exact maxIndepSets_indep hI u (mem_Forced.1 hu I hI) v hvI huv

/-- Two distinct forced vertices never share a bag. -/
lemma bagOf_injOn_Forced : Set.InjOn D.bagOf (D.ForcedV : Set V) := by
  obtain ⟨I, hI⟩ := maxIndepSets_nonempty D
  refine Set.InjOn.mono ?_ (bagOf_injOn_indep (maxIndepSets_indep hI))
  intro v hv
  exact_mod_cast mem_Forced.1 (by exact_mod_cast hv) I hI

/-- A bag containing a forced vertex contains no other allowed vertex. -/
lemma bagSet_eq_of_forced {b : D.Bag} {f : V} (hf : f ∈ D.bagSet b) (hfF : f ∈ D.ForcedV)
    {v : V} (hv : v ∈ D.bagSet b) (hvne : v ≠ f) : v ∈ D.ForbiddenV :=
  forbidden_of_adj_forced hfF (bagSet_adj hf hv (Ne.symm hvne))

/-- If a bag contains no forced vertex, then it consists of two flexible vertices. -/
lemma bagSet_flex_of_no_forced {b : D.Bag} (hb : ∀ v ∈ D.bagSet b, v ∉ D.ForcedV) :
    (D.bagSet b).card = 2 ∧ ∀ v ∈ D.bagSet b, v ∈ D.FlexV := by
  obtain ⟨v, hv, hvA⟩ := exists_allowed_mem_bagSet (D := D) b
  have hvF : v ∉ D.ForcedV := hb v hv
  obtain ⟨I, hI, hvI⟩ : ∃ I ∈ D.maxIndepSets, v ∉ I := by
    by_contra hc
    push_neg at hc
    exact hvF (mem_Forced.2 hc)
  have h1 := card_inter_bagSet hI b
  obtain ⟨u, hu⟩ := Finset.card_eq_one.1 h1
  have huI : u ∈ I := (Finset.mem_inter.1 (hu ▸ Finset.mem_singleton_self u)).1
  have hub : u ∈ D.bagSet b := (Finset.mem_inter.1 (hu ▸ Finset.mem_singleton_self u)).2
  have huv : u ≠ v := by rintro rfl; exact hvI huI
  have hsub : ({u, v} : Finset V) ⊆ D.bagSet b := by
    intro w hw
    rcases Finset.mem_insert.1 hw with rfl | hw
    · exact hub
    · rw [Finset.mem_singleton] at hw; subst hw; exact hv
  have h2 : ({u, v} : Finset V).card = 2 := by
    rw [Finset.card_insert_of_notMem (by simp [huv]), Finset.card_singleton]
  have hle := card_bagSet_le (D := D) b
  have hcard2 : (D.bagSet b).card = 2 := by
    have := Finset.card_le_card hsub
    omega
  have heq : ({u, v} : Finset V) = D.bagSet b := Finset.eq_of_subset_of_card_le hsub (by omega)
  refine ⟨hcard2, fun w hw => ?_⟩
  refine mem_Flex.2 ⟨?_, hb w hw⟩
  rw [← heq] at hw
  rcases Finset.mem_insert.1 hw with rfl | hw
  · exact mem_Allowed.2 ⟨I, hI, huI⟩
  · rw [Finset.mem_singleton] at hw; subst hw; exact hvA

/-- A flexible vertex lies in a bag with two flexible vertices. -/
lemma bag_of_flex {v : V} (hv : v ∈ D.FlexV) :
    (D.bagSet (D.bagOf v)).card = 2 ∧ ∀ w ∈ D.bagSet (D.bagOf v), w ∈ D.FlexV := by
  apply bagSet_flex_of_no_forced
  intro w hw hwF
  by_cases hwv : w = v
  · subst hwv; exact (mem_Flex.1 hv).2 hwF
  · exact Finset.disjoint_left.1 Allowed_disjoint_Forbidden (mem_Flex.1 hv).1
      (bagSet_eq_of_forced hw hwF (mem_bagSet_bagOf v) (Ne.symm hwv))

/-- A forbidden vertex lies in a bag whose other vertex is forced. -/
lemma bag_of_forbidden {v : V} (hv : v ∈ D.ForbiddenV) :
    ∃ f ∈ D.bagSet (D.bagOf v), f ∈ D.ForcedV ∧ f ≠ v := by
  obtain ⟨w, hw, hwA⟩ := exists_allowed_mem_bagSet (D := D) (D.bagOf v)
  have hwv : w ≠ v := by
    rintro rfl
    exact Finset.disjoint_left.1 Allowed_disjoint_Forbidden hwA hv
  refine ⟨w, hw, ?_, hwv⟩
  by_contra hwF
  have hflex : w ∈ D.FlexV := mem_Flex.2 ⟨hwA, hwF⟩
  have hb := (bag_of_flex hflex).2
  have hvb : v ∈ D.bagSet (D.bagOf w) := by
    rw [bagOf_eq_of_mem hw]; exact mem_bagSet_bagOf v
  exact Finset.disjoint_left.1 Allowed_disjoint_Forbidden (mem_Flex.1 (hb v hvb)).1 hv

end TreeMatching

end MatchingBag
