import RequestProject.ClanAudit.TreeCriterion
import RequestProject.ClanAudit.TreeC2

/-!
# The explicit models are trees with at most two branch vertices

The recognition theorems of `SpiderRecognition.lean` and `ConnectorRecognition.lean` show
that every finite tree with at most two branch vertices is one of the models.  This file
proves the converse: `spider n len` and `connGraph t m a n b` really are trees, and their
only possible branch vertices are the hubs.  Together this characterizes the `C₂` family:

`IsC2Model G ↔ G.IsTree ∧ (branchVerts G).card ≤ 2`.
-/

namespace ClanAudit

open SimpleGraph Finset

/-! ### The spider is a tree -/

variable {n : ℕ} {len : Fin n → ℕ}

/-- The neighbour of a spider vertex one step closer to the hub. -/
def spiderParent : SpiderV n len → SpiderV n len
  | none => none
  | some ⟨i, j⟩ =>
      if j.val = 0 then none
      else some ⟨i, ⟨j.val - 1, lt_of_le_of_lt (Nat.sub_le _ _) j.isLt⟩⟩

theorem spiderLevel_spiderParent_lt {v : SpiderV n len} (hv : v ≠ none) :
    spiderLevel (spiderParent v) < spiderLevel v := by
  match v with
  | none => exact absurd rfl hv
  | some ⟨i, j⟩ =>
      by_cases h : j.val = 0
      · simp only [spiderParent, if_pos h, spiderLevel]
        omega
      · simp only [spiderParent, if_neg h, spiderLevel]
        omega

theorem spider_adj_spiderParent {v : SpiderV n len} (hv : v ≠ none) :
    (spider n len).Adj (spiderParent v) v := by
  match v with
  | none => exact absurd rfl hv
  | some ⟨i, j⟩ =>
      by_cases h : j.val = 0
      · rw [spiderParent, if_pos h]
        exact spider_adj_hub (q := ⟨i, j⟩) h
      · rw [spiderParent, if_neg h]
        exact spider_adj_arm (i := i) (j := ⟨j.val - 1, _⟩) (j' := j) (by simp; omega)

theorem spiderLevel_ne_of_adj {x y : SpiderV n len} (h : (spider n len).Adj x y) :
    spiderLevel x ≠ spiderLevel y := by
  have := spiderLevel_adj h
  omega

theorem spiderLevel_hub_le (v : SpiderV n len) :
    spiderLevel (none : SpiderV n len) ≤ spiderLevel v := Nat.zero_le _

theorem spiderParent_unique {u v : SpiderV n len} (h : (spider n len).Adj u v)
    (hlev : spiderLevel u < spiderLevel v) : u = spiderParent v := by
  match u, v with
  | none, none => exact absurd h (by simp [spider])
  | none, some ⟨i, j⟩ =>
      have hj : j.val = 0 := h
      rw [spiderParent, if_pos hj]
  | some ⟨i, j⟩, none =>
      have hlev' : j.val + 1 < 0 := hlev
      omega
  | some ⟨i, j⟩, some ⟨i', j'⟩ =>
      have hi : i = i' := h.1
      have hlev' : j.val + 1 < j'.val + 1 := hlev
      have hj : j.val + 1 = j'.val := by
        rcases h.2 with h' | h' <;> simp only at h' <;> omega
      have hj0 : ¬ j'.val = 0 := by omega
      rw [spiderParent, if_neg hj0]
      exact congrArg some (sigma_fin_ext hi (show j.val = j'.val - 1 by omega))

theorem spider_isTree (n : ℕ) (len : Fin n → ℕ) : (spider n len).IsTree := by
  classical
  exact isTree_of_rank_parent (spider_connected n len) none spiderLevel spiderParent
    (fun v hv => spiderLevel_spiderParent_lt hv)
    (fun v hv => spider_adj_spiderParent hv)
    (fun _ _ h => spiderLevel_ne_of_adj h)
    spiderLevel_hub_le
    (fun _ _ h hlev => spiderParent_unique h hlev)

/-! ### Only the hub of a spider can be a branch vertex -/

theorem spider_eq_of_adj_of_level_eq {x y z : SpiderV n len} (hx : x ≠ none)
    (hy : (spider n len).Adj x y) (hz : (spider n len).Adj x z)
    (hlev : spiderLevel y = spiderLevel z) : y = z := by
  match x with
  | none => exact absurd rfl hx
  | some ⟨i, j⟩ =>
      match y, z with
      | none, none => rfl
      | none, some ⟨i', j'⟩ => exact absurd hlev (by simp [spiderLevel])
      | some ⟨i', j'⟩, none => exact absurd hlev (by simp [spiderLevel])
      | some ⟨i', j'⟩, some ⟨i'', j''⟩ =>
          have h1 : i = i' := hy.1
          have h2 : i = i'' := hz.1
          have hval : j'.val = j''.val := by
            simp only [spiderLevel] at hlev
            omega
          exact congrArg some (sigma_fin_ext (h1 ▸ h2) hval)

theorem spider_degree_le_two [DecidableRel (spider n len).Adj] (x : SpiderV n len)
    (hx : x ≠ none) : (spider n len).degree x ≤ 2 := by
  classical
  have hmap : Set.MapsTo (fun y => decide (spiderLevel y = spiderLevel x + 1))
      (((spider n len).neighborFinset x : Finset (SpiderV n len)) : Set (SpiderV n len))
      ((Finset.univ : Finset Bool) : Set Bool) := by
    intro y _
    exact Finset.mem_coe.mpr (Finset.mem_univ _)
  have hinj : Set.InjOn (fun y => decide (spiderLevel y = spiderLevel x + 1))
      (((spider n len).neighborFinset x : Finset (SpiderV n len)) : Set (SpiderV n len)) := by
    intro y hy z hz hyz
    rw [Finset.mem_coe, SimpleGraph.mem_neighborFinset] at hy hz
    have hly := spiderLevel_adj hy
    have hlz := spiderLevel_adj hz
    refine spider_eq_of_adj_of_level_eq hx hy hz ?_
    by_cases hy1 : spiderLevel y = spiderLevel x + 1
    · have hz1 : spiderLevel z = spiderLevel x + 1 := by
        by_contra hc
        simp [hy1, hc] at hyz
      omega
    · have hz1 : ¬ spiderLevel z = spiderLevel x + 1 := by
        by_contra hc
        simp [hy1, hc] at hyz
      omega
  have hcard := Finset.card_le_card_of_injOn _ hmap hinj
  simpa using hcard

theorem spider_branchVerts_card_le_one [DecidableRel (spider n len).Adj] :
    (branchVerts (spider n len)).card ≤ 1 := by
  classical
  have hsub : branchVerts (spider n len) ⊆ {none} := by
    intro x hx
    rw [mem_branchVerts] at hx
    rw [Finset.mem_singleton]
    by_contra hc
    have := spider_degree_le_two x hc
    omega
  calc (branchVerts (spider n len)).card ≤ ({none} : Finset (SpiderV n len)).card :=
        Finset.card_le_card hsub
    _ = 1 := Finset.card_singleton _

/-! ### The connector tree is a tree -/

section Connector

variable {t m : ℕ} {a : Fin m → ℕ} {b : Fin n → ℕ}

/-- The neighbour of a connector-tree vertex one step closer to the left hub. -/
def connParent : ConnV t m a n b → ConnV t m a n b
  | Sum.inl p => Sum.inl (spiderParent p)
  | Sum.inr (Sum.inl i) =>
      if h : i.val = 0 then Sum.inl none
      else Sum.inr (Sum.inl ⟨i.val - 1, lt_of_le_of_lt (Nat.sub_le _ _) i.isLt⟩)
  | Sum.inr (Sum.inr none) =>
      if h : t = 0 then Sum.inl none else Sum.inr (Sum.inl ⟨t - 1, by omega⟩)
  | Sum.inr (Sum.inr (some q)) => Sum.inr (Sum.inr (spiderParent (some q)))

theorem clevel_connParent_lt {v : ConnV t m a n b} (hv : v ≠ chubL t m a n b) :
    clevel (connParent v) < clevel v := by
  match v with
  | Sum.inl none => exact absurd rfl hv
  | Sum.inl (some p) =>
      exact spiderLevel_spiderParent_lt (v := some p) (by simp)
  | Sum.inr (Sum.inl i) =>
      by_cases h : i.val = 0
      · rw [connParent, dif_pos h]
        show spiderLevel (none : SpiderV m a) < i.val + 1
        simp [spiderLevel]
      · rw [connParent, dif_neg h]
        show i.val - 1 + 1 < i.val + 1
        omega
  | Sum.inr (Sum.inr none) =>
      by_cases h : t = 0
      · rw [connParent, dif_pos h]
        show spiderLevel (none : SpiderV m a) < spiderLevel (none : SpiderV n b) + t + 1
        simp [spiderLevel]
      · rw [connParent, dif_neg h]
        show t - 1 + 1 < spiderLevel (none : SpiderV n b) + t + 1
        simp [spiderLevel]
        omega
  | Sum.inr (Sum.inr (some q)) =>
      have := spiderLevel_spiderParent_lt (v := some q) (by simp)
      show spiderLevel (spiderParent (some q)) + t + 1 < spiderLevel (some q) + t + 1
      omega

theorem connGraph_adj_connParent {v : ConnV t m a n b} (hv : v ≠ chubL t m a n b) :
    (connGraph t m a n b).Adj (connParent v) v := by
  match v with
  | Sum.inl none => exact absurd rfl hv
  | Sum.inl (some p) =>
      exact connGraph_adj_ll.mpr (spider_adj_spiderParent (v := some p) (by simp))
  | Sum.inr (Sum.inl i) =>
      by_cases h : i.val = 0
      · rw [connParent, dif_pos h]
        exact connGraph_adj_lc.mpr ⟨rfl, h⟩
      · rw [connParent, dif_neg h]
        exact connGraph_adj_cc.mpr (Or.inl (by simp; omega))
  | Sum.inr (Sum.inr none) =>
      by_cases h : t = 0
      · rw [connParent, dif_pos h]
        exact connGraph_adj_lr.mpr ⟨rfl, rfl, h⟩
      · rw [connParent, dif_neg h]
        exact connGraph_adj_cr.mpr ⟨rfl, by simp; omega⟩
  | Sum.inr (Sum.inr (some q)) =>
      exact connGraph_adj_rr.mpr (spider_adj_spiderParent (v := some q) (by simp))

theorem clevel_ne_of_adj {x y : ConnV t m a n b} (h : (connGraph t m a n b).Adj x y) :
    clevel x ≠ clevel y := by
  have := clevel_adj h
  omega

theorem clevel_hubL_le (v : ConnV t m a n b) : clevel (chubL t m a n b) ≤ clevel v :=
  Nat.zero_le _

theorem connParent_unique {u v : ConnV t m a n b} (h : (connGraph t m a n b).Adj u v)
    (hlev : clevel u < clevel v) : u = connParent v := by
  match v with
  | Sum.inl p =>
      match u with
      | Sum.inl p' =>
          exact congrArg Sum.inl (spiderParent_unique (connGraph_adj_ll.mp h) hlev)
      | Sum.inr (Sum.inl j) =>
          exfalso
          obtain ⟨hp, -⟩ := h
          subst hp
          have : j.val + 1 < spiderLevel (none : SpiderV m a) := hlev
          simp [spiderLevel] at this
      | Sum.inr (Sum.inr q) =>
          exfalso
          obtain ⟨hp, hq, ht⟩ := h
          subst hp
          have : spiderLevel q + t + 1 < spiderLevel (none : SpiderV m a) := hlev
          simp [spiderLevel] at this
  | Sum.inr (Sum.inl i) =>
      match u with
      | Sum.inl p =>
          obtain ⟨hp, hi⟩ := connGraph_adj_lc.mp h
          subst hp
          rw [connParent, dif_pos hi]
      | Sum.inr (Sum.inl j) =>
          have hcc := connGraph_adj_cc.mp h
          have hlev' : j.val + 1 < i.val + 1 := hlev
          have hji : j.val + 1 = i.val := by omega
          have hi0 : ¬ i.val = 0 := by omega
          rw [connParent, dif_neg hi0]
          exact congrArg (fun z => Sum.inr (Sum.inl z)) (Fin.ext (show j.val = i.val - 1 by omega))
      | Sum.inr (Sum.inr q) =>
          exfalso
          obtain ⟨hq, hi⟩ := h
          subst hq
          have hlev' : spiderLevel (none : SpiderV n b) + t + 1 < i.val + 1 := hlev
          simp [spiderLevel] at hlev'
          omega
  | Sum.inr (Sum.inr none) =>
      match u with
      | Sum.inl p =>
          obtain ⟨hp, -, ht⟩ := connGraph_adj_lr.mp h
          subst hp
          rw [connParent, dif_pos ht]
      | Sum.inr (Sum.inl j) =>
          obtain ⟨-, hj⟩ := connGraph_adj_cr.mp h
          have ht : ¬ t = 0 := by omega
          rw [connParent, dif_neg ht]
          exact congrArg (fun z => Sum.inr (Sum.inl z)) (Fin.ext (show j.val = t - 1 by omega))
      | Sum.inr (Sum.inr q) =>
          exfalso
          have hq : (spider n b).Adj q none := connGraph_adj_rr.mp h
          have hlev' : spiderLevel q + t + 1 < spiderLevel (none : SpiderV n b) + t + 1 := hlev
          have := spiderLevel_adj hq
          simp [spiderLevel] at hlev'
  | Sum.inr (Sum.inr (some q)) =>
      match u with
      | Sum.inl p =>
          exact absurd h.2.1 (by simp)
      | Sum.inr (Sum.inl j) =>
          exact absurd h.1 (by simp)
      | Sum.inr (Sum.inr q') =>
          have hq : (spider n b).Adj q' (some q) := connGraph_adj_rr.mp h
          have hlev' : spiderLevel q' + t + 1 < spiderLevel (some q) + t + 1 := hlev
          exact congrArg (fun z => Sum.inr (Sum.inr z)) (spiderParent_unique hq (by omega))

theorem connGraph_isTree (t m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) :
    (connGraph t m a n b).IsTree := by
  classical
  exact isTree_of_rank_parent (connGraph_connected t m a n b) (chubL t m a n b) clevel connParent
    (fun v hv => clevel_connParent_lt hv)
    (fun v hv => connGraph_adj_connParent hv)
    (fun _ _ h => clevel_ne_of_adj h)
    clevel_hubL_le
    (fun _ _ h hlev => connParent_unique h hlev)

theorem connGraph_branchVerts_card_le_two [DecidableRel (connGraph t m a n b).Adj] :
    (branchVerts (connGraph t m a n b)).card ≤ 2 := by
  classical
  have hsub : branchVerts (connGraph t m a n b) ⊆ {chubL t m a n b, chubR t m a n b} := by
    intro x hx
    rw [mem_branchVerts] at hx
    have := connGraph_branch_vertex x (by convert hx using 2)
    rcases this with h | h
    · exact Finset.mem_insert.mpr (Or.inl h)
    · exact Finset.mem_insert.mpr (Or.inr (Finset.mem_singleton.mpr h))
  refine le_trans (Finset.card_le_card hsub) ?_
  exact le_trans (Finset.card_insert_le _ _) (by simp)

end Connector

/-! ### The characterization of the `C₂` family -/

theorem isTree_of_isC2Model {V : Type*} [Fintype V] [DecidableEq V] {G : SimpleGraph V}
    (hG : IsC2Model G) : G.IsTree := by
  classical
  rcases hG with ⟨k, l, ⟨e⟩⟩ | ⟨t, m, a, k, c, ⟨e⟩⟩
  · exact (Iso.isTree_iff e).mpr (spider_isTree k l)
  · exact (Iso.isTree_iff e).mpr (connGraph_isTree t m a k c)

/-- The number of branch vertices is an isomorphism invariant. -/
theorem branchVerts_card_iso {V W : Type*} [Fintype V] [DecidableEq V] [Fintype W]
    [DecidableEq W] {G : SimpleGraph V} {H : SimpleGraph W} [DecidableRel G.Adj]
    [DecidableRel H.Adj] (e : G ≃g H) : (branchVerts G).card = (branchVerts H).card := by
  refine Finset.card_bij (fun v _ => e v) ?_ ?_ ?_
  · intro v hv
    rw [mem_branchVerts] at hv
    rw [mem_branchVerts, Iso.degree_eq e v]
    exact hv
  · intro u _ v _ huv
    exact e.injective huv
  · intro w hw
    refine ⟨e.symm w, ?_, by simp⟩
    rw [mem_branchVerts] at hw ⊢
    rw [← Iso.degree_eq e (e.symm w)]
    simpa using hw

theorem branchVerts_card_le_two_of_isC2Model {V : Type*} [Fintype V] [DecidableEq V]
    {G : SimpleGraph V} [DecidableRel G.Adj] (hG : IsC2Model G) :
    (branchVerts G).card ≤ 2 := by
  classical
  rcases hG with ⟨k, l, ⟨e⟩⟩ | ⟨t, m, a, k, c, ⟨e⟩⟩
  · rw [branchVerts_card_iso e]
    exact le_trans spider_branchVerts_card_le_one (by norm_num)
  · rw [branchVerts_card_iso e]
    exact connGraph_branchVerts_card_le_two

/-- **Characterization of the `C₂` family.**  A finite graph is a `C₂` model exactly when
it is a tree with at most two branch vertices. -/
theorem isC2Model_iff_isTree_branchVerts {V : Type*} [Fintype V] [DecidableEq V]
    {G : SimpleGraph V} [DecidableRel G.Adj] :
    IsC2Model G ↔ G.IsTree ∧ (branchVerts G).card ≤ 2 :=
  ⟨fun h => ⟨isTree_of_isC2Model h, branchVerts_card_le_two_of_isC2Model h⟩,
    fun h => isC2Model_of_isTree_of_branchVerts_card_le_two h.1 h.2⟩

end ClanAudit
