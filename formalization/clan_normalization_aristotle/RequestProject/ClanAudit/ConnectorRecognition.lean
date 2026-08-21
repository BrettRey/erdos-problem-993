import RequestProject.ClanAudit.TreeTwoHub
import RequestProject.ClanAudit.Connector

/-!
# Recognising a tree with two branch vertices as the connector model

Let `G` be a finite tree with two distinct vertices `h1 ≠ h2` such that every *other*
vertex has degree at most two.  This file proves that `G` is isomorphic to the explicit
connector model `connGraph t m a n b`, where

* `t + 1` is the distance from `h1` to `h2` (the number of connector edges),
* `m` is the number of arms of `h1` other than the one containing `h2`, with `a i` the
  size of the `i`-th such arm,
* `n` and `b j` are the same data at `h2`.

The vertex map is the "polar coordinate" map on each side: on the `h1` side a vertex is
recorded by its arm at `h1` and its depth from `h1`, on the `h2` side by its arm at `h2`
and its depth from `h2`, and the vertices strictly between the hubs become the interior
connector vertices.

The structural facts about the two rootings are collected in the first section.
-/

namespace ClanAudit

open SimpleGraph Finset

variable {V : Type*} [Fintype V] [DecidableEq V] {G : SimpleGraph V} [DecidableRel G.Adj]

/-! ### The two sides of a two-hub tree -/

variable (hG : G.IsTree) (h1 h2 : V)

/-- `v` is on the far side of `h2`: its root path from `h1` passes through `h2`. -/
def IsRightSide (v : V) : Prop := h2 ∈ (treePath hG h1 v).support

variable {hG h1 h2}

omit [Fintype V] [DecidableEq V] [DecidableRel G.Adj] in
theorem depth_hub_pos (hne : h1 ≠ h2) : 0 < treeDepth hG h1 h2 :=
  treeDepth_pos (hG := hG) (h := h1) (Ne.symm hne)

omit [Fintype V] [DecidableEq V] [DecidableRel G.Adj] in
theorem isRightSide_self : IsRightSide hG h1 h2 h2 := Walk.end_mem_support _

omit [Fintype V] [DecidableRel G.Adj] in
/-- Every vertex of the far side lies in the arm of `h1` that contains `h2`. -/
theorem arm1_of_right {v : V} (hne : h1 ≠ h2) (hr : IsRightSide hG h1 h2 v) :
    treeArm hG h1 v = treeArm hG h1 h2 :=
  treeArm_of_mem_support (hG := hG) (h := h1) hr (Ne.symm hne)

omit [Fintype V] [DecidableRel G.Adj] in
/-- Depths add up across `h2`. -/
theorem depth1_of_right {v : V} (hr : IsRightSide hG h1 h2 v) :
    treeDepth hG h1 v = treeDepth hG h1 h2 + treeDepth hG h2 v :=
  treeDepth_add (hG := hG) (h := h1) hr

omit [Fintype V] [DecidableRel G.Adj] in
theorem not_right_of_depth_lt {v : V} (hlt : treeDepth hG h1 v < treeDepth hG h1 h2) :
    ¬ IsRightSide hG h1 h2 v := by
  intro hr
  have := depth1_of_right hr
  omega

omit [DecidableEq V] in
/-- Inside the arm of `h1` containing `h2`, every vertex strictly shallower than `h2` has
degree at most two. -/
theorem deg_le_two_of_lt (hdeg : ∀ w : V, w ≠ h1 → w ≠ h2 → G.degree w ≤ 2) {w : V}
    (hw1 : w ≠ h1) (hlt : treeDepth hG h1 w < treeDepth hG h1 h2) : G.degree w ≤ 2 := by
  refine hdeg w hw1 ?_
  intro hc
  rw [hc] at hlt
  omega

/-- Depth-injectivity inside the arm of `h1` containing `h2`, up to the depth of `h2`. -/
theorem uniq_arm_hub (hdeg : ∀ w : V, w ≠ h1 → w ≠ h2 → G.degree w ≤ 2) :
    ∀ k : ℕ, k ≤ treeDepth hG h1 h2 → ∀ u v : V, u ≠ h1 → v ≠ h1 →
      treeArm hG h1 u = treeArm hG h1 h2 → treeArm hG h1 v = treeArm hG h1 h2 →
      treeDepth hG h1 u = k → treeDepth hG h1 v = k → u = v :=
  treeArm_depth_inj_lt (hG := hG) (h := h1) (x := treeArm hG h1 h2)
    (D := treeDepth hG h1 h2) (fun _ hw _ hlt => deg_le_two_of_lt hdeg hw hlt)

/-- A vertex of the arm of `h1` containing `h2` which is not beyond `h2` is strictly
shallower than `h2`. -/
theorem depth1_lt_of_not_right (hne : h1 ≠ h2)
    (hdeg : ∀ w : V, w ≠ h1 → w ≠ h2 → G.degree w ≤ 2) {v : V}
    (hnr : ¬ IsRightSide hG h1 h2 v) (hc : treeArm hG h1 v = treeArm hG h1 h2) :
    treeDepth hG h1 v < treeDepth hG h1 h2 := by
  by_contra hge
  push_neg at hge
  obtain ⟨w, hw, hwd⟩ := exists_ancestor (hG := hG) (h := h1) v (treeDepth hG h1 h2) hge
  have hw1 : w ≠ h1 := by
    intro hcw
    rw [hcw, treeDepth_root (hG := hG)] at hwd
    exact absurd hwd.symm (by have := depth_hub_pos (hG := hG) hne; omega)
  have hwc : treeArm hG h1 w = treeArm hG h1 h2 := by
    rw [← (treeArm_of_mem_support (hG := hG) (h := h1) hw hw1), hc]
  have hwh2 : w = h2 :=
    uniq_arm_hub hdeg (treeDepth hG h1 h2) le_rfl w h2 hw1 (Ne.symm hne) hwc rfl hwd rfl
  exact hnr (hwh2 ▸ hw)

/-! ### Comparing the two rootings -/

omit [Fintype V] [DecidableRel G.Adj] in
/-- A vertex beyond `h2` is not the left hub. -/
theorem ne_h1_of_right (hne : h1 ≠ h2) {v : V} (hr : IsRightSide hG h1 h2 v) : v ≠ h1 := by
  intro hc
  have h2' := depth1_of_right hr
  have h3' := depth_hub_pos (hG := hG) (h1 := h1) (h2 := h2) hne
  rw [hc, treeDepth_root (hG := hG)] at h2'
  omega

omit [Fintype V] [DecidableRel G.Adj] in
/-- A vertex beyond `h2` does not lie in the arm of `h2` pointing back at `h1`. -/
theorem arm2_ne_of_right (hne : h1 ≠ h2) {v : V} (hr : IsRightSide hG h1 h2 v)
    (hv2 : v ≠ h2) : treeArm hG h2 v ≠ treeArm hG h2 h1 := by
  intro harm
  have hdec : treePath hG h1 v = (treePath hG h1 h2).append (treePath hG h2 v) :=
    treePath_append (hG := hG) (h := h1) hr
  -- the arm of `h2` towards `h1` lies on both pieces
  have hc'ne : treeArm hG h2 h1 ≠ h2 := treeArm_ne_root (hG := hG) (h := h2) hne
  have hmem1 : treeArm hG h2 h1 ∈ (treePath hG h1 h2).support := by
    refine mem_support_treePath_comm (hG := hG) (h := h2) ?_
    exact Walk.getVert_mem_support _ _
  have hmem2 : treeArm hG h2 h1 ∈ (treePath hG h2 v).support := by
    obtain ⟨w, hw, hwd⟩ := exists_ancestor (hG := hG) (h := h2) v 1
      (treeDepth_pos (hG := hG) (h := h2) hv2)
    have hw2 : w ≠ h2 := by
      intro hcw
      rw [hcw, treeDepth_root (hG := hG)] at hwd
      omega
    have : treeArm hG h2 v = w := by
      rw [treeArm_of_mem_support (hG := hG) (h := h2) hw hw2,
        treeArm_of_depth_one (hG := hG) (h := h2) hwd]
    rw [← harm, this]
    exact hw
  -- but the concatenation is a path
  have hnd : (treePath hG h1 v).support.Nodup :=
    (treePath_isPath (hG := hG) (h := h1) v).support_nodup
  rw [hdec, Walk.support_append] at hnd
  have hdisj := (List.nodup_append.mp hnd).2.2
  have htail : treeArm hG h2 h1 ∈ (treePath hG h2 v).support.tail := by
    have hcons : (treePath hG h2 v).support = h2 :: (treePath hG h2 v).support.tail :=
      Walk.support_eq_cons _
    rw [hcons, List.mem_cons] at hmem2
    exact hmem2.resolve_left hc'ne
  exact hdisj _ hmem1 _ htail rfl

omit [Fintype V] [DecidableRel G.Adj] in
/-- Conversely, a vertex outside the arm of `h2` pointing back at `h1` is beyond `h2`. -/
theorem right_of_arm2_ne {v : V}
    (harm : treeArm hG h2 v ≠ treeArm hG h2 h1) : IsRightSide hG h1 h2 v := by
  have hpath : ((treePath hG h1 h2).append (treePath hG h2 v)).IsPath := by
    refine walk_isPath_append (treePath_isPath (hG := hG) (h := h1) h2)
      (treePath_isPath (hG := hG) (h := h2) v) ?_
    intro x hx hx'
    by_contra hxh2
    have hx1 : x ∈ (treePath hG h2 h1).support :=
      mem_support_treePath_comm (hG := hG) (h := h1) hx
    have e1 : treeArm hG h2 h1 = treeArm hG h2 x :=
      treeArm_of_mem_support (hG := hG) (h := h2) hx1 hxh2
    have e2 : treeArm hG h2 v = treeArm hG h2 x :=
      treeArm_of_mem_support (hG := hG) (h := h2) hx' hxh2
    exact harm (by rw [e1, e2])
  have := treePath_eq (hG := hG) (h := h1) (v := v) hpath
  show h2 ∈ (treePath hG h1 v).support
  rw [← this, Walk.support_append]
  exact List.mem_append.mpr (Or.inl (Walk.end_mem_support _))

omit [Fintype V] [DecidableRel G.Adj] in
/-- On the near side of `h2`, every vertex lies in the arm of `h2` pointing at `h1`. -/
theorem arm2_of_not_right {v : V} (hnr : ¬ IsRightSide hG h1 h2 v) :
    treeArm hG h2 v = treeArm hG h2 h1 := by
  by_contra hc
  exact hnr (right_of_arm2_ne hc)

omit [DecidableEq V] in
/-- Every vertex of an arm of `h1` other than the one containing `h2` has degree at most
two. -/
theorem deg_le_two_left_arm (hdeg : ∀ w : V, w ≠ h1 → w ≠ h2 → G.degree w ≤ 2) {x : V}
    (hx : x ≠ treeArm hG h1 h2) {w : V} (hw1 : w ≠ h1) (hwx : treeArm hG h1 w = x) :
    G.degree w ≤ 2 := by
  refine hdeg w hw1 ?_
  intro hc
  rw [hc] at hwx
  exact hx hwx.symm

omit [DecidableEq V] in
/-- Every vertex of an arm of `h2` other than the one pointing at `h1` has degree at most
two. -/
theorem deg_le_two_right_arm (hdeg : ∀ w : V, w ≠ h1 → w ≠ h2 → G.degree w ≤ 2) {y : V}
    (hy : y ≠ treeArm hG h2 h1) {w : V} (hw2 : w ≠ h2) (hwy : treeArm hG h2 w = y) :
    G.degree w ≤ 2 := by
  refine hdeg w ?_ hw2
  intro hc
  rw [hc] at hwy
  exact hy hwy.symm

omit [Fintype V] [DecidableRel G.Adj] in
/-- Every vertex of an arm of `h1` other than the one containing `h2` is on the near
side. -/
theorem not_right_of_arm1_ne (hne : h1 ≠ h2) {v : V}
    (hvc : treeArm hG h1 v ≠ treeArm hG h1 h2) : ¬ IsRightSide hG h1 h2 v := by
  intro hr
  exact hvc (arm1_of_right hne hr)

omit [Fintype V] [DecidableEq V] [DecidableRel G.Adj] in
/-- Every vertex of an arm of `h2` other than the one pointing at `h1` is different from
`h1`. -/
theorem ne_h1_of_arm2_ne {v : V} (hvc : treeArm hG h2 v ≠ treeArm hG h2 h1) : v ≠ h1 := by
  intro hc
  exact hvc (by rw [hc])

/-! ### The parameters of the connector model -/

variable (hG h1 h2)

/-- The arms of `h1` other than the one containing `h2`. -/
noncomputable def leftArms : Finset V := (G.neighborFinset h1).erase (treeArm hG h1 h2)

/-- The arms of `h2` other than the one pointing at `h1`. -/
noncomputable def rightArms : Finset V := (G.neighborFinset h2).erase (treeArm hG h2 h1)

/-- The `i`-th arm of `h1`. -/
noncomputable def leftArm (i : Fin (leftArms hG h1 h2).card) : V :=
  ((leftArms hG h1 h2).equivFin.symm i : V)

/-- The `j`-th arm of `h2`. -/
noncomputable def rightArm (j : Fin (rightArms hG h1 h2).card) : V :=
  ((rightArms hG h1 h2).equivFin.symm j : V)

/-- The number of vertices of the `i`-th arm of `h1`. -/
noncomputable def leftLen (i : Fin (leftArms hG h1 h2).card) : ℕ :=
  (armFinset hG h1 (leftArm hG h1 h2 i)).card

/-- The number of vertices of the `j`-th arm of `h2`. -/
noncomputable def rightLen (j : Fin (rightArms hG h1 h2).card) : ℕ :=
  (armFinset hG h2 (rightArm hG h1 h2 j)).card

variable {hG h1 h2}

theorem leftArm_mem (i : Fin (leftArms hG h1 h2).card) : leftArm hG h1 h2 i ∈ leftArms hG h1 h2 :=
  ((leftArms hG h1 h2).equivFin.symm i).2

theorem rightArm_mem (j : Fin (rightArms hG h1 h2).card) :
    rightArm hG h1 h2 j ∈ rightArms hG h1 h2 :=
  ((rightArms hG h1 h2).equivFin.symm j).2

theorem leftArm_injective : Function.Injective (leftArm hG h1 h2) := fun _ _ hij =>
  (leftArms hG h1 h2).equivFin.symm.injective (Subtype.coe_injective hij)

theorem rightArm_injective : Function.Injective (rightArm hG h1 h2) := fun _ _ hij =>
  (rightArms hG h1 h2).equivFin.symm.injective (Subtype.coe_injective hij)

theorem leftArm_ne (i : Fin (leftArms hG h1 h2).card) :
    leftArm hG h1 h2 i ≠ treeArm hG h1 h2 :=
  (Finset.mem_erase.mp (leftArm_mem i)).1

theorem rightArm_ne (j : Fin (rightArms hG h1 h2).card) :
    rightArm hG h1 h2 j ≠ treeArm hG h2 h1 :=
  (Finset.mem_erase.mp (rightArm_mem j)).1

/-- The index of the arm of `h1` containing a near-side vertex. -/
noncomputable def leftIdx {v : V} (hv1 : v ≠ h1) (hvc : treeArm hG h1 v ≠ treeArm hG h1 h2) :
    Fin (leftArms hG h1 h2).card :=
  (leftArms hG h1 h2).equivFin ⟨treeArm hG h1 v, by
    rw [leftArms, Finset.mem_erase]
    exact ⟨hvc, by rw [mem_neighborFinset]; exact treeArm_adj_root (hG := hG) (h := h1) hv1⟩⟩

/-- The index of the arm of `h2` containing a far-side vertex. -/
noncomputable def rightIdx {v : V} (hv2 : v ≠ h2) (hvc : treeArm hG h2 v ≠ treeArm hG h2 h1) :
    Fin (rightArms hG h1 h2).card :=
  (rightArms hG h1 h2).equivFin ⟨treeArm hG h2 v, by
    rw [rightArms, Finset.mem_erase]
    exact ⟨hvc, by rw [mem_neighborFinset]; exact treeArm_adj_root (hG := hG) (h := h2) hv2⟩⟩

theorem leftArm_leftIdx {v : V} (hv1 : v ≠ h1) (hvc : treeArm hG h1 v ≠ treeArm hG h1 h2) :
    leftArm hG h1 h2 (leftIdx hv1 hvc) = treeArm hG h1 v := by
  simp [leftArm, leftIdx]

theorem rightArm_rightIdx {v : V} (hv2 : v ≠ h2) (hvc : treeArm hG h2 v ≠ treeArm hG h2 h1) :
    rightArm hG h1 h2 (rightIdx hv2 hvc) = treeArm hG h2 v := by
  simp [rightArm, rightIdx]

theorem leftIdx_eq {v : V} (hv1 : v ≠ h1) (hvc : treeArm hG h1 v ≠ treeArm hG h1 h2)
    (i : Fin (leftArms hG h1 h2).card) (hi : treeArm hG h1 v = leftArm hG h1 h2 i) :
    leftIdx hv1 hvc = i := by
  refine leftArm_injective ?_
  rw [leftArm_leftIdx hv1 hvc, hi]

theorem rightIdx_eq {v : V} (hv2 : v ≠ h2) (hvc : treeArm hG h2 v ≠ treeArm hG h2 h1)
    (j : Fin (rightArms hG h1 h2).card) (hj : treeArm hG h2 v = rightArm hG h1 h2 j) :
    rightIdx hv2 hvc = j := by
  refine rightArm_injective ?_
  rw [rightArm_rightIdx hv2 hvc, hj]

/-! ### The vertex map -/

variable (hG)

/-- The model recognized: a connector tree with `treeDepth hG h1 h2 - 1` interior
connector vertices and the arms of the two hubs. -/
noncomputable abbrev connModel : SimpleGraph (ConnV (treeDepth hG h1 h2 - 1)
    (leftArms hG h1 h2).card (leftLen hG h1 h2) (rightArms hG h1 h2).card (rightLen hG h1 h2)) :=
  connGraph (treeDepth hG h1 h2 - 1) (leftArms hG h1 h2).card (leftLen hG h1 h2)
    (rightArms hG h1 h2).card (rightLen hG h1 h2)

open scoped Classical in
/-- The polar-coordinate map of a two-hub tree into the connector model. -/
noncomputable def toConn (hne : h1 ≠ h2) (hdeg : ∀ w : V, w ≠ h1 → w ≠ h2 → G.degree w ≤ 2) :
    V → ConnV (treeDepth hG h1 h2 - 1) (leftArms hG h1 h2).card (leftLen hG h1 h2)
      (rightArms hG h1 h2).card (rightLen hG h1 h2) := fun v =>
  if hv1 : v = h1 then Sum.inl none
  else if hv2 : v = h2 then Sum.inr (Sum.inr none)
  else if hr : IsRightSide hG h1 h2 v then
    Sum.inr (Sum.inr (some ⟨rightIdx hv2 (arm2_ne_of_right hne hr hv2),
      ⟨treeDepth hG h2 v - 1, by
        have h1' : rightArm hG h1 h2 (rightIdx hv2 (arm2_ne_of_right hne hr hv2))
            = treeArm hG h2 v := rightArm_rightIdx _ _
        have h2' : treeDepth hG h2 v ≤ (armFinset hG h2 (treeArm hG h2 v)).card :=
          treeDepth_le_card_arm (hG := hG) (h := h2)
        have h3' : 0 < treeDepth hG h2 v := treeDepth_pos (hG := hG) (h := h2) hv2
        simp only [rightLen, h1']
        omega⟩⟩))
  else if hc : treeArm hG h1 v = treeArm hG h1 h2 then
    Sum.inr (Sum.inl ⟨treeDepth hG h1 v - 1, by
      have hlt : treeDepth hG h1 v < treeDepth hG h1 h2 :=
        depth1_lt_of_not_right hne hdeg hr hc
      have h3' : 0 < treeDepth hG h1 v := treeDepth_pos (hG := hG) (h := h1) hv1
      omega⟩)
  else
    Sum.inl (some ⟨leftIdx hv1 hc, ⟨treeDepth hG h1 v - 1, by
      have h1' : leftArm hG h1 h2 (leftIdx hv1 hc) = treeArm hG h1 v := leftArm_leftIdx _ _
      have h2' : treeDepth hG h1 v ≤ (armFinset hG h1 (treeArm hG h1 v)).card :=
        treeDepth_le_card_arm (hG := hG) (h := h1)
      have h3' : 0 < treeDepth hG h1 v := treeDepth_pos (hG := hG) (h := h1) hv1
      simp only [leftLen, h1']
      omega⟩⟩)

variable (hne : h1 ≠ h2) (hdeg : ∀ w : V, w ≠ h1 → w ≠ h2 → G.degree w ≤ 2)

theorem toConn_h1 : toConn hG hne hdeg h1 = Sum.inl none := by
  simp [toConn]

theorem toConn_h2 : toConn hG hne hdeg h2 = Sum.inr (Sum.inr none) := by
  simp [toConn, Ne.symm hne]

theorem toConn_right {v : V} (hv2 : v ≠ h2) (hr : IsRightSide hG h1 h2 v) :
    ∃ hlt : treeDepth hG h2 v - 1 <
        rightLen hG h1 h2 (rightIdx hv2 (arm2_ne_of_right hne hr hv2)),
      toConn hG hne hdeg v =
        Sum.inr (Sum.inr (some ⟨rightIdx hv2 (arm2_ne_of_right hne hr hv2),
          ⟨treeDepth hG h2 v - 1, hlt⟩⟩)) := by
  have hv1 : v ≠ h1 := by
    intro hc
    subst hc
    rw [IsRightSide, treePath_root (hG := hG)] at hr
    simp only [Walk.support_nil, List.mem_singleton] at hr
    exact hne hr.symm
  refine ⟨?_, ?_⟩
  · have h1' : rightArm hG h1 h2 (rightIdx hv2 (arm2_ne_of_right hne hr hv2))
        = treeArm hG h2 v := rightArm_rightIdx _ _
    have h2' : treeDepth hG h2 v ≤ (armFinset hG h2 (treeArm hG h2 v)).card :=
      treeDepth_le_card_arm (hG := hG) (h := h2)
    have h3' : 0 < treeDepth hG h2 v := treeDepth_pos (hG := hG) (h := h2) hv2
    simp only [rightLen, h1']
    omega
  · simp only [toConn, dif_neg hv1, dif_neg hv2, dif_pos hr]

theorem toConn_interior {v : V} (hv1 : v ≠ h1) (hv2 : v ≠ h2)
    (hnr : ¬ IsRightSide hG h1 h2 v) (hc : treeArm hG h1 v = treeArm hG h1 h2) :
    ∃ hlt : treeDepth hG h1 v - 1 < treeDepth hG h1 h2 - 1,
      toConn hG hne hdeg v = Sum.inr (Sum.inl ⟨treeDepth hG h1 v - 1, hlt⟩) := by
  refine ⟨?_, ?_⟩
  · have hlt : treeDepth hG h1 v < treeDepth hG h1 h2 := depth1_lt_of_not_right hne hdeg hnr hc
    have h3' : 0 < treeDepth hG h1 v := treeDepth_pos (hG := hG) (h := h1) hv1
    omega
  · simp only [toConn, dif_neg hv1, dif_neg hv2, dif_neg hnr, dif_pos hc]

theorem toConn_left {v : V} (hv1 : v ≠ h1) (hc : treeArm hG h1 v ≠ treeArm hG h1 h2) :
    ∃ hlt : treeDepth hG h1 v - 1 < leftLen hG h1 h2 (leftIdx hv1 hc),
      toConn hG hne hdeg v = Sum.inl (some ⟨leftIdx hv1 hc, ⟨treeDepth hG h1 v - 1, hlt⟩⟩) := by
  have hv2 : v ≠ h2 := by
    intro hcc
    subst hcc
    exact hc rfl
  have hnr : ¬ IsRightSide hG h1 h2 v := not_right_of_arm1_ne hne hc
  refine ⟨?_, ?_⟩
  · have h1' : leftArm hG h1 h2 (leftIdx hv1 hc) = treeArm hG h1 v := leftArm_leftIdx _ _
    have h2' : treeDepth hG h1 v ≤ (armFinset hG h1 (treeArm hG h1 v)).card :=
      treeDepth_le_card_arm (hG := hG) (h := h1)
    have h3' : 0 < treeDepth hG h1 v := treeDepth_pos (hG := hG) (h := h1) hv1
    simp only [leftLen, h1']
    omega
  · simp only [toConn, dif_neg hv1, dif_neg hv2, dif_neg hnr, dif_neg hc]

/-! ### The level of a vertex of the model -/

/-- The distance from the left hub, read off from a vertex of the connector model
(measured from the right hub on the right-hand side). -/
def connLevel {t m n : ℕ} {a : Fin m → ℕ} {b : Fin n → ℕ} (x : ConnV t m a n b) : ℕ :=
  match x with
  | Sum.inl none => 0
  | Sum.inl (some p) => p.2.val + 1
  | Sum.inr (Sum.inl j) => j.val + 1
  | Sum.inr (Sum.inr none) => 0
  | Sum.inr (Sum.inr (some q)) => q.2.val + 1

theorem connLevel_toConn_right {v : V} (hv2 : v ≠ h2) (hr : IsRightSide hG h1 h2 v) :
    connLevel (toConn hG hne hdeg v) = treeDepth hG h2 v := by
  obtain ⟨_, hev⟩ := toConn_right hG hne hdeg hv2 hr
  rw [hev]
  have := treeDepth_pos (hG := hG) (h := h2) hv2
  show treeDepth hG h2 v - 1 + 1 = treeDepth hG h2 v
  omega

theorem connLevel_toConn_interior {v : V} (hv1 : v ≠ h1) (hv2 : v ≠ h2)
    (hnr : ¬ IsRightSide hG h1 h2 v) (hc : treeArm hG h1 v = treeArm hG h1 h2) :
    connLevel (toConn hG hne hdeg v) = treeDepth hG h1 v := by
  obtain ⟨_, hev⟩ := toConn_interior hG hne hdeg hv1 hv2 hnr hc
  rw [hev]
  have := treeDepth_pos (hG := hG) (h := h1) hv1
  show treeDepth hG h1 v - 1 + 1 = treeDepth hG h1 v
  omega

theorem connLevel_toConn_left {v : V} (hv1 : v ≠ h1)
    (hc : treeArm hG h1 v ≠ treeArm hG h1 h2) :
    connLevel (toConn hG hne hdeg v) = treeDepth hG h1 v := by
  obtain ⟨_, hev⟩ := toConn_left hG hne hdeg hv1 hc
  rw [hev]
  have := treeDepth_pos (hG := hG) (h := h1) hv1
  show treeDepth hG h1 v - 1 + 1 = treeDepth hG h1 v
  omega

omit [Fintype V] [DecidableRel G.Adj] in
/-- The five kinds of vertex of a two-hub tree. -/
theorem toConn_classify (v : V) :
    v = h1 ∨ v = h2 ∨ (v ≠ h2 ∧ IsRightSide hG h1 h2 v) ∨
      (v ≠ h1 ∧ v ≠ h2 ∧ ¬ IsRightSide hG h1 h2 v ∧ treeArm hG h1 v = treeArm hG h1 h2) ∨
      (v ≠ h1 ∧ treeArm hG h1 v ≠ treeArm hG h1 h2) := by
  by_cases hv1 : v = h1
  · exact Or.inl hv1
  by_cases hv2 : v = h2
  · exact Or.inr (Or.inl hv2)
  by_cases hr : IsRightSide hG h1 h2 v
  · exact Or.inr (Or.inr (Or.inl ⟨hv2, hr⟩))
  by_cases hc : treeArm hG h1 v = treeArm hG h1 h2
  · exact Or.inr (Or.inr (Or.inr (Or.inl ⟨hv1, hv2, hr, hc⟩)))
  · exact Or.inr (Or.inr (Or.inr (Or.inr ⟨hv1, hc⟩)))

/-! ### The isomorphism -/

theorem toConn_injective : Function.Injective (toConn hG hne hdeg) := by
  intro u v huv
  rcases toConn_classify (hG := hG) (h1 := h1) (h2 := h2) u with
      hu | hu | ⟨hu2, hur⟩ | ⟨hu1, hu2, hunr, huc⟩ | ⟨hu1, huc⟩ <;>
    rcases toConn_classify (hG := hG) (h1 := h1) (h2 := h2) v with
      hv | hv | ⟨hv2, hvr⟩ | ⟨hv1, hv2, hvnr, hvc⟩ | ⟨hv1, hvc⟩
  -- u = h1
  · rw [hu, hv]
  · exfalso
    rw [hu, hv, toConn_h1, toConn_h2] at huv
    exact absurd huv (by simp)
  · exfalso
    obtain ⟨_, hev⟩ := toConn_right hG hne hdeg hv2 hvr
    rw [hu, toConn_h1, hev] at huv
    exact absurd huv (by simp)
  · exfalso
    obtain ⟨_, hev⟩ := toConn_interior hG hne hdeg hv1 hv2 hvnr hvc
    rw [hu, toConn_h1, hev] at huv
    exact absurd huv (by simp)
  · exfalso
    obtain ⟨_, hev⟩ := toConn_left hG hne hdeg hv1 hvc
    rw [hu, toConn_h1, hev] at huv
    exact absurd huv (by simp)
  -- u = h2
  · exfalso
    rw [hu, hv, toConn_h1, toConn_h2] at huv
    exact absurd huv (by simp)
  · rw [hu, hv]
  · exfalso
    obtain ⟨_, hev⟩ := toConn_right hG hne hdeg hv2 hvr
    rw [hu, toConn_h2, hev] at huv
    exact absurd huv (by simp)
  · exfalso
    obtain ⟨_, hev⟩ := toConn_interior hG hne hdeg hv1 hv2 hvnr hvc
    rw [hu, toConn_h2, hev] at huv
    exact absurd huv (by simp)
  · exfalso
    obtain ⟨_, hev⟩ := toConn_left hG hne hdeg hv1 hvc
    rw [hu, toConn_h2, hev] at huv
    exact absurd huv (by simp)
  -- u on the right
  · exfalso
    obtain ⟨_, heu⟩ := toConn_right hG hne hdeg hu2 hur
    rw [hv, toConn_h1, heu] at huv
    exact absurd huv (by simp)
  · exfalso
    obtain ⟨_, heu⟩ := toConn_right hG hne hdeg hu2 hur
    rw [hv, toConn_h2, heu] at huv
    exact absurd huv (by simp)
  · have hlev := congrArg connLevel huv
    rw [connLevel_toConn_right hG hne hdeg hu2 hur,
      connLevel_toConn_right hG hne hdeg hv2 hvr] at hlev
    obtain ⟨_, heu⟩ := toConn_right hG hne hdeg hu2 hur
    obtain ⟨_, hev⟩ := toConn_right hG hne hdeg hv2 hvr
    rw [heu, hev] at huv
    simp only [Sum.inr.injEq, Option.some.injEq, Sigma.mk.injEq] at huv
    have harm : treeArm hG h2 u = treeArm hG h2 v := by
      rw [← rightArm_rightIdx hu2 (arm2_ne_of_right hne hur hu2),
        ← rightArm_rightIdx hv2 (arm2_ne_of_right hne hvr hv2), huv.1]
    have hxne : treeArm hG h2 u ≠ treeArm hG h2 h1 := arm2_ne_of_right hne hur hu2
    exact treeArm_depth_inj (hG := hG) (h := h2) (x := treeArm hG h2 u)
      (fun w hw hwx => deg_le_two_right_arm hdeg hxne hw hwx) (treeDepth hG h2 u) u v hu2 hv2
      rfl harm.symm rfl hlev.symm
  · exfalso
    obtain ⟨_, heu⟩ := toConn_right hG hne hdeg hu2 hur
    obtain ⟨_, hev⟩ := toConn_interior hG hne hdeg hv1 hv2 hvnr hvc
    rw [heu, hev] at huv
    exact absurd huv (by simp)
  · exfalso
    obtain ⟨_, heu⟩ := toConn_right hG hne hdeg hu2 hur
    obtain ⟨_, hev⟩ := toConn_left hG hne hdeg hv1 hvc
    rw [heu, hev] at huv
    exact absurd huv (by simp)
  -- u interior
  · exfalso
    obtain ⟨_, heu⟩ := toConn_interior hG hne hdeg hu1 hu2 hunr huc
    rw [hv, toConn_h1, heu] at huv
    exact absurd huv (by simp)
  · exfalso
    obtain ⟨_, heu⟩ := toConn_interior hG hne hdeg hu1 hu2 hunr huc
    rw [hv, toConn_h2, heu] at huv
    exact absurd huv (by simp)
  · exfalso
    obtain ⟨_, heu⟩ := toConn_interior hG hne hdeg hu1 hu2 hunr huc
    obtain ⟨_, hev⟩ := toConn_right hG hne hdeg hv2 hvr
    rw [heu, hev] at huv
    exact absurd huv (by simp)
  · have hlev := congrArg connLevel huv
    rw [connLevel_toConn_interior hG hne hdeg hu1 hu2 hunr huc,
      connLevel_toConn_interior hG hne hdeg hv1 hv2 hvnr hvc] at hlev
    have hlt : treeDepth hG h1 u < treeDepth hG h1 h2 :=
      depth1_lt_of_not_right hne hdeg hunr huc
    exact uniq_arm_hub hdeg (treeDepth hG h1 u) (le_of_lt hlt) u v hu1 hv1 huc hvc rfl hlev.symm
  · exfalso
    obtain ⟨_, heu⟩ := toConn_interior hG hne hdeg hu1 hu2 hunr huc
    obtain ⟨_, hev⟩ := toConn_left hG hne hdeg hv1 hvc
    rw [heu, hev] at huv
    exact absurd huv (by simp)
  -- u on a left arm
  · exfalso
    obtain ⟨_, heu⟩ := toConn_left hG hne hdeg hu1 huc
    rw [hv, toConn_h1, heu] at huv
    exact absurd huv (by simp)
  · exfalso
    obtain ⟨_, heu⟩ := toConn_left hG hne hdeg hu1 huc
    rw [hv, toConn_h2, heu] at huv
    exact absurd huv (by simp)
  · exfalso
    obtain ⟨_, heu⟩ := toConn_left hG hne hdeg hu1 huc
    obtain ⟨_, hev⟩ := toConn_right hG hne hdeg hv2 hvr
    rw [heu, hev] at huv
    exact absurd huv (by simp)
  · exfalso
    obtain ⟨_, heu⟩ := toConn_left hG hne hdeg hu1 huc
    obtain ⟨_, hev⟩ := toConn_interior hG hne hdeg hv1 hv2 hvnr hvc
    rw [heu, hev] at huv
    exact absurd huv (by simp)
  · have hlev := congrArg connLevel huv
    rw [connLevel_toConn_left hG hne hdeg hu1 huc,
      connLevel_toConn_left hG hne hdeg hv1 hvc] at hlev
    obtain ⟨_, heu⟩ := toConn_left hG hne hdeg hu1 huc
    obtain ⟨_, hev⟩ := toConn_left hG hne hdeg hv1 hvc
    rw [heu, hev] at huv
    simp only [Sum.inl.injEq, Option.some.injEq, Sigma.mk.injEq] at huv
    have harm : treeArm hG h1 u = treeArm hG h1 v := by
      rw [← leftArm_leftIdx hu1 huc, ← leftArm_leftIdx hv1 hvc, huv.1]
    exact treeArm_depth_inj (hG := hG) (h := h1) (x := treeArm hG h1 u)
      (fun w hw hwx => deg_le_two_left_arm hdeg huc hw hwx) (treeDepth hG h1 u) u v hu1 hv1
      rfl harm.symm rfl hlev.symm

/-- Equality of dependent pairs indexing the arms. -/
theorem sigma_fin_ext {N : ℕ} {f : Fin N → ℕ} {i j : Fin N} {x : Fin (f i)} {y : Fin (f j)}
    (hij : i = j) (hxy : x.val = y.val) : (⟨i, x⟩ : Σ i : Fin N, Fin (f i)) = ⟨j, y⟩ := by
  subst hij
  simp only [Sigma.mk.injEq, heq_eq_eq, true_and]
  exact Fin.ext hxy

theorem toConn_surjective : Function.Surjective (toConn hG hne hdeg) := by
  rintro (x | (j | y))
  · match x with
    | none => exact ⟨h1, toConn_h1 hG hne hdeg⟩
    | some ⟨i, k⟩ =>
        have hinj : ∀ u v : V, u ∈ armFinset hG h1 (leftArm hG h1 h2 i) →
            v ∈ armFinset hG h1 (leftArm hG h1 h2 i) →
            treeDepth hG h1 u = treeDepth hG h1 v → u = v := by
          intro u v hu hv hd
          rw [mem_armFinset] at hu hv
          exact treeArm_depth_inj (hG := hG) (h := h1) (x := leftArm hG h1 h2 i)
            (fun w hw hwx => deg_le_two_left_arm hdeg (leftArm_ne i) hw hwx)
            (treeDepth hG h1 u) u v hu.1 hv.1 hu.2 hv.2 rfl hd.symm
        obtain ⟨v, hv, hvd⟩ := exists_mem_armFinset_depth (hG := hG) (h := h1) hinj (k.val + 1)
          (by omega) (by simp [leftLen])
        rw [mem_armFinset] at hv
        have hc : treeArm hG h1 v ≠ treeArm hG h1 h2 := by
          rw [hv.2]; exact leftArm_ne i
        obtain ⟨_, hev⟩ := toConn_left hG hne hdeg hv.1 hc
        have hidx : leftIdx hv.1 hc = i := leftIdx_eq hv.1 hc i hv.2
        refine ⟨v, ?_⟩
        rw [hev]
        exact congrArg (fun z => Sum.inl (some z))
          (sigma_fin_ext hidx (show treeDepth hG h1 v - 1 = (k : ℕ) by omega))
  · obtain ⟨w, hw, hwd⟩ := exists_ancestor (hG := hG) (h := h1) h2 (j.val + 1)
      (by have := j.isLt; omega)
    have hw1 : w ≠ h1 := by
      intro hc
      rw [hc, treeDepth_root (hG := hG)] at hwd
      omega
    have hc : treeArm hG h1 w = treeArm hG h1 h2 :=
      (treeArm_of_mem_support (hG := hG) (h := h1) hw hw1).symm
    have hw2 : w ≠ h2 := by
      intro hcc
      rw [hcc] at hwd
      have := j.isLt
      omega
    have hnr : ¬ IsRightSide hG h1 h2 w := by
      refine not_right_of_depth_lt ?_
      have := j.isLt
      omega
    obtain ⟨_, hew⟩ := toConn_interior hG hne hdeg hw1 hw2 hnr hc
    refine ⟨w, ?_⟩
    rw [hew]
    simp only [Sum.inr.injEq, Sum.inl.injEq, Fin.ext_iff]
    omega
  · match y with
    | none => exact ⟨h2, toConn_h2 hG hne hdeg⟩
    | some ⟨i, k⟩ =>
        have hinj : ∀ u v : V, u ∈ armFinset hG h2 (rightArm hG h1 h2 i) →
            v ∈ armFinset hG h2 (rightArm hG h1 h2 i) →
            treeDepth hG h2 u = treeDepth hG h2 v → u = v := by
          intro u v hu hv hd
          rw [mem_armFinset] at hu hv
          exact treeArm_depth_inj (hG := hG) (h := h2) (x := rightArm hG h1 h2 i)
            (fun w hw hwx => deg_le_two_right_arm hdeg (rightArm_ne i) hw hwx)
            (treeDepth hG h2 u) u v hu.1 hv.1 hu.2 hv.2 rfl hd.symm
        obtain ⟨v, hv, hvd⟩ := exists_mem_armFinset_depth (hG := hG) (h := h2) hinj (k.val + 1)
          (by omega) (by simp [rightLen])
        rw [mem_armFinset] at hv
        have hc : treeArm hG h2 v ≠ treeArm hG h2 h1 := by
          rw [hv.2]; exact rightArm_ne i
        have hr : IsRightSide hG h1 h2 v := right_of_arm2_ne hc
        obtain ⟨_, hev⟩ := toConn_right hG hne hdeg hv.1 hr
        have hidx : rightIdx hv.1 (arm2_ne_of_right hne hr hv.1) = i :=
          rightIdx_eq hv.1 _ i hv.2
        refine ⟨v, ?_⟩
        rw [hev]
        exact congrArg (fun z => Sum.inr (Sum.inr (some z)))
          (sigma_fin_ext hidx (show treeDepth hG h2 v - 1 = (k : ℕ) by omega))

/-! ### Adjacency, case by case -/

omit [Fintype V] [DecidableEq V] [DecidableRel G.Adj] in
theorem adj_iff_symm {W W' : Type*} {A : SimpleGraph W} {B : SimpleGraph W'} {x y : W}
    {p q : W'} (h : A.Adj x y ↔ B.Adj p q) : A.Adj y x ↔ B.Adj q p :=
  ⟨fun hh => (h.mp hh.symm).symm, fun hh => (h.mpr hh.symm).symm⟩

theorem conn_adj_h1_h2 :
    (connModel hG (h1 := h1) (h2 := h2)).Adj (toConn hG hne hdeg h1) (toConn hG hne hdeg h2)
      ↔ G.Adj h1 h2 := by
  rw [toConn_h1, toConn_h2, connGraph_adj_lr, adj_root_iff (hG := hG) (h := h1) (Ne.symm hne)]
  have hpos := depth_hub_pos (hG := hG) (h1 := h1) (h2 := h2) hne
  constructor
  · rintro ⟨-, -, hz⟩
    omega
  · intro hd
    exact ⟨rfl, rfl, by omega⟩

theorem conn_adj_h1_right {v : V} (hv2 : v ≠ h2) (hr : IsRightSide hG h1 h2 v) :
    (connModel hG (h1 := h1) (h2 := h2)).Adj (toConn hG hne hdeg h1) (toConn hG hne hdeg v)
      ↔ G.Adj h1 v := by
  obtain ⟨_, hev⟩ := toConn_right hG hne hdeg hv2 hr
  rw [toConn_h1, hev, connGraph_adj_lr]
  have hv1 : v ≠ h1 := ne_h1_of_right hne hr
  constructor
  · rintro ⟨-, hq, -⟩
    exact absurd hq (by simp)
  · intro hadj
    exfalso
    have e1 := (adj_root_iff (hG := hG) (h := h1) hv1).mp hadj
    have e2 := depth1_of_right hr
    have e3 := depth_hub_pos (hG := hG) (h1 := h1) (h2 := h2) hne
    have e4 := treeDepth_pos (hG := hG) (h := h2) hv2
    omega

theorem conn_adj_h1_interior {w : V} (hw1 : w ≠ h1) (hw2 : w ≠ h2)
    (hnr : ¬ IsRightSide hG h1 h2 w) (hc : treeArm hG h1 w = treeArm hG h1 h2) :
    (connModel hG (h1 := h1) (h2 := h2)).Adj (toConn hG hne hdeg h1) (toConn hG hne hdeg w)
      ↔ G.Adj h1 w := by
  obtain ⟨_, hew⟩ := toConn_interior hG hne hdeg hw1 hw2 hnr hc
  rw [toConn_h1, hew, connGraph_adj_lc, adj_root_iff (hG := hG) (h := h1) hw1]
  have hpos := treeDepth_pos (hG := hG) (h := h1) hw1
  constructor
  · rintro ⟨-, hz⟩
    have hz' : treeDepth hG h1 w - 1 = 0 := hz
    omega
  · intro hd
    exact ⟨rfl, show treeDepth hG h1 w - 1 = 0 by omega⟩

theorem conn_adj_h1_left {v : V} (hv1 : v ≠ h1) (hvc : treeArm hG h1 v ≠ treeArm hG h1 h2) :
    (connModel hG (h1 := h1) (h2 := h2)).Adj (toConn hG hne hdeg h1) (toConn hG hne hdeg v)
      ↔ G.Adj h1 v := by
  obtain ⟨_, hev⟩ := toConn_left hG hne hdeg hv1 hvc
  rw [toConn_h1, hev, connGraph_adj_ll, adj_root_iff (hG := hG) (h := h1) hv1]
  have hpos := treeDepth_pos (hG := hG) (h := h1) hv1
  constructor
  · intro hh
    have hh' : treeDepth hG h1 v - 1 = 0 := hh
    omega
  · intro hd
    show treeDepth hG h1 v - 1 = 0
    omega

theorem conn_adj_h2_right {v : V} (hv2 : v ≠ h2) (hr : IsRightSide hG h1 h2 v) :
    (connModel hG (h1 := h1) (h2 := h2)).Adj (toConn hG hne hdeg h2) (toConn hG hne hdeg v)
      ↔ G.Adj h2 v := by
  obtain ⟨_, hev⟩ := toConn_right hG hne hdeg hv2 hr
  rw [toConn_h2, hev, connGraph_adj_rr, adj_root_iff (hG := hG) (h := h2) hv2]
  have hpos := treeDepth_pos (hG := hG) (h := h2) hv2
  constructor
  · intro hh
    have hh' : treeDepth hG h2 v - 1 = 0 := hh
    omega
  · intro hd
    show treeDepth hG h2 v - 1 = 0
    omega

theorem conn_adj_interior_h2 {w : V} (hw1 : w ≠ h1) (hw2 : w ≠ h2)
    (hnr : ¬ IsRightSide hG h1 h2 w) (hc : treeArm hG h1 w = treeArm hG h1 h2) :
    (connModel hG (h1 := h1) (h2 := h2)).Adj (toConn hG hne hdeg w) (toConn hG hne hdeg h2)
      ↔ G.Adj w h2 := by
  obtain ⟨_, hew⟩ := toConn_interior hG hne hdeg hw1 hw2 hnr hc
  rw [hew, toConn_h2, connGraph_adj_cr]
  have hlt : treeDepth hG h1 w < treeDepth hG h1 h2 := depth1_lt_of_not_right hne hdeg hnr hc
  have hpos := treeDepth_pos (hG := hG) (h := h1) hw1
  constructor
  · rintro ⟨-, hval⟩
    have hval' : treeDepth hG h1 w - 1 + 1 = treeDepth hG h1 h2 - 1 := hval
    refine adj_of_depth_succ_of_uniq (hG := hG) (h := h1) hw1 (Ne.symm hne) hc (by omega) ?_
    intro z hz1 hzarm hzd
    exact uniq_arm_hub hdeg (treeDepth hG h1 w) (le_of_lt hlt) z w hz1 hw1
      (by rw [hzarm, hc]) hc hzd rfl
  · intro hadj
    refine ⟨rfl, ?_⟩
    show treeDepth hG h1 w - 1 + 1 = treeDepth hG h1 h2 - 1
    rcases treeDepth_adj (hG := hG) (h := h1) hadj with e | e <;> omega

theorem conn_adj_left_h2 {u : V} (hu1 : u ≠ h1) (huc : treeArm hG h1 u ≠ treeArm hG h1 h2) :
    (connModel hG (h1 := h1) (h2 := h2)).Adj (toConn hG hne hdeg u) (toConn hG hne hdeg h2)
      ↔ G.Adj u h2 := by
  obtain ⟨_, heu⟩ := toConn_left hG hne hdeg hu1 huc
  rw [heu, toConn_h2, connGraph_adj_lr]
  constructor
  · rintro ⟨hp, -, -⟩
    exact absurd hp (by simp)
  · intro hadj
    exact absurd (treeArm_eq_of_adj (hG := hG) (h := h1) hu1 (Ne.symm hne) hadj) huc

theorem conn_adj_right_right {u v : V} (hu2 : u ≠ h2) (hur : IsRightSide hG h1 h2 u)
    (hv2 : v ≠ h2) (hvr : IsRightSide hG h1 h2 v) :
    (connModel hG (h1 := h1) (h2 := h2)).Adj (toConn hG hne hdeg u) (toConn hG hne hdeg v)
      ↔ G.Adj u v := by
  obtain ⟨_, heu⟩ := toConn_right hG hne hdeg hu2 hur
  obtain ⟨_, hev⟩ := toConn_right hG hne hdeg hv2 hvr
  rw [heu, hev, connGraph_adj_rr]
  have hpu := treeDepth_pos (hG := hG) (h := h2) hu2
  have hpv := treeDepth_pos (hG := hG) (h := h2) hv2
  constructor
  · rintro ⟨hidx, hd⟩
    have hidx' : rightIdx hu2 (arm2_ne_of_right hne hur hu2)
        = rightIdx hv2 (arm2_ne_of_right hne hvr hv2) := hidx
    have harm : treeArm hG h2 u = treeArm hG h2 v := by
      rw [← rightArm_rightIdx hu2 (arm2_ne_of_right hne hur hu2),
        ← rightArm_rightIdx hv2 (arm2_ne_of_right hne hvr hv2), hidx']
    have hd' : treeDepth hG h2 u - 1 + 1 = treeDepth hG h2 v - 1 ∨
        treeDepth hG h2 v - 1 + 1 = treeDepth hG h2 u - 1 := hd
    rcases hd' with e | e
    · refine adj_of_depth_succ_of_uniq (hG := hG) (h := h2) hu2 hv2 harm (by omega) ?_
      intro z hz hzarm hzd
      exact treeArm_depth_inj (hG := hG) (h := h2) (x := treeArm hG h2 u)
        (fun w hw hwx => deg_le_two_right_arm hdeg (arm2_ne_of_right hne hur hu2) hw hwx)
        (treeDepth hG h2 u) z u hz hu2 hzarm rfl hzd rfl
    · refine (adj_of_depth_succ_of_uniq (hG := hG) (h := h2) hv2 hu2 harm.symm (by omega) ?_).symm
      intro z hz hzarm hzd
      exact treeArm_depth_inj (hG := hG) (h := h2) (x := treeArm hG h2 v)
        (fun w hw hwx => deg_le_two_right_arm hdeg (arm2_ne_of_right hne hvr hv2) hw hwx)
        (treeDepth hG h2 v) z v hz hv2 hzarm rfl hzd rfl
  · intro hadj
    have harm : treeArm hG h2 u = treeArm hG h2 v :=
      treeArm_eq_of_adj (hG := hG) (h := h2) hu2 hv2 hadj
    have hidx : rightIdx hu2 (arm2_ne_of_right hne hur hu2)
        = rightIdx hv2 (arm2_ne_of_right hne hvr hv2) := by
      refine rightIdx_eq hu2 (arm2_ne_of_right hne hur hu2) _ ?_
      rw [harm, rightArm_rightIdx hv2 (arm2_ne_of_right hne hvr hv2)]
    refine ⟨hidx, ?_⟩
    show treeDepth hG h2 u - 1 + 1 = treeDepth hG h2 v - 1 ∨
      treeDepth hG h2 v - 1 + 1 = treeDepth hG h2 u - 1
    rcases treeDepth_adj (hG := hG) (h := h2) hadj with e | e
    · left; omega
    · right; omega

theorem conn_adj_interior_right {u v : V} (hu1 : u ≠ h1) (hu2 : u ≠ h2)
    (hunr : ¬ IsRightSide hG h1 h2 u) (huc : treeArm hG h1 u = treeArm hG h1 h2)
    (hv2 : v ≠ h2) (hvr : IsRightSide hG h1 h2 v) :
    (connModel hG (h1 := h1) (h2 := h2)).Adj (toConn hG hne hdeg u) (toConn hG hne hdeg v)
      ↔ G.Adj u v := by
  obtain ⟨_, heu⟩ := toConn_interior hG hne hdeg hu1 hu2 hunr huc
  obtain ⟨_, hev⟩ := toConn_right hG hne hdeg hv2 hvr
  rw [heu, hev, connGraph_adj_cr]
  constructor
  · rintro ⟨hq, -⟩
    exact absurd hq (by simp)
  · intro hadj
    exfalso
    have e1 := depth1_of_right hvr
    have e2 : treeDepth hG h1 u < treeDepth hG h1 h2 := depth1_lt_of_not_right hne hdeg hunr huc
    have e3 := treeDepth_pos (hG := hG) (h := h2) hv2
    rcases treeDepth_adj (hG := hG) (h := h1) hadj with e | e <;> omega

theorem conn_adj_left_right {u v : V} (hu1 : u ≠ h1) (huc : treeArm hG h1 u ≠ treeArm hG h1 h2)
    (hv2 : v ≠ h2) (hvr : IsRightSide hG h1 h2 v) :
    (connModel hG (h1 := h1) (h2 := h2)).Adj (toConn hG hne hdeg u) (toConn hG hne hdeg v)
      ↔ G.Adj u v := by
  obtain ⟨_, heu⟩ := toConn_left hG hne hdeg hu1 huc
  obtain ⟨_, hev⟩ := toConn_right hG hne hdeg hv2 hvr
  rw [heu, hev, connGraph_adj_lr]
  constructor
  · rintro ⟨hp, -, -⟩
    exact absurd hp (by simp)
  · intro hadj
    exfalso
    have hv1 : v ≠ h1 := ne_h1_of_right hne hvr
    have harm := treeArm_eq_of_adj (hG := hG) (h := h1) hu1 hv1 hadj
    exact huc (by rw [harm, arm1_of_right hne hvr])

theorem conn_adj_interior_interior {u v : V} (hu1 : u ≠ h1) (hu2 : u ≠ h2)
    (hunr : ¬ IsRightSide hG h1 h2 u) (huc : treeArm hG h1 u = treeArm hG h1 h2)
    (hv1 : v ≠ h1) (hv2 : v ≠ h2) (hvnr : ¬ IsRightSide hG h1 h2 v)
    (hvc : treeArm hG h1 v = treeArm hG h1 h2) :
    (connModel hG (h1 := h1) (h2 := h2)).Adj (toConn hG hne hdeg u) (toConn hG hne hdeg v)
      ↔ G.Adj u v := by
  obtain ⟨_, heu⟩ := toConn_interior hG hne hdeg hu1 hu2 hunr huc
  obtain ⟨_, hev⟩ := toConn_interior hG hne hdeg hv1 hv2 hvnr hvc
  rw [heu, hev, connGraph_adj_cc]
  have hpu := treeDepth_pos (hG := hG) (h := h1) hu1
  have hpv := treeDepth_pos (hG := hG) (h := h1) hv1
  have hltu : treeDepth hG h1 u < treeDepth hG h1 h2 := depth1_lt_of_not_right hne hdeg hunr huc
  have hltv : treeDepth hG h1 v < treeDepth hG h1 h2 := depth1_lt_of_not_right hne hdeg hvnr hvc
  constructor
  · intro hh
    have hh' : treeDepth hG h1 u - 1 + 1 = treeDepth hG h1 v - 1 ∨
        treeDepth hG h1 v - 1 + 1 = treeDepth hG h1 u - 1 := hh
    rcases hh' with e | e
    · refine adj_of_depth_succ_of_uniq (hG := hG) (h := h1) hu1 hv1 (by rw [huc, hvc])
        (by omega) ?_
      intro z hz hzarm hzd
      exact uniq_arm_hub hdeg (treeDepth hG h1 u) (le_of_lt hltu) z u hz hu1
        (by rw [hzarm, huc]) huc hzd rfl
    · refine (adj_of_depth_succ_of_uniq (hG := hG) (h := h1) hv1 hu1 (by rw [hvc, huc])
        (by omega) ?_).symm
      intro z hz hzarm hzd
      exact uniq_arm_hub hdeg (treeDepth hG h1 v) (le_of_lt hltv) z v hz hv1
        (by rw [hzarm, hvc]) hvc hzd rfl
  · intro hadj
    show treeDepth hG h1 u - 1 + 1 = treeDepth hG h1 v - 1 ∨
      treeDepth hG h1 v - 1 + 1 = treeDepth hG h1 u - 1
    rcases treeDepth_adj (hG := hG) (h := h1) hadj with e | e
    · left; omega
    · right; omega

theorem conn_adj_left_interior {u v : V} (hu1 : u ≠ h1)
    (huc : treeArm hG h1 u ≠ treeArm hG h1 h2) (hv1 : v ≠ h1) (hv2 : v ≠ h2)
    (hvnr : ¬ IsRightSide hG h1 h2 v) (hvc : treeArm hG h1 v = treeArm hG h1 h2) :
    (connModel hG (h1 := h1) (h2 := h2)).Adj (toConn hG hne hdeg u) (toConn hG hne hdeg v)
      ↔ G.Adj u v := by
  obtain ⟨_, heu⟩ := toConn_left hG hne hdeg hu1 huc
  obtain ⟨_, hev⟩ := toConn_interior hG hne hdeg hv1 hv2 hvnr hvc
  rw [heu, hev, connGraph_adj_lc]
  constructor
  · rintro ⟨hp, -⟩
    exact absurd hp (by simp)
  · intro hadj
    exfalso
    have harm := treeArm_eq_of_adj (hG := hG) (h := h1) hu1 hv1 hadj
    exact huc (by rw [harm, hvc])

theorem conn_adj_left_left {u v : V} (hu1 : u ≠ h1) (huc : treeArm hG h1 u ≠ treeArm hG h1 h2)
    (hv1 : v ≠ h1) (hvc : treeArm hG h1 v ≠ treeArm hG h1 h2) :
    (connModel hG (h1 := h1) (h2 := h2)).Adj (toConn hG hne hdeg u) (toConn hG hne hdeg v)
      ↔ G.Adj u v := by
  obtain ⟨_, heu⟩ := toConn_left hG hne hdeg hu1 huc
  obtain ⟨_, hev⟩ := toConn_left hG hne hdeg hv1 hvc
  rw [heu, hev, connGraph_adj_ll]
  have hpu := treeDepth_pos (hG := hG) (h := h1) hu1
  have hpv := treeDepth_pos (hG := hG) (h := h1) hv1
  constructor
  · rintro ⟨hidx, hd⟩
    have hidx' : leftIdx hu1 huc = leftIdx hv1 hvc := hidx
    have harm : treeArm hG h1 u = treeArm hG h1 v := by
      rw [← leftArm_leftIdx hu1 huc, ← leftArm_leftIdx hv1 hvc, hidx']
    have hd' : treeDepth hG h1 u - 1 + 1 = treeDepth hG h1 v - 1 ∨
        treeDepth hG h1 v - 1 + 1 = treeDepth hG h1 u - 1 := hd
    rcases hd' with e | e
    · refine adj_of_depth_succ_of_uniq (hG := hG) (h := h1) hu1 hv1 harm (by omega) ?_
      intro z hz hzarm hzd
      exact treeArm_depth_inj (hG := hG) (h := h1) (x := treeArm hG h1 u)
        (fun w hw hwx => deg_le_two_left_arm hdeg huc hw hwx)
        (treeDepth hG h1 u) z u hz hu1 hzarm rfl hzd rfl
    · refine (adj_of_depth_succ_of_uniq (hG := hG) (h := h1) hv1 hu1 harm.symm (by omega) ?_).symm
      intro z hz hzarm hzd
      exact treeArm_depth_inj (hG := hG) (h := h1) (x := treeArm hG h1 v)
        (fun w hw hwx => deg_le_two_left_arm hdeg hvc hw hwx)
        (treeDepth hG h1 v) z v hz hv1 hzarm rfl hzd rfl
  · intro hadj
    have harm : treeArm hG h1 u = treeArm hG h1 v :=
      treeArm_eq_of_adj (hG := hG) (h := h1) hu1 hv1 hadj
    have hidx : leftIdx hu1 huc = leftIdx hv1 hvc := by
      refine leftIdx_eq hu1 huc _ ?_
      rw [harm, leftArm_leftIdx hv1 hvc]
    refine ⟨hidx, ?_⟩
    show treeDepth hG h1 u - 1 + 1 = treeDepth hG h1 v - 1 ∨
      treeDepth hG h1 v - 1 + 1 = treeDepth hG h1 u - 1
    rcases treeDepth_adj (hG := hG) (h := h1) hadj with e | e
    · left; omega
    · right; omega

theorem toConn_adj (u v : V) :
    (connModel (hG := hG) (h1 := h1) (h2 := h2)).Adj (toConn hG hne hdeg u) (toConn hG hne hdeg v)
      ↔ G.Adj u v := by
  rcases toConn_classify (hG := hG) (h1 := h1) (h2 := h2) u with
      hu | hu | ⟨hu2, hur⟩ | ⟨hu1, hu2, hunr, huc⟩ | ⟨hu1, huc⟩ <;>
    rcases toConn_classify (hG := hG) (h1 := h1) (h2 := h2) v with
      hv | hv | ⟨hv2, hvr⟩ | ⟨hv1, hv2, hvnr, hvc⟩ | ⟨hv1, hvc⟩
  -- `u = h1`
  · subst hu; subst hv; simp
  · subst hu; subst hv; exact conn_adj_h1_h2 hG hne hdeg
  · subst hu; exact conn_adj_h1_right hG hne hdeg hv2 hvr
  · subst hu; exact conn_adj_h1_interior hG hne hdeg hv1 hv2 hvnr hvc
  · subst hu; exact conn_adj_h1_left hG hne hdeg hv1 hvc
  -- `u = h2`
  · subst hu; subst hv; exact adj_iff_symm (conn_adj_h1_h2 hG hne hdeg)
  · subst hu; subst hv; simp
  · subst hu; exact conn_adj_h2_right hG hne hdeg hv2 hvr
  · subst hu; exact adj_iff_symm (conn_adj_interior_h2 hG hne hdeg hv1 hv2 hvnr hvc)
  · subst hu; exact adj_iff_symm (conn_adj_left_h2 hG hne hdeg hv1 hvc)
  -- `u` on the right side
  · subst hv; exact adj_iff_symm (conn_adj_h1_right hG hne hdeg hu2 hur)
  · subst hv; exact adj_iff_symm (conn_adj_h2_right hG hne hdeg hu2 hur)
  · exact conn_adj_right_right hG hne hdeg hu2 hur hv2 hvr
  · exact adj_iff_symm (conn_adj_interior_right hG hne hdeg hv1 hv2 hvnr hvc hu2 hur)
  · exact adj_iff_symm (conn_adj_left_right hG hne hdeg hv1 hvc hu2 hur)
  -- `u` in the connecting path
  · subst hv; exact adj_iff_symm (conn_adj_h1_interior hG hne hdeg hu1 hu2 hunr huc)
  · subst hv; exact conn_adj_interior_h2 hG hne hdeg hu1 hu2 hunr huc
  · exact conn_adj_interior_right hG hne hdeg hu1 hu2 hunr huc hv2 hvr
  · exact conn_adj_interior_interior hG hne hdeg hu1 hu2 hunr huc hv1 hv2 hvnr hvc
  · exact adj_iff_symm (conn_adj_left_interior hG hne hdeg hv1 hvc hu1 hu2 hunr huc)
  -- `u` on the left side
  · subst hv; exact adj_iff_symm (conn_adj_h1_left hG hne hdeg hu1 huc)
  · subst hv; exact conn_adj_left_h2 hG hne hdeg hu1 huc
  · exact conn_adj_left_right hG hne hdeg hu1 huc hv2 hvr
  · exact conn_adj_left_interior hG hne hdeg hu1 huc hv1 hv2 hvnr hvc
  · exact conn_adj_left_left hG hne hdeg hu1 huc hv1 hvc

/-- **Connector recognition.**  A finite tree with two distinct vertices `h1`, `h2` such
that every other vertex has degree at most two is isomorphic to the connector model. -/
noncomputable def connIsoOfTree : G ≃g connModel (hG := hG) (h1 := h1) (h2 := h2) where
  toEquiv := Equiv.ofBijective _ ⟨toConn_injective hG hne hdeg, toConn_surjective hG hne hdeg⟩
  map_rel_iff' := by
    intro u v
    exact toConn_adj hG hne hdeg u v

end ClanAudit
