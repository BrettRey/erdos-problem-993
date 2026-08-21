import RequestProject.ClanAudit.TreeRoot

/-!
# Arms of a rooted tree

Continuing `TreeRoot`, this file proves the two facts that make an arm of a rooted tree
behave like a path:

* `treeArm_depth_inj` — inside an arm all of whose vertices have degree at most two, a
  vertex is determined by its depth;
* `treeDepth_le_card_arm` — the depth of a vertex is at most the number of vertices of its
  arm.

Together they identify an arm with an initial segment of the natural numbers.
-/

namespace ClanAudit

open SimpleGraph Finset

variable {V : Type*} [Fintype V] [DecidableEq V] {G : SimpleGraph V} [DecidableRel G.Adj]
variable {hG : G.IsTree} {h : V}

theorem three_le_degree_of_three_adj {p u v w : V} (hpu : G.Adj p u) (hpv : G.Adj p v)
    (hpw : G.Adj p w) (huv : u ≠ v) (huw : u ≠ w) (hvw : v ≠ w) : 3 ≤ G.degree p := by
  have hsub : ({u, v, w} : Finset V) ⊆ G.neighborFinset p := by
    intro z hz
    simp only [Finset.mem_insert, Finset.mem_singleton] at hz
    rcases hz with rfl | rfl | rfl <;> simp [mem_neighborFinset, hpu, hpv, hpw]
  have := Finset.card_le_card hsub
  rwa [Finset.card_insert_of_notMem (by simp [huv, huw]),
    Finset.card_insert_of_notMem (by simp [hvw]), Finset.card_singleton] at this

/-- Inside an arm whose vertices all have degree at most two, the depth determines the
vertex. -/
theorem treeArm_depth_inj {x : V}
    (hdeg : ∀ w : V, w ≠ h → treeArm hG h w = x → G.degree w ≤ 2) :
    ∀ k : ℕ, ∀ u v : V, u ≠ h → v ≠ h → treeArm hG h u = x → treeArm hG h v = x →
      treeDepth hG h u = k → treeDepth hG h v = k → u = v := by
  intro k
  induction k using Nat.strong_induction_on with
  | _ k IH =>
    intro u v hu hv hau hav hdu hdv
    rcases Nat.lt_or_ge k 2 with hk | hk
    · interval_cases k
      · exact absurd hdu (by have := treeDepth_pos (hG := hG) (h := h) hu; omega)
      · rw [← treeArm_of_depth_one (hG := hG) (h := h) hdu,
          ← treeArm_of_depth_one (hG := hG) (h := h) hdv, hau, hav]
    · set p := treeParent hG h u with hp
      set q := treeParent hG h v with hq
      have hdp : treeDepth hG h p + 1 = k := by
        rw [hp, treeDepth_parent (hG := hG) (h := h) hu, hdu]
      have hdq : treeDepth hG h q + 1 = k := by
        rw [hq, treeDepth_parent (hG := hG) (h := h) hv, hdv]
      have hph : p ≠ h := by
        intro hc
        rw [hc, treeDepth_root (hG := hG)] at hdp
        omega
      have hqh : q ≠ h := by
        intro hc
        rw [hc, treeDepth_root (hG := hG)] at hdq
        omega
      have hap : treeArm hG h p = x := by
        rw [hp, treeArm_parent (hG := hG) (h := h) (by omega), hau]
      have haq : treeArm hG h q = x := by
        rw [hq, treeArm_parent (hG := hG) (h := h) (by omega), hav]
      have hpq : p = q :=
        IH (k - 1) (by omega) p q hph hqh hap haq (by omega) (by omega)
      by_contra hne
      -- `p` has three distinct neighbours: its own parent, `u` and `v`.
      have hgp : treeDepth hG h (treeParent hG h p) + 1 = treeDepth hG h p :=
        treeDepth_parent (hG := hG) (h := h) hph
      have hnu : treeParent hG h p ≠ u := by
        intro hc; rw [hc] at hgp; omega
      have hnv : treeParent hG h p ≠ v := by
        intro hc; rw [hc] at hgp; omega
      have h3 : 3 ≤ G.degree p :=
        three_le_degree_of_three_adj
          (treeParent_adj (hG := hG) (h := h) hph).symm
          (treeParent_adj (hG := hG) (h := h) hu)
          (by rw [hpq]; exact treeParent_adj (hG := hG) (h := h) hv)
          hnu hnv hne
      exact absurd (hdeg p hph hap) (by omega)

/-- The set of vertices of the arm `x`. -/
noncomputable def armFinset (hG : G.IsTree) (h x : V) : Finset V :=
  Finset.univ.filter (fun v => v ≠ h ∧ treeArm hG h v = x)

omit [DecidableRel G.Adj] in
theorem mem_armFinset {x v : V} :
    v ∈ armFinset hG h x ↔ v ≠ h ∧ treeArm hG h v = x := by
  simp [armFinset]

omit [DecidableRel G.Adj] in
/-- Every non-root vertex on the root path of `v` belongs to the arm of `v`. -/
theorem support_sub_armFinset {v : V} :
    ((treePath hG h v).support.toFinset.erase h) ⊆ armFinset hG h (treeArm hG h v) := by
  intro w hw
  rw [Finset.mem_erase, List.mem_toFinset] at hw
  exact mem_armFinset.mpr ⟨hw.1, (treeArm_of_mem_support (hG := hG) (h := h) hw.2 hw.1).symm⟩

omit [DecidableRel G.Adj] in
/-- The depth of a vertex is at most the size of its arm. -/
theorem treeDepth_le_card_arm {v : V} :
    treeDepth hG h v ≤ (armFinset hG h (treeArm hG h v)).card := by
  have hnd : (treePath hG h v).support.Nodup := (treePath_isPath (hG := hG) (h := h) v).support_nodup
  have hcard : ((treePath hG h v).support.toFinset).card = treeDepth hG h v + 1 := by
    rw [List.toFinset_card_of_nodup hnd, Walk.length_support]
    rfl
  have hmem : h ∈ (treePath hG h v).support.toFinset := by
    simp [List.mem_toFinset]
  have herase : ((treePath hG h v).support.toFinset.erase h).card = treeDepth hG h v := by
    rw [Finset.card_erase_of_mem hmem, hcard]
    omega
  calc treeDepth hG h v = ((treePath hG h v).support.toFinset.erase h).card := herase.symm
    _ ≤ (armFinset hG h (treeArm hG h v)).card :=
        Finset.card_le_card (support_sub_armFinset (hG := hG) (h := h) (v := v))

end ClanAudit
