import Mathlib.Combinatorics.SimpleGraph.Acyclic
import Mathlib.Combinatorics.SimpleGraph.Metric
import Mathlib.Tactic

/-!
# Rooted-tree toolkit

For a tree `G` and a chosen root `h`, every vertex `v` is joined to `h` by a unique path
`treePath hG h v`.  This file develops the basic API of that path: the depth of a vertex,
its parent, the arm (the depth-one vertex on its path to the root), and the exact relation
between adjacency and depth.

These are the tools used to recognise abstract finite trees as the explicit spider and
connector models of this development.
-/

namespace ClanAudit

open SimpleGraph

variable {V : Type*} [DecidableEq V] {G : SimpleGraph V}

/-- The unique path from the root `h` to `v` in a tree. -/
noncomputable def treePath (hG : G.IsTree) (h v : V) : G.Walk h v :=
  (hG.existsUnique_path h v).choose

variable {hG : G.IsTree} {h : V}

omit [DecidableEq V] in
theorem treePath_isPath (v : V) : (treePath hG h v).IsPath :=
  (hG.existsUnique_path h v).choose_spec.1

omit [DecidableEq V] in
theorem treePath_eq {v : V} {p : G.Walk h v} (hp : p.IsPath) : p = treePath hG h v :=
  (hG.existsUnique_path h v).choose_spec.2 p hp

/-- The depth of `v`: the length of the unique path from the root to `v`. -/
noncomputable def treeDepth (hG : G.IsTree) (h v : V) : ℕ := (treePath hG h v).length

omit [DecidableEq V] in
theorem treePath_root : treePath hG h h = Walk.nil := (treePath_eq Walk.IsPath.nil).symm

omit [DecidableEq V] in
theorem treeDepth_root : treeDepth hG h h = 0 := by
  rw [treeDepth, treePath_root]; rfl

omit [DecidableEq V] in
theorem treeDepth_eq_zero_iff {v : V} : treeDepth hG h v = 0 ↔ v = h := by
  constructor
  · intro hv
    exact (Walk.eq_of_length_eq_zero (p := treePath hG h v) hv).symm
  · rintro rfl; exact treeDepth_root

omit [DecidableEq V] in
theorem treeDepth_pos {v : V} (hv : v ≠ h) : 0 < treeDepth hG h v := by
  rcases Nat.eq_zero_or_pos (treeDepth hG h v) with hz | hp
  · exact absurd (treeDepth_eq_zero_iff.mp hz) hv
  · exact hp

/-- Any vertex on the unique root path of `v` has the corresponding prefix as its own
root path. -/
theorem treePath_takeUntil {v w : V} (hw : w ∈ (treePath hG h v).support) :
    (treePath hG h v).takeUntil w hw = treePath hG h w :=
  treePath_eq ((treePath_isPath v).takeUntil hw)

theorem treeDepth_le_of_mem_support {v w : V} (hw : w ∈ (treePath hG h v).support) :
    treeDepth hG h w ≤ treeDepth hG h v := by
  show (treePath hG h w).length ≤ (treePath hG h v).length
  rw [← treePath_takeUntil hw]
  exact Walk.length_takeUntil_le _ _

/-- If `v` lies on the root path of `u` and `u` and `v` are adjacent, the root path of `u`
is the root path of `v` extended by the edge. -/
theorem treeDepth_succ_of_mem_support {u v : V} (huv : G.Adj u v)
    (hv : v ∈ (treePath hG h u).support) :
    treeDepth hG h v + 1 = treeDepth hG h u := by
  have hspec := Walk.take_spec (treePath hG h u) hv
  have hdrop : (treePath hG h u).dropUntil v hv = Walk.cons huv.symm Walk.nil := by
    have h1 : ((treePath hG h u).dropUntil v hv).IsPath := (treePath_isPath u).dropUntil hv
    have h2 : (Walk.cons huv.symm Walk.nil : G.Walk v u).IsPath := by
      simp [Walk.isPath_def, huv.ne']
    have := hG.IsAcyclic.path_unique ⟨_, h1⟩ ⟨_, h2⟩
    exact congrArg Subtype.val this
  have hlen := congrArg Walk.length hspec
  rw [Walk.length_append, hdrop, treePath_takeUntil hv] at hlen
  simpa [treeDepth] using hlen

/-- Adjacency changes the depth by exactly one. -/
theorem treeDepth_adj {u v : V} (huv : G.Adj u v) :
    treeDepth hG h u + 1 = treeDepth hG h v ∨ treeDepth hG h v + 1 = treeDepth hG h u := by
  by_cases hv : v ∈ (treePath hG h u).support
  · exact Or.inr (treeDepth_succ_of_mem_support huv hv)
  · left
    have : (treePath hG h u).concat huv = treePath hG h v :=
      treePath_eq ((treePath_isPath u).concat hv huv)
    rw [treeDepth, treeDepth, ← this, Walk.length_concat]

/-- The root path of a deeper adjacent vertex is obtained by extending. -/
theorem treePath_concat {u v : V} (huv : G.Adj u v)
    (hd : treeDepth hG h u + 1 = treeDepth hG h v) :
    (treePath hG h u).concat huv = treePath hG h v := by
  by_cases hv : v ∈ (treePath hG h u).support
  · exact absurd (treeDepth_succ_of_mem_support huv hv) (by omega)
  · exact treePath_eq ((treePath_isPath u).concat hv huv)

/-! ### Parents -/

/-- The parent of `v`: the vertex just before `v` on its root path. -/
noncomputable def treeParent (hG : G.IsTree) (h v : V) : V :=
  (treePath hG h v).getVert (treeDepth hG h v - 1)

omit [DecidableEq V] in
theorem treeParent_mem_support (v : V) :
    treeParent hG h v ∈ (treePath hG h v).support :=
  Walk.getVert_mem_support _ _

omit [DecidableEq V] in
theorem treeParent_adj {v : V} (hv : v ≠ h) : G.Adj (treeParent hG h v) v := by
  have hpos : 0 < treeDepth hG h v := treeDepth_pos hv
  have := Walk.adj_getVert_succ (treePath hG h v)
    (i := treeDepth hG h v - 1) (by simpa [treeDepth] using Nat.sub_lt hpos Nat.one_pos)
  have he : treeDepth hG h v - 1 + 1 = (treePath hG h v).length := by
    simp only [treeDepth] at hpos ⊢; omega
  rw [he, Walk.getVert_length] at this
  exact this

theorem treeDepth_parent {v : V} (hv : v ≠ h) :
    treeDepth hG h (treeParent hG h v) + 1 = treeDepth hG h v := by
  have hle := treeDepth_le_of_mem_support (hG := hG) (h := h) (treeParent_mem_support (hG := hG)
    (h := h) v)
  rcases treeDepth_adj (treeParent_adj (hG := hG) (h := h) hv) with h1 | h2
  · exact h1
  · omega

/-! ### Arms -/

/-- The arm of `v`: the depth-one vertex on the root path of `v`. -/
noncomputable def treeArm (hG : G.IsTree) (h v : V) : V := (treePath hG h v).getVert 1

omit [DecidableEq V] in
theorem treeArm_of_depth_one {v : V} (hd : treeDepth hG h v = 1) : treeArm hG h v = v := by
  have : (treePath hG h v).length = 1 := hd
  rw [treeArm, ← this, Walk.getVert_length]

omit [DecidableEq V] in
theorem treeArm_adj_root {v : V} (hv : v ≠ h) : G.Adj h (treeArm hG h v) := by
  have hpos : 0 < (treePath hG h v).length := treeDepth_pos hv
  have := Walk.adj_getVert_succ (treePath hG h v) (i := 0) hpos
  rwa [Walk.getVert_zero, zero_add] at this

theorem treeDepth_arm {v : V} (hv : v ≠ h) : treeDepth hG h (treeArm hG h v) = 1 := by
  rcases treeDepth_adj (treeArm_adj_root (hG := hG) (h := h) hv) with h1 | h2
  · rw [treeDepth_root (hG := hG)] at h1; omega
  · rw [treeDepth_root (hG := hG)] at h2; omega

theorem treeArm_ne_root {v : V} (hv : v ≠ h) : treeArm hG h v ≠ h := by
  intro hc
  have := treeDepth_arm (hG := hG) (h := h) hv
  rw [hc, treeDepth_root (hG := hG)] at this
  exact absurd this (by omega)

/-- Every non-root vertex on the root path of `v` has the same arm as `v`. -/
theorem treeArm_of_mem_support {v w : V} (hw : w ∈ (treePath hG h v).support) (hwh : w ≠ h) :
    treeArm hG h v = treeArm hG h w := by
  have hspec := Walk.take_spec (treePath hG h v) hw
  rw [treePath_takeUntil hw] at hspec
  have hpos : 0 < (treePath hG h w).length := treeDepth_pos hwh
  rw [treeArm, ← hspec, Walk.getVert_append]
  by_cases hlt : 1 < (treePath hG h w).length
  · rw [if_pos hlt]; rfl
  · rw [if_neg hlt]
    have h1 : (treePath hG h w).length = 1 := by omega
    rw [h1]
    simp only [Nat.sub_self, Walk.getVert_zero]
    exact (treeArm_of_depth_one (hG := hG) (h := h) h1).symm

theorem treeArm_parent {v : V} (hv : 2 ≤ treeDepth hG h v) :
    treeArm hG h (treeParent hG h v) = treeArm hG h v := by
  have hvh : v ≠ h := by
    intro hc; rw [hc, treeDepth_root] at hv; omega
  have hp : treeParent hG h v ≠ h := by
    intro hc
    have := treeDepth_parent (hG := hG) (h := h) hvh
    rw [hc, treeDepth_root] at this
    omega
  exact (treeArm_of_mem_support (hG := hG) (h := h)
    (treeParent_mem_support (hG := hG) (h := h) v) hp).symm

end ClanAudit
