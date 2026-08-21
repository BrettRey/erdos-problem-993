import RequestProject.ClanAudit.TreeArm

/-!
# Rooted-tree toolkit, second layer

Tools needed to recognise a finite tree with two branch vertices as the explicit
connector model:

* `walk_isPath_append` — concatenating two paths meeting only at the joint is a path;
* `treePath_reverse` — the root path from `v` to `h` is the reverse of the root path from
  `h` to `v`;
* `exists_ancestor` — every depth below the depth of `v` is realized on the root path
  of `v`;
* `treeArm_depth_inj_lt` — the depth-injectivity of `TreeArm`, localized to the vertices
  of an arm of depth below a bound (so that a *second* branch vertex further out along
  the arm is allowed);
* `exists_mem_armFinset_depth` — an arm on which depth is injective realizes every depth
  between `1` and its size;
* `treePath_append` and `treeDepth_add` — the decomposition of a root path through an
  intermediate vertex.
-/

namespace ClanAudit

open SimpleGraph Finset

variable {V : Type*} [Fintype V] [DecidableEq V] {G : SimpleGraph V} [DecidableRel G.Adj]

/-! ### Paths -/

omit [Fintype V] [DecidableEq V] [DecidableRel G.Adj] in
/-- Two paths meeting only at their joint concatenate to a path. -/
theorem walk_isPath_append {a b c : V} {p : G.Walk a b} {q : G.Walk b c} (hp : p.IsPath)
    (hq : q.IsPath) (hjoint : ∀ x : V, x ∈ p.support → x ∈ q.support → x = b) :
    (p.append q).IsPath := by
  rw [Walk.isPath_def, Walk.support_append]
  refine List.Nodup.append hp.support_nodup (hq.support_nodup.tail) ?_
  intro x hx hx'
  have hxq : x ∈ q.support := List.mem_of_mem_tail hx'
  have hxb : x = b := hjoint x hx hxq
  subst hxb
  -- `b` is the head of `q.support`, and `q.support` has no duplicates
  have hcons : q.support = x :: q.support.tail := (Walk.support_eq_cons q)
  have hnd : (x :: q.support.tail).Nodup := hcons ▸ hq.support_nodup
  exact (List.nodup_cons.mp hnd).1 hx'

variable {hG : G.IsTree} {h : V}

omit [Fintype V] [DecidableEq V] [DecidableRel G.Adj] in
/-- The root path from `v` to `h` is the reverse of the root path from `h` to `v`. -/
theorem treePath_reverse (v : V) : (treePath hG h v).reverse = treePath hG v h :=
  treePath_eq (treePath_isPath (hG := hG) (h := h) v).reverse

omit [Fintype V] [DecidableEq V] [DecidableRel G.Adj] in
theorem mem_support_treePath_comm {u v : V} (huv : u ∈ (treePath hG h v).support) :
    u ∈ (treePath hG v h).support := by
  rw [← treePath_reverse, Walk.support_reverse]
  simpa using huv

/-! ### Ancestors -/

omit [Fintype V] [DecidableRel G.Adj] in
/-- Every depth at most the depth of `v` occurs on the root path of `v`. -/
theorem exists_ancestor (v : V) (k : ℕ) (hk : k ≤ treeDepth hG h v) :
    ∃ w : V, w ∈ (treePath hG h v).support ∧ treeDepth hG h w = k := by
  generalize hj : treeDepth hG h v = j at hk
  induction j using Nat.strong_induction_on generalizing v with
  | _ j IH =>
    rcases eq_or_lt_of_le hk with rfl | hlt
    · exact ⟨v, Walk.end_mem_support _, hj⟩
    · have hvh : v ≠ h := by
        intro hc
        rw [hc, treeDepth_root (hG := hG)] at hj
        omega
      have hp : treeDepth hG h (treeParent hG h v) + 1 = treeDepth hG h v :=
        treeDepth_parent (hG := hG) (h := h) hvh
      obtain ⟨w, hw, hwk⟩ := IH (j - 1) (by omega) (treeParent hG h v) (by omega) (by omega)
      refine ⟨w, ?_, hwk⟩
      have hconcat := treePath_concat (hG := hG) (h := h)
        (treeParent_adj (hG := hG) (h := h) hvh) hp
      rw [← hconcat, Walk.support_concat, List.concat_eq_append]
      exact List.mem_append.mpr (Or.inl hw)

/-! ### Depth-injectivity along an arm, up to a bound -/

/-- The depth determines the vertex inside an arm, as long as all *strictly shallower*
vertices of the arm have degree at most two.  (A branch vertex at the far end of the arm
is therefore allowed.) -/
theorem treeArm_depth_inj_lt {x : V} {D : ℕ}
    (hdeg : ∀ w : V, w ≠ h → treeArm hG h w = x → treeDepth hG h w < D → G.degree w ≤ 2) :
    ∀ k : ℕ, k ≤ D → ∀ u v : V, u ≠ h → v ≠ h → treeArm hG h u = x → treeArm hG h v = x →
      treeDepth hG h u = k → treeDepth hG h v = k → u = v := by
  intro k
  induction k using Nat.strong_induction_on with
  | _ k IH =>
    intro hkD u v hu hv hau hav hdu hdv
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
        IH (k - 1) (by omega) (by omega) p q hph hqh hap haq (by omega) (by omega)
      by_contra hne
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
      exact absurd (hdeg p hph hap (by omega)) (by omega)

/-! ### An arm on which depth is injective realizes every depth -/

omit [DecidableRel G.Adj] in
theorem exists_mem_armFinset_depth {x : V}
    (hinj : ∀ u v : V, u ∈ armFinset hG h x → v ∈ armFinset hG h x →
      treeDepth hG h u = treeDepth hG h v → u = v)
    (k : ℕ) (hk1 : 1 ≤ k) (hk2 : k ≤ (armFinset hG h x).card) :
    ∃ v : V, v ∈ armFinset hG h x ∧ treeDepth hG h v = k := by
  have hmaps : ∀ v : V, ∀ _ : v ∈ armFinset hG h x,
      treeDepth hG h v ∈ Finset.Icc 1 (armFinset hG h x).card := by
    intro v hv
    rw [mem_armFinset] at hv
    refine Finset.mem_Icc.mpr ⟨treeDepth_pos (hG := hG) (h := h) hv.1, ?_⟩
    have := treeDepth_le_card_arm (hG := hG) (h := h) (v := v)
    rwa [hv.2] at this
  have hcard : (Finset.Icc 1 (armFinset hG h x).card).card ≤ (armFinset hG h x).card := by
    rw [Nat.card_Icc]
    omega
  obtain ⟨v, hv, hvk⟩ := Finset.surj_on_of_inj_on_of_card_le
    (s := armFinset hG h x) (t := Finset.Icc 1 (armFinset hG h x).card)
    (fun v _ => treeDepth hG h v) hmaps (fun u v hu hv he => hinj u v hu hv he) hcard k
    (Finset.mem_Icc.mpr ⟨hk1, hk2⟩)
  exact ⟨v, hv, hvk.symm⟩

/-! ### Adjacency in terms of depth and arm -/

omit [Fintype V] [DecidableRel G.Adj] in
/-- The neighbours of the root are exactly the vertices of depth one. -/
theorem adj_root_iff {v : V} (hv : v ≠ h) : G.Adj h v ↔ treeDepth hG h v = 1 := by
  constructor
  · intro hadj
    rcases treeDepth_adj (hG := hG) (h := h) hadj with h1 | h2
    · rw [treeDepth_root (hG := hG)] at h1; omega
    · rw [treeDepth_root (hG := hG)] at h2
      have := treeDepth_pos (hG := hG) (h := h) hv
      omega
  · intro hd
    have := treeArm_adj_root (hG := hG) (h := h) hv
    rwa [treeArm_of_depth_one (hG := hG) (h := h) hd] at this

omit [Fintype V] [DecidableRel G.Adj] in
/-- Two adjacent non-root vertices lie in the same arm. -/
theorem treeArm_eq_of_adj {u v : V} (hu : u ≠ h) (hv : v ≠ h) (hadj : G.Adj u v) :
    treeArm hG h u = treeArm hG h v := by
  rcases treeDepth_adj (hG := hG) (h := h) hadj with h1 | h2
  · have hmem : u ∈ (treePath hG h v).support := by
      rw [← treePath_concat (hG := hG) (h := h) hadj h1, Walk.support_concat,
        List.concat_eq_append]
      exact List.mem_append.mpr (Or.inl (Walk.end_mem_support _))
    exact (treeArm_of_mem_support (hG := hG) (h := h) hmem hu).symm
  · have hmem : v ∈ (treePath hG h u).support := by
      rw [← treePath_concat (hG := hG) (h := h) hadj.symm h2, Walk.support_concat,
        List.concat_eq_append]
      exact List.mem_append.mpr (Or.inl (Walk.end_mem_support _))
    exact treeArm_of_mem_support (hG := hG) (h := h) hmem hv

omit [Fintype V] [DecidableRel G.Adj] in
/-- Inside an arm on which depth is injective, consecutive depths are adjacent. -/
theorem adj_of_depth_succ_of_uniq {u v : V} (hu : u ≠ h) (hv : v ≠ h)
    (harm : treeArm hG h u = treeArm hG h v)
    (hd : treeDepth hG h u + 1 = treeDepth hG h v)
    (huniq : ∀ w : V, w ≠ h → treeArm hG h w = treeArm hG h u →
      treeDepth hG h w = treeDepth hG h u → w = u) : G.Adj u v := by
  have hpos : 0 < treeDepth hG h u := treeDepth_pos (hG := hG) (h := h) hu
  have hdp : treeDepth hG h (treeParent hG h v) + 1 = treeDepth hG h v :=
    treeDepth_parent (hG := hG) (h := h) hv
  have hph : treeParent hG h v ≠ h := by
    intro hc
    rw [hc, treeDepth_root (hG := hG)] at hdp
    omega
  have harmp : treeArm hG h (treeParent hG h v) = treeArm hG h u := by
    rw [treeArm_parent (hG := hG) (h := h) (by omega), harm]
  have : treeParent hG h v = u := huniq _ hph harmp (by omega)
  rw [← this]
  exact treeParent_adj (hG := hG) (h := h) hv

/-! ### Decomposing a root path through an intermediate vertex -/

omit [Fintype V] [DecidableRel G.Adj] in
/-- If `w` lies on the root path of `v`, that path splits at `w`. -/
theorem treePath_append {w v : V} (hw : w ∈ (treePath hG h v).support) :
    treePath hG h v = (treePath hG h w).append (treePath hG w v) := by
  have hspec := Walk.take_spec (treePath hG h v) hw
  have hdrop : (treePath hG h v).dropUntil w hw = treePath hG w v :=
    treePath_eq ((treePath_isPath (hG := hG) (h := h) v).dropUntil hw)
  rw [← hspec, treePath_takeUntil hw, hdrop]

omit [Fintype V] [DecidableRel G.Adj] in
/-- Depths add along a root path through an intermediate vertex. -/
theorem treeDepth_add {w v : V} (hw : w ∈ (treePath hG h v).support) :
    treeDepth hG h v = treeDepth hG h w + treeDepth hG w v := by
  have := congrArg Walk.length (treePath_append (hG := hG) (h := h) hw)
  rw [Walk.length_append] at this
  exact this

end ClanAudit
