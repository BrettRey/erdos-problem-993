import RequestProject.ClanAudit.TreeArm
import RequestProject.ClanAudit.Spider

/-!
# Recognising a rooted tree with no other branch vertex as a spider

If `G` is a finite tree with a vertex `h` such that every *other* vertex has degree at
most two, then `G` is isomorphic to the explicit model `spider n len`, where `n` is the
degree of `h` and `len i` is the number of vertices of the `i`-th arm.

The isomorphism is the "polar coordinate" map: a vertex `v ≠ h` is sent to the pair
consisting of the arm of `v` (the depth-one vertex on its root path) and its depth minus
one; the root goes to the hub.  Injectivity is `treeArm_depth_inj`; surjectivity is a
counting argument (the arms partition `V \ {h}`).
-/

namespace ClanAudit

open SimpleGraph Finset

variable {V : Type*} [Fintype V] [DecidableEq V] {G : SimpleGraph V} [DecidableRel G.Adj]

variable (G) in
/-- The `i`-th arm of the root `h`, i.e. the `i`-th neighbour of `h`. -/
noncomputable def hubArm (h : V) (i : Fin (G.neighborFinset h).card) : V :=
  ((G.neighborFinset h).equivFin.symm i : V)

omit [DecidableEq V] in
theorem hubArm_injective (h : V) : Function.Injective (hubArm G h) := by
  intro i j hij
  exact (G.neighborFinset h).equivFin.symm.injective (Subtype.coe_injective hij)

omit [DecidableEq V] in
theorem hubArm_adj (h : V) (i : Fin (G.neighborFinset h).card) : G.Adj h (hubArm G h i) := by
  have := ((G.neighborFinset h).equivFin.symm i).2
  rwa [mem_neighborFinset] at this

/-- The length of the `i`-th arm: the number of vertices whose root path passes through
the `i`-th neighbour of `h`. -/
noncomputable def hubArmLen (hG : G.IsTree) (h : V) (i : Fin (G.neighborFinset h).card) : ℕ :=
  (armFinset hG h (hubArm G h i)).card

/-- The index of the arm of a non-root vertex. -/
noncomputable def armIdx (hG : G.IsTree) (h : V) (v : V) (hv : v ≠ h) :
    Fin (G.neighborFinset h).card :=
  (G.neighborFinset h).equivFin
    ⟨treeArm hG h v, by rw [mem_neighborFinset]; exact treeArm_adj_root (hG := hG) (h := h) hv⟩

omit [DecidableEq V] in
theorem hubArm_armIdx (hG : G.IsTree) (h : V) {v : V} (hv : v ≠ h) :
    hubArm G h (armIdx hG h v hv) = treeArm hG h v := by
  simp [hubArm, armIdx]

omit [DecidableEq V] in
theorem armIdx_eq_iff (hG : G.IsTree) (h : V) {u v : V} (hu : u ≠ h) (hv : v ≠ h) :
    armIdx hG h u hu = armIdx hG h v hv ↔ treeArm hG h u = treeArm hG h v := by
  constructor
  · intro hij
    rw [← hubArm_armIdx hG h hu, ← hubArm_armIdx hG h hv, hij]
  · intro he
    exact hubArm_injective h (by rw [hubArm_armIdx hG h hu, hubArm_armIdx hG h hv, he])

/-- The polar-coordinate map of a rooted tree into the spider model. -/
noncomputable def toSpider (hG : G.IsTree) (h : V) :
    V → SpiderV (G.neighborFinset h).card (hubArmLen hG h) := fun v =>
  if hv : v = h then none
  else
    some ⟨armIdx hG h v hv, ⟨treeDepth hG h v - 1, by
      have h1 : hubArm G h (armIdx hG h v hv) = treeArm hG h v := hubArm_armIdx hG h hv
      have h2 : treeDepth hG h v ≤ (armFinset hG h (treeArm hG h v)).card :=
        treeDepth_le_card_arm (hG := hG) (h := h)
      have h3 : 0 < treeDepth hG h v := treeDepth_pos (hG := hG) (h := h) hv
      simp only [hubArmLen, h1]
      omega⟩⟩

theorem toSpider_root (hG : G.IsTree) (h : V) : toSpider hG h h = none := by
  simp [toSpider]

theorem toSpider_of_ne (hG : G.IsTree) (h : V) {v : V} (hv : v ≠ h) :
    ∃ hlt : treeDepth hG h v - 1 < hubArmLen hG h (armIdx hG h v hv),
      toSpider hG h v = some ⟨armIdx hG h v hv, ⟨treeDepth hG h v - 1, hlt⟩⟩ := by
  refine ⟨?_, by simp [toSpider, dif_neg hv]⟩
  have h1 : hubArm G h (armIdx hG h v hv) = treeArm hG h v := hubArm_armIdx hG h hv
  have h2 : treeDepth hG h v ≤ (armFinset hG h (treeArm hG h v)).card :=
    treeDepth_le_card_arm (hG := hG) (h := h)
  have h3 : 0 < treeDepth hG h v := treeDepth_pos (hG := hG) (h := h) hv
  simp only [hubArmLen, h1]
  omega

theorem toSpider_ne_none (hG : G.IsTree) (h : V) {v : V} (hv : v ≠ h) :
    toSpider hG h v ≠ none := by
  obtain ⟨_, he⟩ := toSpider_of_ne hG h hv
  rw [he]
  exact Option.some_ne_none _

/-- The level of the image is the depth of the vertex. -/
theorem spiderLevel_toSpider (hG : G.IsTree) (h : V) (v : V) :
    spiderLevel (toSpider hG h v) = treeDepth hG h v := by
  by_cases hv : v = h
  · subst hv; rw [toSpider_root, treeDepth_root (hG := hG)]; rfl
  · obtain ⟨_, he⟩ := toSpider_of_ne hG h hv
    have h3 : 0 < treeDepth hG h v := treeDepth_pos (hG := hG) (h := h) hv
    rw [he]
    show treeDepth hG h v - 1 + 1 = treeDepth hG h v
    omega

/-! ### Injectivity -/

theorem toSpider_injective (hG : G.IsTree) (h : V) (hdeg : ∀ v : V, v ≠ h → G.degree v ≤ 2) :
    Function.Injective (toSpider hG h) := by
  intro u v huv
  by_cases hu : u = h
  · by_cases hv : v = h
    · rw [hu, hv]
    · refine absurd ?_ (toSpider_ne_none hG h hv)
      rw [← huv, hu, toSpider_root]
  · by_cases hv : v = h
    · refine absurd ?_ (toSpider_ne_none hG h hu)
      rw [huv, hv, toSpider_root]
    · obtain ⟨hu', heu⟩ := toSpider_of_ne hG h hu
      obtain ⟨hv', hev⟩ := toSpider_of_ne hG h hv
      have hdepth : treeDepth hG h u = treeDepth hG h v := by
        have := congrArg spiderLevel huv
        rwa [spiderLevel_toSpider, spiderLevel_toSpider] at this
      have hidx : armIdx hG h u hu = armIdx hG h v hv := by
        rw [heu, hev] at huv
        exact congrArg (fun x => x.1) (Option.some.inj huv)
      have harm : treeArm hG h u = treeArm hG h v := (armIdx_eq_iff hG h hu hv).mp hidx
      exact treeArm_depth_inj (hG := hG) (h := h)
        (x := treeArm hG h u) (fun w hw _ => hdeg w hw) (treeDepth hG h u) u v hu hv rfl
        harm.symm rfl hdepth.symm

/-! ### Surjectivity by counting -/

theorem card_eq_sum_hubArmLen (hG : G.IsTree) (h : V) :
    Fintype.card V = (∑ i, hubArmLen hG h i) + 1 := by
  have hfib : (univ.erase h).card
      = ∑ x ∈ G.neighborFinset h, ((univ.erase h).filter (fun v => treeArm hG h v = x)).card := by
    refine Finset.card_eq_sum_card_fiberwise ?_
    intro v hv
    simp only [Finset.mem_coe, Finset.mem_erase] at hv
    simp only [Finset.mem_coe, mem_neighborFinset]
    exact treeArm_adj_root (hG := hG) (h := h) hv.1
  have hfilter : ∀ x : V, (univ.erase h).filter (fun v => treeArm hG h v = x)
      = armFinset hG h x := by
    intro x
    ext v
    simp [Finset.mem_erase, mem_armFinset, and_comm]
  have hsum : ∑ x ∈ G.neighborFinset h, (armFinset hG h x).card = ∑ i, hubArmLen hG h i := by
    rw [← Finset.sum_coe_sort (G.neighborFinset h) (fun x => (armFinset hG h x).card)]
    exact (Equiv.sum_comp (G.neighborFinset h).equivFin.symm
      (fun x : (G.neighborFinset h : Finset V) => (armFinset hG h (x : V)).card)).symm
  have hcard : (univ.erase h).card + 1 = Fintype.card V := by
    rw [Finset.card_erase_of_mem (Finset.mem_univ h)]
    have : 0 < Fintype.card V := Fintype.card_pos_iff.mpr ⟨h⟩
    simp only [Finset.card_univ]
    omega
  rw [← hcard, hfib]
  simp only [hfilter]
  rw [hsum]

theorem card_spiderV (n : ℕ) (len : Fin n → ℕ) :
    Fintype.card (SpiderV n len) = (∑ i, len i) + 1 := by
  simp [SpiderV, Fintype.card_sigma]

theorem toSpider_bijective (hG : G.IsTree) (h : V) (hdeg : ∀ v : V, v ≠ h → G.degree v ≤ 2) :
    Function.Bijective (toSpider hG h) := by
  refine (Fintype.bijective_iff_injective_and_card _).mpr ⟨toSpider_injective hG h hdeg, ?_⟩
  rw [card_eq_sum_hubArmLen hG h, card_spiderV]

/-! ### The isomorphism -/

theorem toSpider_adj (hG : G.IsTree) (h : V) (hdeg : ∀ v : V, v ≠ h → G.degree v ≤ 2)
    (u v : V) :
    (spider (G.neighborFinset h).card (hubArmLen hG h)).Adj (toSpider hG h u) (toSpider hG h v)
      ↔ G.Adj u v := by
  by_cases hu : u = h
  · by_cases hv : v = h
    · rw [hu, hv]
      simp [toSpider_root, spider]
    · rw [hu]
      obtain ⟨hv', hev⟩ := toSpider_of_ne hG h hv
      rw [toSpider_root, hev]
      have hpos : 0 < treeDepth hG h v := treeDepth_pos (hG := hG) (h := h) hv
      show treeDepth hG h v - 1 = 0 ↔ _
      constructor
      · intro hd
        have hd1 : treeDepth hG h v = 1 := by omega
        have := treeArm_adj_root (hG := hG) (h := h) hv
        rwa [treeArm_of_depth_one (hG := hG) (h := h) hd1] at this
      · intro hadj
        rcases treeDepth_adj (hG := hG) (h := h) hadj with h1 | h2
        · rw [treeDepth_root (hG := hG)] at h1; omega
        · rw [treeDepth_root (hG := hG)] at h2; omega
  · by_cases hv : v = h
    · rw [hv]
      obtain ⟨hu', heu⟩ := toSpider_of_ne hG h hu
      rw [toSpider_root, heu]
      have hpos : 0 < treeDepth hG h u := treeDepth_pos (hG := hG) (h := h) hu
      show treeDepth hG h u - 1 = 0 ↔ _
      constructor
      · intro hd
        have hd1 : treeDepth hG h u = 1 := by omega
        have := treeArm_adj_root (hG := hG) (h := h) hu
        rw [treeArm_of_depth_one (hG := hG) (h := h) hd1] at this
        exact this.symm
      · intro hadj
        rcases treeDepth_adj (hG := hG) (h := h) hadj.symm with h1 | h2
        · rw [treeDepth_root (hG := hG)] at h1; omega
        · rw [treeDepth_root (hG := hG)] at h2; omega
    · obtain ⟨hu', heu⟩ := toSpider_of_ne hG h hu
      obtain ⟨hv', hev⟩ := toSpider_of_ne hG h hv
      rw [heu, hev]
      have hpu : 0 < treeDepth hG h u := treeDepth_pos (hG := hG) (h := h) hu
      have hpv : 0 < treeDepth hG h v := treeDepth_pos (hG := hG) (h := h) hv
      show (armIdx hG h u hu = armIdx hG h v hv ∧
          (treeDepth hG h u - 1 + 1 = treeDepth hG h v - 1 ∨
            treeDepth hG h v - 1 + 1 = treeDepth hG h u - 1)) ↔ _
      rw [armIdx_eq_iff hG h hu hv]
      constructor
      · rintro ⟨harm, hd⟩
        -- the shallower of the two is the parent of the deeper one
        rcases hd with hd | hd
        · have hdd : treeDepth hG h u + 1 = treeDepth hG h v := by omega
          have hp : treeParent hG h v ≠ h := by
            intro hc
            have := treeDepth_parent (hG := hG) (h := h) hv
            rw [hc, treeDepth_root (hG := hG)] at this
            omega
          have hap : treeArm hG h (treeParent hG h v) = treeArm hG h u := by
            rw [treeArm_parent (hG := hG) (h := h) (by omega), harm]
          have : u = treeParent hG h v :=
            treeArm_depth_inj (hG := hG) (h := h) (x := treeArm hG h u)
              (fun w hw _ => hdeg w hw) (treeDepth hG h u) u (treeParent hG h v) hu hp rfl
              hap rfl (by have := treeDepth_parent (hG := hG) (h := h) hv; omega)
          rw [this]
          exact treeParent_adj (hG := hG) (h := h) hv
        · have hdd : treeDepth hG h v + 1 = treeDepth hG h u := by omega
          have hp : treeParent hG h u ≠ h := by
            intro hc
            have := treeDepth_parent (hG := hG) (h := h) hu
            rw [hc, treeDepth_root (hG := hG)] at this
            omega
          have hap : treeArm hG h (treeParent hG h u) = treeArm hG h u := by
            rw [treeArm_parent (hG := hG) (h := h) (by omega)]
          have : v = treeParent hG h u :=
            treeArm_depth_inj (hG := hG) (h := h) (x := treeArm hG h u)
              (fun w hw _ => hdeg w hw) (treeDepth hG h v) v (treeParent hG h u) hv hp harm.symm
              hap rfl (by have := treeDepth_parent (hG := hG) (h := h) hu; omega)
          rw [this]
          exact (treeParent_adj (hG := hG) (h := h) hu).symm
      · intro hadj
        rcases treeDepth_adj (hG := hG) (h := h) hadj with h1 | h2
        · have hmem : u ∈ (treePath hG h v).support := by
            rw [← treePath_concat (hG := hG) (h := h) hadj h1]
            simp [Walk.support_concat]
          refine ⟨(treeArm_of_mem_support (hG := hG) (h := h) hmem hu).symm, Or.inl (by omega)⟩
        · have hmem : v ∈ (treePath hG h u).support := by
            rw [← treePath_concat (hG := hG) (h := h) hadj.symm h2]
            simp [Walk.support_concat]
          refine ⟨treeArm_of_mem_support (hG := hG) (h := h) hmem hv, Or.inr (by omega)⟩

/-- **Spider recognition.**  A finite tree in which every vertex other than `h` has degree
at most two is isomorphic to the spider whose arms are the arms of `h`. -/
noncomputable def spiderIsoOfTree (hG : G.IsTree) (h : V)
    (hdeg : ∀ v : V, v ≠ h → G.degree v ≤ 2) :
    G ≃g spider (G.neighborFinset h).card (hubArmLen hG h) where
  toEquiv := Equiv.ofBijective _ (toSpider_bijective hG h hdeg)
  map_rel_iff' := by
    intro u v
    exact toSpider_adj hG h hdeg u v

end ClanAudit
