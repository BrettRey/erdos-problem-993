import RequestProject.B1Zero.Corona

/-!
# The codimension profile of a corona

`CoronaData.vec C X A` is the vector of the first five codimension counts of the independent
subsets of the corona over `X` that avoid `A`.  The two recursions matching
`RootedForestCertificate.profile` are proved here:

* `vec_union`  — a corona splitting with no edges across multiplies profiles;
* `vec_peel`, `vec_peel_root` — peeling a base vertex off a corona.
-/

open Finset SimpleGraph

namespace MatchingBag

namespace B1Zero

open DepthThree.RootedForestCertificate

attribute [local instance 1] Classical.propDecidable

variable {V : Type*} [Fintype V] [DecidableEq V] {G : SimpleGraph V}

/-! ### Elementary facts about `cntRev` -/

lemma cntRev_succ_succ (G : SimpleGraph V) (W : Finset V) (n j : ℕ) :
    cntRev G W (n + 1) (j + 1) = cntRev G W n j := by
  rw [cntRev, cntRev, revSets, revSets]
  congr 1
  refine Finset.filter_congr fun S _ => ?_
  constructor
  · rintro ⟨h1, h2⟩; exact ⟨h1, by omega⟩
  · rintro ⟨h1, h2⟩; exact ⟨h1, by omega⟩

lemma cntRev_succ_zero {G : SimpleGraph V} {W : Finset V} {n : ℕ}
    (hb : ∀ S ⊆ W, IndepOn G S → S.card ≤ n) : cntRev G W (n + 1) 0 = 0 := by
  rw [cntRev, Finset.card_eq_zero, Finset.eq_empty_iff_forall_notMem]
  intro S hS
  rw [mem_revSets] at hS
  have := hb S hS.1 hS.2.1
  omega

namespace CoronaData

variable (C : CoronaData G)

/-- The number of independent subsets of the corona over `X` avoiding `A`, of codimension
`j` below `|X|`. -/
noncomputable def vecAt (X A : Finset V) (j : ℕ) : ℕ := cntRev G (C.cor X \ A) X.card j

/-- The first five codimension counts of the corona over `X` avoiding `A`. -/
noncomputable def vec (X A : Finset V) : TopFive :=
  ⟨C.vecAt X A 0, C.vecAt X A 1, C.vecAt X A 2, C.vecAt X A 3, C.vecAt X A 4⟩

variable {C}

lemma vecAt_empty_avoid (X : Finset V) (j : ℕ) :
    C.vecAt X ∅ j = cntRev G (C.cor X) X.card j := by
  rw [vecAt, Finset.sdiff_empty]

lemma vec_congr {X A A' : Finset V} (h : C.cor X \ A = C.cor X \ A') :
    C.vec X A = C.vec X A' := by
  simp only [vec, vecAt, h]

lemma vec_empty (A : Finset V) : C.vec (∅ : Finset V) A = ⟨1, 0, 0, 0, 0⟩ := by
  have h : ∀ j, C.vecAt (∅ : Finset V) A j = if j = 0 then 1 else 0 := by
    intro j
    rw [vecAt, cor_empty, Finset.empty_sdiff]
    rcases j with _ | j
    · have hset : revSets G (∅ : Finset V) (∅ : Finset V).card 0 = {∅} := by
        ext S
        rw [mem_revSets, Finset.mem_singleton, Finset.subset_empty]
        constructor
        · rintro ⟨rfl, -, -⟩; rfl
        · rintro rfl
          exact ⟨rfl, by intro u hu; simp at hu, by simp⟩
      rw [cntRev, hset, Finset.card_singleton, if_pos rfl]
    · rw [if_neg (by omega), cntRev, Finset.card_eq_zero, Finset.eq_empty_iff_forall_notMem]
      intro S hS
      rw [mem_revSets, Finset.subset_empty] at hS
      obtain ⟨rfl, -, hc⟩ := hS
      simp at hc
  rw [vec]
  simp [h]

/-- Disjoint base sets have disjoint coronas. -/
lemma cor_disjoint {X Y : Finset V} (hX : X ⊆ C.base) (hY : Y ⊆ C.base)
    (hdisj : Disjoint X Y) : Disjoint (C.cor X) (C.cor Y) := by
  rw [Finset.disjoint_left]
  intro z hz hz'
  rcases mem_cor.1 hz with h | ⟨b, hb, rfl⟩
  · rcases mem_cor.1 hz' with h' | ⟨b', hb', hbz⟩
    · exact Finset.disjoint_left.1 hdisj h h'
    · exact C.pend_notMem b' (hY hb') (hbz ▸ hX h)
  · rcases mem_cor.1 hz' with h' | ⟨b', hb', hbz⟩
    · exact C.pend_notMem b (hX hb) (hY h')
    · exact Finset.disjoint_left.1 hdisj hb (C.pend_inj b' (hY hb') b (hX hb) hbz ▸ hb')

/-- **The product rule.** -/
theorem vec_union {X Y : Finset V} (hX : X ⊆ C.base) (hY : Y ⊆ C.base)
    (hdisj : Disjoint X Y) (hcross : ∀ u ∈ X, ∀ w ∈ Y, ¬ G.Adj u w) (A : Finset V) :
    C.vec (X ∪ Y) A = mul (C.vec X A) (C.vec Y A) := by
  have key : ∀ j, C.vecAt (X ∪ Y) A j
      = ∑ ij ∈ Finset.antidiagonal j, C.vecAt X A ij.1 * C.vecAt Y A ij.2 := by
    intro j
    have hW : C.cor (X ∪ Y) \ A = (C.cor X \ A) ∪ (C.cor Y \ A) := by
      rw [cor_union, Finset.union_sdiff_distrib]
    have hcd : (X ∪ Y).card = X.card + Y.card := Finset.card_union_of_disjoint hdisj
    rw [vecAt, hW, hcd]
    refine cntRev_union ((cor_disjoint hX hY hdisj).mono Finset.sdiff_subset Finset.sdiff_subset)
      ?_ ?_ ?_ j
    · intro u hu w hw
      exact cor_no_cross hX hY hdisj hcross u (Finset.mem_sdiff.1 hu).1 w
        (Finset.mem_sdiff.1 hw).1
    · intro S hS hind
      exact card_le_of_indep hX (hS.trans Finset.sdiff_subset) hind
    · intro S hS hind
      exact card_le_of_indep hY (hS.trans Finset.sdiff_subset) hind
  have e0 := key 0
  have e1 := key 1
  have e2 := key 2
  have e3 := key 3
  have e4 := key 4
  simp only [Finset.Nat.sum_antidiagonal_eq_sum_range_succ_mk, Finset.sum_range_succ,
    Finset.sum_range_zero, Nat.zero_add] at e0 e1 e2 e3 e4
  simp only [vec, mul, TopFive.mk.injEq]
  refine ⟨?_, ?_, ?_, ?_, ?_⟩
  · rw [e0]
  · rw [e1]
  · rw [e2]
  · rw [e3]
  · rw [e4]

/-! ### Peeling a base vertex -/

lemma not_adj_root {Ch : Finset V} (hCh : Ch ⊆ C.base) {v : V} (hv : v ∈ C.base)
    (hvCh : v ∉ Ch) {w : V} (hw : w ∈ C.cor Ch \ Ch.filter (fun c => G.Adj v c)) :
    ¬ G.Adj v w := by
  rw [Finset.mem_sdiff, Finset.mem_filter] at hw
  rcases mem_cor.1 hw.1 with h | ⟨b, hb, rfl⟩
  · exact fun hadj => hw.2 ⟨h, hadj⟩
  · intro hadj
    have := C.pend_nbr b (hCh hb) v (cor_mono (Finset.Subset.refl _) (subset_cor hv)) hadj.symm
    exact hvCh (this ▸ hb)

/-- Sets containing the root. -/
lemma peel_root_fibre {Ch : Finset V} (hCh : Ch ⊆ C.base) {v : V} (hv : v ∈ C.base)
    (hvCh : v ∉ Ch) (j : ℕ) :
    (revSets G (C.cor (insert v Ch)) (Ch.card + 1) j).filter (fun S => v ∈ S)
      = (revSets G (C.cor Ch \ Ch.filter (fun c => G.Adj v c)) Ch.card j).image (insert v) := by
  have hvC : v ∉ C.cor Ch := notMem_cor_of_base hCh hv hvCh
  have hpadj : G.Adj v (C.pend v) := C.adj_pend v hv
  ext S
  simp only [Finset.mem_filter, Finset.mem_image, mem_revSets]
  constructor
  · rintro ⟨⟨hsub, hind, hcd⟩, hvS⟩
    refine ⟨S.erase v, ⟨?_, hind.subset (Finset.erase_subset _ _), ?_⟩, ?_⟩
    · intro x hx
      rw [Finset.mem_erase] at hx
      have hxS := hsub hx.2
      rw [cor_insert, Finset.mem_insert, Finset.mem_insert] at hxS
      rcases hxS with rfl | rfl | hxc
      · exact absurd rfl hx.1
      · exact absurd hpadj (hind v hvS _ hx.2)
      · rw [Finset.mem_sdiff, Finset.mem_filter]
        exact ⟨hxc, fun hxn => hind v hvS x hx.2 hxn.2⟩
    · rw [Finset.card_erase_of_mem hvS]
      have : 1 ≤ S.card := Finset.card_pos.2 ⟨v, hvS⟩
      omega
    · exact Finset.insert_erase hvS
  · rintro ⟨T, ⟨hsub, hind, hcd⟩, rfl⟩
    have hvT : v ∉ T := fun h => hvC (Finset.mem_sdiff.1 (hsub h)).1
    refine ⟨⟨?_, ?_, ?_⟩, Finset.mem_insert_self _ _⟩
    · intro x hx
      rw [Finset.mem_insert] at hx
      rw [cor_insert]
      rcases hx with rfl | hx
      · exact Finset.mem_insert_self _ _
      · exact Finset.mem_insert_of_mem
          (Finset.mem_insert_of_mem (Finset.mem_sdiff.1 (hsub hx)).1)
    · intro a ha b hb
      rw [Finset.mem_insert] at ha hb
      rcases ha with rfl | ha <;> rcases hb with rfl | hb
      · exact G.irrefl
      · exact not_adj_root hCh hv hvCh (hsub hb)
      · exact fun hadj => not_adj_root hCh hv hvCh (hsub ha) hadj.symm
      · exact hind a ha b hb
    · rw [Finset.card_insert_of_notMem hvT]
      omega

/-- Sets containing the pendant of the root but not the root. -/
lemma peel_pend_fibre {Ch : Finset V} (hCh : Ch ⊆ C.base) {v : V} (hv : v ∈ C.base)
    (hvCh : v ∉ Ch) (j : ℕ) :
    (revSets G (C.cor (insert v Ch)) (Ch.card + 1) j).filter
        (fun S => v ∉ S ∧ C.pend v ∈ S)
      = (revSets G (C.cor Ch) Ch.card j).image (insert (C.pend v)) := by
  have hvC : v ∉ C.cor Ch := notMem_cor_of_base hCh hv hvCh
  have hpC : C.pend v ∉ C.cor Ch := pend_notMem_cor hCh hv hvCh
  have hpv : C.pend v ≠ v := (C.adj_pend v hv).ne'
  ext S
  simp only [Finset.mem_filter, Finset.mem_image, mem_revSets]
  constructor
  · rintro ⟨⟨hsub, hind, hcd⟩, hvS, hpS⟩
    refine ⟨S.erase (C.pend v), ⟨?_, hind.subset (Finset.erase_subset _ _), ?_⟩, ?_⟩
    · intro x hx
      rw [Finset.mem_erase] at hx
      have hxS := hsub hx.2
      rw [cor_insert, Finset.mem_insert, Finset.mem_insert] at hxS
      rcases hxS with rfl | rfl | hxc
      · exact absurd hx.2 hvS
      · exact absurd rfl hx.1
      · exact hxc
    · rw [Finset.card_erase_of_mem hpS]
      have : 1 ≤ S.card := Finset.card_pos.2 ⟨_, hpS⟩
      omega
    · exact Finset.insert_erase hpS
  · rintro ⟨T, ⟨hsub, hind, hcd⟩, rfl⟩
    have hpT : C.pend v ∉ T := fun h => hpC (hsub h)
    refine ⟨⟨?_, ?_, ?_⟩, ?_, Finset.mem_insert_self _ _⟩
    · intro x hx
      rw [Finset.mem_insert] at hx
      rw [cor_insert]
      rcases hx with rfl | hx
      · exact Finset.mem_insert_of_mem (Finset.mem_insert_self _ _)
      · exact Finset.mem_insert_of_mem (Finset.mem_insert_of_mem (hsub hx))
    · intro a ha b hb
      rw [Finset.mem_insert] at ha hb
      rcases ha with rfl | ha <;> rcases hb with rfl | hb
      · exact G.irrefl
      · exact not_adj_pend_of_notMem hCh hv hvCh b (hsub hb)
      · exact fun hadj => not_adj_pend_of_notMem hCh hv hvCh a (hsub ha) hadj.symm
      · exact hind a ha b hb
    · rw [Finset.card_insert_of_notMem hpT]
      omega
    · rw [Finset.mem_insert]
      rintro (h | h)
      · exact hpv h.symm
      · exact hvC (hsub h)

/-- Sets avoiding the root and its pendant. -/
lemma peel_rest_fibre {Ch : Finset V} (hCh : Ch ⊆ C.base) {v : V} (hv : v ∈ C.base)
    (hvCh : v ∉ Ch) (j : ℕ) :
    (revSets G (C.cor (insert v Ch)) (Ch.card + 1) j).filter
        (fun S => v ∉ S ∧ C.pend v ∉ S)
      = revSets G (C.cor Ch) (Ch.card + 1) j := by
  have hvC : v ∉ C.cor Ch := notMem_cor_of_base hCh hv hvCh
  have hpC : C.pend v ∉ C.cor Ch := pend_notMem_cor hCh hv hvCh
  ext S
  simp only [Finset.mem_filter, mem_revSets]
  constructor
  · rintro ⟨⟨hsub, hind, hcd⟩, hvS, hpS⟩
    refine ⟨fun x hx => ?_, hind, hcd⟩
    have hxS := hsub hx
    rw [cor_insert, Finset.mem_insert, Finset.mem_insert] at hxS
    rcases hxS with rfl | rfl | hxc
    · exact absurd hx hvS
    · exact absurd hx hpS
    · exact hxc
  · rintro ⟨hsub, hind, hcd⟩
    refine ⟨⟨fun x hx => ?_, hind, hcd⟩, fun h => hvC (hsub h), fun h => hpC (hsub h)⟩
    rw [cor_insert]
    exact Finset.mem_insert_of_mem (Finset.mem_insert_of_mem (hsub hx))

/-- **Peeling the root.** -/
theorem vecAt_peel {Ch : Finset V} (hCh : Ch ⊆ C.base) {v : V} (hv : v ∈ C.base)
    (hvCh : v ∉ Ch) (j : ℕ) :
    C.vecAt (insert v Ch) ∅ j
      = C.vecAt Ch (Ch.filter (fun c => G.Adj v c)) j + C.vecAt Ch ∅ j
        + cntRev G (C.cor Ch) (Ch.card + 1) j := by
  classical
  have hcard : (insert v Ch).card = Ch.card + 1 := Finset.card_insert_of_notMem hvCh
  rw [vecAt_empty_avoid, hcard, cntRev]
  set 𝒜 := revSets G (C.cor (insert v Ch)) (Ch.card + 1) j with h𝒜
  have h1 := Finset.card_filter_add_card_filter_not (s := 𝒜) (p := fun S => v ∈ S)
  have h2 := Finset.card_filter_add_card_filter_not
    (s := 𝒜.filter (fun S => v ∉ S)) (p := fun S => C.pend v ∈ S)
  rw [Finset.filter_filter, Finset.filter_filter] at h2
  have ha : (𝒜.filter (fun S => v ∈ S)).card
      = C.vecAt Ch (Ch.filter (fun c => G.Adj v c)) j := by
    have hinj : Set.InjOn (insert v)
        ((revSets G (C.cor Ch \ Ch.filter (fun c => G.Adj v c)) Ch.card j :
          Finset (Finset V)) : Set (Finset V)) := by
      intro T hT T' hT' heq
      rw [Finset.mem_coe] at hT hT'
      have hvT : v ∉ T := fun h =>
        notMem_cor_of_base hCh hv hvCh (Finset.mem_sdiff.1 ((mem_revSets.1 hT).1 h)).1
      have hvT' : v ∉ T' := fun h =>
        notMem_cor_of_base hCh hv hvCh (Finset.mem_sdiff.1 ((mem_revSets.1 hT').1 h)).1
      rw [← Finset.erase_insert hvT, ← Finset.erase_insert hvT', heq]
    rw [h𝒜, peel_root_fibre hCh hv hvCh j, vecAt, Finset.card_image_of_injOn hinj, cntRev]
  have hb : (𝒜.filter (fun S => v ∉ S ∧ C.pend v ∈ S)).card = C.vecAt Ch ∅ j := by
    have hinj : Set.InjOn (insert (C.pend v))
        ((revSets G (C.cor Ch) Ch.card j : Finset (Finset V)) : Set (Finset V)) := by
      intro T hT T' hT' heq
      rw [Finset.mem_coe] at hT hT'
      have hpT : C.pend v ∉ T := fun h =>
        pend_notMem_cor hCh hv hvCh ((mem_revSets.1 hT).1 h)
      have hpT' : C.pend v ∉ T' := fun h =>
        pend_notMem_cor hCh hv hvCh ((mem_revSets.1 hT').1 h)
      rw [← Finset.erase_insert hpT, ← Finset.erase_insert hpT', heq]
    rw [h𝒜, peel_pend_fibre hCh hv hvCh j, vecAt_empty_avoid, cntRev,
      Finset.card_image_of_injOn hinj]
  have hc : (𝒜.filter (fun S => v ∉ S ∧ C.pend v ∉ S)).card
      = cntRev G (C.cor Ch) (Ch.card + 1) j := by
    rw [h𝒜, peel_rest_fibre hCh hv hvCh j, cntRev]
  omega

/-- **Peeling the root when the root itself is avoided.** -/
theorem vecAt_peel_root {Ch : Finset V} (hCh : Ch ⊆ C.base) {v : V} (hv : v ∈ C.base)
    (hvCh : v ∉ Ch) (j : ℕ) :
    C.vecAt (insert v Ch) {v} j
      = C.vecAt Ch ∅ j + cntRev G (C.cor Ch) (Ch.card + 1) j := by
  classical
  have hvC : v ∉ C.cor Ch := notMem_cor_of_base hCh hv hvCh
  have hpv : C.pend v ≠ v := (C.adj_pend v hv).ne'
  have hcard : (insert v Ch).card = Ch.card + 1 := Finset.card_insert_of_notMem hvCh
  have hWeq : C.cor (insert v Ch) \ {v} = insert (C.pend v) (C.cor Ch) := by
    ext x
    simp only [cor_insert, Finset.mem_sdiff, Finset.mem_insert, Finset.mem_singleton]
    constructor
    · rintro ⟨rfl | rfl | h, hx⟩
      · exact absurd rfl hx
      · exact Or.inl rfl
      · exact Or.inr h
    · rintro (rfl | hx)
      · exact ⟨Or.inr (Or.inl rfl), hpv⟩
      · exact ⟨Or.inr (Or.inr hx), fun h => hvC (h ▸ hx)⟩
  rw [vecAt, hWeq, hcard, cntRev]
  set ℬ := revSets G (insert (C.pend v) (C.cor Ch)) (Ch.card + 1) j with hℬ
  have h1 := Finset.card_filter_add_card_filter_not (s := ℬ)
    (p := fun S => C.pend v ∈ S)
  have hb : (ℬ.filter (fun S => C.pend v ∈ S)).card = C.vecAt Ch ∅ j := by
    have himg : ℬ.filter (fun S => C.pend v ∈ S)
        = (revSets G (C.cor Ch) Ch.card j).image (insert (C.pend v)) := by
      ext S
      simp only [hℬ, Finset.mem_filter, Finset.mem_image, mem_revSets]
      constructor
      · rintro ⟨⟨hsub, hind, hcd⟩, hpS⟩
        refine ⟨S.erase (C.pend v), ⟨?_, hind.subset (Finset.erase_subset _ _), ?_⟩, ?_⟩
        · intro x hx
          rw [Finset.mem_erase] at hx
          have := hsub hx.2
          rw [Finset.mem_insert] at this
          rcases this with rfl | h
          · exact absurd rfl hx.1
          · exact h
        · rw [Finset.card_erase_of_mem hpS]
          have : 1 ≤ S.card := Finset.card_pos.2 ⟨_, hpS⟩
          omega
        · exact Finset.insert_erase hpS
      · rintro ⟨T, ⟨hsub, hind, hcd⟩, rfl⟩
        have hpT : C.pend v ∉ T := fun h => pend_notMem_cor hCh hv hvCh (hsub h)
        refine ⟨⟨?_, ?_, ?_⟩, Finset.mem_insert_self _ _⟩
        · intro x hx
          rw [Finset.mem_insert] at hx ⊢
          rcases hx with rfl | hx
          · exact Or.inl rfl
          · exact Or.inr (hsub hx)
        · intro a ha b hb
          rw [Finset.mem_insert] at ha hb
          rcases ha with rfl | ha <;> rcases hb with rfl | hb
          · exact G.irrefl
          · exact not_adj_pend_of_notMem hCh hv hvCh b (hsub hb)
          · exact fun hadj => not_adj_pend_of_notMem hCh hv hvCh a (hsub ha) hadj.symm
          · exact hind a ha b hb
        · rw [Finset.card_insert_of_notMem hpT]
          omega
    have hinj : Set.InjOn (insert (C.pend v))
        ((revSets G (C.cor Ch) Ch.card j : Finset (Finset V)) : Set (Finset V)) := by
      intro T hT T' hT' heq
      rw [Finset.mem_coe] at hT hT'
      have hpT : C.pend v ∉ T := fun h => pend_notMem_cor hCh hv hvCh ((mem_revSets.1 hT).1 h)
      have hpT' : C.pend v ∉ T' := fun h => pend_notMem_cor hCh hv hvCh ((mem_revSets.1 hT').1 h)
      rw [← Finset.erase_insert hpT, ← Finset.erase_insert hpT', heq]
    rw [himg, vecAt_empty_avoid, cntRev, Finset.card_image_of_injOn hinj]
  have hc : (ℬ.filter (fun S => ¬ C.pend v ∈ S)).card
      = cntRev G (C.cor Ch) (Ch.card + 1) j := by
    have hset : ℬ.filter (fun S => ¬ C.pend v ∈ S) = revSets G (C.cor Ch) (Ch.card + 1) j := by
      ext S
      simp only [hℬ, Finset.mem_filter, mem_revSets]
      constructor
      · rintro ⟨⟨hsub, hind, hcd⟩, hpS⟩
        refine ⟨fun x hx => ?_, hind, hcd⟩
        have := hsub hx
        rw [Finset.mem_insert] at this
        rcases this with rfl | h
        · exact absurd hx hpS
        · exact h
      · rintro ⟨hsub, hind, hcd⟩
        exact ⟨⟨fun x hx => Finset.mem_insert_of_mem (hsub hx), hind, hcd⟩,
          fun h => pend_notMem_cor hCh hv hvCh (hsub h)⟩
    rw [hset, cntRev]
  omega

/-- The `TopFive` form of the peeling identity. -/
theorem vec_peel {Ch : Finset V} (hCh : Ch ⊆ C.base) {v : V} (hv : v ∈ C.base)
    (hvCh : v ∉ Ch) :
    C.vec (insert v Ch) ∅
      = add (add (C.vec Ch ∅) (shift (C.vec Ch ∅))) (C.vec Ch (Ch.filter (fun c => G.Adj v c))) := by
  have hzero : cntRev G (C.cor Ch) (Ch.card + 1) 0 = 0 :=
    cntRev_succ_zero (fun S hS hind => card_le_of_indep hCh hS hind)
  have hsucc : ∀ j, cntRev G (C.cor Ch) (Ch.card + 1) (j + 1) = C.vecAt Ch ∅ j := by
    intro j
    rw [cntRev_succ_succ, vecAt_empty_avoid]
  have e0 := vecAt_peel hCh hv hvCh 0
  have e1 := vecAt_peel hCh hv hvCh 1
  have e2 := vecAt_peel hCh hv hvCh 2
  have e3 := vecAt_peel hCh hv hvCh 3
  have e4 := vecAt_peel hCh hv hvCh 4
  rw [hzero] at e0
  rw [hsucc 0] at e1
  rw [hsucc 1] at e2
  rw [hsucc 2] at e3
  rw [hsucc 3] at e4
  simp only [vec, add, shift, TopFive.mk.injEq]
  refine ⟨by omega, by omega, by omega, by omega, by omega⟩

/-- The `TopFive` form of the peeling identity when the root is avoided. -/
theorem vec_peel_root {Ch : Finset V} (hCh : Ch ⊆ C.base) {v : V} (hv : v ∈ C.base)
    (hvCh : v ∉ Ch) :
    C.vec (insert v Ch) {v} = add (C.vec Ch ∅) (shift (C.vec Ch ∅)) := by
  have hzero : cntRev G (C.cor Ch) (Ch.card + 1) 0 = 0 :=
    cntRev_succ_zero (fun S hS hind => card_le_of_indep hCh hS hind)
  have hsucc : ∀ j, cntRev G (C.cor Ch) (Ch.card + 1) (j + 1) = C.vecAt Ch ∅ j := by
    intro j
    rw [cntRev_succ_succ, vecAt_empty_avoid]
  have e0 := vecAt_peel_root hCh hv hvCh 0
  have e1 := vecAt_peel_root hCh hv hvCh 1
  have e2 := vecAt_peel_root hCh hv hvCh 2
  have e3 := vecAt_peel_root hCh hv hvCh 3
  have e4 := vecAt_peel_root hCh hv hvCh 4
  rw [hzero] at e0
  rw [hsucc 0] at e1
  rw [hsucc 1] at e2
  rw [hsucc 2] at e3
  rw [hsucc 3] at e4
  simp only [vec, add, shift, TopFive.mk.injEq]
  refine ⟨by omega, by omega, by omega, by omega, by omega⟩

end CoronaData

end B1Zero

end MatchingBag
