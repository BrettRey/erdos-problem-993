import RequestProject.ClanAudit.BlockWeight

/-!
# The explicit global partition of the degree-`N` clan-map space

This is the global bookkeeping obligation of the adjacent two-hub target.  For the tree
`dgraph m a n b` — two adjacent hubs carrying arbitrary finite ordered families of
positive-length pendant paths — we build an *explicit* partition of the finite set of all
clan maps of total order `N` into

1. four-state blocks, obtained by applying the canonical local transformation at both
   active hubs;
2. two-state blocks, obtained by applying it at exactly one active hub;
3. singletons.

The partition is realised as the fibres of an explicit idempotent, order-preserving
representative map `rep`.  A map or a partner therefore occurs in *exactly one* block by
construction, which settles the collision audit: the single-hub and double-hub strata are
separated because `repSide` is computed side by side and the transformation lands in the
hub-inactive stratum, disjoint from its own source stratum.

* `fiber_eq_image` — each block is exactly the explicit set of 4/2/1 maps;
* `card_block_cases` — the blocks have exactly 4, 2 or 1 elements;
* `blocks_pairwise_disjoint`, `blocks_cover` — disjointness and exhaustion;
* `sum_repSide`, `sum_rep` — the total order is preserved, so each block lies inside a
  single degree;
* `isCU_totalW` — the total normalized weight of a degree is centrally unimodal;
* `sum_normalizedTwoRowCoeff_nonneg` — **the final degreewise theorem**.
-/

namespace ClanAudit

open Finset LaurentPolynomial SimpleGraph
open scoped Classical

/-! ### The representative of one side -/

variable {k : ℕ} {len : Fin k → ℕ}

/-- The image stratum of the canonical transformation. -/
def ImgSide (s : SpiderV k len → ℕ) : Prop := ∃ t, ActiveSide t ∧ transf t = s

/-- **The source and the image strata are disjoint**: the transformation switches the hub
off, while an active side has an active hub. -/
theorem not_activeSide_of_imgSide {s : SpiderV k len → ℕ} (h : ImgSide s) :
    ¬ ActiveSide s := by
  obtain ⟨t, ht, rfl⟩ := h
  intro hc
  have h1 := hc.hub
  rw [transf_none] at h1
  omega

theorem transf_ne_self {s : SpiderV k len → ℕ} (h : ActiveSide s) : transf s ≠ s := by
  intro hc
  have h1 : transf s none = 0 := transf_none s
  rw [hc, h.hub] at h1
  omega

/-- The canonical representative of the block of one side: an active side represents
itself, an image is represented by its (unique) source, and anything else represents
itself. -/
noncomputable def repSide (s : SpiderV k len → ℕ) : SpiderV k len → ℕ :=
  if ActiveSide s then s else if h : ImgSide s then h.choose else s

theorem repSide_of_active {s : SpiderV k len → ℕ} (h : ActiveSide s) : repSide s = s := by
  rw [repSide, if_pos h]

theorem repSide_of_imgSide {s : SpiderV k len → ℕ} (h1 : ¬ ActiveSide s) (h2 : ImgSide s) :
    ActiveSide (repSide s) ∧ transf (repSide s) = s := by
  rw [repSide, if_neg h1, dif_pos h2]
  exact h2.choose_spec

theorem repSide_of_neither {s : SpiderV k len → ℕ} (h1 : ¬ ActiveSide s) (h2 : ¬ ImgSide s) :
    repSide s = s := by
  rw [repSide, if_neg h1, dif_neg h2]

theorem repSide_spec (s : SpiderV k len → ℕ) :
    (ActiveSide (repSide s) ∧ transf (repSide s) = s) ∨ repSide s = s := by
  by_cases h1 : ActiveSide s
  · exact Or.inr (repSide_of_active h1)
  by_cases h2 : ImgSide s
  · exact Or.inl (repSide_of_imgSide h1 h2)
  · exact Or.inr (repSide_of_neither h1 h2)

/-- **The image of an active side is represented by that side.**  This is where the
injectivity of the published local map is consumed. -/
theorem repSide_transf {s : SpiderV k len → ℕ} (h : ActiveSide s) : repSide (transf s) = s := by
  have h1 : ¬ ActiveSide (transf s) := not_activeSide_of_imgSide ⟨s, h, rfl⟩
  obtain ⟨hact, heq⟩ := repSide_of_imgSide h1 ⟨s, h, rfl⟩
  exact transf_injective hact h heq

theorem repSide_idem (s : SpiderV k len → ℕ) : repSide (repSide s) = repSide s := by
  rcases repSide_spec s with ⟨h, _⟩ | h
  · exact repSide_of_active h
  · rw [h]
    exact h

/-- **The representative has the same total order.**  This is where the order preservation
of the published local map is consumed. -/
theorem sum_repSide (s : SpiderV k len → ℕ) : ∑ x, repSide s x = ∑ x, s x := by
  rcases repSide_spec s with ⟨h, he⟩ | h
  · calc ∑ x, repSide s x = ∑ x, transf (repSide s) x := (sum_transf h).symm
      _ = ∑ x, s x := by rw [he]
  · rw [h]

/-! ### The block of one side -/

/-- The block of one side: the pair `{s, transf s}` at an active side, the singleton
`{s}` otherwise. -/
noncomputable def sideBlock (s : SpiderV k len → ℕ) : Finset (SpiderV k len → ℕ) :=
  if ActiveSide s then {s, transf s} else {s}

theorem mem_sideBlock {s t : SpiderV k len → ℕ} (hs : repSide s = s) :
    t ∈ sideBlock s ↔ repSide t = s := by
  by_cases h : ActiveSide s
  · rw [sideBlock, if_pos h]
    simp only [Finset.mem_insert, Finset.mem_singleton]
    constructor
    · rintro (rfl | rfl)
      · exact hs
      · exact repSide_transf h
    · intro ht
      rcases repSide_spec t with ⟨ha, he⟩ | he
      · rw [ht] at he
        exact Or.inr he.symm
      · exact Or.inl (he.symm.trans ht)
  · rw [sideBlock, if_neg h]
    simp only [Finset.mem_singleton]
    constructor
    · rintro rfl
      exact hs
    · intro ht
      rcases repSide_spec t with ⟨ha, he⟩ | he
      · exact absurd (ht ▸ ha) h
      · exact he.symm.trans ht

theorem sum_sideBlock_active {s : SpiderV k len → ℕ} (h : ActiveSide s)
    (f : (SpiderV k len → ℕ) → LaurentPolynomial ℚ) :
    ∑ t ∈ sideBlock s, f t = f s + f (transf s) := by
  rw [sideBlock, if_pos h]
  exact Finset.sum_pair (Ne.symm (transf_ne_self h))

theorem sum_sideBlock_not_active {s : SpiderV k len → ℕ} (h : ¬ ActiveSide s)
    (f : (SpiderV k len → ℕ) → LaurentPolynomial ℚ) :
    ∑ t ∈ sideBlock s, f t = f s := by
  rw [sideBlock, if_neg h]
  exact Finset.sum_singleton _ _

theorem card_sideBlock_active {s : SpiderV k len → ℕ} (h : ActiveSide s) :
    (sideBlock s).card = 2 := by
  rw [sideBlock, if_pos h]
  exact Finset.card_pair (Ne.symm (transf_ne_self h))

theorem card_sideBlock_not_active {s : SpiderV k len → ℕ} (h : ¬ ActiveSide s) :
    (sideBlock s).card = 1 := by
  rw [sideBlock, if_neg h]
  exact Finset.card_singleton _

/-! ### The global representative and the degree-`N` map space -/

variable {m n : ℕ} {a : Fin m → ℕ} {b : Fin n → ℕ}

/-- The canonical representative of the block of a clan map of the two-hub tree: the two
sides are normalized independently. -/
noncomputable def rep (α : DV m a n b → ℕ) : DV m a n b → ℕ :=
  Sum.elim (repSide (fun x => α (Sum.inl x))) (repSide (fun y => α (Sum.inr y)))

theorem sum_elim_comp (α : DV m a n b → ℕ) :
    Sum.elim (fun x => α (Sum.inl x)) (fun y => α (Sum.inr y)) = α := by
  funext z
  rcases z with z | z <;> rfl

theorem rep_eq_iff {α β : DV m a n b → ℕ} :
    rep α = β ↔ repSide (fun x => α (Sum.inl x)) = (fun x => β (Sum.inl x))
      ∧ repSide (fun y => α (Sum.inr y)) = (fun y => β (Sum.inr y)) := by
  constructor
  · intro h
    exact ⟨funext fun x => congrFun h (Sum.inl x), funext fun y => congrFun h (Sum.inr y)⟩
  · rintro ⟨h1, h2⟩
    funext z
    rcases z with z | z
    · exact congrFun h1 z
    · exact congrFun h2 z

theorem rep_idem (α : DV m a n b → ℕ) : rep (rep α) = rep α := by
  show Sum.elim (repSide (repSide (fun x => α (Sum.inl x))))
      (repSide (repSide (fun y => α (Sum.inr y)))) = rep α
  rw [repSide_idem, repSide_idem]
  rfl

theorem sum_rep (α : DV m a n b → ℕ) : ∑ x, rep α x = ∑ x, α x := by
  rw [Fintype.sum_sum_type, Fintype.sum_sum_type]
  exact congrArg₂ (· + ·) (sum_repSide _) (sum_repSide _)

/-- The finite set of all clan maps of the adjacent two-hub tree of total order `N`. -/
noncomputable def mapsOfOrder (m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) (N : ℕ) :
    Finset (DV m a n b → ℕ) :=
  (Fintype.piFinset fun _ => Finset.range (N + 1)).filter (fun α => ∑ x, α x = N)

theorem mem_mapsOfOrder {N : ℕ} {α : DV m a n b → ℕ} :
    α ∈ mapsOfOrder m a n b N ↔ ∑ x, α x = N := by
  rw [mapsOfOrder, Finset.mem_filter, Fintype.mem_piFinset]
  refine ⟨fun h => h.2, fun h => ⟨fun x => Finset.mem_range.mpr ?_, h⟩⟩
  have hx : α x ≤ ∑ y, α y :=
    Finset.single_le_sum (f := α) (fun i _ => Nat.zero_le _) (Finset.mem_univ x)
  omega

theorem rep_mem_mapsOfOrder {N : ℕ} {α : DV m a n b → ℕ} (h : α ∈ mapsOfOrder m a n b N) :
    rep α ∈ mapsOfOrder m a n b N :=
  mem_mapsOfOrder.mpr ((sum_rep α).trans (mem_mapsOfOrder.mp h))

/-! ### The blocks -/

/-- The block of the representative `β`: the fibre of `rep` over `β`. -/
noncomputable def block (m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) (N : ℕ)
    (β : DV m a n b → ℕ) : Finset (DV m a n b → ℕ) :=
  (mapsOfOrder m a n b N).filter (fun α => rep α = β)

/-- **Disjointness.**  Distinct representatives have disjoint blocks. -/
theorem blocks_pairwise_disjoint (N : ℕ) {β γ : DV m a n b → ℕ} (h : β ≠ γ) :
    Disjoint (block m a n b N β) (block m a n b N γ) := by
  rw [Finset.disjoint_left]
  intro α hβ hγ
  rw [block, Finset.mem_filter] at hβ hγ
  exact h (hβ.2.symm.trans hγ.2)

/-- **Exhaustion.**  Every map of order `N` lies in the block of its representative, and
that representative is itself a map of order `N`. -/
theorem blocks_cover {N : ℕ} {α : DV m a n b → ℕ} (h : α ∈ mapsOfOrder m a n b N) :
    α ∈ block m a n b N (rep α) ∧ rep α ∈ mapsOfOrder m a n b N :=
  ⟨Finset.mem_filter.mpr ⟨h, rfl⟩, rep_mem_mapsOfOrder h⟩

/-- **A map lies in exactly one block.**  The blocks are the fibres of a function, so
this is automatic — which is precisely the point of building the partition as a fibre
decomposition rather than as a list of postulated pairings. -/
theorem block_unique {N : ℕ} {α β γ : DV m a n b → ℕ} (hβ : α ∈ block m a n b N β)
    (hγ : α ∈ block m a n b N γ) : β = γ := by
  rw [block, Finset.mem_filter] at hβ hγ
  exact hβ.2.symm.trans hγ.2

/-- **The collision audit.**  The image of the transformation at an active side is *not*
itself a source: the two strata are separated because the transformation switches the hub
off (`not_activeSide_of_imgSide`).  Consequently a map whose left side is the image of an
active side lands in the *same* block as the untransformed map — it is not the start of a
second, competing block.  The same holds on the right. -/
theorem collision_resolved_left {u : SpiderV m a → ℕ} (h : ActiveSide u) (v : SpiderV n b → ℕ) :
    rep (Sum.elim (transf u) v) = rep (Sum.elim u v) := by
  show Sum.elim (repSide (transf u)) (repSide v) = Sum.elim (repSide u) (repSide v)
  rw [repSide_transf h, repSide_of_active h]

theorem collision_resolved_right (u : SpiderV m a → ℕ) {v : SpiderV n b → ℕ}
    (h : ActiveSide v) : rep (Sum.elim u (transf v)) = rep (Sum.elim u v) := by
  show Sum.elim (repSide u) (repSide (transf v)) = Sum.elim (repSide u) (repSide v)
  rw [repSide_transf h, repSide_of_active h]

/-- **A block is empty unless its index is a representative.** -/
theorem block_eq_empty_of_not_fixed {N : ℕ} {β : DV m a n b → ℕ} (h : rep β ≠ β) :
    block m a n b N β = ∅ := by
  rw [block, Finset.filter_eq_empty_iff]
  intro α _ hc
  exact h (by rw [← hc, rep_idem])

theorem repSide_left_of_fixed {β : DV m a n b → ℕ} (hβ : rep β = β) :
    repSide (fun x => β (Sum.inl x)) = fun x => β (Sum.inl x) := (rep_eq_iff.mp hβ).1

theorem repSide_right_of_fixed {β : DV m a n b → ℕ} (hβ : rep β = β) :
    repSide (fun y => β (Sum.inr y)) = fun y => β (Sum.inr y) := (rep_eq_iff.mp hβ).2

theorem sum_elim_injective :
    Function.Injective (fun p : (SpiderV m a → ℕ) × (SpiderV n b → ℕ) => Sum.elim p.1 p.2) := by
  rintro ⟨p1, p2⟩ ⟨q1, q2⟩ h
  have h1 : p1 = q1 := funext fun x => congrFun h (Sum.inl x)
  have h2 : p2 = q2 := funext fun y => congrFun h (Sum.inr y)
  rw [h1, h2]

/-- **The explicit description of a block.**  The block of a representative `β` is exactly
the set of maps obtained by independently keeping or transforming each of its two sides:
four maps if both hubs are active, two if exactly one is, one otherwise. -/
theorem fiber_eq_image {N : ℕ} {β : DV m a n b → ℕ} (hβ : rep β = β)
    (hmem : β ∈ mapsOfOrder m a n b N) :
    block m a n b N β
      = ((sideBlock (fun x => β (Sum.inl x))) ×ˢ (sideBlock (fun y => β (Sum.inr y)))).image
          (fun p => Sum.elim p.1 p.2) := by
  have hu := repSide_left_of_fixed hβ
  have hv := repSide_right_of_fixed hβ
  have hsum : ∑ x, β x = N := mem_mapsOfOrder.mp hmem
  ext α
  simp only [block, Finset.mem_filter, Finset.mem_image, Finset.mem_product]
  constructor
  · rintro ⟨hα, hrep⟩
    obtain ⟨h1, h2⟩ := rep_eq_iff.mp hrep
    exact ⟨(fun x => α (Sum.inl x), fun y => α (Sum.inr y)),
      ⟨(mem_sideBlock hu).mpr h1, (mem_sideBlock hv).mpr h2⟩, sum_elim_comp α⟩
  · rintro ⟨⟨p1, p2⟩, ⟨hp1, hp2⟩, rfl⟩
    have e1 : repSide p1 = fun x => β (Sum.inl x) := (mem_sideBlock hu).mp hp1
    have e2 : repSide p2 = fun y => β (Sum.inr y) := (mem_sideBlock hv).mp hp2
    refine ⟨mem_mapsOfOrder.mpr ?_, rep_eq_iff.mpr ⟨e1, e2⟩⟩
    have s1 : ∑ x, p1 x = ∑ x, β (Sum.inl x) := by rw [← sum_repSide p1, e1]
    have s2 : ∑ y, p2 y = ∑ y, β (Sum.inr y) := by rw [← sum_repSide p2, e2]
    rw [Fintype.sum_sum_type (f := Sum.elim p1 p2)]
    rw [Fintype.sum_sum_type (f := β)] at hsum
    simp only [Sum.elim_inl, Sum.elim_inr]
    rw [s1, s2]
    exact hsum

/-- **The block sizes.**  Each block has exactly `4`, `2` or `1` element. -/
theorem card_block {N : ℕ} {β : DV m a n b → ℕ} (hβ : rep β = β)
    (hmem : β ∈ mapsOfOrder m a n b N) :
    (block m a n b N β).card
      = (sideBlock (fun x => β (Sum.inl x))).card
        * (sideBlock (fun y => β (Sum.inr y))).card := by
  rw [fiber_eq_image hβ hmem,
    Finset.card_image_of_injective _ sum_elim_injective, Finset.card_product]

theorem card_block_cases {N : ℕ} {β : DV m a n b → ℕ} (hβ : rep β = β)
    (hmem : β ∈ mapsOfOrder m a n b N) :
    (block m a n b N β).card
      = if ActiveSide (fun x => β (Sum.inl x)) then
          (if ActiveSide (fun y => β (Sum.inr y)) then 4 else 2)
        else (if ActiveSide (fun y => β (Sum.inr y)) then 2 else 1) := by
  rw [card_block hβ hmem]
  by_cases hu : ActiveSide (fun x => β (Sum.inl x))
  · rw [card_sideBlock_active hu, if_pos hu]
    by_cases hv : ActiveSide (fun y => β (Sum.inr y))
    · rw [card_sideBlock_active hv, if_pos hv]
    · rw [card_sideBlock_not_active hv, if_neg hv]
  · rw [card_sideBlock_not_active hu, if_neg hu]
    by_cases hv : ActiveSide (fun y => β (Sum.inr y))
    · rw [card_sideBlock_active hv, if_pos hv]
    · rw [card_sideBlock_not_active hv, if_neg hv]

/-! ### The weight of a block -/

theorem sum_block {N : ℕ} {β : DV m a n b → ℕ} (hβ : rep β = β)
    (hmem : β ∈ mapsOfOrder m a n b N) (f : (DV m a n b → ℕ) → LaurentPolynomial ℚ) :
    ∑ α ∈ block m a n b N β, f α
      = ∑ s ∈ sideBlock (fun x => β (Sum.inl x)),
          ∑ t ∈ sideBlock (fun y => β (Sum.inr y)), f (Sum.elim s t) := by
  rw [fiber_eq_image hβ hmem,
    Finset.sum_image (fun p _ q _ h => sum_elim_injective h), Finset.sum_product]

/-- **Blockwise central unimodality.**  Every block of the partition has a centrally
unimodal normalized weight: the four-state blocks give `F(r,s;c,d)` with `r,s ≥ 1` and
`c,d ≥ 1`, the two-state blocks give `A(r,c)` with `r ≥ 1` and `c ≥ 1`, and the singletons
give products of good factors. -/
theorem isCU_sum_block {N : ℕ} {β : DV m a n b → ℕ} (hβ : rep β = β)
    (hmem : β ∈ mapsOfOrder m a n b N) :
    IsCU (∑ α ∈ block m a n b N β, Wpoly (dgraph m a n b) α) := by
  rw [sum_block hβ hmem]
  by_cases hu : ActiveSide (fun x => β (Sum.inl x))
  · by_cases hv : ActiveSide (fun y => β (Sum.inr y))
    · rw [sum_sideBlock_active hu, sum_sideBlock_active hv, sum_sideBlock_active hv]
      exact blockSum_both_active hu hv
    · rw [sum_sideBlock_active hu, sum_sideBlock_not_active hv, sum_sideBlock_not_active hv]
      exact blockSum_left_active hu hv
  · by_cases hv : ActiveSide (fun y => β (Sum.inr y))
    · rw [sum_sideBlock_not_active hu, sum_sideBlock_active hv]
      exact blockSum_right_active hu hv
    · rw [sum_sideBlock_not_active hu, sum_sideBlock_not_active hv]
      exact (blockSum_neither_active hu hv).isCU

/-! ### The final degreewise theorem -/

/-- The total normalized weight of all clan maps of total order `N`. -/
noncomputable def totalW (m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) (N : ℕ) :
    LaurentPolynomial ℚ :=
  ∑ α ∈ mapsOfOrder m a n b N, Wpoly (dgraph m a n b) α

/-- **The total weight of a degree is centrally unimodal.**  This is the partition
theorem combined with the blockwise weight theorem. -/
theorem isCU_totalW (m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) (N : ℕ) :
    IsCU (totalW m a n b N) := by
  have hmaps : ∀ α ∈ mapsOfOrder m a n b N, rep α ∈ mapsOfOrder m a n b N :=
    fun _ h => rep_mem_mapsOfOrder h
  rw [totalW, ← Finset.sum_fiberwise_of_maps_to hmaps]
  refine Finset.sum_induction _ IsCU (fun _ _ => IsCU.add) isCU_zero ?_
  intro β hβmem
  by_cases hfix : rep β = β
  · exact isCU_sum_block hfix hβmem
  · rw [show (mapsOfOrder m a n b N).filter (fun α => rep α = β) = ∅ from
      block_eq_empty_of_not_fixed hfix, Finset.sum_empty]
    exact isCU_zero

/-- **The final degreewise theorem.**  Summed over *all* clan maps of the adjacent
two-hub tree of total order `N = k + l`, every normalized two-row coefficient is
nonnegative. -/
theorem sum_normalizedTwoRowCoeff_nonneg (m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ)
    (N k l : ℕ) (hl : 1 ≤ l) (hkl : l ≤ k) (hN : k + l = N) :
    0 ≤ ∑ α ∈ mapsOfOrder m a n b N, normalizedTwoRowCoeff (dgraph m a n b) α k l := by
  have hCU := isCU_totalW m a n b N
  have hval : ∑ α ∈ mapsOfOrder m a n b N, normalizedTwoRowCoeff (dgraph m a n b) α k l
      = coeffL (totalW m a n b N) ((k : ℤ) - l)
        - coeffL (totalW m a n b N) ((k : ℤ) - l + 2) := by
    rw [Finset.sum_congr rfl (fun α hα => normalizedTwoRowCoeff_eq (dgraph m a n b) α k l hl
      (by rw [mem_mapsOfOrder.mp hα]; omega)), Finset.sum_sub_distrib, totalW,
      coeffL_sum, coeffL_sum]
  rw [hval, sub_nonneg]
  exact hCU.decr _ (by omega)

end ClanAudit
