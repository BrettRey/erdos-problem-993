import RequestProject.ClanAudit.ConnectorZero

/-!
# The explicit global partition of the degree-`N` clan-map space of the connector tree

Exactly as for the adjacent two-hub tree, the finite set of clan maps of the connector
tree `connGraph t m a n b` of total order `N` is partitioned into the fibres of an
explicit idempotent, order-preserving representative map `connRep`, which normalizes the
two hub sides independently and leaves the connector untouched.

* `conn_blocks_pairwise_disjoint`, `conn_blocks_cover`, `conn_block_unique` —
  disjointness, exhaustion and uniqueness;
* `conn_fiber_eq_image`, `card_conn_block` — each block is the explicit set of 4/2/1
  maps;
* `sum_connRep` — the total order is preserved, so a block lies inside a single degree;
* `isCU_connTotalW` — the total normalized weight of a degree is centrally unimodal;
* `conn_sum_normalizedTwoRowCoeff_nonneg` — **the connector degreewise theorem**.
-/

namespace ClanAudit

open Finset LaurentPolynomial SimpleGraph
open scoped Classical

variable {t m n : ℕ} {a : Fin m → ℕ} {b : Fin n → ℕ}

/-- The canonical representative of the block of a clan map of the connector tree: the
two hub sides are normalized independently, the connector is untouched. -/
noncomputable def connRep (α : ConnV t m a n b → ℕ) : ConnV t m a n b → ℕ :=
  Sum.elim (repSide (fun x => α (Sum.inl x)))
    (Sum.elim (fun j => α (Sum.inr (Sum.inl j)))
      (repSide (fun y => α (Sum.inr (Sum.inr y)))))

theorem connRep_eq_iff {α β : ConnV t m a n b → ℕ} :
    connRep α = β ↔ (repSide (fun x => α (Sum.inl x)) = fun x => β (Sum.inl x))
      ∧ (∀ j : Fin t, α (Sum.inr (Sum.inl j)) = β (Sum.inr (Sum.inl j)))
      ∧ (repSide (fun y => α (Sum.inr (Sum.inr y))) = fun y => β (Sum.inr (Sum.inr y))) := by
  constructor
  · intro h
    exact ⟨funext fun x => congrFun h (Sum.inl x),
      fun j => congrFun h (Sum.inr (Sum.inl j)),
      funext fun y => congrFun h (Sum.inr (Sum.inr y))⟩
  · rintro ⟨h1, h2, h3⟩
    funext z
    rcases z with z | z | z
    · exact congrFun h1 z
    · exact h2 z
    · exact congrFun h3 z

theorem connRep_idem (α : ConnV t m a n b → ℕ) : connRep (connRep α) = connRep α := by
  show Sum.elim (repSide (repSide (fun x => α (Sum.inl x))))
      (Sum.elim (fun j => α (Sum.inr (Sum.inl j)))
        (repSide (repSide (fun y => α (Sum.inr (Sum.inr y)))))) = connRep α
  rw [repSide_idem, repSide_idem]
  rfl

theorem sum_connRep (α : ConnV t m a n b → ℕ) : ∑ x, connRep α x = ∑ x, α x := by
  rw [Fintype.sum_sum_type, Fintype.sum_sum_type, Fintype.sum_sum_type, Fintype.sum_sum_type]
  exact congrArg₂ (· + ·) (sum_repSide _) (congrArg₂ (· + ·) rfl (sum_repSide _))

/-- The finite set of all clan maps of the connector tree of total order `N`. -/
noncomputable def connMapsOfOrder (t m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) (N : ℕ) :
    Finset (ConnV t m a n b → ℕ) :=
  (Fintype.piFinset fun _ => Finset.range (N + 1)).filter (fun α => ∑ x, α x = N)

theorem mem_connMapsOfOrder {N : ℕ} {α : ConnV t m a n b → ℕ} :
    α ∈ connMapsOfOrder t m a n b N ↔ ∑ x, α x = N := by
  rw [connMapsOfOrder, Finset.mem_filter, Fintype.mem_piFinset]
  refine ⟨fun h => h.2, fun h => ⟨fun x => Finset.mem_range.mpr ?_, h⟩⟩
  have hx : α x ≤ ∑ y, α y :=
    Finset.single_le_sum (f := α) (fun i _ => Nat.zero_le _) (Finset.mem_univ x)
  omega

theorem connRep_mem_connMapsOfOrder {N : ℕ} {α : ConnV t m a n b → ℕ}
    (h : α ∈ connMapsOfOrder t m a n b N) : connRep α ∈ connMapsOfOrder t m a n b N :=
  mem_connMapsOfOrder.mpr ((sum_connRep α).trans (mem_connMapsOfOrder.mp h))

/-- The block of the representative `β`: the fibre of `connRep` over `β`. -/
noncomputable def connBlock (t m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) (N : ℕ)
    (β : ConnV t m a n b → ℕ) : Finset (ConnV t m a n b → ℕ) :=
  (connMapsOfOrder t m a n b N).filter (fun α => connRep α = β)

/-- **Disjointness.** -/
theorem conn_blocks_pairwise_disjoint (N : ℕ) {β γ : ConnV t m a n b → ℕ} (h : β ≠ γ) :
    Disjoint (connBlock t m a n b N β) (connBlock t m a n b N γ) := by
  rw [Finset.disjoint_left]
  intro α hβ hγ
  rw [connBlock, Finset.mem_filter] at hβ hγ
  exact h (hβ.2.symm.trans hγ.2)

/-- **Exhaustion.** -/
theorem conn_blocks_cover {N : ℕ} {α : ConnV t m a n b → ℕ}
    (h : α ∈ connMapsOfOrder t m a n b N) :
    α ∈ connBlock t m a n b N (connRep α) ∧ connRep α ∈ connMapsOfOrder t m a n b N :=
  ⟨Finset.mem_filter.mpr ⟨h, rfl⟩, connRep_mem_connMapsOfOrder h⟩

/-- **A map lies in exactly one block.** -/
theorem conn_block_unique {N : ℕ} {α β γ : ConnV t m a n b → ℕ}
    (hβ : α ∈ connBlock t m a n b N β) (hγ : α ∈ connBlock t m a n b N γ) : β = γ := by
  rw [connBlock, Finset.mem_filter] at hβ hγ
  exact hβ.2.symm.trans hγ.2

theorem conn_block_eq_empty_of_not_fixed {N : ℕ} {β : ConnV t m a n b → ℕ}
    (h : connRep β ≠ β) : connBlock t m a n b N β = ∅ := by
  rw [connBlock, Finset.filter_eq_empty_iff]
  intro α _ hc
  exact h (by rw [← hc, connRep_idem])

theorem conn_sum_elim_injective (w : Fin t → ℕ) :
    Function.Injective
      (fun p : (SpiderV m a → ℕ) × (SpiderV n b → ℕ) => Sum.elim p.1 (Sum.elim w p.2)) := by
  rintro ⟨p1, p2⟩ ⟨q1, q2⟩ h
  have h1 : p1 = q1 := funext fun x => congrFun h (Sum.inl x)
  have h2 : p2 = q2 := funext fun y => congrFun h (Sum.inr (Sum.inr y))
  rw [h1, h2]

/-- **The explicit description of a block.**  The block of a representative `β` is
exactly the set of maps obtained by independently keeping or transforming each of its
two hub sides; the connector is copied unchanged. -/
theorem conn_fiber_eq_image {N : ℕ} {β : ConnV t m a n b → ℕ} (hβ : connRep β = β)
    (hmem : β ∈ connMapsOfOrder t m a n b N) :
    connBlock t m a n b N β
      = ((sideBlock (fun x => β (Sum.inl x)))
          ×ˢ (sideBlock (fun y => β (Sum.inr (Sum.inr y))))).image
          (fun p => Sum.elim p.1 (Sum.elim (fun j => β (Sum.inr (Sum.inl j))) p.2)) := by
  obtain ⟨hu, hc, hv⟩ := connRep_eq_iff.mp hβ
  have hsum : ∑ x, β x = N := mem_connMapsOfOrder.mp hmem
  ext α
  simp only [connBlock, Finset.mem_filter, Finset.mem_image, Finset.mem_product]
  constructor
  · rintro ⟨hα, hrep⟩
    obtain ⟨h1, h2, h3⟩ := connRep_eq_iff.mp hrep
    refine ⟨(fun x => α (Sum.inl x), fun y => α (Sum.inr (Sum.inr y))),
      ⟨(mem_sideBlock hu).mpr h1, (mem_sideBlock hv).mpr h3⟩, ?_⟩
    funext z
    rcases z with z | z | z
    · rfl
    · exact (h2 z).symm
    · rfl
  · rintro ⟨⟨p1, p2⟩, ⟨hp1, hp2⟩, rfl⟩
    have e1 : repSide p1 = fun x => β (Sum.inl x) := (mem_sideBlock hu).mp hp1
    have e2 : repSide p2 = fun y => β (Sum.inr (Sum.inr y)) := (mem_sideBlock hv).mp hp2
    refine ⟨mem_connMapsOfOrder.mpr ?_, connRep_eq_iff.mpr ⟨e1, fun j => rfl, e2⟩⟩
    have s1 : ∑ x, p1 x = ∑ x, β (Sum.inl x) := by rw [← sum_repSide p1, e1]
    have s2 : ∑ y, p2 y = ∑ y, β (Sum.inr (Sum.inr y)) := by rw [← sum_repSide p2, e2]
    rw [Fintype.sum_sum_type, Fintype.sum_sum_type]
    rw [Fintype.sum_sum_type (f := β), Fintype.sum_sum_type
      (f := fun z => β (Sum.inr z))] at hsum
    simp only [Sum.elim_inl, Sum.elim_inr]
    rw [s1, s2]
    exact hsum

/-- **The block sizes.**  Each block has exactly `4`, `2` or `1` element. -/
theorem card_conn_block {N : ℕ} {β : ConnV t m a n b → ℕ} (hβ : connRep β = β)
    (hmem : β ∈ connMapsOfOrder t m a n b N) :
    (connBlock t m a n b N β).card
      = (sideBlock (fun x => β (Sum.inl x))).card
        * (sideBlock (fun y => β (Sum.inr (Sum.inr y)))).card := by
  rw [conn_fiber_eq_image hβ hmem,
    Finset.card_image_of_injective _ (conn_sum_elim_injective _), Finset.card_product]

theorem card_conn_block_cases {N : ℕ} {β : ConnV t m a n b → ℕ} (hβ : connRep β = β)
    (hmem : β ∈ connMapsOfOrder t m a n b N) :
    (connBlock t m a n b N β).card
      = if ActiveSide (fun x => β (Sum.inl x)) then
          (if ActiveSide (fun y => β (Sum.inr (Sum.inr y))) then 4 else 2)
        else (if ActiveSide (fun y => β (Sum.inr (Sum.inr y))) then 2 else 1) := by
  rw [card_conn_block hβ hmem]
  by_cases hu : ActiveSide (fun x => β (Sum.inl x))
  · rw [card_sideBlock_active hu, if_pos hu]
    by_cases hv : ActiveSide (fun y => β (Sum.inr (Sum.inr y)))
    · rw [card_sideBlock_active hv, if_pos hv]
    · rw [card_sideBlock_not_active hv, if_neg hv]
  · rw [card_sideBlock_not_active hu, if_neg hu]
    by_cases hv : ActiveSide (fun y => β (Sum.inr (Sum.inr y)))
    · rw [card_sideBlock_active hv, if_pos hv]
    · rw [card_sideBlock_not_active hv, if_neg hv]

theorem conn_sum_block {N : ℕ} {β : ConnV t m a n b → ℕ} (hβ : connRep β = β)
    (hmem : β ∈ connMapsOfOrder t m a n b N)
    (f : (ConnV t m a n b → ℕ) → LaurentPolynomial ℚ) :
    ∑ α ∈ connBlock t m a n b N β, f α
      = ∑ s ∈ sideBlock (fun x => β (Sum.inl x)),
          ∑ r ∈ sideBlock (fun y => β (Sum.inr (Sum.inr y))),
            f (Sum.elim s (Sum.elim (fun j => β (Sum.inr (Sum.inl j))) r)) := by
  rw [conn_fiber_eq_image hβ hmem,
    Finset.sum_image (fun p _ q _ h => conn_sum_elim_injective _ h), Finset.sum_product]

/-- **Blockwise central unimodality for the connector tree.** -/
theorem isCU_conn_sum_block {N : ℕ} {β : ConnV t m a n b → ℕ} (hβ : connRep β = β)
    (hmem : β ∈ connMapsOfOrder t m a n b N) :
    IsCU (∑ α ∈ connBlock t m a n b N β, Wpoly (connGraph t m a n b) α) := by
  rw [conn_sum_block hβ hmem]
  exact connBlockSum_isCU_all _ _ _

/-- The total normalized weight of all clan maps of total order `N`. -/
noncomputable def connTotalW (t m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) (N : ℕ) :
    LaurentPolynomial ℚ :=
  ∑ α ∈ connMapsOfOrder t m a n b N, Wpoly (connGraph t m a n b) α

/-- **The total weight of a degree is centrally unimodal.** -/
theorem isCU_connTotalW (t m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) (N : ℕ) :
    IsCU (connTotalW t m a n b N) := by
  have hmaps : ∀ α ∈ connMapsOfOrder t m a n b N, connRep α ∈ connMapsOfOrder t m a n b N :=
    fun _ h => connRep_mem_connMapsOfOrder h
  rw [connTotalW, ← Finset.sum_fiberwise_of_maps_to hmaps]
  refine Finset.sum_induction _ IsCU (fun _ _ => IsCU.add) isCU_zero ?_
  intro β hβmem
  by_cases hfix : connRep β = β
  · exact isCU_conn_sum_block hfix hβmem
  · rw [show (connMapsOfOrder t m a n b N).filter (fun α => connRep α = β) = ∅ from
      conn_block_eq_empty_of_not_fixed hfix, Finset.sum_empty]
    exact isCU_zero

/-- **The connector degreewise theorem.**  Summed over *all* clan maps of the connector
tree of total order `N = k + l`, every normalized two-row coefficient is nonnegative. -/
theorem conn_sum_normalizedTwoRowCoeff_nonneg (t m : ℕ) (a : Fin m → ℕ) (n : ℕ)
    (b : Fin n → ℕ) (N k l : ℕ) (hl : 1 ≤ l) (hkl : l ≤ k) (hN : k + l = N) :
    0 ≤ ∑ α ∈ connMapsOfOrder t m a n b N,
      normalizedTwoRowCoeff (connGraph t m a n b) α k l := by
  have hCU := isCU_connTotalW t m a n b N
  have hval : ∑ α ∈ connMapsOfOrder t m a n b N,
      normalizedTwoRowCoeff (connGraph t m a n b) α k l
      = coeffL (connTotalW t m a n b N) ((k : ℤ) - l)
        - coeffL (connTotalW t m a n b N) ((k : ℤ) - l + 2) := by
    rw [Finset.sum_congr rfl (fun α hα => normalizedTwoRowCoeff_eq (connGraph t m a n b) α k l
      hl (by rw [mem_connMapsOfOrder.mp hα]; omega)), Finset.sum_sub_distrib, connTotalW,
      coeffL_sum, coeffL_sum]
  rw [hval, sub_nonneg]
  exact hCU.decr _ (by omega)

end ClanAudit
