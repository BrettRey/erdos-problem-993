import RequestProject.ClanAudit.Bridge
import RequestProject.ClanAudit.GlobalPartition

/-!
# Log-concavity of the independence polynomial of a spider (and of a path)

The one-hub case of the clan argument.  The degree-`N` clan-map space of a spider is
partitioned into the fibres of the idempotent, order-preserving representative map
`repSide`; each fibre is a pair `{s, transf s}` at an active side and a singleton
otherwise, with normalized weight `A(p-1, c) * ∏ tails` (respectively a good weight).

Since a path with `N ≥ 1` vertices is the spider with one arm of length `N - 1`, this
also covers all paths, that is, all trees with no branch vertex.
-/

namespace ClanAudit

open Finset LaurentPolynomial SimpleGraph
open scoped Classical

variable {k : ℕ} {len : Fin k → ℕ}

/-! ### The block weight -/

/-- **The block of one side of a spider is centrally unimodal**: at an active side its
weight is `A(p-1, c)` times a product of good factors, with `p - 1 ≥ 1` and a derived
scalar `c = 2^e ≥ 1`; at any other side it is good. -/
theorem spiderBlockSum_isCU (s : SpiderV k len → ℕ) :
    IsCU (∑ s' ∈ sideBlock s, Wpoly (spider k len) s') := by
  by_cases hact : ActiveSide s
  · rw [sum_sideBlock_active hact]
    obtain ⟨c, hc, he⟩ := Wpoly_side_image hact
    have hp : 2 ≤ pNum s := hact.two_odd
    have key : Wpoly (spider k len) s + Wpoly (spider k len) (transf s)
        = Ablock (pNum s - 1) c * ∏ i : Fin k, tailW s i := by
      rw [Wpoly_spider_active hact, he, Ablock]
      simp only [smul_eq_C_mul]
      ring
    rw [key]
    exact (Ablock_isCU (by omega) hc).mul (isGoodW_prod_tailW s).isCU
  · rw [sum_sideBlock_not_active hact]
    exact (isGoodW_Wpoly_spider_of_not_active hact).isCU

/-! ### The partition -/

/-- The finite set of all clan maps of a spider of total order `N`. -/
noncomputable def spiderMapsOfOrder (k : ℕ) (len : Fin k → ℕ) (N : ℕ) :
    Finset (SpiderV k len → ℕ) :=
  (Fintype.piFinset fun _ => Finset.range (N + 1)).filter (fun s => ∑ x, s x = N)

theorem mem_spiderMapsOfOrder {N : ℕ} {s : SpiderV k len → ℕ} :
    s ∈ spiderMapsOfOrder k len N ↔ ∑ x, s x = N := by
  rw [spiderMapsOfOrder, Finset.mem_filter, Fintype.mem_piFinset]
  refine ⟨fun h => h.2, fun h => ⟨fun x => Finset.mem_range.mpr ?_, h⟩⟩
  have hx : s x ≤ ∑ y, s y :=
    Finset.single_le_sum (f := s) (fun i _ => Nat.zero_le _) (Finset.mem_univ x)
  omega

theorem repSide_mem_spiderMapsOfOrder {N : ℕ} {s : SpiderV k len → ℕ}
    (h : s ∈ spiderMapsOfOrder k len N) : repSide s ∈ spiderMapsOfOrder k len N :=
  mem_spiderMapsOfOrder.mpr ((sum_repSide s).trans (mem_spiderMapsOfOrder.mp h))

/-- The block of the representative `β`: the fibre of `repSide` over `β`. -/
noncomputable def spiderBlock (k : ℕ) (len : Fin k → ℕ) (N : ℕ) (β : SpiderV k len → ℕ) :
    Finset (SpiderV k len → ℕ) :=
  (spiderMapsOfOrder k len N).filter (fun s => repSide s = β)

/-- **Disjointness.** -/
theorem spider_blocks_pairwise_disjoint (N : ℕ) {β γ : SpiderV k len → ℕ} (h : β ≠ γ) :
    Disjoint (spiderBlock k len N β) (spiderBlock k len N γ) := by
  rw [Finset.disjoint_left]
  intro s hβ hγ
  rw [spiderBlock, Finset.mem_filter] at hβ hγ
  exact h (hβ.2.symm.trans hγ.2)

/-- **Exhaustion.** -/
theorem spider_blocks_cover {N : ℕ} {s : SpiderV k len → ℕ}
    (h : s ∈ spiderMapsOfOrder k len N) :
    s ∈ spiderBlock k len N (repSide s) ∧ repSide s ∈ spiderMapsOfOrder k len N :=
  ⟨Finset.mem_filter.mpr ⟨h, rfl⟩, repSide_mem_spiderMapsOfOrder h⟩

theorem spider_block_eq_empty_of_not_fixed {N : ℕ} {β : SpiderV k len → ℕ}
    (h : repSide β ≠ β) : spiderBlock k len N β = ∅ := by
  rw [spiderBlock, Finset.filter_eq_empty_iff]
  intro s _ hc
  exact h (by rw [← hc, repSide_idem])

/-- **Each block is exactly the two-element set `{β, transf β}`, or a singleton.** -/
theorem spider_fiber_eq_sideBlock {N : ℕ} {β : SpiderV k len → ℕ} (hβ : repSide β = β)
    (hmem : β ∈ spiderMapsOfOrder k len N) :
    spiderBlock k len N β = sideBlock β := by
  ext s
  simp only [spiderBlock, Finset.mem_filter]
  rw [mem_sideBlock hβ]
  refine ⟨fun h => h.2, fun h => ⟨mem_spiderMapsOfOrder.mpr ?_, h⟩⟩
  rw [← sum_repSide s, h]
  exact mem_spiderMapsOfOrder.mp hmem

/-- The total normalized weight of all clan maps of total order `N`. -/
noncomputable def spiderTotalW (k : ℕ) (len : Fin k → ℕ) (N : ℕ) : LaurentPolynomial ℚ :=
  ∑ s ∈ spiderMapsOfOrder k len N, Wpoly (spider k len) s

/-- **The total weight of a degree is centrally unimodal.** -/
theorem isCU_spiderTotalW (k : ℕ) (len : Fin k → ℕ) (N : ℕ) : IsCU (spiderTotalW k len N) := by
  have hmaps : ∀ s ∈ spiderMapsOfOrder k len N, repSide s ∈ spiderMapsOfOrder k len N :=
    fun _ h => repSide_mem_spiderMapsOfOrder h
  rw [spiderTotalW, ← Finset.sum_fiberwise_of_maps_to hmaps]
  refine Finset.sum_induction _ IsCU (fun _ _ => IsCU.add) isCU_zero ?_
  intro β hβmem
  by_cases hfix : repSide β = β
  · rw [show (spiderMapsOfOrder k len N).filter (fun s => repSide s = β)
      = spiderBlock k len N β from rfl, spider_fiber_eq_sideBlock hfix hβmem]
    exact spiderBlockSum_isCU β
  · rw [show (spiderMapsOfOrder k len N).filter (fun s => repSide s = β) = ∅ from
      spider_block_eq_empty_of_not_fixed hfix, Finset.sum_empty]
    exact isCU_zero

/-- **The spider degreewise theorem.** -/
theorem spider_sum_normalizedTwoRowCoeff_nonneg (k : ℕ) (len : Fin k → ℕ) (N j l : ℕ)
    (hl : 1 ≤ l) (hjl : l ≤ j) (hN : j + l = N) :
    0 ≤ ∑ s ∈ spiderMapsOfOrder k len N, normalizedTwoRowCoeff (spider k len) s j l := by
  have hCU := isCU_spiderTotalW k len N
  have hval : ∑ s ∈ spiderMapsOfOrder k len N, normalizedTwoRowCoeff (spider k len) s j l
      = coeffL (spiderTotalW k len N) ((j : ℤ) - l)
        - coeffL (spiderTotalW k len N) ((j : ℤ) - l + 2) := by
    rw [Finset.sum_congr rfl (fun s hs => normalizedTwoRowCoeff_eq (spider k len) s j l hl
      (by rw [mem_spiderMapsOfOrder.mp hs]; omega)), Finset.sum_sub_distrib, spiderTotalW,
      coeffL_sum, coeffL_sum]
  rw [hval, sub_nonneg]
  exact hCU.decr _ (by omega)

/-! ### Log-concavity -/

theorem spiderMapsOfOrder_eq_multMaps (k : ℕ) (len : Fin k → ℕ) (N : ℕ) :
    spiderMapsOfOrder k len N = multMaps (SpiderV k len) N := by
  ext s
  rw [mem_spiderMapsOfOrder, mem_multMaps]

theorem spider_degreewise_nonneg (k : ℕ) (len : Fin k → ℕ) (N j l : ℕ) (hl : 1 ≤ l)
    (hjl : l ≤ j) (hN : j + l = N) :
    0 ≤ ∑ s ∈ multMaps (SpiderV k len) N, normalizedTwoRowCoeff (spider k len) s j l := by
  rw [← spiderMapsOfOrder_eq_multMaps]
  exact spider_sum_normalizedTwoRowCoeff_nonneg k len N j l hl hjl hN

/-- **Log-concavity of the independence-count sequence of a spider.** -/
theorem spider_indepCount_logConcave (k : ℕ) (len : Fin k → ℕ) (j : ℕ) :
    indepCount (spider k len) j * indepCount (spider k len) (j + 2)
      ≤ indepCount (spider k len) (j + 1) * indepCount (spider k len) (j + 1) :=
  indepCount_logConcave_of_degreewise_nonneg (spider k len)
    (fun N j' l hl hjl hN => spider_degreewise_nonneg k len N j' l hl hjl hN) j

/-- **Log-concavity of the independence polynomial of a spider.** -/
theorem spider_indepPoly_logConcave (k : ℕ) (len : Fin k → ℕ) (j : ℕ) :
    (indepPoly (spider k len)).coeff j * (indepPoly (spider k len)).coeff (j + 2)
      ≤ (indepPoly (spider k len)).coeff (j + 1) * (indepPoly (spider k len)).coeff (j + 1) :=
  indepPoly_logConcave_of_degreewise_nonneg (spider k len)
    (fun N j' l hl hjl hN => spider_degreewise_nonneg k len N j' l hl hjl hN) j

end ClanAudit
