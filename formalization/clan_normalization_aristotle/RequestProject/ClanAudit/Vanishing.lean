import RequestProject.ClanAudit.TwoHubWeight

/-!
# Vanishing and good-shape classification of the non-active clan maps

The global partition needs an exact classification of *every* clan map, not only of the
maps to which the published local transformation applies.  This file supplies the two
missing facts.

* `Wpoly_eq_zero_of_mult_two` — a multiplicity `≥ 2` next to a positive vertex makes the
  clan graph contain a triangle, so the normalized weight vanishes;
* `exists_spider_triangle_of_not_admRun` — a side whose hub is positive and whose initial
  positive run is *not* made of ones always contains such a pair;
* `isGoodW_Wpoly_spider_of_hub_ne_one` — a side whose hub does not carry multiplicity one
  has a weight of the good shape `2^j (z + z⁻¹)^i` (possibly zero).

Together these prove that every clan map is either weight-zero, of good shape, or falls
into the hub-active stratum handled by `Prefix.lean` / `SideTransform.lean`.
-/

namespace ClanAudit

open Finset LaurentPolynomial SimpleGraph

/-! ### Triangles kill the weight -/

variable {V : Type*} [Fintype V] [DecidableEq V]

/-- **A multiplicity `≥ 2` next to a positive vertex kills the normalized weight.** -/
theorem Wpoly_eq_zero_of_mult_two {G : SimpleGraph V} {alpha : V → ℕ} {u v : V}
    (huv : G.Adj u v) (hu : 2 ≤ alpha u) (hv : 1 ≤ alpha v) : Wpoly G alpha = 0 := by
  rw [Wpoly, clan_imbalanceGF_eq_zero_of_mult_two huv hu hv, smul_zero]

/-! ### Sides whose initial run is not made of ones -/

variable {k : ℕ} {len : Fin k → ℕ}

/-- **A side with a positive hub violating `AdmRun` contains a clan triangle.** -/
theorem exists_spider_triangle_of_not_admRun {s : SpiderV k len → ℕ}
    (hhub : 0 < s none) (h : ¬ AdmRun s) :
    ∃ x y : SpiderV k len, (spider k len).Adj x y ∧ 2 ≤ s x ∧ 1 ≤ s y := by
  rw [AdmRun] at h
  push_neg at h
  obtain ⟨i, j, hrun, hne⟩ := h
  have hpos : 0 < armVal s i j := hrun j le_rfl
  have hlt : j < len i := by
    by_contra hc
    rw [armVal_of_ge s i (by omega)] at hpos
    omega
  have h2 : 2 ≤ armVal s i j := by omega
  rw [armVal, dif_pos hlt] at h2
  match j with
  | 0 => exact ⟨some ⟨i, ⟨0, hlt⟩⟩, none, rfl, h2, hhub⟩
  | (j + 1) =>
      have hlt' : j < len i := by omega
      have hp : 0 < armVal s i j := hrun j (by omega)
      rw [armVal, dif_pos hlt'] at hp
      exact ⟨some ⟨i, ⟨j + 1, hlt⟩⟩, some ⟨i, ⟨j, hlt'⟩⟩, ⟨rfl, Or.inr rfl⟩, h2, hp⟩

/-- A side whose hub is positive and whose initial positive run is not made of ones has
weight zero. -/
theorem Wpoly_spider_eq_zero_of_not_admRun {s : SpiderV k len → ℕ}
    (hhub : 0 < s none) (h : ¬ AdmRun s) : Wpoly (spider k len) s = 0 := by
  obtain ⟨x, y, hadj, hx, hy⟩ := exists_spider_triangle_of_not_admRun hhub h
  exact Wpoly_eq_zero_of_mult_two hadj hx hy

/-! ### Sides whose hub does not carry multiplicity one -/

/-- **A side whose hub does not carry multiplicity one has a good weight.**  Either the
hub is switched off, and the side is a product of arm weights; or the hub carries
multiplicity `≥ 2`, and then either the first vertex of some arm is positive — producing
a triangle — or again the side splits. -/
theorem isGoodW_Wpoly_spider_of_hub_ne_one {s : SpiderV k len → ℕ} (h : s none ≠ 1) :
    IsGoodW (Wpoly (spider k len) s) := by
  classical
  by_cases hcross : ∀ (i : Fin k), 0 < len i → s none = 0 ∨ armVal s i 0 = 0
  · rw [Wpoly_spider_hub s hcross]
    refine IsGoodW.mul ?_ (isGoodW_dotW _)
    exact Finset.prod_induction _ IsGoodW (fun _ _ => IsGoodW.mul) IsGoodW.one
      (fun i _ => Wpoly_pathGraph_isGood _ _)
  · push_neg at hcross
    obtain ⟨i, hlen, hhub, harm⟩ := hcross
    have hh2 : 2 ≤ s none := by omega
    have harm' : 1 ≤ s (some ⟨i, ⟨0, hlen⟩⟩) := by
      rw [armVal, dif_pos hlen] at harm
      omega
    refine Or.inl (Wpoly_eq_zero_of_mult_two (G := spider k len) (u := none)
      (v := some ⟨i, ⟨0, hlen⟩⟩) rfl hh2 harm')

/-! ### Good weights of the whole non-active stratum -/

theorem isGoodW_Pw_le_one {t : ℕ} (h : t ≤ 1) : IsGoodW (Pw t) := by
  have ht : t = t % 2 := by omega
  rw [ht]
  exact IsGoodW.Pw_mod_two t

theorem isGoodW_prod_tailW (s : SpiderV k len → ℕ) : IsGoodW (∏ i : Fin k, tailW s i) :=
  Finset.prod_induction _ IsGoodW (fun _ _ => IsGoodW.mul) IsGoodW.one
    (fun i _ => isGoodW_tailW s i)

/-- **Every side that is not active has a good weight.**  Either the hub does not carry
multiplicity one, or the initial run is not made of ones (weight zero), or fewer than
two arms have odd active prefix, in which case the hub component has imbalance at most
one. -/
theorem isGoodW_Wpoly_spider_of_not_active {s : SpiderV k len → ℕ} (h : ¬ ActiveSide s) :
    IsGoodW (Wpoly (spider k len) s) := by
  by_cases h1 : s none = 1
  · by_cases h2 : AdmRun s
    · have hp : pNum s ≤ 1 := by
        by_contra hc
        exact h ⟨h1, h2, by omega⟩
      rw [Wpoly_side_active s h1 (fun i hi => h2.stop i hi)]
      exact IsGoodW.mul (isGoodW_Pw_le_one (by omega)) (isGoodW_prod_tailW s)
    · exact Or.inl (Wpoly_spider_eq_zero_of_not_admRun (by omega) h2)
  · exact isGoodW_Wpoly_spider_of_hub_ne_one h1

/-! ### The two-hub versions -/

variable {m n : ℕ} {a : Fin m → ℕ} {b : Fin n → ℕ}

/-- If both hubs are positive and one of them carries multiplicity `≥ 2`, the weight of
the two-hub map vanishes. -/
theorem Wpoly_dgraph_eq_zero_of_hub_two_left (su : SpiderV m a → ℕ) (sv : SpiderV n b → ℕ)
    (hu : 2 ≤ su none) (hv : 1 ≤ sv none) :
    Wpoly (dgraph m a n b) (Sum.elim su sv) = 0 :=
  Wpoly_eq_zero_of_mult_two (G := dgraph m a n b) (u := Sum.inl none) (v := Sum.inr none)
    (dgraph_adj_lr.mpr ⟨rfl, rfl⟩) hu hv

theorem Wpoly_dgraph_eq_zero_of_hub_two_right (su : SpiderV m a → ℕ) (sv : SpiderV n b → ℕ)
    (hu : 1 ≤ su none) (hv : 2 ≤ sv none) :
    Wpoly (dgraph m a n b) (Sum.elim su sv) = 0 :=
  Wpoly_eq_zero_of_mult_two (G := dgraph m a n b) (u := Sum.inr none) (v := Sum.inl none)
    (dgraph_adj_rl.mpr ⟨rfl, rfl⟩) hv hu

/-- A two-hub map whose left side has a positive hub and a non-admissible run vanishes. -/
theorem Wpoly_dgraph_eq_zero_of_not_admRun_left (su : SpiderV m a → ℕ) (sv : SpiderV n b → ℕ)
    (hhub : 0 < su none) (h : ¬ AdmRun su) :
    Wpoly (dgraph m a n b) (Sum.elim su sv) = 0 := by
  obtain ⟨x, y, hadj, hx, hy⟩ := exists_spider_triangle_of_not_admRun hhub h
  exact Wpoly_eq_zero_of_mult_two (G := dgraph m a n b) (u := Sum.inl x) (v := Sum.inl y)
    (dgraph_adj_ll.mpr hadj) hx hy

/-- A two-hub map whose right side has a positive hub and a non-admissible run vanishes. -/
theorem Wpoly_dgraph_eq_zero_of_not_admRun_right (su : SpiderV m a → ℕ) (sv : SpiderV n b → ℕ)
    (hhub : 0 < sv none) (h : ¬ AdmRun sv) :
    Wpoly (dgraph m a n b) (Sum.elim su sv) = 0 := by
  obtain ⟨x, y, hadj, hx, hy⟩ := exists_spider_triangle_of_not_admRun hhub h
  exact Wpoly_eq_zero_of_mult_two (G := dgraph m a n b) (u := Sum.inr x) (v := Sum.inr y)
    (dgraph_adj_rr.mpr hadj) hx hy

end ClanAudit
