import RequestProject.ClanAudit.ConnectorBlock

/-!
# The connector tree with no interior vertex

For `t = 0` — that is, for a single connector edge — the connector tree *is* the
adjacent two-hub tree `dgraph m a n b`, and the block sums are exactly the ones computed
in `BlockWeight.lean`.  This file records the isomorphism and extends
`connBlockSum_isCU` to all `t`.
-/

namespace ClanAudit

open Finset LaurentPolynomial SimpleGraph

variable {t m n : ℕ} {a : Fin m → ℕ} {b : Fin n → ℕ}

/-- With no interior connector vertex the vertex type is that of the adjacent two-hub
tree. -/
def conn0Equiv (m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) :
    ConnV 0 m a n b ≃ DV m a n b where
  toFun x :=
    match x with
    | Sum.inl p => Sum.inl p
    | Sum.inr (Sum.inl j) => j.elim0
    | Sum.inr (Sum.inr q) => Sum.inr q
  invFun x :=
    match x with
    | Sum.inl p => Sum.inl p
    | Sum.inr q => Sum.inr (Sum.inr q)
  left_inv := by
    rintro (p | j | q)
    · rfl
    · exact j.elim0
    · rfl
  right_inv := by rintro (p | q) <;> rfl

/-- With no interior connector vertex the connector tree is the adjacent two-hub tree. -/
def conn0Iso (m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) :
    connGraph 0 m a n b ≃g dgraph m a n b where
  toEquiv := conn0Equiv m a n b
  map_rel_iff' := by
    rintro (p | j | q) (p' | j' | q')
    · exact Iff.rfl
    · exact j'.elim0
    · show (dgraph m a n b).Adj (Sum.inl p) (Sum.inr q') ↔ _
      rw [dgraph_adj_lr, connGraph_adj_lr]
      exact ⟨fun h => ⟨h.1, h.2, rfl⟩, fun h => ⟨h.1, h.2.1⟩⟩
    · exact j.elim0
    · exact j.elim0
    · exact j.elim0
    · show (dgraph m a n b).Adj (Sum.inr q) (Sum.inl p') ↔ _
      rw [dgraph_adj_rl]
      show _ ↔ (p' = none ∧ q = none ∧ (0 : ℕ) = 0)
      exact ⟨fun h => ⟨h.2, h.1, rfl⟩, fun h => ⟨h.2.1, h.1⟩⟩
    · exact j'.elim0
    · exact Iff.rfl

theorem Wpoly_conn0 (x : SpiderV m a → ℕ) (w : Fin 0 → ℕ) (y : SpiderV n b → ℕ) :
    Wpoly (connGraph 0 m a n b) (Sum.elim x (Sum.elim w y))
      = Wpoly (dgraph m a n b) (Sum.elim x y) := by
  rw [← Wpoly_of_iso (conn0Iso m a n b) (Sum.elim x y)]
  congr 1
  funext z
  rcases z with p | j | q
  · rfl
  · exact j.elim0
  · rfl

/-- **The block sum for a single connector edge** is the adjacent two-hub block sum. -/
theorem connBlockSum_zero_isCU (u : SpiderV m a → ℕ) (w : Fin 0 → ℕ)
    (v : SpiderV n b → ℕ) :
    IsCU (∑ u' ∈ sideBlock u, ∑ v' ∈ sideBlock v,
      Wpoly (connGraph 0 m a n b) (Sum.elim u' (Sum.elim w v'))) := by
  classical
  rw [Finset.sum_congr rfl (fun u' _ => Finset.sum_congr rfl (fun v' _ =>
    Wpoly_conn0 u' w v'))]
  by_cases hu : ActiveSide u
  · by_cases hv : ActiveSide v
    · rw [sum_sideBlock_active hu, sum_sideBlock_active hv, sum_sideBlock_active hv]
      exact blockSum_both_active hu hv
    · rw [sum_sideBlock_active hu, sum_sideBlock_not_active hv, sum_sideBlock_not_active hv]
      exact blockSum_left_active hu hv
  · by_cases hv : ActiveSide v
    · rw [sum_sideBlock_not_active hu, sum_sideBlock_active hv]
      exact blockSum_right_active hu hv
    · rw [sum_sideBlock_not_active hu, sum_sideBlock_not_active hv]
      exact (blockSum_neither_active hu hv).isCU

/-- **Every block of the connector partition is centrally unimodal**, for every number
`t ≥ 0` of interior connector vertices, that is, for every number `ell = t + 1 ≥ 1` of
connector edges. -/
theorem connBlockSum_isCU_all (u : SpiderV m a → ℕ) (w : Fin t → ℕ) (v : SpiderV n b → ℕ) :
    IsCU (∑ u' ∈ sideBlock u, ∑ v' ∈ sideBlock v,
      Wpoly (connGraph t m a n b) (Sum.elim u' (Sum.elim w v'))) := by
  rcases Nat.eq_zero_or_pos t with ht | ht
  · subst ht
    exact connBlockSum_zero_isCU u w v
  · exact connBlockSum_isCU ht u w v

end ClanAudit
