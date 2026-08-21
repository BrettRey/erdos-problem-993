import RequestProject.ClanAudit.EvenArms

/-!
# Cutting a clan graph along inactive vertices

The global bookkeeping of the adjacent two-hub target needs to split a clan graph into
independent pieces, where the pieces are *not* separated in the underlying graph `G`:
they are separated only after the vertices of multiplicity zero have been deleted.

`Wpoly_split_sum` is the basic tool: if the vertex type of `G` is a sum type, all edges
inside each summand are recorded by graphs `G₁`, `G₂`, and every edge between the two
summands has an endpoint of multiplicity zero, then the normalized weight is the product
of the normalized weights of the two sides.  `Wpoly_split_equiv` is the same statement
transported along an arbitrary bijection of the vertex type with a sum type.
-/

namespace ClanAudit

open Finset LaurentPolynomial

variable {V V₁ V₂ : Type*} [Fintype V] [DecidableEq V] [Fintype V₁] [DecidableEq V₁]
  [Fintype V₂] [DecidableEq V₂]

/-- **Cutting a clan graph.**  If every edge of `G` joining the two summands of the
vertex type has an endpoint of multiplicity zero, the normalized weight factors. -/
theorem Wpoly_split_sum (G : SimpleGraph (V₁ ⊕ V₂)) (alpha : V₁ ⊕ V₂ → ℕ)
    (G₁ : SimpleGraph V₁) (G₂ : SimpleGraph V₂)
    (hl : ∀ x y, G.Adj (Sum.inl x) (Sum.inl y) ↔ G₁.Adj x y)
    (hr : ∀ x y, G.Adj (Sum.inr x) (Sum.inr y) ↔ G₂.Adj x y)
    (hc : ∀ x y, G.Adj (Sum.inl x) (Sum.inr y) →
      alpha (Sum.inl x) = 0 ∨ alpha (Sum.inr y) = 0) :
    Wpoly G alpha
      = Wpoly G₁ (fun x => alpha (Sum.inl x)) * Wpoly G₂ (fun y => alpha (Sum.inr y)) := by
  classical
  have hclan : clan G alpha = clan (G₁.sum G₂) alpha := by
    ext x y
    have hval : ∀ z : Σ v : V₁ ⊕ V₂, Fin (alpha v), 0 < alpha z.1 := by
      rintro ⟨v, i⟩
      exact i.pos
    have key : G.Adj x.1 y.1 ↔ (G₁.sum G₂).Adj x.1 y.1 := by
      obtain ⟨v, i⟩ := x
      obtain ⟨w, j⟩ := y
      have hv := hval ⟨v, i⟩
      have hw := hval ⟨w, j⟩
      simp only at hv hw ⊢
      match v, w with
      | Sum.inl v, Sum.inl w => rw [SimpleGraph.sum_adj]; exact hl v w
      | Sum.inr v, Sum.inr w => rw [SimpleGraph.sum_adj]; exact hr v w
      | Sum.inl v, Sum.inr w =>
          constructor
          · intro h; rcases hc v w h with h' | h' <;> omega
          · intro h; rw [SimpleGraph.sum_adj] at h; exact h.elim
      | Sum.inr v, Sum.inl w =>
          constructor
          · intro h; rcases hc w v h.symm with h' | h' <;> omega
          · intro h; rw [SimpleGraph.sum_adj] at h; exact h.elim
    simp only [clan]
    rw [key]
  have halpha : alpha
      = Sum.elim (fun x => alpha (Sum.inl x)) (fun y => alpha (Sum.inr y)) := by
    funext x; cases x <;> rfl
  calc Wpoly G alpha = Wpoly (G₁.sum G₂) alpha := by rw [Wpoly, hclan, ← Wpoly]
    _ = Wpoly (G₁.sum G₂)
          (Sum.elim (fun x => alpha (Sum.inl x)) (fun y => alpha (Sum.inr y))) := by
        rw [← halpha]
    _ = _ := Wpoly_sum_elim _ _ _ _

/-- **Cutting a clan graph along an arbitrary splitting of the vertex type.** -/
theorem Wpoly_split_equiv (G : SimpleGraph V) (alpha : V → ℕ)
    (G₁ : SimpleGraph V₁) (G₂ : SimpleGraph V₂) (f : V ≃ V₁ ⊕ V₂)
    (hl : ∀ x y, G.Adj (f.symm (Sum.inl x)) (f.symm (Sum.inl y)) ↔ G₁.Adj x y)
    (hr : ∀ x y, G.Adj (f.symm (Sum.inr x)) (f.symm (Sum.inr y)) ↔ G₂.Adj x y)
    (hc : ∀ x y, G.Adj (f.symm (Sum.inl x)) (f.symm (Sum.inr y)) →
      alpha (f.symm (Sum.inl x)) = 0 ∨ alpha (f.symm (Sum.inr y)) = 0) :
    Wpoly G alpha
      = Wpoly G₁ (fun x => alpha (f.symm (Sum.inl x)))
        * Wpoly G₂ (fun y => alpha (f.symm (Sum.inr y))) := by
  classical
  set H : SimpleGraph (V₁ ⊕ V₂) := G.comap f.symm with hH
  have hWH : Wpoly H (fun p => alpha (f.symm p)) = Wpoly G alpha :=
    Wpoly_of_iso (G₁ := H) (G₂ := G)
      { toEquiv := f.symm, map_rel_iff' := Iff.rfl } alpha
  rw [← hWH]
  exact Wpoly_split_sum H (fun p => alpha (f.symm p)) G₁ G₂ hl hr hc

end ClanAudit
