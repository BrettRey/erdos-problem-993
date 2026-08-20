import RequestProject.ClanAudit.Normalized

/-!
# Normalized weights of clan graphs under disjoint unions and cloned `K₂`s

The published spider transformation replaces a run of multiplicity-one vertices by an
alternating run `2,0,2,0,…`.  Each multiplicity-two vertex that becomes isolated in the
clan graph produces a *cloned `K₂`* component.  The request emphasises that the
transformation "does not preserve the raw factorial denominator term-by-term; the cloned
`K2` orientation counts are supposed to cancel the new factors of `2!`", and that this
cancellation "is one of the main facts to prove, not an assumption".

This file proves exactly that fact:

* `imbalanceGF_iso` — the imbalance generating function is an isomorphism invariant;
* `Wpoly_sum_elim` — `W` is multiplicative over disjoint unions of clan graphs;
* `Wpoly_isolated_two` — `W(single vertex, alpha = 2) = 1`: the cloned `K₂` contributes
  the orientation count `2`, which is exactly cancelled by the new factor `2! = 2`;
* `Wpoly_add_cloned_K2` — adding an isolated multiplicity-two vertex does not change the
  normalized imbalance polynomial, hence changes no normalized two-row coefficient.
-/

namespace ClanAudit

open Finset LaurentPolynomial

variable {W₁ W₂ : Type*} [Fintype W₁] [DecidableEq W₁] [Fintype W₂] [DecidableEq W₂]

/-! ### Isomorphism invariance -/

theorem isStablePair_image {H₁ : SimpleGraph W₁} {H₂ : SimpleGraph W₂} (e : H₁ ≃g H₂)
    {A B : Finset W₁} (h : IsStablePair H₁ A B) :
    IsStablePair H₂ (A.image e) (B.image e) := by
  obtain ⟨hA, hB, hd, hu⟩ := h
  have hinj : Function.Injective (e : W₁ → W₂) := e.toEquiv.injective
  refine ⟨?_, ?_, ?_, ?_⟩
  · intro x hx y hy hadj
    simp only [Finset.mem_image] at hx hy
    obtain ⟨a, ha, rfl⟩ := hx
    obtain ⟨b, hb, rfl⟩ := hy
    exact hA a ha b hb (e.map_rel_iff.mp hadj)
  · intro x hx y hy hadj
    simp only [Finset.mem_image] at hx hy
    obtain ⟨a, ha, rfl⟩ := hx
    obtain ⟨b, hb, rfl⟩ := hy
    exact hB a ha b hb (e.map_rel_iff.mp hadj)
  · rw [Finset.disjoint_left]
    intro x hx hx'
    simp only [Finset.mem_image] at hx hx'
    obtain ⟨a, ha, rfl⟩ := hx
    obtain ⟨b, hb, hab⟩ := hx'
    exact (Finset.disjoint_left.mp hd) ha (hinj hab ▸ hb)
  · rw [← Finset.image_union, hu]
    ext x
    simp only [Finset.mem_image, Finset.mem_univ, iff_true]
    exact ⟨e.symm x, by simp⟩

/-- **The imbalance generating function is an isomorphism invariant.** -/
theorem imbalanceGF_iso {H₁ : SimpleGraph W₁} {H₂ : SimpleGraph W₂} (e : H₁ ≃g H₂) :
    imbalanceGF H₁ = imbalanceGF H₂ := by
  classical
  have hinj : Function.Injective (e : W₁ → W₂) := e.toEquiv.injective
  have hinj' : Function.Injective ((e.symm : W₂ → W₁)) := e.symm.toEquiv.injective
  refine Finset.sum_nbij' (i := fun p => (p.1.image e, p.2.image e))
    (j := fun q => (q.1.image e.symm, q.2.image e.symm)) ?_ ?_ ?_ ?_ ?_
  · intro p hp
    rw [mem_stablePairs] at hp ⊢
    exact isStablePair_image e hp
  · intro q hq
    rw [mem_stablePairs] at hq ⊢
    exact isStablePair_image e.symm hq
  · intro p _
    have h1 : ∀ A : Finset W₁, (A.image (e : W₁ → W₂)).image (e.symm : W₂ → W₁) = A := by
      intro A
      rw [Finset.image_image]
      simp
    exact Prod.ext_iff.mpr ⟨h1 p.1, h1 p.2⟩
  · intro q _
    have h1 : ∀ A : Finset W₂, (A.image (e.symm : W₂ → W₁)).image (e : W₁ → W₂) = A := by
      intro A
      rw [Finset.image_image]
      simp
    exact Prod.ext_iff.mpr ⟨h1 q.1, h1 q.2⟩
  · intro p _
    rw [Finset.card_image_of_injective _ hinj, Finset.card_image_of_injective _ hinj]

/-! ### Clan graphs of disjoint unions -/

variable {V₁ V₂ : Type*} [Fintype V₁] [DecidableEq V₁] [Fintype V₂] [DecidableEq V₂]

/-- The clan graph of a disjoint union is the disjoint union of the clan graphs. -/
def clanSumIso (G₁ : SimpleGraph V₁) (G₂ : SimpleGraph V₂) (a₁ : V₁ → ℕ) (a₂ : V₂ → ℕ) :
    clan (G₁.sum G₂) (Sum.elim a₁ a₂) ≃g (clan G₁ a₁).sum (clan G₂ a₂) where
  toFun x := match x with
    | ⟨Sum.inl v, k⟩ => Sum.inl ⟨v, k⟩
    | ⟨Sum.inr v, k⟩ => Sum.inr ⟨v, k⟩
  invFun y := match y with
    | Sum.inl ⟨v, k⟩ => ⟨Sum.inl v, k⟩
    | Sum.inr ⟨v, k⟩ => ⟨Sum.inr v, k⟩
  left_inv := by rintro ⟨v | v, k⟩ <;> rfl
  right_inv := by rintro (⟨v, k⟩ | ⟨v, k⟩) <;> rfl
  map_rel_iff' := by
    rintro ⟨u | u, i⟩ ⟨v | v, j⟩ <;>
      simp [clan, SimpleGraph.sum_adj, Sigma.ext_iff]

theorem imbalanceGF_clan_sum (G₁ : SimpleGraph V₁) (G₂ : SimpleGraph V₂) (a₁ : V₁ → ℕ)
    (a₂ : V₂ → ℕ) :
    imbalanceGF (clan (G₁.sum G₂) (Sum.elim a₁ a₂))
      = imbalanceGF (clan G₁ a₁) * imbalanceGF (clan G₂ a₂) := by
  rw [imbalanceGF_iso (clanSumIso G₁ G₂ a₁ a₂), imbalanceGF_sum]

omit [DecidableEq V₁] [DecidableEq V₂] in
theorem clanFactorial_sum_elim (a₁ : V₁ → ℕ) (a₂ : V₂ → ℕ) :
    clanFactorial (Sum.elim a₁ a₂) = clanFactorial a₁ * clanFactorial a₂ := by
  rw [clanFactorial, clanFactorial, clanFactorial, Fintype.prod_sum_type]
  rfl

/-- **`W` is multiplicative over disjoint unions**, normalization included. -/
theorem Wpoly_sum_elim (G₁ : SimpleGraph V₁) (G₂ : SimpleGraph V₂) (a₁ : V₁ → ℕ)
    (a₂ : V₂ → ℕ) :
    Wpoly (G₁.sum G₂) (Sum.elim a₁ a₂) = Wpoly G₁ a₁ * Wpoly G₂ a₂ := by
  rw [Wpoly, Wpoly, Wpoly, imbalanceGF_clan_sum, clanFactorial_sum_elim, smul_mul_smul_comm,
    mul_inv]

/-! ### A cloned `K₂` normalizes to `1` -/

/-- The clan graph of a single vertex of multiplicity two is `K₂`, which is connected. -/
theorem clan_dot_two_connected :
    (clan (⊥ : SimpleGraph (Fin 1)) (fun _ => 2)).Connected := by
  rw [SimpleGraph.connected_iff]
  refine ⟨fun x y => ?_, ⟨⟨0, ⟨0, by norm_num⟩⟩⟩⟩
  by_cases h : x = y
  · rw [h]
  · exact SimpleGraph.Adj.reachable (Or.inl ⟨Subsingleton.elim _ _, h⟩)

/-- The cloned `K₂` has orientation count `2 = z⁰ + z⁰`. -/
theorem imbalanceGF_clan_dot_two :
    imbalanceGF (clan (⊥ : SimpleGraph (Fin 1)) (fun _ => 2)) = 2 := by
  classical
  set x₀ : (Σ _ : Fin 1, Fin 2) := ⟨0, ⟨0, by norm_num⟩⟩ with hx₀
  set x₁ : (Σ _ : Fin 1, Fin 2) := ⟨0, ⟨1, by norm_num⟩⟩ with hx₁
  have hne : x₀ ≠ x₁ := by decide
  have hst : IsStablePair (clan (⊥ : SimpleGraph (Fin 1)) (fun _ => 2)) {x₀} {x₁} := by
    refine ⟨?_, ?_, ?_, ?_⟩
    · intro a ha b hb hadj
      rw [Finset.mem_singleton] at ha hb
      subst ha; subst hb
      exact (clan (⊥ : SimpleGraph (Fin 1)) (fun _ => 2)).irrefl hadj
    · intro a ha b hb hadj
      rw [Finset.mem_singleton] at ha hb
      subst ha; subst hb
      exact (clan (⊥ : SimpleGraph (Fin 1)) (fun _ => 2)).irrefl hadj
    · simpa using hne
    · decide
  rw [imbalanceGF_connected clan_dot_two_connected hst]
  simp only [Finset.card_singleton, Nat.cast_one, sub_self, T_zero]
  norm_num

/-- **The cloned `K₂` normalizes to `1`.**  Its orientation count `2` is exactly
cancelled by the factorial `2! = 2` it introduces in the denominator. -/
theorem Wpoly_isolated_two : Wpoly (⊥ : SimpleGraph (Fin 1)) (fun _ => 2) = 1 := by
  have hf : clanFactorial (fun _ : Fin 1 => (2 : ℕ)) = 2 := by
    simp [clanFactorial, Nat.factorial]
  rw [Wpoly, imbalanceGF_clan_dot_two, hf, smul_eq_C_mul,
    show (2 : LaurentPolynomial ℚ) = C (2 : ℚ) from
      (LaurentPolynomial.ext_iff.mpr (congrFun rfl)).symm, ← map_mul]
  norm_num

/-- **The `K₂` cancellation.**  Adjoining an isolated vertex of multiplicity two — i.e.
a cloned `K₂` component of the clan graph — leaves the normalized imbalance polynomial
unchanged, hence changes no normalized two-row coefficient. -/
theorem Wpoly_add_cloned_K2 (G : SimpleGraph V₁) (alpha : V₁ → ℕ) :
    Wpoly (G.sum (⊥ : SimpleGraph (Fin 1))) (Sum.elim alpha (fun _ => 2)) = Wpoly G alpha := by
  rw [Wpoly_sum_elim, Wpoly_isolated_two, mul_one]

end ClanAudit
