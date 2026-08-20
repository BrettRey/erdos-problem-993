import RequestProject.ClanAudit.Clan
import RequestProject.ClanAudit.Stable

/-!
# Normalized two-row coefficients of clan graphs

Combining the clan normalization (`RequestProject.ClanAudit.Clan`) with the stable
two-block partition bookkeeping (`RequestProject.ClanAudit.Stable`) we define

* `clanFactorial alpha = ∏ v, (alpha v)!`, the normalization denominator;
* `Wpoly G alpha`, the normalized imbalance Laurent polynomial
  `W(G, alpha; z) = (1 / ∏ (alpha v)!) * (imbalance polynomial of the clan graph)`;
* `normalizedTwoRowCoeff G alpha k l`, the normalized two-row coefficient.

The main statement `normalizedTwoRowCoeff_eq` is exactly the interior formula asked for
in the request:

```text
normalizedTwoRowCoeff alpha k l = coeff W (k-l) - coeff W (k-l+2).
```

Together with `imbalanceGF_connected` (each connected bipartite component contributes
`z^delta + z^(-delta)`, including the balanced case `delta = 0`, which contributes
`2`, and isolated vertices, which contribute `z + z⁻¹`) and `imbalanceGF_sum`
(multiplicativity over disjoint unions), this *derives* the candidate Laurent
polynomial `W` of the request from `stableCount`, rather than assuming it.
-/

namespace ClanAudit

open Finset LaurentPolynomial

variable {V : Type*} [Fintype V] [DecidableEq V]

/-- The clan normalization denominator `∏ v (alpha v)!`. -/
noncomputable def clanFactorial (alpha : V → ℕ) : ℚ :=
  ∏ v : V, (Nat.factorial (alpha v) : ℚ)

omit [DecidableEq V] in
theorem clanFactorial_pos (alpha : V → ℕ) : 0 < clanFactorial alpha := by
  refine Finset.prod_pos fun v _ => ?_
  exact_mod_cast Nat.factorial_pos (alpha v)

omit [DecidableEq V] in
theorem clanFactorial_ne_zero (alpha : V → ℕ) : clanFactorial alpha ≠ 0 :=
  ne_of_gt (clanFactorial_pos alpha)

omit [DecidableEq V] in
/-- The order of the clan graph is `∑ v, alpha v`. -/
theorem clan_order (alpha : V → ℕ) :
    Fintype.card (Σ v : V, Fin (alpha v)) = ∑ v : V, alpha v := by
  simp [Fintype.card_sigma]

/-- The normalized imbalance Laurent polynomial `W(G, alpha; z)`. -/
noncomputable def Wpoly (G : SimpleGraph V) (alpha : V → ℕ) : LaurentPolynomial ℚ :=
  (clanFactorial alpha)⁻¹ • imbalanceGF (clan G alpha)

/-- The normalized two-row coefficient of a clan graph. -/
noncomputable def normalizedTwoRowCoeff (G : SimpleGraph V) (alpha : V → ℕ) (k l : ℕ) : ℚ :=
  (twoRowCoeff (clan G alpha) k l : ℚ) / clanFactorial alpha

/-- **The normalized two-row coefficient is a difference of two `W`-coefficients.**
This is the interior formula `normalizedTwoRowCoeff alpha k l = coeff W (k-l) - coeff W (k-l+2)`
of the request, derived from `stableCount`. -/
theorem normalizedTwoRowCoeff_eq (G : SimpleGraph V) (alpha : V → ℕ) (k l : ℕ) (hl : 1 ≤ l)
    (h : k + l = ∑ v : V, alpha v) :
    normalizedTwoRowCoeff G alpha k l
      = coeffL (Wpoly G alpha) ((k : ℤ) - l) - coeffL (Wpoly G alpha) ((k : ℤ) - l + 2) := by
  have hcard : k + l = Fintype.card (Σ v : V, Fin (alpha v)) := by rw [clan_order]; exact h
  rw [normalizedTwoRowCoeff, twoRowCoeff_eq_coeff (clan G alpha) k l hl hcard, Wpoly,
    coeffL_smul, coeffL_smul]
  field_simp

/-! ### Vanishing from triangles in the clan graph -/

/-- If a neighbour of an active vertex has multiplicity at least two, the clan graph
contains a triangle, so its imbalance polynomial – and hence every two-row coefficient –
vanishes. -/
theorem clan_imbalanceGF_eq_zero_of_mult_two {G : SimpleGraph V} {alpha : V → ℕ} {u v : V}
    (huv : G.Adj u v) (hu : 2 ≤ alpha u) (hv : 1 ≤ alpha v) :
    imbalanceGF (clan G alpha) = 0 := by
  refine imbalanceGF_eq_zero_of_triangle
    (x := ⟨u, ⟨0, by omega⟩⟩) (y := ⟨u, ⟨1, by omega⟩⟩) (z := ⟨v, ⟨0, by omega⟩⟩) ?_ ?_ ?_
  · refine Or.inl ⟨rfl, ?_⟩
    simp only [ne_eq, Sigma.mk.injEq, heq_eq_eq, true_and]
    intro hcon
    exact absurd (congrArg Fin.val hcon) (by simp)
  · exact Or.inr huv
  · exact Or.inr huv

theorem clan_normalizedTwoRowCoeff_eq_zero_of_mult_two {G : SimpleGraph V} {alpha : V → ℕ}
    {u v : V} (huv : G.Adj u v) (hu : 2 ≤ alpha u) (hv : 1 ≤ alpha v) (k l : ℕ) (hl : 1 ≤ l)
    (h : k + l = ∑ v : V, alpha v) : normalizedTwoRowCoeff G alpha k l = 0 := by
  rw [normalizedTwoRowCoeff_eq G alpha k l hl h, Wpoly,
    clan_imbalanceGF_eq_zero_of_mult_two huv hu hv]
  simp

/-! ### The derived product form of `W` -/

variable {W₁ W₂ : Type*} [Fintype W₁] [DecidableEq W₁] [Fintype W₂] [DecidableEq W₂]

/-- **Derived weighted bipartition formula.**  For a graph that is the disjoint union of
two connected bipartite graphs with colour-class sizes `(a₁,b₁)` and `(a₂,b₂)`, the
imbalance polynomial is the product of the two orientation weights
`z^(a-b) + z^(b-a)`.  Balanced components (`a = b`) contribute the factor `2`, and an
isolated vertex contributes `z + z⁻¹`; no separate boundary convention is needed. -/
theorem imbalanceGF_sum_connected {H₁ : SimpleGraph W₁} {H₂ : SimpleGraph W₂}
    (hc₁ : H₁.Connected) (hc₂ : H₂.Connected) {A₁ B₁ : Finset W₁} {A₂ B₂ : Finset W₂}
    (h₁ : IsStablePair H₁ A₁ B₁) (h₂ : IsStablePair H₂ A₂ B₂) :
    imbalanceGF (H₁.sum H₂)
      = (T ((A₁.card : ℤ) - B₁.card) + T ((B₁.card : ℤ) - A₁.card))
        * (T ((A₂.card : ℤ) - B₂.card) + T ((B₂.card : ℤ) - A₂.card)) := by
  rw [imbalanceGF_sum, imbalanceGF_connected hc₁ h₁, imbalanceGF_connected hc₂ h₂]

end ClanAudit
