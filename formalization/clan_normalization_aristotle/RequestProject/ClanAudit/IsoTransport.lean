import RequestProject.ClanAudit.Bridge

/-!
# Isomorphism invariance of the independence counts

The independence-count sequence `indepCount` and hence the log-concavity statement are
invariant under graph isomorphism.  This is what makes the explicit models of this
development usable for arbitrary graphs of the same shape.
-/

namespace ClanAudit

open Finset

variable {V W : Type*} [Fintype V] [DecidableEq V] [Fintype W] [DecidableEq W]

open scoped Classical

omit [Fintype V] [DecidableEq V] [Fintype W] in
/-- Independent sets are transported by a graph isomorphism. -/
theorem isIndep_map {G : SimpleGraph V} {H : SimpleGraph W} (e : G ≃g H) {S : Finset V}
    (hS : IsIndep G S) : IsIndep H (S.image e) := by
  intro x hx y hy hxy
  obtain ⟨u, hu, rfl⟩ := Finset.mem_image.mp hx
  obtain ⟨v, hv, rfl⟩ := Finset.mem_image.mp hy
  exact hS u hu v hv ((e.map_adj_iff).mp hxy)

/-- The independence counts are isomorphism invariants. -/
theorem indepCount_iso {G : SimpleGraph V} {H : SimpleGraph W} (e : G ≃g H) (j : ℕ) :
    indepCount G j = indepCount H j := by
  have hinj : Function.Injective (e : V → W) := e.injective
  have hinj' : Function.Injective (e.symm : W → V) := e.symm.injective
  refine Finset.card_bij (fun S _ => S.image e) ?_ ?_ ?_
  · intro S hS
    obtain ⟨h1, h2⟩ := mem_indepSets.mp hS
    exact mem_indepSets.mpr ⟨isIndep_map e h1,
      by rw [Finset.card_image_of_injective _ hinj, h2]⟩
  · intro S hS T hT hST
    have hST' : S.image (e : V → W) = T.image (e : V → W) := hST
    have : (S.image (e : V → W)).image (e.symm : W → V)
        = (T.image (e : V → W)).image (e.symm : W → V) := by rw [hST']
    simpa [Finset.image_image, Function.comp_def] using this
  · intro T hT
    refine ⟨T.image (e.symm : W → V), ?_, ?_⟩
    · obtain ⟨h1, h2⟩ := mem_indepSets.mp hT
      exact mem_indepSets.mpr ⟨isIndep_map e.symm h1,
        by rw [Finset.card_image_of_injective _ hinj', h2]⟩
    · simp [Finset.image_image, Function.comp_def]

/-- Log-concavity of the independence counts is transported along an isomorphism. -/
theorem indepCount_logConcave_of_iso {G : SimpleGraph V} {H : SimpleGraph W} (e : G ≃g H)
    (hG : ∀ j : ℕ, indepCount G j * indepCount G (j + 2)
      ≤ indepCount G (j + 1) * indepCount G (j + 1)) (j : ℕ) :
    indepCount H j * indepCount H (j + 2) ≤ indepCount H (j + 1) * indepCount H (j + 1) := by
  rw [← indepCount_iso e, ← indepCount_iso e, ← indepCount_iso e]
  exact hG j

/-- The same, phrased for the coefficients of the independence polynomial. -/
theorem indepPoly_logConcave_of_iso {G : SimpleGraph V} {H : SimpleGraph W} (e : G ≃g H)
    (hG : ∀ j : ℕ, indepCount G j * indepCount G (j + 2)
      ≤ indepCount G (j + 1) * indepCount G (j + 1)) (j : ℕ) :
    (indepPoly H).coeff j * (indepPoly H).coeff (j + 2)
      ≤ (indepPoly H).coeff (j + 1) * (indepPoly H).coeff (j + 1) := by
  simpa [coeff_indepPoly] using indepCount_logConcave_of_iso e hG j

end ClanAudit
