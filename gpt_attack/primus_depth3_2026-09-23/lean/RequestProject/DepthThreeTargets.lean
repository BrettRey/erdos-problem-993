import RequestProject.LowDensity

/-!
# The two target propositions of the frozen specification

* `DepthThree.blocked_shadow_of_forest : DepthThree.BlockedShadowTarget`
* `DepthThree.low_density_of_forest : DepthThree.LowDensityTarget`

Both are proved for an arbitrary finite forest (disconnected forests and isolated vertices
included), with no extra hypotheses beyond the ones in the frozen statements.
-/

open MatchingBag

namespace DepthThree

/-- **Target 1.**  Every finite forest with `α ≥ 4` satisfies `(α - 4) b₂ ≤ 6 b₃`. -/
theorem blocked_shadow_of_forest : DepthThree.BlockedShadowTarget := by
  intro V _ G hG ha
  classical
  obtain ⟨D, rfl⟩ := MatchingBag.exists_treeMatching_of_isAcyclic G hG
  rw [TreeMatching.indepNum_eq_card_Bag] at ha ⊢
  rw [TreeMatching.b_eq_card_blkSets D (show 2 ≤ Fintype.card D.Bag by omega),
    TreeMatching.b_eq_card_blkSets D (show 3 ≤ Fintype.card D.Bag by omega)]
  exact TreeMatching.blocked_shadow D ha

/-- **Target 2.**  Every finite forest with `α ∈ {17, 18, 19}` and
`3 (α - 3) b₄ ≤ (α - 7) e₄` satisfies `s₂ s₄ < s₃²`. -/
theorem low_density_of_forest : DepthThree.LowDensityTarget := by
  intro V _ G hG h17 h19 hden
  classical
  obtain ⟨D, rfl⟩ := MatchingBag.exists_treeMatching_of_isAcyclic G hG
  rw [TreeMatching.indepNum_eq_card_Bag] at h17 h19 hden
  exact TreeMatching.depthThreeStrict_of_lowDensity D h17 h19 hden

#print axioms blocked_shadow_of_forest
#print axioms low_density_of_forest

end DepthThree
