import RequestProject.B1Zero.EBounds
import RequestProject.B1Zero.Numeric

/-!
# The low-density hypothesis for at most three forbidden vertices

For a forest in the window with `b₁ = 0` and at most three forbidden vertices, we prove the
graph density inequality `3 (α - 3) b₄ ≤ (α - 7) e₄` required by
`DepthThree.low_density_of_forest`.
-/

open Finset SimpleGraph

namespace MatchingBag

attribute [local instance 1] Classical.propDecidable

variable {V : Type*} [Fintype V] [DecidableEq V]

namespace TreeMatching

variable {D : TreeMatching V}

lemma EE_eq_sum_range (cc k : ℕ) :
    D.EE cc k = ∑ i ∈ Finset.range (k + 1), Nat.choose cc i * D.pflex (k - i) := by
  rw [EE, Finset.Nat.sum_antidiagonal_eq_sum_range_succ_mk]

/-- The number of forbidden vertices is at most `δ - 1`, hence at most four in the window. -/
lemma two_mul_card_forbidden_lt (hb1 : D.B1Zero) (ha : 3 ≤ Fintype.card D.Bag)
    (hr : 1 ≤ D.ForbiddenV.card) : 2 * D.ForbiddenV.card + 1 ≤ D.ForcedV.card := by
  have hne : D.ForbiddenV.Nonempty := Finset.card_pos.1 hr
  have h := card_forcedNbrsSet_ge hb1 ha (subset_refl _) hne
  have h2 : (D.forcedNbrsSet D.ForbiddenV).card ≤ D.ForcedV.card :=
    Finset.card_le_card (Finset.filter_subset _ _)
  omega

/-- **The density bound for at most three forbidden vertices.** -/
theorem density_small (hb1 : D.B1Zero)
    (h17 : 17 ≤ Fintype.card D.Bag) (h19 : Fintype.card D.Bag ≤ 19)
    (hn1 : 33 ≤ Fintype.card V) (hn2 : Fintype.card V ≤ 38)
    (hr3 : D.ForbiddenV.card ≤ 3) :
    3 * (Fintype.card D.Bag - 3) * (D.blkSets (Fintype.card D.Bag - 4)).card
      ≤ (Fintype.card D.Bag - 7) * (D.extSets (Fintype.card D.Bag - 4)).card := by
  classical
  set a := Fintype.card D.Bag with ha
  set c := D.ForcedV.card with hc
  set r := D.ForbiddenV.card with hr
  have hnc : Fintype.card V + c = 2 * a + r := card_verts_add_card_forced D
  have hcm : D.flexBags.card + c = a := card_flexBags D
  have hca : c ≤ a := card_forced_le_bags D
  have hrc0 : r ≤ c := card_forbidden_le_forced D
  have hrc : 1 ≤ r → 2 * r + 1 ≤ c := two_mul_card_forbidden_lt hb1 (by omega)
  have hm : a - c = D.flexBags.card := by omega
  have hz : ∀ j, a - c < j → D.pflex j = 0 := by
    intro j hj
    exact pflex_eq_zero_of_gt (by omega)
  have hcone : ∀ l, (l + 1) * D.pflex (l + 1) ≤ 2 * (a - c - l) * D.pflex l := by
    intro l
    rw [hm]
    exact pflex_cone D l
  -- the extendable count
  have he : (D.extSets (a - 4)).card
      = ∑ i ∈ Finset.range (a - 4 + 1), Nat.choose c i * D.pflex (a - 4 - i) := by
    rw [card_extSets_eq_EE hb1, EE_eq_sum_range]
  -- the blocked bound
  have hb : (D.blkSets (a - 4)).card
      ≤ ∑ t ∈ Finset.Icc 1 3, Nat.choose r t *
          ∑ i ∈ Finset.range (a - 4 - t + 1), Nat.choose (c - (2 * t + 1)) i *
            D.pflex (a - 4 - t - i) := by
    refine le_trans (card_blkSets_le hb1 (by omega) (a - 4)) ?_
    have hsub : Finset.Icc 1 r ⊆ Finset.Icc 1 3 := by
      intro t ht
      rw [Finset.mem_Icc] at ht ⊢
      omega
    refine le_trans (le_of_eq (Finset.sum_congr rfl fun t _ => ?_))
      (Finset.sum_le_sum_of_subset_of_nonneg hsub (fun _ _ _ => Nat.zero_le _))
    rw [EE_eq_sum_range]
  rw [he]
  refine le_trans (Nat.mul_le_mul_left _ hb) ?_
  exact DepthThree.B1Zero.density_numeric a c r D.pflex h17 h19 hca hrc0
    (by omega) (by omega) hrc hr3 hz hcone

end TreeMatching

end MatchingBag
