import RequestProject.B1Zero.DensitySmall
import RequestProject.B1Zero.FourStructure

/-!
# The `b₁ = 0` depth-three window theorem

`DepthThree.b1_zero_window : DepthThree.B1ZeroWindowTarget`.

The proof runs the structural reduction of the written argument: with `b₁ = 0` a forest's
allowed subgraph is well covered, forbidden vertices have at least three forced neighbours,
and the number `r` of forbidden vertices satisfies `2r + 1 ≤ |Forced| = r + δ ≤ r + 5`, hence
`r ≤ 4`.  For `r ≤ 3` the density bound `3 (α-3) b₄ ≤ (α-7) e₄` follows from the
coefficientwise blocked bound and the matching-bag deletion inequalities
(`density_small`); the remaining case `r = 4` forces `(|V|, α, |Forced|, #flexBags)
= (33, 19, 9, 10)` and is settled by `density_four_aux`, which represents the flexible part
as the corona of a ten-vertex forest and applies the finite enumeration certificate.  Feeding the density bound
into the already proved `DepthThree.low_density_of_forest` gives the target.
-/

open Finset SimpleGraph

namespace MatchingBag

attribute [local instance 1] Classical.propDecidable

variable {V : Type*} [Fintype V] [DecidableEq V]

namespace TreeMatching

variable {D : TreeMatching V}

/-- **The density bound for exactly four forbidden vertices** (the case
`(n, α, δ, r) = (33, 19, 5, 4)`). -/
theorem density_four (hb1 : D.B1Zero) (hconn : D.G.Connected)
    (h17 : 17 ≤ Fintype.card D.Bag) (h19 : Fintype.card D.Bag ≤ 19)
    (hn1 : 33 ≤ Fintype.card V) (hn2 : Fintype.card V ≤ 38)
    (hr4 : D.ForbiddenV.card = 4) :
    3 * (Fintype.card D.Bag - 3) * (D.blkSets (Fintype.card D.Bag - 4)).card
      ≤ (Fintype.card D.Bag - 7) * (D.extSets (Fintype.card D.Bag - 4)).card := by
  classical
  have hnc := card_verts_add_card_forced D
  have hU : D.ForbiddenV.Nonempty := Finset.card_pos.1 (by omega)
  have hexp := card_forcedNbrsSet_ge hb1 (by omega) (Finset.Subset.refl D.ForbiddenV) hU
  have hNsub : D.forcedNbrsSet D.ForbiddenV ⊆ D.ForcedV := Finset.filter_subset _ _
  have hNcard := Finset.card_le_card hNsub
  have ha : Fintype.card D.Bag = 19 := by omega
  have hc : D.ForcedV.card = 9 := by omega
  have hmm := card_flexBags D
  have hm : D.flexBags.card = 10 := by omega
  have hmain := density_four_aux hb1 hconn ha hc hr4 hm
  rw [ha]
  norm_num
  omega

/-- The density bound in the whole `b₁ = 0` window. -/
theorem density_window (hb1 : D.B1Zero) (hconn : D.G.Connected)
    (h17 : 17 ≤ Fintype.card D.Bag) (h19 : Fintype.card D.Bag ≤ 19)
    (hn1 : 33 ≤ Fintype.card V) (hn2 : Fintype.card V ≤ 38)
    (hdelta : 2 * Fintype.card D.Bag ≤ Fintype.card V + 5) :
    3 * (Fintype.card D.Bag - 3) * (D.blkSets (Fintype.card D.Bag - 4)).card
      ≤ (Fintype.card D.Bag - 7) * (D.extSets (Fintype.card D.Bag - 4)).card := by
  classical
  have hnc : Fintype.card V + D.ForcedV.card = 2 * Fintype.card D.Bag + D.ForbiddenV.card :=
    card_verts_add_card_forced D
  have hrc : 1 ≤ D.ForbiddenV.card → 2 * D.ForbiddenV.card + 1 ≤ D.ForcedV.card :=
    two_mul_card_forbidden_lt hb1 (by omega)
  by_cases hr3 : D.ForbiddenV.card ≤ 3
  · exact density_small hb1 h17 h19 hn1 hn2 hr3
  · have hr4 : D.ForbiddenV.card = 4 := by
      have := hrc (by omega)
      omega
    exact density_four hb1 hconn h17 h19 hn1 hn2 hr4

end TreeMatching

end MatchingBag

open MatchingBag

namespace DepthThree

/-- **The frozen `b₁ = 0` window target.** -/
theorem b1_zero_window : DepthThree.B1ZeroWindowTarget := by
  intro V _ G hG hn1 hn2 h17 h19 hdelta hb1
  classical
  obtain ⟨D, rfl⟩ := MatchingBag.exists_treeMatching_of_isAcyclic G hG.IsAcyclic
  rw [TreeMatching.indepNum_eq_card_Bag] at h17 h19 hdelta
  -- `b₁ = 0` in the specification's sense gives `D.B1Zero`
  have hb1' : D.B1Zero := by
    intro S hcard hind
    rw [TreeMatching.b_eq_card_blkSets D (show 1 ≤ Fintype.card D.Bag by omega),
      Finset.card_eq_zero] at hb1
    by_contra hcon
    push_neg at hcon
    have : S ∈ D.blkSets (Fintype.card D.Bag - 1) := by
      rw [TreeMatching.mem_blkSets]
      exact ⟨hcard, hind, fun ⟨I, hI, hSI⟩ => hcon I hI hSI⟩
    rw [hb1] at this
    exact absurd this (Finset.notMem_empty S)
  have hden : 3 * (Fintype.card D.Bag - 3) * DepthThree.b D.G 4
      ≤ (Fintype.card D.Bag - 7) * DepthThree.e D.G 4 := by
    rw [TreeMatching.b_eq_card_blkSets D (show 4 ≤ Fintype.card D.Bag by omega),
      TreeMatching.e_eq_card_extSets D (show 4 ≤ Fintype.card D.Bag by omega)]
    exact TreeMatching.density_window hb1' hG.isConnected h17 h19 hn1 hn2 hdelta
  exact TreeMatching.depthThreeStrict_of_lowDensity D h17 h19 hden

#print axioms b1_zero_window

end DepthThree
