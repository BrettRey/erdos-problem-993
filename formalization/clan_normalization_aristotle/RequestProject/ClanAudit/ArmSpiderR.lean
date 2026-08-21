import RequestProject.ClanAudit.ArmSpider

/-!
# The mirrored `armSpider`, and the all-ones extra arm

Two complements to `RequestProject/ClanAudit/ArmSpider.lean`:

* `Wpoly_armSpider_active_ones` — the weight of an `armSpider` whose extra arm carries
  all multiplicities one and whose hub is active: the extra arm simply contributes its
  parity to the number of odd arms;
* `armSpiderR` — the same graph with the extra path written on the *left* of the vertex
  sum, which is the shape in which the right half of a cut connector appears.  It is
  isomorphic to `armSpider` by reversing the path, so all results transfer.
-/

namespace ClanAudit

open Finset LaurentPolynomial SimpleGraph

variable {m L : ℕ} {a : Fin m → ℕ}

/-! ### The extra arm with all multiplicities one -/

theorem gVal_ones (j : ℕ) (hj : j < L) : gVal (fun _ : Fin L => 1) j = 1 := by
  rw [gVal, dif_pos hj]

theorem gPref_ones (u : SpiderV m a → ℕ) : gPref u (fun _ : Fin L => 1) = L := by
  rw [gPref_eq]
  refine pref_eq (combMap u (fun _ : Fin L => 1)) (Fin.last m) L (le_of_eq aExt_last.symm)
    (fun j hj => ?_) (Or.inl aExt_last.symm)
  rw [armVal_combMap_last, gVal_ones j hj]

theorem tailW_ones_last (u : SpiderV m a → ℕ) :
    tailW (combMap u (fun _ : Fin L => 1)) (Fin.last m) = 1 := by
  have hlen : aExt m a L (Fin.last m) - pref (combMap u (fun _ : Fin L => 1)) (Fin.last m) = 0 := by
    rw [aExt_last, ← gPref_eq, gPref_ones]
    omega
  rw [tailW]
  rw [Wpoly_pathGraph_len_congr hlen
    (fun x => armVal (combMap u (fun _ : Fin L => 1)) (Fin.last m)
      (pref (combMap u (fun _ : Fin L => 1)) (Fin.last m) + x))]
  exact Wpoly_pathGraph_zero _

theorem admRun_combMap_ones {u : SpiderV m a → ℕ} (h : AdmRun u) :
    AdmRun (combMap u (fun _ : Fin L => 1)) := by
  intro i j
  refine Fin.lastCases ?_ ?_ i
  · intro hpos
    by_cases hj : j < L
    · rw [armVal_combMap_last, gVal_ones j hj]
    · have := hpos j le_rfl
      rw [armVal_combMap_last, gVal, dif_neg hj] at this
      omega
  · intro i' hpos
    rw [armVal_combMap_castSucc]
    refine h i' j (fun j' hj' => ?_)
    have := hpos j' hj'
    rwa [armVal_combMap_castSucc] at this

/-- **The weight of an `armSpider` with an all-ones extra arm and an active hub.** -/
theorem Wpoly_armSpider_active_ones (u : SpiderV m a → ℕ) (hhub : u none = 1)
    (hrun : AdmRun u) :
    Wpoly (armSpider m a L) (Sum.elim u (fun _ => 1))
      = Pw ((1 - ((pNum u : ℤ) + ((L % 2 : ℕ) : ℤ))).natAbs) * ∏ i : Fin m, tailW u i := by
  have hadm : AdmRun (combMap u (fun _ : Fin L => 1)) := admRun_combMap_ones hrun
  rw [Wpoly_armSpider_eq,
    Wpoly_side_active _ (by rw [combMap_none]; exact hhub) (fun i hi => hadm.stop i hi),
    prod_tailW_combMap, tailW_ones_last, mul_one, pNum_combMap, gPref_ones]
  congr 2
  rcases Nat.mod_two_eq_zero_or_one L with h | h <;> rw [h] <;> simp

/-! ### The mirrored graph -/

/-- A spider with one extra pendant path of `L` vertices attached to the hub, with the
path written to the left of the hub: vertex `L - 1` of the path is the neighbour of the
hub. -/
def armSpiderR (m : ℕ) (a : Fin m → ℕ) (L : ℕ) : SimpleGraph (Fin L ⊕ SpiderV m a) where
  Adj x y :=
    match x, y with
    | Sum.inl i, Sum.inl j => i.val + 1 = j.val ∨ j.val + 1 = i.val
    | Sum.inl i, Sum.inr q => q = none ∧ i.val + 1 = L
    | Sum.inr q, Sum.inl i => q = none ∧ i.val + 1 = L
    | Sum.inr p, Sum.inr q => (spider m a).Adj p q
  symm := by
    rintro (i | p) (j | q) h <;> first | exact h | exact h.symm
  loopless := ⟨by
    rintro (i | p) h
    · rcases h with h | h <;> omega
    · exact (spider m a).irrefl h⟩

/-- Reversing a path. -/
def finRev (L : ℕ) : Fin L ≃ Fin L where
  toFun i := ⟨L - 1 - i.val, by have := i.isLt; omega⟩
  invFun i := ⟨L - 1 - i.val, by have := i.isLt; omega⟩
  left_inv i := Fin.ext (by have := i.isLt; simp; omega)
  right_inv i := Fin.ext (by have := i.isLt; simp; omega)

/-- The mirrored graph is isomorphic to `armSpider` by reversing the extra path. -/
def armSpiderRIso (m : ℕ) (a : Fin m → ℕ) (L : ℕ) :
    armSpiderR m a L ≃g armSpider m a L where
  toEquiv := (Equiv.sumCongr (finRev L) (Equiv.refl _)).trans (Equiv.sumComm _ _)
  map_rel_iff' := by
    rintro (i | p) (j | q)
    · show ((armSpider m a L).Adj (Sum.inr (finRev L i)) (Sum.inr (finRev L j)))
        ↔ (armSpiderR m a L).Adj (Sum.inl i) (Sum.inl j)
      rw [armSpider_adj_rr]
      have hi := i.isLt
      have hj := j.isLt
      show ((L - 1 - i.val) + 1 = L - 1 - j.val ∨ (L - 1 - j.val) + 1 = L - 1 - i.val)
        ↔ (i.val + 1 = j.val ∨ j.val + 1 = i.val)
      omega
    · show ((armSpider m a L).Adj (Sum.inr (finRev L i)) (Sum.inl q))
        ↔ (armSpiderR m a L).Adj (Sum.inl i) (Sum.inr q)
      rw [armSpider_adj_rl]
      have hi := i.isLt
      constructor
      · rintro ⟨hq, hval⟩
        exact ⟨hq, by simp only [finRev, Equiv.coe_fn_mk] at hval; omega⟩
      · rintro ⟨hq, hval⟩
        exact ⟨hq, by simp only [finRev, Equiv.coe_fn_mk]; omega⟩
    · show ((armSpider m a L).Adj (Sum.inl p) (Sum.inr (finRev L j)))
        ↔ (armSpiderR m a L).Adj (Sum.inr p) (Sum.inl j)
      rw [armSpider_adj_lr]
      have hj := j.isLt
      constructor
      · rintro ⟨hp, hval⟩
        exact ⟨hp, by simp only [finRev, Equiv.coe_fn_mk] at hval; omega⟩
      · rintro ⟨hp, hval⟩
        exact ⟨hp, by simp only [finRev, Equiv.coe_fn_mk]; omega⟩
    · exact Iff.rfl

theorem Wpoly_armSpiderR (g : Fin L → ℕ) (u : SpiderV m a → ℕ) :
    Wpoly (armSpiderR m a L) (Sum.elim g u)
      = Wpoly (armSpider m a L) (Sum.elim u (fun i => g ((finRev L).symm i))) := by
  rw [← Wpoly_of_iso (armSpiderRIso m a L) (Sum.elim u (fun i => g ((finRev L).symm i)))]
  congr 1
  funext x
  rcases x with i | p
  · show g i = g ((finRev L).symm (finRev L i))
    rw [Equiv.symm_apply_apply]
  · rfl

/-- **The two-map block of a mirrored `armSpider`.** -/
theorem armSpiderR_sideBlock_isCU (g : Fin L → ℕ) (u : SpiderV m a → ℕ) :
    IsCU (∑ u' ∈ sideBlock u, Wpoly (armSpiderR m a L) (Sum.elim g u')) := by
  rw [Finset.sum_congr rfl (fun u' _ => Wpoly_armSpiderR g u')]
  exact armSpider_sideBlock_isCU u (fun i => g ((finRev L).symm i))

/-- **The weight of a mirrored `armSpider` with an all-ones extra arm and an active
hub.** -/
theorem Wpoly_armSpiderR_active_ones (u : SpiderV m a → ℕ) (hhub : u none = 1)
    (hrun : AdmRun u) :
    Wpoly (armSpiderR m a L) (Sum.elim (fun _ => 1) u)
      = Pw ((1 - ((pNum u : ℤ) + ((L % 2 : ℕ) : ℤ))).natAbs) * ∏ i : Fin m, tailW u i := by
  rw [Wpoly_armSpiderR]
  exact Wpoly_armSpider_active_ones u hhub hrun

end ClanAudit
