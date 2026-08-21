import RequestProject.ClanAudit.GlobalPartition

/-!
# A spider with one extra pendant path attached to its hub

`armSpider m a L` is the graph consisting of a spider (hub with `m` pendant paths of
lengths `a i`) together with one extra pendant path of `L` vertices attached to the hub.

This is exactly the shape of every piece into which the connector tree falls apart once
the connector is cut: the connector prefix becomes an *ordinary ordered arm* of the hub
it is attached to.

The point of this file is that `armSpider m a L` is isomorphic to the spider
`spider (m+1) (aExt m a L)` whose last arm has length `L`, so that the whole existing
spider machinery (prefixes, active sides, weights) applies verbatim.  What is *not*
inherited is the canonical transformation: on `armSpider` it must act on the `m` original
arms only, leaving the extra arm — which belongs to the connector — untouched.
-/

namespace ClanAudit

open Finset LaurentPolynomial SimpleGraph

/-! ### The graph -/

/-- A spider with one extra pendant path of `L` vertices attached to the hub. -/
def armSpider (m : ℕ) (a : Fin m → ℕ) (L : ℕ) : SimpleGraph (SpiderV m a ⊕ Fin L) where
  Adj x y :=
    match x, y with
    | Sum.inl p, Sum.inl q => (spider m a).Adj p q
    | Sum.inl p, Sum.inr j => p = none ∧ j.val = 0
    | Sum.inr j, Sum.inl p => p = none ∧ j.val = 0
    | Sum.inr i, Sum.inr j => i.val + 1 = j.val ∨ j.val + 1 = i.val
  symm := by
    rintro (p | i) (q | j) h <;> first | exact h | exact h.symm
  loopless := ⟨by
    rintro (p | i) h
    · exact (spider m a).irrefl h
    · rcases h with h | h <;> omega⟩

variable {m L : ℕ} {a : Fin m → ℕ}

theorem armSpider_adj_ll {p q : SpiderV m a} :
    (armSpider m a L).Adj (Sum.inl p) (Sum.inl q) ↔ (spider m a).Adj p q := Iff.rfl

theorem armSpider_adj_lr {p : SpiderV m a} {j : Fin L} :
    (armSpider m a L).Adj (Sum.inl p) (Sum.inr j) ↔ (p = none ∧ j.val = 0) := Iff.rfl

theorem armSpider_adj_rl {j : Fin L} {p : SpiderV m a} :
    (armSpider m a L).Adj (Sum.inr j) (Sum.inl p) ↔ (p = none ∧ j.val = 0) := Iff.rfl

theorem armSpider_adj_rr {i j : Fin L} :
    (armSpider m a L).Adj (Sum.inr i) (Sum.inr j) ↔ (i.val + 1 = j.val ∨ j.val + 1 = i.val) :=
  Iff.rfl

/-! ### The extended arm family -/

/-- The arm family of `m` arms of lengths `a i` extended by one arm of length `L`. -/
def aExt (m : ℕ) (a : Fin m → ℕ) (L : ℕ) : Fin (m + 1) → ℕ :=
  fun i => if h : i.val < m then a ⟨i.val, h⟩ else L

theorem aExt_castSucc (i : Fin m) : aExt m a L i.castSucc = a i := by
  rw [aExt, dif_pos (show i.castSucc.val < m from i.isLt)]
  rfl

theorem aExt_last : aExt m a L (Fin.last m) = L := by
  rw [aExt, dif_neg (by simp)]

/-- The vertex bijection between `armSpider m a L` and the spider with one more arm. -/
def armEquiv (m : ℕ) (a : Fin m → ℕ) (L : ℕ) :
    (SpiderV m a ⊕ Fin L) ≃ SpiderV (m + 1) (aExt m a L) where
  toFun x :=
    match x with
    | Sum.inl none => none
    | Sum.inl (some ⟨i, j⟩) =>
        some ⟨i.castSucc, ⟨j.val, Nat.lt_of_lt_of_eq j.isLt (aExt_castSucc i).symm⟩⟩
    | Sum.inr k => some ⟨Fin.last m, ⟨k.val, Nat.lt_of_lt_of_eq k.isLt aExt_last.symm⟩⟩
  invFun y :=
    match y with
    | none => Sum.inl none
    | some ⟨i, j⟩ =>
        if h : i.val < m then
          Sum.inl (some ⟨⟨i.val, h⟩, ⟨j.val, Nat.lt_of_lt_of_eq j.isLt (by rw [aExt, dif_pos h])⟩⟩)
        else
          Sum.inr ⟨j.val, Nat.lt_of_lt_of_eq j.isLt (by
            have hi : i = Fin.last m := Fin.ext (by have := i.isLt; simp at h ⊢; omega)
            rw [hi, aExt_last])⟩
  left_inv := by
    rintro ((_ | ⟨i, j⟩) | k)
    · rfl
    · simp only
      rw [dif_pos (show (i.castSucc).val < m from i.isLt)]
      rfl
    · simp only
      rw [dif_neg (show ¬ (Fin.last m).val < m from by simp)]
  right_inv := by
    rintro (_ | ⟨i, j⟩)
    · rfl
    · by_cases h : i.val < m
      · simp only [dif_pos h]
        rfl
      · have hi : i = Fin.last m := Fin.ext (by have := i.isLt; simp at h ⊢; omega)
        subst hi
        simp only [dif_neg h]

/-- The graph isomorphism between `armSpider m a L` and the spider with one more arm. -/
def armIso (m : ℕ) (a : Fin m → ℕ) (L : ℕ) :
    armSpider m a L ≃g spider (m + 1) (aExt m a L) where
  toEquiv := armEquiv m a L
  map_rel_iff' := by
    rintro ((_ | ⟨i, j⟩) | k) ((_ | ⟨i', j'⟩) | k')
    · exact Iff.rfl
    · exact Iff.rfl
    · exact ⟨fun h => ⟨rfl, h⟩, fun h => h.2⟩
    · exact Iff.rfl
    · exact ⟨fun h => ⟨Fin.castSucc_injective m h.1, h.2⟩,
        fun h => ⟨congrArg Fin.castSucc h.1, h.2⟩⟩
    · exact ⟨fun h => absurd h.1 (Fin.castSucc_lt_last i).ne,
        fun h => absurd h.1 (by simp)⟩
    · exact ⟨fun h => ⟨rfl, h⟩, fun h => h.2⟩
    · exact ⟨fun h => absurd h.1.symm (Fin.castSucc_lt_last i').ne,
        fun h => absurd h.1 (by simp)⟩
    · exact ⟨fun h => h.2, fun h => ⟨rfl, h⟩⟩

/-- Two paths of equal length carry the same weight for a value function given by index. -/
theorem Wpoly_pathGraph_len_congr {N N' : ℕ} (h : N = N') (f : ℕ → ℕ) :
    Wpoly (pathGraph N) (fun j : Fin N => f j.val)
      = Wpoly (pathGraph N') (fun j : Fin N' => f j.val) := by
  subst h; rfl

/-! ### The clan map of the extended spider -/

/-- A clan map of `armSpider m a L`, read as a clan map of the extended spider. -/
def combMap (u : SpiderV m a → ℕ) (g : Fin L → ℕ) : SpiderV (m + 1) (aExt m a L) → ℕ :=
  fun x => Sum.elim u g ((armEquiv m a L).symm x)

theorem Wpoly_armSpider_eq (u : SpiderV m a → ℕ) (g : Fin L → ℕ) :
    Wpoly (armSpider m a L) (Sum.elim u g)
      = Wpoly (spider (m + 1) (aExt m a L)) (combMap u g) := by
  rw [← Wpoly_of_iso (armIso m a L) (combMap u g)]
  congr 1
  funext x
  show Sum.elim u g x = Sum.elim u g ((armEquiv m a L).symm ((armEquiv m a L) x))
  rw [Equiv.symm_apply_apply]

theorem combMap_none (u : SpiderV m a → ℕ) (g : Fin L → ℕ) : combMap u g none = u none := rfl

theorem armEquiv_symm_some_castSucc (i : Fin m) (k : Fin (aExt m a L i.castSucc)) :
    (armEquiv m a L).symm (some ⟨i.castSucc, k⟩)
      = Sum.inl (some ⟨i, ⟨k.val, Nat.lt_of_lt_of_eq k.isLt (aExt_castSucc i)⟩⟩) := by
  simp only [armEquiv, Equiv.coe_fn_symm_mk]
  rw [dif_pos (show (i.castSucc).val < m from i.isLt)]
  rfl

theorem armEquiv_symm_some_last (k : Fin (aExt m a L (Fin.last m))) :
    (armEquiv m a L).symm (some ⟨Fin.last m, k⟩)
      = Sum.inr ⟨k.val, Nat.lt_of_lt_of_eq k.isLt aExt_last⟩ := by
  simp only [armEquiv, Equiv.coe_fn_symm_mk]
  rw [dif_neg (show ¬ (Fin.last m).val < m from by simp)]

theorem armVal_combMap_castSucc (u : SpiderV m a → ℕ) (g : Fin L → ℕ) (i : Fin m) (j : ℕ) :
    armVal (combMap u g) i.castSucc j = armVal u i j := by
  by_cases h : j < a i
  · rw [armVal, dif_pos (show j < aExt m a L i.castSucc from by rw [aExt_castSucc]; exact h),
      armVal, dif_pos h]
    show Sum.elim u g ((armEquiv m a L).symm (some ⟨i.castSucc, _⟩)) = _
    rw [armEquiv_symm_some_castSucc]
    rfl
  · rw [armVal_of_ge _ _ (show aExt m a L i.castSucc ≤ j from by rw [aExt_castSucc]; omega),
      armVal_of_ge _ _ (by omega)]

/-- The multiplicity of position `j` of the extra arm, extended by `0` past its end. -/
def gVal (g : Fin L → ℕ) (j : ℕ) : ℕ := if h : j < L then g ⟨j, h⟩ else 0

theorem armVal_combMap_last (u : SpiderV m a → ℕ) (g : Fin L → ℕ) (j : ℕ) :
    armVal (combMap u g) (Fin.last m) j = gVal g j := by
  by_cases h : j < L
  · rw [armVal, dif_pos (show j < aExt m a L (Fin.last m) from by rw [aExt_last]; exact h),
      gVal, dif_pos h]
    show Sum.elim u g ((armEquiv m a L).symm (some ⟨Fin.last m, _⟩)) = _
    rw [armEquiv_symm_some_last]
    rfl
  · rw [armVal_of_ge _ _ (show aExt m a L (Fin.last m) ≤ j from by rw [aExt_last]; omega),
      gVal, dif_neg h]

theorem pref_combMap_castSucc (u : SpiderV m a → ℕ) (g : Fin L → ℕ) (i : Fin m) :
    pref (combMap u g) i.castSucc = pref u i := by
  rw [pref, pref, aExt_castSucc]
  congr 1
  funext j
  rw [armVal_combMap_castSucc]

/-- The active prefix of the extra arm. -/
def gPref (u : SpiderV m a → ℕ) (g : Fin L → ℕ) : ℕ :=
  pref (combMap u g) (Fin.last m)

theorem gPref_eq (u : SpiderV m a → ℕ) (g : Fin L → ℕ) :
    gPref u g = pref (combMap u g) (Fin.last m) := rfl

theorem gPref_le (u : SpiderV m a → ℕ) (g : Fin L → ℕ) : gPref u g ≤ L := by
  have := pref_le (combMap u g) (Fin.last m)
  rwa [aExt_last] at this

/-- **The number of odd arms of the extended spider**: the arms of the spider, plus the
extra arm if its active prefix is odd. -/
theorem pNum_combMap (u : SpiderV m a → ℕ) (g : Fin L → ℕ) :
    pNum (combMap u g) = pNum u + (if gPref u g % 2 = 1 then 1 else 0) := by
  classical
  have hsum : ∀ i : Fin m, (if pref (combMap u g) i.castSucc % 2 = 1 then 1 else 0)
      = (if pref u i % 2 = 1 then 1 else 0) := by
    intro i; rw [pref_combMap_castSucc]
  rw [pNum, pNum, Finset.card_filter, Finset.card_filter, Fin.sum_univ_castSucc,
    Finset.sum_congr rfl (fun i (_ : i ∈ Finset.univ) => hsum i), gPref_eq]

theorem tailW_combMap_castSucc (u : SpiderV m a → ℕ) (g : Fin L → ℕ) (i : Fin m) :
    tailW (combMap u g) i.castSucc = tailW u i := by
  have hlen : aExt m a L i.castSucc - pref (combMap u g) i.castSucc = a i - pref u i := by
    rw [aExt_castSucc, pref_combMap_castSucc]
  rw [tailW, tailW]
  refine (Wpoly_pathGraph_len_congr hlen
    (fun x => armVal (combMap u g) i.castSucc (pref (combMap u g) i.castSucc + x))).trans ?_
  refine congrArg _ (funext fun x => ?_)
  rw [armVal_combMap_castSucc, pref_combMap_castSucc]

theorem prod_tailW_combMap (u : SpiderV m a → ℕ) (g : Fin L → ℕ) :
    ∏ i : Fin (m + 1), tailW (combMap u g) i
      = (∏ i : Fin m, tailW u i) * tailW (combMap u g) (Fin.last m) := by
  rw [Fin.prod_univ_castSucc]
  congr 1
  exact Finset.prod_congr rfl fun i _ => tailW_combMap_castSucc u g i

/-! ### Splitting off the extra arm at an inactive hub -/

theorem Wpoly_armSpider_hub_zero (u : SpiderV m a → ℕ) (g : Fin L → ℕ) (h : u none = 0) :
    Wpoly (armSpider m a L) (Sum.elim u g)
      = Wpoly (spider m a) u * Wpoly (pathGraph L) g := by
  refine Wpoly_split_sum (armSpider m a L) (Sum.elim u g) (spider m a) (pathGraph L)
    (fun _ _ => Iff.rfl) (fun x y => by rw [pathGraph_adj]; exact Iff.rfl) ?_
  intro x y hadj
  exact Or.inl (by rw [Sum.elim_inl, (armSpider_adj_lr.mp hadj).1]; exact h)

/-! ### The weight of the extra arm -/

theorem Wpoly_extra_arm (u : SpiderV m a → ℕ) (g : Fin L → ℕ)
    (hstop : gPref u g < L → gVal g (gPref u g) = 0) :

    Wpoly (pathGraph L) g = prefW (gPref u g) * tailW (combMap u g) (Fin.last m) := by
  have key := Wpoly_arm (combMap u g) (Fin.last m) (by
    intro h
    rw [armVal_combMap_last]
    refine hstop ?_
    rwa [aExt_last] at h)
  rw [show (fun j : Fin (aExt m a L (Fin.last m)) => armVal (combMap u g) (Fin.last m) j.val)
      = (fun j : Fin (aExt m a L (Fin.last m)) => gVal g j.val) from
    funext fun j => armVal_combMap_last u g j.val] at key
  have hg : Wpoly (pathGraph L) g = Wpoly (pathGraph L) (fun j : Fin L => gVal g j.val) :=
    congrArg _ (funext fun j => by rw [gVal, dif_pos j.isLt])
  rw [gPref_eq, ← key, hg]
  exact (Wpoly_pathGraph_len_congr (aExt_last (m := m) (a := a) (L := L)) (gVal g)).symm

/-! ### The block of an `armSpider` -/

theorem admRun_of_admRun_combMap {u : SpiderV m a → ℕ} {g : Fin L → ℕ}
    (h : AdmRun (combMap u g)) : AdmRun u := by
  intro i j hpos
  have := h i.castSucc j (fun j' hj' => by
    rw [armVal_combMap_castSucc]
    exact hpos j' hj')
  rwa [armVal_combMap_castSucc] at this

/-- **The two-map block of an `armSpider`.**  Summed over the canonical block of the hub
side — the side itself and, if it is active, its canonical transform — the normalized
weight is centrally unimodal, for an *arbitrary* multiplicity pattern `g` on the extra
arm.  The extra arm is never touched by the transformation. -/
theorem armSpider_sideBlock_isCU (u : SpiderV m a → ℕ) (g : Fin L → ℕ) :
    IsCU (∑ u' ∈ sideBlock u, Wpoly (armSpider m a L) (Sum.elim u' g)) := by
  classical
  by_cases hu : ActiveSide u
  · rw [sum_sideBlock_active hu]
    obtain ⟨c, hc, hcu⟩ := Wpoly_side_image hu
    have hp : 2 ≤ pNum u := hu.two_odd
    have himg : Wpoly (armSpider m a L) (Sum.elim (transf u) g)
        = (c • (zz ^ (pNum u - 1) * ∏ i : Fin m, tailW u i)) * Wpoly (pathGraph L) g := by
      rw [Wpoly_armSpider_hub_zero _ _ (transf_none u), hcu]
    by_cases hadm : AdmRun (combMap u g)
    · have hact : Wpoly (armSpider m a L) (Sum.elim u g)
          = Pw ((1 - (pNum (combMap u g) : ℤ)).natAbs)
            * ∏ i : Fin (m + 1), tailW (combMap u g) i := by
        rw [Wpoly_armSpider_eq]
        exact Wpoly_side_active _ (by rw [combMap_none]; exact hu.hub)
          (fun i hi => hadm.stop i hi)
      have hextra : Wpoly (pathGraph L) g
          = prefW (gPref u g) * tailW (combMap u g) (Fin.last m) := by
        refine Wpoly_extra_arm u g (fun hlt => ?_)
        have := hadm.stop (Fin.last m) (by rw [aExt_last]; exact hlt)
        rwa [armVal_combMap_last] at this
      set T := ∏ i : Fin m, tailW u i with hT
      set tE := tailW (combMap u g) (Fin.last m) with htE
      have hgood : IsGoodW (T * tE) :=
        IsGoodW.mul (isGoodW_prod_tailW u) (isGoodW_tailW (combMap u g) (Fin.last m))
      rw [hact, prod_tailW_combMap, himg, hextra, pNum_combMap]
      rcases Nat.eq_zero_or_pos (gPref u g) with hrho | hrho
      · have key : Pw ((1 - ((pNum u + (if gPref u g % 2 = 1 then 1 else 0) : ℕ) : ℤ)).natAbs)
              * (T * tE) + c • (zz ^ (pNum u - 1) * T) * (prefW (gPref u g) * tE)
            = Ablock (pNum u - 1) c * (T * tE) := by
          rw [hrho, prefW_zero, if_neg (by omega), Ablock,
            show (1 - ((pNum u + 0 : ℕ) : ℤ)).natAbs = pNum u - 1 from by omega]
          simp only [smul_eq_C_mul]
          ring
        rw [key]
        exact (Ablock_isCU (by omega) hc).mul hgood.isCU
      · rcases Nat.mod_two_eq_zero_or_one (gPref u g) with hpar | hpar
        · have key : Pw ((1 - ((pNum u + (if gPref u g % 2 = 1 then 1 else 0) : ℕ) : ℤ)).natAbs)
                * (T * tE) + c • (zz ^ (pNum u - 1) * T) * (prefW (gPref u g) * tE)
              = Ablock (pNum u - 1) (2 * c) * (T * tE) := by
            rw [prefW_even hpar (by omega), if_neg (by omega), Ablock,
              show (1 - ((pNum u + 0 : ℕ) : ℤ)).natAbs = pNum u - 1 from by omega]
            simp only [smul_eq_C_mul, map_mul, mul_one]
            ring
          rw [key]
          exact (Ablock_isCU (by omega) (by linarith)).mul hgood.isCU
        · have key : Pw ((1 - ((pNum u + (if gPref u g % 2 = 1 then 1 else 0) : ℕ) : ℤ)).natAbs)
                * (T * tE) + c • (zz ^ (pNum u - 1) * T) * (prefW (gPref u g) * tE)
              = Ablock (pNum u) c * (T * tE) := by
            have hzz : zz ^ (pNum u - 1) * zz = zz ^ pNum u := by
              rw [← pow_succ]
              congr 1
              omega
            rw [prefW_odd hpar, if_pos hpar, Ablock,
              show (1 - ((pNum u + 1 : ℕ) : ℤ)).natAbs = pNum u from by omega]
            simp only [smul_eq_C_mul]
            rw [← hzz]
            ring
          rw [key]
          exact (Ablock_isCU (by omega) hc).mul hgood.isCU
    · have h0 : Wpoly (armSpider m a L) (Sum.elim u g) = 0 := by
        rw [Wpoly_armSpider_eq]
        exact Wpoly_spider_eq_zero_of_not_admRun
          (by rw [combMap_none, hu.hub]; norm_num) hadm
      have key : Wpoly (armSpider m a L) (Sum.elim u g)
            + Wpoly (armSpider m a L) (Sum.elim (transf u) g)
          = c • (zz ^ (pNum u - 1)
              * ((∏ i : Fin m, tailW u i) * Wpoly (pathGraph L) g)) := by
        rw [h0, zero_add, himg]
        simp only [smul_eq_C_mul]
        ring
      rw [key]
      exact ((isCU_zz_pow _).mul
        (IsGoodW.mul (isGoodW_prod_tailW u) (Wpoly_pathGraph_isGood L g)).isCU).smul
        (by linarith)
  · rw [sum_sideBlock_not_active hu, Wpoly_armSpider_eq]
    by_cases hact : ActiveSide (combMap u g)
    · have h1 : pNum u ≤ 1 := by
        by_contra hcon
        exact hu ⟨by rw [← combMap_none u g]; exact hact.hub,
          admRun_of_admRun_combMap hact.run, by omega⟩
      have h2 : pNum (combMap u g) ≤ 2 := by
        rw [pNum_combMap]
        split <;> omega
      rw [Wpoly_side_active _ hact.hub (fun i hi => hact.run.stop i hi)]
      exact (IsGoodW.mul (isGoodW_Pw_le_one (by omega))
        (isGoodW_prod_tailW (combMap u g))).isCU
    · exact (isGoodW_Wpoly_spider_of_not_active hact).isCU

end ClanAudit
