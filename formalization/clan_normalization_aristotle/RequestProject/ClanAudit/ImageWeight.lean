import RequestProject.ClanAudit.SideTransform

/-!
# The normalized weight of a transformed side

`Wpoly_side_image`: if the side `s` is active with `p` arms of odd active prefix, the
image `transf s` has normalized weight `c * (z+z⁻¹)^(p-1) * ∏ tails` with a *derived*
scalar `c = 2^e ≥ 1`, `e` the number of arms with even positive active prefix.

Together with `Wpoly_side_active` (weight `(z^(p-1) + z^(-(p-1))) * ∏ tails` for the
active side itself) this produces the one-hub block `A(p-1, 2^e; z)`.
-/

namespace ClanAudit

open Finset LaurentPolynomial SimpleGraph

variable {k : ℕ} {len : Fin k → ℕ}

/-- Splitting an arm at a position whose value, or whose predecessor's value, is zero. -/
theorem Wpoly_arm_split (g : ℕ → ℕ) (N r t : ℕ) (hN : N = r + t)
    (hcut : 0 < r → 0 < t → (g (r - 1) = 0 ∨ g r = 0)) :
    Wpoly (pathGraph N) (fun j : Fin N => g j.val)
      = Wpoly (pathGraph r) (fun j : Fin r => g j.val)
        * Wpoly (pathGraph t) (fun j : Fin t => g (r + j.val)) := by
  rw [Wpoly_pathGraph_split' r t hN _ ?_]
  · rfl
  · intro x y hxy
    have hx : x.val = r - 1 := by have := x.isLt; omega
    have hy : y.val = 0 := by have := x.isLt; omega
    have i1 : (finCongr hN.symm (Fin.castAdd t x)).val = r - 1 := hx
    have i2 : (finCongr hN.symm (Fin.natAdd r y)).val = r := by
      show r + y.val = r
      omega
    rw [i1, i2]
    exact hcut (by have := x.isLt; omega) (by have := y.isLt; omega)

/-- The arm weight only depends on the arm. -/
theorem Wpoly_arm_congr {s s' : SpiderV k len → ℕ} {i : Fin k}
    (heq : ∀ j, armVal s i j = armVal s' i j) :
    Wpoly (pathGraph (len i)) (fun j : Fin (len i) => armVal s i j.val)
      = Wpoly (pathGraph (len i)) (fun j : Fin (len i) => armVal s' i j.val) :=
  congrArg _ (funext fun j => heq j.val)

/-- **The normalized weight of a transformed side.**  The derived scalar is a power of
two, hence at least one — the hypothesis under which `Ablock_isCU` and `Fblock_decr`
apply. -/
theorem Wpoly_side_image {s : SpiderV k len → ℕ} (h : ActiveSide s) :
    ∃ c : ℚ, 1 ≤ c ∧
      Wpoly (spider k len) (transf s) = c • (zz ^ (pNum s - 1) * ∏ i : Fin k, tailW s i) := by
  classical
  obtain ⟨i₀, i₁, hidx0, hidx1, hplen, hodd0, hodd1, hne, hle⟩ := canonical_spec h
  have hL0 : pref s i₀ ≤ len i₀ := pref_le s i₀
  have hL1 : pref s i₀ ≤ pref s i₁ := by rw [hplen]; exact hle
  have hP1 : pref s i₁ ≤ len i₁ := pref_le s i₁
  have hne' : i₀.val ≠ i₁.val := fun hc => hne (Fin.ext hc)
  have hnei0 : i₀.val ≠ idx1 s := by rw [← hidx1]; exact hne'
  have hnei1 : i₁.val ≠ idx0 s := by rw [← hidx0]; exact fun hc => hne' hc.symm
  -- the transformed arm values
  have hval0 : ∀ j, j < pref s i₀ → armVal (transf s) i₀ j = altVal j := by
    intro j hj
    rw [transf, armVal_transfSide s (idx0 s) (idx1 s) (plen s) i₀ j (by omega),
      if_pos ⟨hidx0, by omega⟩]
  have hval0' : ∀ j, pref s i₀ ≤ j → armVal (transf s) i₀ j = armVal s i₀ j := by
    intro j hj
    by_cases hlt : j < len i₀
    · rw [transf, armVal_transfSide s (idx0 s) (idx1 s) (plen s) i₀ j hlt,
        if_neg (by omega), if_neg (by tauto)]
    · rw [armVal_of_ge _ _ (by omega), armVal_of_ge _ _ (by omega)]
  have hval1 : ∀ j, j < pref s i₀ - 1 → armVal (transf s) i₁ j = altVal j := by
    intro j hj
    rw [transf, armVal_transfSide s (idx0 s) (idx1 s) (plen s) i₁ j (by omega),
      if_neg (by tauto), if_pos ⟨hidx1, by omega⟩]
  have hval1' : ∀ j, pref s i₀ - 1 ≤ j → armVal (transf s) i₁ j = armVal s i₁ j := by
    intro j hj
    by_cases hlt : j < len i₁
    · rw [transf, armVal_transfSide s (idx0 s) (idx1 s) (plen s) i₁ j hlt,
        if_neg (by tauto), if_neg (by omega)]
    · rw [armVal_of_ge _ _ (by omega), armVal_of_ge _ _ (by omega)]
  -- the factor of the arm A: the rewritten prefix contributes 1
  have hfac0 : Wpoly (pathGraph (len i₀)) (fun j : Fin (len i₀) => armVal (transf s) i₀ j.val)
      = tailW s i₀ := by
    rw [Wpoly_arm_split _ (len i₀) (pref s i₀) (len i₀ - pref s i₀) (by omega) ?_]
    · have e1 : (fun j : Fin (pref s i₀) => armVal (transf s) i₀ j.val)
          = fun j : Fin (pref s i₀) => altVal j.val :=
        funext fun j => hval0 j.val j.isLt
      rw [e1, Wpoly_pathGraph_alt, one_mul]
      exact congrArg _ (funext fun j => hval0' _ (by omega))
    · intro _ h2
      right
      rw [hval0' _ le_rfl]
      exact h.run.stop i₀ (by omega)
  -- the factor of the arm B: the remaining odd run of ones contributes `z + z⁻¹`
  have hfac1 : Wpoly (pathGraph (len i₁)) (fun j : Fin (len i₁) => armVal (transf s) i₁ j.val)
      = zz * tailW s i₁ := by
    rw [Wpoly_arm_split _ (len i₁) (pref s i₀ - 1) (len i₁ - pref s i₀ + 1) (by omega) ?_]
    · have e1 : (fun j : Fin (pref s i₀ - 1) => armVal (transf s) i₁ j.val)
          = fun j : Fin (pref s i₀ - 1) => altVal j.val :=
        funext fun j => hval1 j.val j.isLt
      rw [e1, Wpoly_pathGraph_alt, one_mul]
      rw [Wpoly_arm_split (fun x => armVal (transf s) i₁ (pref s i₀ - 1 + x))
        (len i₁ - pref s i₀ + 1) (pref s i₁ - pref s i₀ + 1)
        (len i₁ - pref s i₁) (by omega) ?_]
      · have e2 : (fun j : Fin (pref s i₁ - pref s i₀ + 1) =>
              armVal (transf s) i₁ (pref s i₀ - 1 + j.val)) = fun _ => 1 := by
          funext j
          rw [hval1' _ (by omega)]
          exact pref_ones s i₁ (by have := j.isLt; omega)
        have e3 : (fun j : Fin (len i₁ - pref s i₁) =>
              armVal (transf s) i₁ (pref s i₀ - 1 + (pref s i₁ - pref s i₀ + 1 + j.val)))
            = fun j : Fin (len i₁ - pref s i₁) => armVal s i₁ (pref s i₁ + j.val) := by
          funext j
          rw [hval1' _ (by omega)]
          congr 1
          omega
        rw [e2, e3, Wpoly_prefix_ones _ _ (fun _ => rfl),
          prefW_odd (show (pref s i₁ - pref s i₀ + 1) % 2 = 1 by omega)]
        rfl
      · intro _ h2
        right
        show armVal (transf s) i₁ (pref s i₀ - 1 + (pref s i₁ - pref s i₀ + 1)) = 0
        rw [hval1' _ (by omega),
          show pref s i₀ - 1 + (pref s i₁ - pref s i₀ + 1) = pref s i₁ from by omega]
        exact h.run.stop i₁ (by omega)
    · intro h1 _
      left
      rw [hval1 _ (by omega), altVal, if_neg (by omega)]
  -- the untouched arms
  have hfacO : ∀ i : Fin k, i ≠ i₀ → i ≠ i₁ →
      Wpoly (pathGraph (len i)) (fun j : Fin (len i) => armVal (transf s) i j.val)
        = prefW (pref s i) * tailW s i := by
    intro i h0 h1
    have hne0 : i.val ≠ idx0 s := fun hc => h0 (Fin.ext (hc.trans hidx0.symm))
    have hne1 : i.val ≠ idx1 s := fun hc => h1 (Fin.ext (hc.trans hidx1.symm))
    rw [Wpoly_arm_congr (s := transf s) (s' := s)
      (fun j => armVal_transfSide_other s (i₀ := idx0 s) (i₁ := idx1 s) (L := plen s)
        hne0 hne1 j)]
    exact Wpoly_arm s i (h.run.stop i)
  -- assemble the product over all arms
  set rest : Finset (Fin k) := (Finset.univ.erase i₀).erase i₁ with hrest
  have hi₁mem : i₁ ∈ Finset.univ.erase i₀ :=
    Finset.mem_erase.mpr ⟨fun hc => hne hc.symm, mem_univ _⟩
  have hsplitP : ∀ f : Fin k → LaurentPolynomial ℚ,
      ∏ i : Fin k, f i = f i₀ * (f i₁ * ∏ i ∈ rest, f i) := by
    intro f
    rw [hrest, Finset.mul_prod_erase _ f hi₁mem, Finset.mul_prod_erase _ f (mem_univ i₀)]
  obtain ⟨e, he⟩ := prod_prefW s rest
  have hcard : (rest.filter (fun i => pref s i % 2 = 1)).card = pNum s - 2 := by
    have hm0 : i₀ ∈ Finset.univ.filter (fun i : Fin k => pref s i % 2 = 1) :=
      Finset.mem_filter.mpr ⟨mem_univ _, hodd0⟩
    have hm1 : i₁ ∈ (Finset.univ.filter (fun i : Fin k => pref s i % 2 = 1)).erase i₀ :=
      Finset.mem_erase.mpr ⟨fun hc => hne hc.symm, Finset.mem_filter.mpr ⟨mem_univ _, hodd1⟩⟩
    rw [hrest, Finset.filter_erase, Finset.filter_erase, Finset.card_erase_of_mem hm1,
      Finset.card_erase_of_mem hm0, pNum]
    omega
  have hrestprod : ∏ i ∈ rest, Wpoly (pathGraph (len i))
      (fun j : Fin (len i) => armVal (transf s) i j.val)
      = (∏ i ∈ rest, prefW (pref s i)) * ∏ i ∈ rest, tailW s i := by
    rw [← Finset.prod_mul_distrib]
    exact Finset.prod_congr rfl fun i hi => hfacO i
      (Finset.mem_erase.mp (Finset.mem_erase.mp hi).2).1 (Finset.mem_erase.mp hi).1
  have hp2 : 2 ≤ pNum s := h.two_odd
  have hzz : zz ^ (pNum s - 1) = zz * zz ^ (pNum s - 2) := by
    rw [show pNum s - 1 = (pNum s - 2) + 1 from by omega, pow_succ, mul_comm]
  refine ⟨(2 : ℚ) ^ e, one_le_pow₀ (by norm_num), ?_⟩
  rw [Wpoly_spider_hub (transf s) (fun i _ => Or.inl rfl)]
  have hdot : dotW (transf s none) = 1 := rfl
  rw [hdot, mul_one, hsplitP (fun i => Wpoly (pathGraph (len i))
    (fun j : Fin (len i) => armVal (transf s) i j.val)), hfac0, hfac1, hrestprod, he, hcard,
    hsplitP (fun i => tailW s i), hzz, smul_eq_C_mul, smul_eq_C_mul]
  ring

end ClanAudit
