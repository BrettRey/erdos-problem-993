import RequestProject.ClanAudit.Vanishing

/-!
# The exact normalized weight of each block of the partition

The global partition of the clan maps of the adjacent two-hub tree (built in
`RequestProject/ClanAudit/GlobalPartition.lean`) has three kinds of blocks.  Here we
compute the exact normalized Laurent sum of each of them and show it is centrally
unimodal.

* `blockSum_both_active` — the four-map block is `F(r,s;c,d)` times a product of good
  factors, with `r, s ≥ 1` and *derived* scalars `c, d = 2^e ≥ 1`;
* `blockSum_left_active`, `blockSum_right_active` — the two-map block is `A(r,c)` times a
  product of good factors, or a good multiple of `(z+z⁻¹)^r`, or zero;
* `blockSum_neither_active` — the singleton weight is a product of good factors, or zero.

All outside component factors (the arm tails and the whole inactive side) are carried
along explicitly; no factor is postulated.
-/

namespace ClanAudit

open Finset LaurentPolynomial SimpleGraph

/-! ### Preliminaries -/

variable {k : ℕ} {len : Fin k → ℕ}

theorem transf_none (s : SpiderV k len → ℕ) : transf s none = 0 := rfl

/-- The weight of an active side, with the imbalance written as `p - 1`. -/
theorem Wpoly_spider_active {s : SpiderV k len → ℕ} (h : ActiveSide s) :
    Wpoly (spider k len) s = Pw (pNum s - 1) * ∏ i : Fin k, tailW s i := by
  rw [Wpoly_side_active s h.hub (fun i hi => h.run.stop i hi)]
  have h2 := h.two_odd
  rw [show (1 - (pNum s : ℤ)).natAbs = pNum s - 1 from by omega]

/-- The four-map identity, with the last term written with an absolute value so that no
ordering of `r` and `s` is assumed. -/
theorem Fblock_expand_abs (r s : ℕ) (c d : ℚ) :
    Fblock r s c d = (c * d) • zz ^ (r + s) + c • (zz ^ r * Pw s) + d • (zz ^ s * Pw r)
      + Pw (((r : ℤ) - s).natAbs) := by
  rcases le_total s r with h | h
  · rw [Fblock_expand h, show ((r : ℤ) - s).natAbs = r - s from by omega]
  · rw [Fblock_comm, Fblock_expand h, show ((r : ℤ) - s).natAbs = s - r from by omega,
      Nat.add_comm s r, mul_comm d c]
    abel

/-! ### The four-map block -/

variable {m n : ℕ} {a : Fin m → ℕ} {b : Fin n → ℕ}

/-- Both hubs active: the two hubs and all active prefixes form one component of
imbalance `|r - s|`, `r = p - 1`, `s = q - 1`. -/
theorem Wpoly_dgraph_both_active {u : SpiderV m a → ℕ} {v : SpiderV n b → ℕ}
    (hu : ActiveSide u) (hv : ActiveSide v) :
    Wpoly (dgraph m a n b) (Sum.elim u v)
      = Pw ((((pNum u - 1 : ℕ) : ℤ) - ((pNum v - 1 : ℕ) : ℤ)).natAbs)
        * ((∏ i : Fin m, tailW u i) * ∏ j : Fin n, tailW v j) := by
  rw [Wpoly_dgraph_joint u v hu.hub hv.hub (fun i hi => hu.run.stop i hi)
    (fun j hj => hv.run.stop j hj)]
  have h1 := hu.two_odd
  have h2 := hv.two_odd
  rw [show ((pNum u : ℤ) - (pNum v : ℤ)).natAbs
      = (((pNum u - 1 : ℕ) : ℤ) - ((pNum v - 1 : ℕ) : ℤ)).natAbs from by omega]

/-- **The four-map block.**  Its sum is `F(r,s;c,d)` times a product of good factors,
with `r, s ≥ 1` and derived scalars `c, d ≥ 1`; hence it is centrally unimodal. -/
theorem blockSum_both_active {u : SpiderV m a → ℕ} {v : SpiderV n b → ℕ}
    (hu : ActiveSide u) (hv : ActiveSide v) :
    IsCU (Wpoly (dgraph m a n b) (Sum.elim u v)
      + Wpoly (dgraph m a n b) (Sum.elim u (transf v))
      + (Wpoly (dgraph m a n b) (Sum.elim (transf u) v)
        + Wpoly (dgraph m a n b) (Sum.elim (transf u) (transf v)))) := by
  obtain ⟨cu, hcu, heu⟩ := Wpoly_side_image hu
  obtain ⟨cv, hcv, hev⟩ := Wpoly_side_image hv
  have hpu : 2 ≤ pNum u := hu.two_odd
  have hpv : 2 ≤ pNum v := hv.two_odd
  have h1 := Wpoly_dgraph_both_active hu hv
  have h2 : Wpoly (dgraph m a n b) (Sum.elim u (transf v))
      = (Pw (pNum u - 1) * ∏ i : Fin m, tailW u i)
        * (cv • (zz ^ (pNum v - 1) * ∏ j : Fin n, tailW v j)) := by
    rw [Wpoly_dgraph_split u (transf v) (Or.inr (transf_none v)), Wpoly_spider_active hu, hev]
  have h3 : Wpoly (dgraph m a n b) (Sum.elim (transf u) v)
      = (cu • (zz ^ (pNum u - 1) * ∏ i : Fin m, tailW u i))
        * (Pw (pNum v - 1) * ∏ j : Fin n, tailW v j) := by
    rw [Wpoly_dgraph_split (transf u) v (Or.inl (transf_none u)), Wpoly_spider_active hv, heu]
  have h4 : Wpoly (dgraph m a n b) (Sum.elim (transf u) (transf v))
      = (cu • (zz ^ (pNum u - 1) * ∏ i : Fin m, tailW u i))
        * (cv • (zz ^ (pNum v - 1) * ∏ j : Fin n, tailW v j)) := by
    rw [Wpoly_dgraph_split (transf u) (transf v) (Or.inl (transf_none u)), heu, hev]
  have key : Wpoly (dgraph m a n b) (Sum.elim u v)
      + Wpoly (dgraph m a n b) (Sum.elim u (transf v))
      + (Wpoly (dgraph m a n b) (Sum.elim (transf u) v)
        + Wpoly (dgraph m a n b) (Sum.elim (transf u) (transf v)))
      = Fblock (pNum u - 1) (pNum v - 1) cu cv
        * ((∏ i : Fin m, tailW u i) * ∏ j : Fin n, tailW v j) := by
    rw [h1, h2, h3, h4, Fblock_expand_abs]
    simp only [smul_eq_C_mul, map_mul]
    ring
  rw [key]
  exact (Fblock_isCU (by omega) (by omega) hcu hcv).mul
    ((IsGoodW.mul (isGoodW_prod_tailW u) (isGoodW_prod_tailW v)).isCU)

/-! ### The two-map blocks -/

/-- **The two-map block at the left hub.** -/
theorem blockSum_left_active {u : SpiderV m a → ℕ} {v : SpiderV n b → ℕ}
    (hu : ActiveSide u) (hv : ¬ ActiveSide v) :
    IsCU (Wpoly (dgraph m a n b) (Sum.elim u v)
      + Wpoly (dgraph m a n b) (Sum.elim (transf u) v)) := by
  obtain ⟨cu, hcu, heu⟩ := Wpoly_side_image hu
  have hpu : 2 ≤ pNum u := hu.two_odd
  have h3 : Wpoly (dgraph m a n b) (Sum.elim (transf u) v)
      = (cu • (zz ^ (pNum u - 1) * ∏ i : Fin m, tailW u i)) * Wpoly (spider n b) v := by
    rw [Wpoly_dgraph_split (transf u) v (Or.inl (transf_none u)), heu]
  by_cases hv1 : v none = 1
  · by_cases hrun : AdmRun v
    · have hq : pNum v ≤ 1 := by
        by_contra hc
        exact hv ⟨hv1, hrun, by omega⟩
      have h1 : Wpoly (dgraph m a n b) (Sum.elim u v)
          = Pw (((pNum u : ℤ) - (pNum v : ℤ)).natAbs)
            * ((∏ i : Fin m, tailW u i) * ∏ j : Fin n, tailW v j) :=
        Wpoly_dgraph_joint u v hu.hub hv1 (fun i hi => hu.run.stop i hi)
          (fun j hj => hrun.stop j hj)
      have hvw : Wpoly (spider n b) v
          = Pw ((1 - (pNum v : ℤ)).natAbs) * ∏ j : Fin n, tailW v j :=
        Wpoly_side_active v hv1 (fun j hj => hrun.stop j hj)
      rcases (show pNum v = 0 ∨ pNum v = 1 from by omega) with hq0 | hq1
      · have key : Wpoly (dgraph m a n b) (Sum.elim u v)
            + Wpoly (dgraph m a n b) (Sum.elim (transf u) v)
            = Ablock (pNum u - 1 + 1) cu
              * ((∏ i : Fin m, tailW u i) * ∏ j : Fin n, tailW v j) := by
          rw [h1, h3, hvw, hq0, Ablock]
          rw [show ((pNum u : ℤ) - ((0 : ℕ) : ℤ)).natAbs = pNum u - 1 + 1 from by omega,
            show (1 - ((0 : ℕ) : ℤ)).natAbs = 1 from by omega, Pw_one, pow_succ]
          simp only [smul_eq_C_mul]
          ring
        rw [key]
        exact (Ablock_isCU (by omega) hcu).mul
          ((IsGoodW.mul (isGoodW_prod_tailW u) (isGoodW_prod_tailW v)).isCU)
      · have key : Wpoly (dgraph m a n b) (Sum.elim u v)
            + Wpoly (dgraph m a n b) (Sum.elim (transf u) v)
            = Ablock (pNum u - 1) (2 * cu)
              * ((∏ i : Fin m, tailW u i) * ∏ j : Fin n, tailW v j) := by
          rw [h1, h3, hvw, hq1, Ablock]
          rw [show ((pNum u : ℤ) - ((1 : ℕ) : ℤ)).natAbs = pNum u - 1 from by omega,
            show (1 - ((1 : ℕ) : ℤ)).natAbs = 0 from by omega, Pw_zero]
          simp only [smul_eq_C_mul, map_mul]
          norm_num [C_two_eq]
          ring
        rw [key]
        exact (Ablock_isCU (by omega) (by linarith)).mul
          ((IsGoodW.mul (isGoodW_prod_tailW u) (isGoodW_prod_tailW v)).isCU)
    · have z1 : Wpoly (dgraph m a n b) (Sum.elim u v) = 0 :=
        Wpoly_dgraph_eq_zero_of_not_admRun_right u v (by omega) hrun
      have z2 : Wpoly (dgraph m a n b) (Sum.elim (transf u) v) = 0 :=
        Wpoly_dgraph_eq_zero_of_not_admRun_right (transf u) v (by omega) hrun
      rw [z1, z2, add_zero]
      exact isCU_zero
  · have hgood : IsGoodW (Wpoly (spider n b) v) := isGoodW_Wpoly_spider_of_hub_ne_one hv1
    by_cases hv0 : v none = 0
    · have h1 : Wpoly (dgraph m a n b) (Sum.elim u v)
          = (Pw (pNum u - 1) * ∏ i : Fin m, tailW u i) * Wpoly (spider n b) v := by
        rw [Wpoly_dgraph_split u v (Or.inr hv0), Wpoly_spider_active hu]
      have key : Wpoly (dgraph m a n b) (Sum.elim u v)
          + Wpoly (dgraph m a n b) (Sum.elim (transf u) v)
          = Ablock (pNum u - 1) cu
            * ((∏ i : Fin m, tailW u i) * Wpoly (spider n b) v) := by
        rw [h1, h3, Ablock]
        simp only [smul_eq_C_mul]
        ring
      rw [key]
      exact (Ablock_isCU (by omega) hcu).mul
        ((IsGoodW.mul (isGoodW_prod_tailW u) hgood).isCU)
    · have h1 : Wpoly (dgraph m a n b) (Sum.elim u v) = 0 :=
        Wpoly_dgraph_eq_zero_of_hub_two_right u v (by rw [hu.hub]) (by omega)
      have key : Wpoly (dgraph m a n b) (Sum.elim u v)
          + Wpoly (dgraph m a n b) (Sum.elim (transf u) v)
          = cu • (zz ^ (pNum u - 1)
              * ((∏ i : Fin m, tailW u i) * Wpoly (spider n b) v)) := by
        rw [h1, h3, zero_add]
        simp only [smul_eq_C_mul]
        ring
      rw [key]
      exact (((isCU_zz_pow _).mul
        ((IsGoodW.mul (isGoodW_prod_tailW u) hgood).isCU)).smul (by linarith))

/-- **The two-map block at the right hub.** -/
theorem blockSum_right_active {u : SpiderV m a → ℕ} {v : SpiderV n b → ℕ}
    (hu : ¬ ActiveSide u) (hv : ActiveSide v) :
    IsCU (Wpoly (dgraph m a n b) (Sum.elim u v)
      + Wpoly (dgraph m a n b) (Sum.elim u (transf v))) := by
  obtain ⟨cv, hcv, hev⟩ := Wpoly_side_image hv
  have hpv : 2 ≤ pNum v := hv.two_odd
  have h3 : Wpoly (dgraph m a n b) (Sum.elim u (transf v))
      = Wpoly (spider m a) u * (cv • (zz ^ (pNum v - 1) * ∏ j : Fin n, tailW v j)) := by
    rw [Wpoly_dgraph_split u (transf v) (Or.inr (transf_none v)), hev]
  by_cases hu1 : u none = 1
  · by_cases hrun : AdmRun u
    · have hq : pNum u ≤ 1 := by
        by_contra hc
        exact hu ⟨hu1, hrun, by omega⟩
      have h1 : Wpoly (dgraph m a n b) (Sum.elim u v)
          = Pw (((pNum u : ℤ) - (pNum v : ℤ)).natAbs)
            * ((∏ i : Fin m, tailW u i) * ∏ j : Fin n, tailW v j) :=
        Wpoly_dgraph_joint u v hu1 hv.hub (fun i hi => hrun.stop i hi)
          (fun j hj => hv.run.stop j hj)
      have huw : Wpoly (spider m a) u
          = Pw ((1 - (pNum u : ℤ)).natAbs) * ∏ i : Fin m, tailW u i :=
        Wpoly_side_active u hu1 (fun i hi => hrun.stop i hi)
      rcases (show pNum u = 0 ∨ pNum u = 1 from by omega) with hq0 | hq1
      · have key : Wpoly (dgraph m a n b) (Sum.elim u v)
            + Wpoly (dgraph m a n b) (Sum.elim u (transf v))
            = Ablock (pNum v - 1 + 1) cv
              * ((∏ i : Fin m, tailW u i) * ∏ j : Fin n, tailW v j) := by
          rw [h1, h3, huw, hq0, Ablock]
          rw [show (((0 : ℕ) : ℤ) - (pNum v : ℤ)).natAbs = pNum v - 1 + 1 from by omega,
            show (1 - ((0 : ℕ) : ℤ)).natAbs = 1 from by omega, Pw_one, pow_succ]
          simp only [smul_eq_C_mul]
          ring
        rw [key]
        exact (Ablock_isCU (by omega) hcv).mul
          ((IsGoodW.mul (isGoodW_prod_tailW u) (isGoodW_prod_tailW v)).isCU)
      · have key : Wpoly (dgraph m a n b) (Sum.elim u v)
            + Wpoly (dgraph m a n b) (Sum.elim u (transf v))
            = Ablock (pNum v - 1) (2 * cv)
              * ((∏ i : Fin m, tailW u i) * ∏ j : Fin n, tailW v j) := by
          rw [h1, h3, huw, hq1, Ablock]
          rw [show (((1 : ℕ) : ℤ) - (pNum v : ℤ)).natAbs = pNum v - 1 from by omega,
            show (1 - ((1 : ℕ) : ℤ)).natAbs = 0 from by omega, Pw_zero]
          simp only [smul_eq_C_mul, map_mul]
          norm_num [C_two_eq]
          ring
        rw [key]
        exact (Ablock_isCU (by omega) (by linarith)).mul
          ((IsGoodW.mul (isGoodW_prod_tailW u) (isGoodW_prod_tailW v)).isCU)
    · have z1 : Wpoly (dgraph m a n b) (Sum.elim u v) = 0 :=
        Wpoly_dgraph_eq_zero_of_not_admRun_left u v (by omega) hrun
      have z2 : Wpoly (dgraph m a n b) (Sum.elim u (transf v)) = 0 :=
        Wpoly_dgraph_eq_zero_of_not_admRun_left u (transf v) (by omega) hrun
      rw [z1, z2, add_zero]
      exact isCU_zero
  · have hgood : IsGoodW (Wpoly (spider m a) u) := isGoodW_Wpoly_spider_of_hub_ne_one hu1
    by_cases hu0 : u none = 0
    · have h1 : Wpoly (dgraph m a n b) (Sum.elim u v)
          = Wpoly (spider m a) u * (Pw (pNum v - 1) * ∏ j : Fin n, tailW v j) := by
        rw [Wpoly_dgraph_split u v (Or.inl hu0), Wpoly_spider_active hv]
      have key : Wpoly (dgraph m a n b) (Sum.elim u v)
          + Wpoly (dgraph m a n b) (Sum.elim u (transf v))
          = Ablock (pNum v - 1) cv
            * (Wpoly (spider m a) u * ∏ j : Fin n, tailW v j) := by
        rw [h1, h3, Ablock]
        simp only [smul_eq_C_mul]
        ring
      rw [key]
      exact (Ablock_isCU (by omega) hcv).mul
        ((IsGoodW.mul hgood (isGoodW_prod_tailW v)).isCU)
    · have h1 : Wpoly (dgraph m a n b) (Sum.elim u v) = 0 :=
        Wpoly_dgraph_eq_zero_of_hub_two_left u v (by omega) (by rw [hv.hub])
      have key : Wpoly (dgraph m a n b) (Sum.elim u v)
          + Wpoly (dgraph m a n b) (Sum.elim u (transf v))
          = cv • (zz ^ (pNum v - 1)
              * (Wpoly (spider m a) u * ∏ j : Fin n, tailW v j)) := by
        rw [h1, h3, zero_add]
        simp only [smul_eq_C_mul]
        ring
      rw [key]
      exact (((isCU_zz_pow _).mul
        ((IsGoodW.mul hgood (isGoodW_prod_tailW v)).isCU)).smul (by linarith))

/-! ### The singletons -/

/-- **A singleton block.**  Neither hub is active, and the weight is a product of good
factors (or zero): centrally unimodal, and in fact of the good shape. -/
theorem blockSum_neither_active {u : SpiderV m a → ℕ} {v : SpiderV n b → ℕ}
    (hu : ¬ ActiveSide u) (hv : ¬ ActiveSide v) :
    IsGoodW (Wpoly (dgraph m a n b) (Sum.elim u v)) := by
  by_cases hu0 : u none = 0
  · rw [Wpoly_dgraph_split u v (Or.inl hu0)]
    exact IsGoodW.mul (isGoodW_Wpoly_spider_of_not_active hu)
      (isGoodW_Wpoly_spider_of_not_active hv)
  by_cases hv0 : v none = 0
  · rw [Wpoly_dgraph_split u v (Or.inr hv0)]
    exact IsGoodW.mul (isGoodW_Wpoly_spider_of_not_active hu)
      (isGoodW_Wpoly_spider_of_not_active hv)
  by_cases hu2 : 2 ≤ u none
  · exact Or.inl (Wpoly_dgraph_eq_zero_of_hub_two_left u v hu2 (by omega))
  by_cases hv2 : 2 ≤ v none
  · exact Or.inl (Wpoly_dgraph_eq_zero_of_hub_two_right u v (by omega) hv2)
  have hu1 : u none = 1 := by omega
  have hv1 : v none = 1 := by omega
  by_cases hru : AdmRun u
  · by_cases hrv : AdmRun v
    · have hpu : pNum u ≤ 1 := by
        by_contra hc
        exact hu ⟨hu1, hru, by omega⟩
      have hpv : pNum v ≤ 1 := by
        by_contra hc
        exact hv ⟨hv1, hrv, by omega⟩
      rw [Wpoly_dgraph_joint u v hu1 hv1 (fun i hi => hru.stop i hi) (fun j hj => hrv.stop j hj)]
      exact IsGoodW.mul (isGoodW_Pw_le_one (by omega))
        (IsGoodW.mul (isGoodW_prod_tailW u) (isGoodW_prod_tailW v))
    · exact Or.inl (Wpoly_dgraph_eq_zero_of_not_admRun_right u v (by omega) hrv)
  · exact Or.inl (Wpoly_dgraph_eq_zero_of_not_admRun_left u v (by omega) hru)

end ClanAudit
