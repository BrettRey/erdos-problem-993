import RequestProject.ClanAudit.ConnectorJoint
import RequestProject.ClanAudit.EvenBlock

/-!
# The exact normalized weight of each block of the connector partition

The blocks of the partition of the clan maps of the connector tree are indexed by the
two hub sides; the connector itself is never touched.  This file computes the exact
normalized Laurent sum of every block and shows it is centrally unimodal.

The classification of the connector states is *proved*, not assumed:

1. if some connector vertex is switched off, the tree splits into an `armSpider` and a
   mirrored `armSpider`, and the block sum is the product of the two arm-spider block
   sums (`armSpider_sideBlock_isCU`, `armSpiderR_sideBlock_isCU`);
2. if every connector vertex is positive but some carries multiplicity at least two,
   every map with a positive hub gives a triangle and hence weight zero, and the
   surviving maps have both hubs switched off;
3. if the connector is unit-positive end to end, the four-state analysis applies:
   for an odd number of connector edges the four-map block is `Fblock` plus an explicit
   good summand, and for an even number of connector edges it is exactly the
   even-connector block `Gblock = EvenConnector.Gpoly`.
-/

namespace ClanAudit

open Finset LaurentPolynomial SimpleGraph

variable {t m n : ℕ} {a : Fin m → ℕ} {b : Fin n → ℕ}

/-! ### The hub-switched-off part of a side block -/

/-- The part of a side block consisting of the maps with a switched-off hub is
centrally unimodal. -/
theorem sideBlock_hubZero_isCU {k : ℕ} {len : Fin k → ℕ} (s : SpiderV k len → ℕ) :
    IsCU (∑ s' ∈ sideBlock s, if s' none = 0 then Wpoly (spider k len) s' else 0) := by
  classical
  by_cases hact : ActiveSide s
  · rw [sum_sideBlock_active hact, if_neg (by rw [hact.hub]; omega),
      if_pos (transf_none s), zero_add]
    obtain ⟨c, hc, he⟩ := Wpoly_side_image hact
    rw [he]
    exact ((isCU_zz_pow _).mul (isGoodW_prod_tailW s).isCU).smul (by linarith)
  · rw [sum_sideBlock_not_active hact]
    by_cases h0 : s none = 0
    · rw [if_pos h0]
      exact (isGoodW_Wpoly_spider_of_not_active hact).isCU
    · rw [if_neg h0]
      exact isCU_zero

/-! ### The four-map block -/

/-- **The four-map block of the connector tree.**  For an odd number of connector edges
it is `Fblock r s c d` plus a good multiple of `(z+z⁻¹)^(r+s)`, and for an even number
of connector edges it is exactly the even-connector block `Gblock r s c d`, in both
cases with `r, s ≥ 1` and derived scalars `c, d = 2^e ≥ 1`. -/
theorem connBlockSum_both_active (h0 : 0 < t) {u : SpiderV m a → ℕ} {v : SpiderV n b → ℕ}
    (hu : ActiveSide u) (hv : ActiveSide v) :
    IsCU (Wpoly (connGraph t m a n b) (Sum.elim u (Sum.elim (fun _ => 1) v))
      + Wpoly (connGraph t m a n b) (Sum.elim u (Sum.elim (fun _ => 1) (transf v)))
      + (Wpoly (connGraph t m a n b) (Sum.elim (transf u) (Sum.elim (fun _ => 1) v))
        + Wpoly (connGraph t m a n b)
            (Sum.elim (transf u) (Sum.elim (fun _ => 1) (transf v))))) := by
  obtain ⟨cu, hcu, heu⟩ := Wpoly_side_image hu
  obtain ⟨cv, hcv, hev⟩ := Wpoly_side_image hv
  have hpu : 2 ≤ pNum u := hu.two_odd
  have hpv : 2 ≤ pNum v := hv.two_odd
  have hgood : IsGoodW ((∏ i : Fin m, tailW u i) * ∏ j : Fin n, tailW v j) :=
    IsGoodW.mul (isGoodW_prod_tailW u) (isGoodW_prod_tailW v)
  have h2 : Wpoly (connGraph t m a n b) (Sum.elim u (Sum.elim (fun _ => 1) (transf v)))
      = (Pw ((1 - ((pNum u : ℤ) + ((t % 2 : ℕ) : ℤ))).natAbs) * ∏ i : Fin m, tailW u i)
        * (cv • (zz ^ (pNum v - 1) * ∏ j : Fin n, tailW v j)) := by
    rw [Wpoly_conn_hubR_zero u (fun _ => 1) (transf v) (transf_none v),
      Wpoly_armSpider_active_ones u hu.hub hu.run, hev]
  have h3 : Wpoly (connGraph t m a n b) (Sum.elim (transf u) (Sum.elim (fun _ => 1) v))
      = (cu • (zz ^ (pNum u - 1) * ∏ i : Fin m, tailW u i))
        * (Pw ((1 - ((pNum v : ℤ) + ((t % 2 : ℕ) : ℤ))).natAbs) * ∏ j : Fin n, tailW v j) := by
    rw [Wpoly_conn_hubL_zero (transf u) (fun _ => 1) v (transf_none u),
      Wpoly_armSpiderR_active_ones v hv.hub hv.run, heu]
  have h4 : Wpoly (connGraph t m a n b)
        (Sum.elim (transf u) (Sum.elim (fun _ => 1) (transf v)))
      = (cu • (zz ^ (pNum u - 1) * ∏ i : Fin m, tailW u i))
        * (Pw (t % 2) * (cv • (zz ^ (pNum v - 1) * ∏ j : Fin n, tailW v j))) := by
    rw [Wpoly_conn_bothHubs_zero (transf u) (fun _ => 1) (transf v) (transf_none u)
      (transf_none v), Wpoly_pathGraph_ones t h0, heu, hev]
  rcases Nat.mod_two_eq_zero_or_one t with ht | ht
  · have h1 : Wpoly (connGraph t m a n b) (Sum.elim u (Sum.elim (fun _ => 1) v))
        = Pw ((((pNum u - 1 : ℕ) : ℤ) - ((pNum v - 1 : ℕ) : ℤ)).natAbs)
          * ((∏ i : Fin m, tailW u i) * ∏ j : Fin n, tailW v j) := by
      rw [Wpoly_conn_joint_even_t u v hu.hub hv.hub hu.run hv.run ht,
        show ((pNum u : ℤ) - (pNum v : ℤ)).natAbs
          = (((pNum u - 1 : ℕ) : ℤ) - ((pNum v - 1 : ℕ) : ℤ)).natAbs from by omega]
    rw [ht] at h2 h3 h4
    rw [show (1 - ((pNum u : ℤ) + (((0 : ℕ)) : ℤ))).natAbs = pNum u - 1 from by omega] at h2
    rw [show (1 - ((pNum v : ℤ) + (((0 : ℕ)) : ℤ))).natAbs = pNum v - 1 from by omega] at h3
    have key : Wpoly (connGraph t m a n b) (Sum.elim u (Sum.elim (fun _ => 1) v))
        + Wpoly (connGraph t m a n b) (Sum.elim u (Sum.elim (fun _ => 1) (transf v)))
        + (Wpoly (connGraph t m a n b) (Sum.elim (transf u) (Sum.elim (fun _ => 1) v))
          + Wpoly (connGraph t m a n b)
              (Sum.elim (transf u) (Sum.elim (fun _ => 1) (transf v))))
        = (Fblock (pNum u - 1) (pNum v - 1) cu cv
            + (cu * cv) • zz ^ ((pNum u - 1) + (pNum v - 1)))
          * ((∏ i : Fin m, tailW u i) * ∏ j : Fin n, tailW v j) := by
      rw [h1, h2, h3, h4, Fblock_expand_abs, Pw_zero]
      simp only [smul_eq_C_mul, map_mul]
      ring
    rw [key]
    exact ((Fblock_isCU (by omega) (by omega) hcu hcv).add
      ((isCU_zz_pow _).smul (by nlinarith))).mul hgood.isCU
  · have h1 : Wpoly (connGraph t m a n b) (Sum.elim u (Sum.elim (fun _ => 1) v))
        = Pw ((pNum u - 1) + (pNum v - 1) + 1)
          * ((∏ i : Fin m, tailW u i) * ∏ j : Fin n, tailW v j) := by
      rw [Wpoly_conn_joint_odd_t u v hu.hub hv.hub hu.run hv.run ht,
        show ((1 : ℤ) - (pNum u : ℤ) - (pNum v : ℤ)).natAbs
          = (pNum u - 1) + (pNum v - 1) + 1 from by omega]
    rw [ht] at h2 h3 h4
    rw [show (1 - ((pNum u : ℤ) + (((1 : ℕ)) : ℤ))).natAbs = (pNum u - 1) + 1 from by
      omega] at h2
    rw [show (1 - ((pNum v : ℤ) + (((1 : ℕ)) : ℤ))).natAbs = (pNum v - 1) + 1 from by
      omega] at h3
    have key : Wpoly (connGraph t m a n b) (Sum.elim u (Sum.elim (fun _ => 1) v))
        + Wpoly (connGraph t m a n b) (Sum.elim u (Sum.elim (fun _ => 1) (transf v)))
        + (Wpoly (connGraph t m a n b) (Sum.elim (transf u) (Sum.elim (fun _ => 1) v))
          + Wpoly (connGraph t m a n b)
              (Sum.elim (transf u) (Sum.elim (fun _ => 1) (transf v))))
        = Gblock (pNum u - 1) (pNum v - 1) cu cv
          * ((∏ i : Fin m, tailW u i) * ∏ j : Fin n, tailW v j) := by
      rw [h1, h2, h3, h4, Gblock, Pw_one]
      simp only [smul_eq_C_mul, map_mul]
      ring
    rw [key]
    exact (Gblock_isCU (by omega) (by omega) hcu hcv).mul hgood.isCU

/-! ### The two-map blocks -/

/-- **The two-map block at the left hub of the connector tree.** -/
theorem connBlockSum_left_active (h0 : 0 < t) {u : SpiderV m a → ℕ} {v : SpiderV n b → ℕ}
    (hu : ActiveSide u) (hv : ¬ ActiveSide v) :
    IsCU (Wpoly (connGraph t m a n b) (Sum.elim u (Sum.elim (fun _ => 1) v))
      + Wpoly (connGraph t m a n b) (Sum.elim (transf u) (Sum.elim (fun _ => 1) v))) := by
  classical
  obtain ⟨cu, hcu, heu⟩ := Wpoly_side_image hu
  have hpu : 2 ≤ pNum u := hu.two_odd
  by_cases hv2 : 2 ≤ v none
  · rw [Wpoly_conn_eq_zero_hubR_two u (fun _ => 1) v h0 hv2 (le_refl 1),
      Wpoly_conn_eq_zero_hubR_two (transf u) (fun _ => 1) v h0 hv2 (le_refl 1), add_zero]
    exact isCU_zero
  by_cases hv0 : v none = 0
  · have hgood : IsGoodW ((∏ i : Fin m, tailW u i) * Wpoly (spider n b) v) :=
      IsGoodW.mul (isGoodW_prod_tailW u) (isGoodW_Wpoly_spider_of_not_active hv)
    have h1 : Wpoly (connGraph t m a n b) (Sum.elim u (Sum.elim (fun _ => 1) v))
        = (Pw ((1 - ((pNum u : ℤ) + ((t % 2 : ℕ) : ℤ))).natAbs) * ∏ i : Fin m, tailW u i)
          * Wpoly (spider n b) v := by
      rw [Wpoly_conn_hubR_zero u (fun _ => 1) v hv0,
        Wpoly_armSpider_active_ones u hu.hub hu.run]
    have h2 : Wpoly (connGraph t m a n b) (Sum.elim (transf u) (Sum.elim (fun _ => 1) v))
        = (cu • (zz ^ (pNum u - 1) * ∏ i : Fin m, tailW u i))
          * (Pw (t % 2) * Wpoly (spider n b) v) := by
      rw [Wpoly_conn_bothHubs_zero (transf u) (fun _ => 1) v (transf_none u) hv0,
        Wpoly_pathGraph_ones t h0, heu]
    rcases Nat.mod_two_eq_zero_or_one t with ht | ht
    · rw [ht] at h1 h2
      rw [show (1 - ((pNum u : ℤ) + (((0 : ℕ)) : ℤ))).natAbs = pNum u - 1 from by omega] at h1
      have key : Wpoly (connGraph t m a n b) (Sum.elim u (Sum.elim (fun _ => 1) v))
          + Wpoly (connGraph t m a n b) (Sum.elim (transf u) (Sum.elim (fun _ => 1) v))
          = Ablock (pNum u - 1) (2 * cu)
            * ((∏ i : Fin m, tailW u i) * Wpoly (spider n b) v) := by
        rw [h1, h2, Pw_zero, Ablock]
        simp only [smul_eq_C_mul, map_mul]
        norm_num [C_two_eq]
        ring
      rw [key]
      exact (Ablock_isCU (by omega) (by linarith)).mul hgood.isCU
    · rw [ht] at h1 h2
      rw [show (1 - ((pNum u : ℤ) + (((1 : ℕ)) : ℤ))).natAbs = (pNum u - 1) + 1 from by
        omega] at h1
      have key : Wpoly (connGraph t m a n b) (Sum.elim u (Sum.elim (fun _ => 1) v))
          + Wpoly (connGraph t m a n b) (Sum.elim (transf u) (Sum.elim (fun _ => 1) v))
          = Ablock ((pNum u - 1) + 1) cu
            * ((∏ i : Fin m, tailW u i) * Wpoly (spider n b) v) := by
        rw [h1, h2, Pw_one, Ablock]
        simp only [smul_eq_C_mul]
        ring
      rw [key]
      exact (Ablock_isCU (by omega) hcu).mul hgood.isCU
  have hv1 : v none = 1 := by omega
  by_cases hrv : AdmRun v
  · have hq : pNum v ≤ 1 := by
      by_contra hc
      exact hv ⟨hv1, hrv, by omega⟩
    have hgood : IsGoodW ((∏ i : Fin m, tailW u i) * ∏ j : Fin n, tailW v j) :=
      IsGoodW.mul (isGoodW_prod_tailW u) (isGoodW_prod_tailW v)
    have h2 : Wpoly (connGraph t m a n b) (Sum.elim (transf u) (Sum.elim (fun _ => 1) v))
        = (cu • (zz ^ (pNum u - 1) * ∏ i : Fin m, tailW u i))
          * (Pw ((1 - ((pNum v : ℤ) + ((t % 2 : ℕ) : ℤ))).natAbs) * ∏ j : Fin n, tailW v j) := by
      rw [Wpoly_conn_hubL_zero (transf u) (fun _ => 1) v (transf_none u),
        Wpoly_armSpiderR_active_ones v hv1 hrv, heu]
    rcases Nat.mod_two_eq_zero_or_one t with ht | ht
    · have h1 : Wpoly (connGraph t m a n b) (Sum.elim u (Sum.elim (fun _ => 1) v))
          = Pw (((pNum u : ℤ) - (pNum v : ℤ)).natAbs)
            * ((∏ i : Fin m, tailW u i) * ∏ j : Fin n, tailW v j) :=
        Wpoly_conn_joint_even_t u v hu.hub hv1 hu.run hrv ht
      rw [ht] at h2
      rcases (show pNum v = 0 ∨ pNum v = 1 from by omega) with hqv | hqv
      · rw [hqv] at h1 h2
        rw [show ((pNum u : ℤ) - ((0 : ℕ) : ℤ)).natAbs = (pNum u - 1) + 1 from by omega] at h1
        rw [show (1 - (((0 : ℕ) : ℤ) + (((0 : ℕ)) : ℤ))).natAbs = 1 from by omega] at h2
        have key : Wpoly (connGraph t m a n b) (Sum.elim u (Sum.elim (fun _ => 1) v))
            + Wpoly (connGraph t m a n b) (Sum.elim (transf u) (Sum.elim (fun _ => 1) v))
            = Ablock ((pNum u - 1) + 1) cu
              * ((∏ i : Fin m, tailW u i) * ∏ j : Fin n, tailW v j) := by
          rw [h1, h2, Pw_one, Ablock]
          simp only [smul_eq_C_mul]
          ring
        rw [key]
        exact (Ablock_isCU (by omega) hcu).mul hgood.isCU
      · rw [hqv] at h1 h2
        rw [show ((pNum u : ℤ) - ((1 : ℕ) : ℤ)).natAbs = pNum u - 1 from by omega] at h1
        rw [show (1 - (((1 : ℕ) : ℤ) + (((0 : ℕ)) : ℤ))).natAbs = 0 from by omega] at h2
        have key : Wpoly (connGraph t m a n b) (Sum.elim u (Sum.elim (fun _ => 1) v))
            + Wpoly (connGraph t m a n b) (Sum.elim (transf u) (Sum.elim (fun _ => 1) v))
            = Ablock (pNum u - 1) (2 * cu)
              * ((∏ i : Fin m, tailW u i) * ∏ j : Fin n, tailW v j) := by
          rw [h1, h2, Pw_zero, Ablock]
          simp only [smul_eq_C_mul, map_mul]
          norm_num [C_two_eq]
          ring
        rw [key]
        exact (Ablock_isCU (by omega) (by linarith)).mul hgood.isCU
    · have h1 : Wpoly (connGraph t m a n b) (Sum.elim u (Sum.elim (fun _ => 1) v))
          = Pw (((1 : ℤ) - (pNum u : ℤ) - (pNum v : ℤ)).natAbs)
            * ((∏ i : Fin m, tailW u i) * ∏ j : Fin n, tailW v j) :=
        Wpoly_conn_joint_odd_t u v hu.hub hv1 hu.run hrv ht
      rw [ht] at h2
      rcases (show pNum v = 0 ∨ pNum v = 1 from by omega) with hqv | hqv
      · rw [hqv] at h1 h2
        rw [show ((1 : ℤ) - (pNum u : ℤ) - ((0 : ℕ) : ℤ)).natAbs = pNum u - 1 from by
          omega] at h1
        rw [show (1 - (((0 : ℕ) : ℤ) + (((1 : ℕ)) : ℤ))).natAbs = 0 from by omega] at h2
        have key : Wpoly (connGraph t m a n b) (Sum.elim u (Sum.elim (fun _ => 1) v))
            + Wpoly (connGraph t m a n b) (Sum.elim (transf u) (Sum.elim (fun _ => 1) v))
            = Ablock (pNum u - 1) (2 * cu)
              * ((∏ i : Fin m, tailW u i) * ∏ j : Fin n, tailW v j) := by
          rw [h1, h2, Pw_zero, Ablock]
          simp only [smul_eq_C_mul, map_mul]
          norm_num [C_two_eq]
          ring
        rw [key]
        exact (Ablock_isCU (by omega) (by linarith)).mul hgood.isCU
      · rw [hqv] at h1 h2
        rw [show ((1 : ℤ) - (pNum u : ℤ) - ((1 : ℕ) : ℤ)).natAbs = (pNum u - 1) + 1 from by
          omega] at h1
        rw [show (1 - (((1 : ℕ) : ℤ) + (((1 : ℕ)) : ℤ))).natAbs = 1 from by omega] at h2
        have key : Wpoly (connGraph t m a n b) (Sum.elim u (Sum.elim (fun _ => 1) v))
            + Wpoly (connGraph t m a n b) (Sum.elim (transf u) (Sum.elim (fun _ => 1) v))
            = Ablock ((pNum u - 1) + 1) cu
              * ((∏ i : Fin m, tailW u i) * ∏ j : Fin n, tailW v j) := by
          rw [h1, h2, Pw_one, Ablock]
          simp only [smul_eq_C_mul]
          ring
        rw [key]
        exact (Ablock_isCU (by omega) hcu).mul hgood.isCU
  · rw [Wpoly_conn_eq_zero_of_not_admRun_right u (fun _ => 1) v (by omega) hrv,
      Wpoly_conn_eq_zero_of_not_admRun_right (transf u) (fun _ => 1) v (by omega) hrv, add_zero]
    exact isCU_zero

/-- **The two-map block at the right hub of the connector tree.** -/
theorem connBlockSum_right_active (h0 : 0 < t) {u : SpiderV m a → ℕ} {v : SpiderV n b → ℕ}
    (hu : ¬ ActiveSide u) (hv : ActiveSide v) :
    IsCU (Wpoly (connGraph t m a n b) (Sum.elim u (Sum.elim (fun _ => 1) v))
      + Wpoly (connGraph t m a n b) (Sum.elim u (Sum.elim (fun _ => 1) (transf v)))) := by
  classical
  obtain ⟨cv, hcv, hev⟩ := Wpoly_side_image hv
  have hpv : 2 ≤ pNum v := hv.two_odd
  by_cases hu2 : 2 ≤ u none
  · rw [Wpoly_conn_eq_zero_hubL_two u (fun _ => 1) v h0 hu2 (le_refl 1),
      Wpoly_conn_eq_zero_hubL_two u (fun _ => 1) (transf v) h0 hu2 (le_refl 1), add_zero]
    exact isCU_zero
  by_cases hu0 : u none = 0
  · have hgood : IsGoodW (Wpoly (spider m a) u * ∏ j : Fin n, tailW v j) :=
      IsGoodW.mul (isGoodW_Wpoly_spider_of_not_active hu) (isGoodW_prod_tailW v)
    have h1 : Wpoly (connGraph t m a n b) (Sum.elim u (Sum.elim (fun _ => 1) v))
        = Wpoly (spider m a) u
          * (Pw ((1 - ((pNum v : ℤ) + ((t % 2 : ℕ) : ℤ))).natAbs) * ∏ j : Fin n, tailW v j) := by
      rw [Wpoly_conn_hubL_zero u (fun _ => 1) v hu0,
        Wpoly_armSpiderR_active_ones v hv.hub hv.run]
    have h2 : Wpoly (connGraph t m a n b) (Sum.elim u (Sum.elim (fun _ => 1) (transf v)))
        = Wpoly (spider m a) u
          * (Pw (t % 2) * (cv • (zz ^ (pNum v - 1) * ∏ j : Fin n, tailW v j))) := by
      rw [Wpoly_conn_bothHubs_zero u (fun _ => 1) (transf v) hu0 (transf_none v),
        Wpoly_pathGraph_ones t h0, hev]
    rcases Nat.mod_two_eq_zero_or_one t with ht | ht
    · rw [ht] at h1 h2
      rw [show (1 - ((pNum v : ℤ) + (((0 : ℕ)) : ℤ))).natAbs = pNum v - 1 from by omega] at h1
      have key : Wpoly (connGraph t m a n b) (Sum.elim u (Sum.elim (fun _ => 1) v))
          + Wpoly (connGraph t m a n b) (Sum.elim u (Sum.elim (fun _ => 1) (transf v)))
          = Ablock (pNum v - 1) (2 * cv)
            * (Wpoly (spider m a) u * ∏ j : Fin n, tailW v j) := by
        rw [h1, h2, Pw_zero, Ablock]
        simp only [smul_eq_C_mul, map_mul]
        norm_num [C_two_eq]
        ring
      rw [key]
      exact (Ablock_isCU (by omega) (by linarith)).mul hgood.isCU
    · rw [ht] at h1 h2
      rw [show (1 - ((pNum v : ℤ) + (((1 : ℕ)) : ℤ))).natAbs = (pNum v - 1) + 1 from by
        omega] at h1
      have key : Wpoly (connGraph t m a n b) (Sum.elim u (Sum.elim (fun _ => 1) v))
          + Wpoly (connGraph t m a n b) (Sum.elim u (Sum.elim (fun _ => 1) (transf v)))
          = Ablock ((pNum v - 1) + 1) cv
            * (Wpoly (spider m a) u * ∏ j : Fin n, tailW v j) := by
        rw [h1, h2, Pw_one, Ablock]
        simp only [smul_eq_C_mul]
        ring
      rw [key]
      exact (Ablock_isCU (by omega) hcv).mul hgood.isCU
  have hu1 : u none = 1 := by omega
  by_cases hru : AdmRun u
  · have hp : pNum u ≤ 1 := by
      by_contra hc
      exact hu ⟨hu1, hru, by omega⟩
    have hgood : IsGoodW ((∏ i : Fin m, tailW u i) * ∏ j : Fin n, tailW v j) :=
      IsGoodW.mul (isGoodW_prod_tailW u) (isGoodW_prod_tailW v)
    have h2 : Wpoly (connGraph t m a n b) (Sum.elim u (Sum.elim (fun _ => 1) (transf v)))
        = (Pw ((1 - ((pNum u : ℤ) + ((t % 2 : ℕ) : ℤ))).natAbs) * ∏ i : Fin m, tailW u i)
          * (cv • (zz ^ (pNum v - 1) * ∏ j : Fin n, tailW v j)) := by
      rw [Wpoly_conn_hubR_zero u (fun _ => 1) (transf v) (transf_none v),
        Wpoly_armSpider_active_ones u hu1 hru, hev]
    rcases Nat.mod_two_eq_zero_or_one t with ht | ht
    · have h1 : Wpoly (connGraph t m a n b) (Sum.elim u (Sum.elim (fun _ => 1) v))
          = Pw (((pNum u : ℤ) - (pNum v : ℤ)).natAbs)
            * ((∏ i : Fin m, tailW u i) * ∏ j : Fin n, tailW v j) :=
        Wpoly_conn_joint_even_t u v hu1 hv.hub hru hv.run ht
      rw [ht] at h2
      rcases (show pNum u = 0 ∨ pNum u = 1 from by omega) with hpu | hpu
      · rw [hpu] at h1 h2
        rw [show (((0 : ℕ) : ℤ) - (pNum v : ℤ)).natAbs = (pNum v - 1) + 1 from by omega] at h1
        rw [show (1 - (((0 : ℕ) : ℤ) + (((0 : ℕ)) : ℤ))).natAbs = 1 from by omega] at h2
        have key : Wpoly (connGraph t m a n b) (Sum.elim u (Sum.elim (fun _ => 1) v))
            + Wpoly (connGraph t m a n b) (Sum.elim u (Sum.elim (fun _ => 1) (transf v)))
            = Ablock ((pNum v - 1) + 1) cv
              * ((∏ i : Fin m, tailW u i) * ∏ j : Fin n, tailW v j) := by
          rw [h1, h2, Pw_one, Ablock]
          simp only [smul_eq_C_mul]
          ring
        rw [key]
        exact (Ablock_isCU (by omega) hcv).mul hgood.isCU
      · rw [hpu] at h1 h2
        rw [show (((1 : ℕ) : ℤ) - (pNum v : ℤ)).natAbs = pNum v - 1 from by omega] at h1
        rw [show (1 - (((1 : ℕ) : ℤ) + (((0 : ℕ)) : ℤ))).natAbs = 0 from by omega] at h2
        have key : Wpoly (connGraph t m a n b) (Sum.elim u (Sum.elim (fun _ => 1) v))
            + Wpoly (connGraph t m a n b) (Sum.elim u (Sum.elim (fun _ => 1) (transf v)))
            = Ablock (pNum v - 1) (2 * cv)
              * ((∏ i : Fin m, tailW u i) * ∏ j : Fin n, tailW v j) := by
          rw [h1, h2, Pw_zero, Ablock]
          simp only [smul_eq_C_mul, map_mul]
          norm_num [C_two_eq]
          ring
        rw [key]
        exact (Ablock_isCU (by omega) (by linarith)).mul hgood.isCU
    · have h1 : Wpoly (connGraph t m a n b) (Sum.elim u (Sum.elim (fun _ => 1) v))
          = Pw (((1 : ℤ) - (pNum u : ℤ) - (pNum v : ℤ)).natAbs)
            * ((∏ i : Fin m, tailW u i) * ∏ j : Fin n, tailW v j) :=
        Wpoly_conn_joint_odd_t u v hu1 hv.hub hru hv.run ht
      rw [ht] at h2
      rcases (show pNum u = 0 ∨ pNum u = 1 from by omega) with hpu | hpu
      · rw [hpu] at h1 h2
        rw [show ((1 : ℤ) - ((0 : ℕ) : ℤ) - (pNum v : ℤ)).natAbs = pNum v - 1 from by
          omega] at h1
        rw [show (1 - (((0 : ℕ) : ℤ) + (((1 : ℕ)) : ℤ))).natAbs = 0 from by omega] at h2
        have key : Wpoly (connGraph t m a n b) (Sum.elim u (Sum.elim (fun _ => 1) v))
            + Wpoly (connGraph t m a n b) (Sum.elim u (Sum.elim (fun _ => 1) (transf v)))
            = Ablock (pNum v - 1) (2 * cv)
              * ((∏ i : Fin m, tailW u i) * ∏ j : Fin n, tailW v j) := by
          rw [h1, h2, Pw_zero, Ablock]
          simp only [smul_eq_C_mul, map_mul]
          norm_num [C_two_eq]
          ring
        rw [key]
        exact (Ablock_isCU (by omega) (by linarith)).mul hgood.isCU
      · rw [hpu] at h1 h2
        rw [show ((1 : ℤ) - ((1 : ℕ) : ℤ) - (pNum v : ℤ)).natAbs = (pNum v - 1) + 1 from by
          omega] at h1
        rw [show (1 - (((1 : ℕ) : ℤ) + (((1 : ℕ)) : ℤ))).natAbs = 1 from by omega] at h2
        have key : Wpoly (connGraph t m a n b) (Sum.elim u (Sum.elim (fun _ => 1) v))
            + Wpoly (connGraph t m a n b) (Sum.elim u (Sum.elim (fun _ => 1) (transf v)))
            = Ablock ((pNum v - 1) + 1) cv
              * ((∏ i : Fin m, tailW u i) * ∏ j : Fin n, tailW v j) := by
          rw [h1, h2, Pw_one, Ablock]
          simp only [smul_eq_C_mul]
          ring
        rw [key]
        exact (Ablock_isCU (by omega) hcv).mul hgood.isCU
  · rw [Wpoly_conn_eq_zero_of_not_admRun_left u (fun _ => 1) v (by omega) hru,
      Wpoly_conn_eq_zero_of_not_admRun_left u (fun _ => 1) (transf v) (by omega) hru, add_zero]
    exact isCU_zero

/-! ### The singletons -/

/-- **A singleton block of the connector tree.** -/
theorem connBlockSum_neither_active (h0 : 0 < t) {u : SpiderV m a → ℕ} {v : SpiderV n b → ℕ}
    (hu : ¬ ActiveSide u) (hv : ¬ ActiveSide v) :
    IsGoodW (Wpoly (connGraph t m a n b) (Sum.elim u (Sum.elim (fun _ => 1) v))) := by
  classical
  by_cases hu2 : 2 ≤ u none
  · exact Or.inl (Wpoly_conn_eq_zero_hubL_two u (fun _ => 1) v h0 hu2 (le_refl 1))
  by_cases hv2 : 2 ≤ v none
  · exact Or.inl (Wpoly_conn_eq_zero_hubR_two u (fun _ => 1) v h0 hv2 (le_refl 1))
  by_cases hu0 : u none = 0
  · by_cases hv0 : v none = 0
    · rw [Wpoly_conn_bothHubs_zero u (fun _ => 1) v hu0 hv0, Wpoly_pathGraph_ones t h0]
      exact IsGoodW.mul (isGoodW_Wpoly_spider_of_not_active hu)
        (IsGoodW.mul (isGoodW_Pw_le_one (by omega))
          (isGoodW_Wpoly_spider_of_not_active hv))
    · have hv1 : v none = 1 := by omega
      by_cases hrv : AdmRun v
      · have hq : pNum v ≤ 1 := by
          by_contra hc
          exact hv ⟨hv1, hrv, by omega⟩
        rw [Wpoly_conn_hubL_zero u (fun _ => 1) v hu0,
          Wpoly_armSpiderR_active_ones v hv1 hrv]
        exact IsGoodW.mul (isGoodW_Wpoly_spider_of_not_active hu)
          (IsGoodW.mul (isGoodW_Pw_le_one (by omega)) (isGoodW_prod_tailW v))
      · exact Or.inl (Wpoly_conn_eq_zero_of_not_admRun_right u (fun _ => 1) v (by omega) hrv)
  have hu1 : u none = 1 := by omega
  by_cases hru : AdmRun u
  · have hp : pNum u ≤ 1 := by
      by_contra hc
      exact hu ⟨hu1, hru, by omega⟩
    by_cases hv0 : v none = 0
    · rw [Wpoly_conn_hubR_zero u (fun _ => 1) v hv0, Wpoly_armSpider_active_ones u hu1 hru]
      exact IsGoodW.mul (IsGoodW.mul (isGoodW_Pw_le_one (by omega)) (isGoodW_prod_tailW u)) (isGoodW_Wpoly_spider_of_not_active hv)
    · have hv1 : v none = 1 := by omega
      by_cases hrv : AdmRun v
      · have hq : pNum v ≤ 1 := by
          by_contra hc
          exact hv ⟨hv1, hrv, by omega⟩
        rcases Nat.mod_two_eq_zero_or_one t with ht | ht
        · rw [Wpoly_conn_joint_even_t u v hu1 hv1 hru hrv ht]
          exact IsGoodW.mul (isGoodW_Pw_le_one (by omega))
            (IsGoodW.mul (isGoodW_prod_tailW u) (isGoodW_prod_tailW v))
        · rw [Wpoly_conn_joint_odd_t u v hu1 hv1 hru hrv ht]
          exact IsGoodW.mul (isGoodW_Pw_le_one (by omega))
            (IsGoodW.mul (isGoodW_prod_tailW u) (isGoodW_prod_tailW v))
      · exact Or.inl (Wpoly_conn_eq_zero_of_not_admRun_right u (fun _ => 1) v (by omega) hrv)
  · exact Or.inl (Wpoly_conn_eq_zero_of_not_admRun_left u (fun _ => 1) v (by omega) hru)

/-! ### The unit-positive connector stratum -/

/-- **The block sum of the unit-positive connector stratum is centrally unimodal.** -/
theorem connBlockSum_ones_isCU (h0 : 0 < t) (u : SpiderV m a → ℕ) (v : SpiderV n b → ℕ) :
    IsCU (∑ u' ∈ sideBlock u, ∑ v' ∈ sideBlock v,
      Wpoly (connGraph t m a n b) (Sum.elim u' (Sum.elim (fun _ => 1) v'))) := by
  classical
  by_cases hu : ActiveSide u
  · by_cases hv : ActiveSide v
    · rw [sum_sideBlock_active hu, sum_sideBlock_active hv, sum_sideBlock_active hv]
      exact connBlockSum_both_active h0 hu hv
    · rw [sum_sideBlock_active hu, sum_sideBlock_not_active hv, sum_sideBlock_not_active hv]
      exact connBlockSum_left_active h0 hu hv
  · by_cases hv : ActiveSide v
    · rw [sum_sideBlock_not_active hu, sum_sideBlock_active hv]
      exact connBlockSum_right_active h0 hu hv
    · rw [sum_sideBlock_not_active hu, sum_sideBlock_not_active hv]
      exact (connBlockSum_neither_active h0 hu hv).isCU

/-! ### The full block sum -/

/-- **Every block of the connector partition has a centrally unimodal weight.**  The
three connector states — a broken connector, a connector with a multiplicity at least
two, and a unit-positive connector — are treated exhaustively. -/
theorem connBlockSum_isCU (h0 : 0 < t) (u : SpiderV m a → ℕ) (w : Fin t → ℕ)
    (v : SpiderV n b → ℕ) :
    IsCU (∑ u' ∈ sideBlock u, ∑ v' ∈ sideBlock v,
      Wpoly (connGraph t m a n b) (Sum.elim u' (Sum.elim w v'))) := by
  classical
  by_cases hbreak : ∃ j : Fin t, w j = 0
  · -- a broken connector: the tree splits into two arm spiders
    obtain ⟨j, hj⟩ := hbreak
    have hj' : w ⟨j.val, j.isLt⟩ = 0 := hj
    have key : (∑ u' ∈ sideBlock u, ∑ v' ∈ sideBlock v,
          Wpoly (connGraph t m a n b) (Sum.elim u' (Sum.elim w v')))
        = (∑ u' ∈ sideBlock u, Wpoly (armSpider m a (j.val + 1))
              (Sum.elim u' (fun k : Fin (j.val + 1) =>
                w ⟨k.val, by have := k.isLt; have := j.isLt; omega⟩)))
          * (∑ v' ∈ sideBlock v, Wpoly (armSpiderR n b (t - j.val - 1))
              (Sum.elim (fun k : Fin (t - j.val - 1) =>
                w ⟨k.val + j.val + 1, by have := k.isLt; have := j.isLt; omega⟩) v')) := by
      rw [Finset.sum_mul]
      refine Finset.sum_congr rfl fun u' _ => ?_
      rw [Finset.mul_sum]
      exact Finset.sum_congr rfl fun v' _ =>
        Wpoly_conn_interior_zero u' w v' j.val j.isLt hj'
    rw [key]
    exact (armSpider_sideBlock_isCU u _).mul (armSpiderR_sideBlock_isCU _ v)
  push_neg at hbreak
  have hpos : ∀ k : Fin t, 1 ≤ w k := fun k => Nat.one_le_iff_ne_zero.mpr (hbreak k)
  by_cases htwo : ∃ j : Fin t, 2 ≤ w j
  · -- a connector multiplicity at least two: only the maps with both hubs off survive
    obtain ⟨j, hj⟩ := htwo
    have hterm : ∀ (u' : SpiderV m a → ℕ) (v' : SpiderV n b → ℕ),
        Wpoly (connGraph t m a n b) (Sum.elim u' (Sum.elim w v'))
          = (if u' none = 0 then Wpoly (spider m a) u' else 0)
            * (Wpoly (pathGraph t) w
              * (if v' none = 0 then Wpoly (spider n b) v' else 0)) := by
      intro u' v'
      by_cases hu0 : u' none = 0
      · by_cases hv0 : v' none = 0
        · rw [if_pos hu0, if_pos hv0, Wpoly_conn_bothHubs_zero u' w v' hu0 hv0]
        · rw [if_pos hu0, if_neg hv0, mul_zero, mul_zero]
          exact Wpoly_conn_eq_zero_of_conn_two_right u' w v' j hj hpos (by omega)
      · rw [if_neg hu0, zero_mul]
        exact Wpoly_conn_eq_zero_of_conn_two_left u' w v' j hj hpos (by omega)
    have key : (∑ u' ∈ sideBlock u, ∑ v' ∈ sideBlock v,
          Wpoly (connGraph t m a n b) (Sum.elim u' (Sum.elim w v')))
        = (∑ u' ∈ sideBlock u, if u' none = 0 then Wpoly (spider m a) u' else 0)
          * (Wpoly (pathGraph t) w
            * (∑ v' ∈ sideBlock v, if v' none = 0 then Wpoly (spider n b) v' else 0)) := by
      rw [Finset.sum_mul]
      refine Finset.sum_congr rfl fun u' _ => ?_
      rw [Finset.mul_sum, Finset.mul_sum]
      exact Finset.sum_congr rfl fun v' _ => hterm u' v'
    rw [key]
    exact (sideBlock_hubZero_isCU u).mul
      ((Wpoly_pathGraph_isGood t w).isCU.mul (sideBlock_hubZero_isCU v))
  · push_neg at htwo
    have hw : w = fun _ => 1 := by
      funext k
      have := hpos k
      have := htwo k
      omega
    subst hw
    exact connBlockSum_ones_isCU h0 u v

end ClanAudit
