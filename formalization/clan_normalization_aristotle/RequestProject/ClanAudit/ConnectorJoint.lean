import RequestProject.ClanAudit.ConnectorSplit

/-!
# The joined stratum of the connector tree

When the connector carries multiplicity one at every vertex and both hubs are active,
the two hubs, the whole connector and all the active arm prefixes form one single
component of the clan graph.  Cutting the tree at the active prefixes (as in
`Wpoly_dgraph_cut`) reduces its weight to the all-ones weight of a smaller connector
tree, which was computed in `ConnectorImbalance.lean`.

The resulting imbalance is `|p - q|` for an odd number of connector edges
(`t` even) and `|1 - p - q|` for an even number of connector edges (`t` odd), where `p`
and `q` are the numbers of odd active prefixes at the two hubs.
-/

namespace ClanAudit

open Finset LaurentPolynomial SimpleGraph

variable {t m n : ℕ} {a : Fin m → ℕ} {b : Fin n → ℕ}

/-- Regrouping a five-fold sum. -/
def connShuffle (A B C D E : Type*) :
    (A ⊕ B) ⊕ (C ⊕ (D ⊕ E)) ≃ (A ⊕ (C ⊕ D)) ⊕ (B ⊕ E) where
  toFun z :=
    match z with
    | Sum.inl (Sum.inl u) => Sum.inl (Sum.inl u)
    | Sum.inl (Sum.inr u) => Sum.inr (Sum.inl u)
    | Sum.inr (Sum.inl u) => Sum.inl (Sum.inr (Sum.inl u))
    | Sum.inr (Sum.inr (Sum.inl u)) => Sum.inl (Sum.inr (Sum.inr u))
    | Sum.inr (Sum.inr (Sum.inr u)) => Sum.inr (Sum.inr u)
  invFun z :=
    match z with
    | Sum.inl (Sum.inl u) => Sum.inl (Sum.inl u)
    | Sum.inl (Sum.inr (Sum.inl u)) => Sum.inr (Sum.inl u)
    | Sum.inl (Sum.inr (Sum.inr u)) => Sum.inr (Sum.inr (Sum.inl u))
    | Sum.inr (Sum.inl u) => Sum.inl (Sum.inr u)
    | Sum.inr (Sum.inr u) => Sum.inr (Sum.inr (Sum.inr u))
  left_inv := by rintro ((u | u) | (u | (u | u))) <;> rfl
  right_inv := by rintro ((u | (u | u)) | (u | u)) <;> rfl

/-- The vertex splitting of the connector tree into the sub-tree spanned by the two
hubs, the connector and the active prefixes of all arms, and the family of arm tails. -/
def connCut (t m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) (pa : Fin m → ℕ)
    (qb : Fin n → ℕ) (hpa : ∀ i, pa i ≤ a i) (hqb : ∀ j, qb j ≤ b j) :
    ConnV t m a n b ≃
      ConnV t m pa n qb ⊕ ((Σ i : Fin m, Fin (a i - pa i)) ⊕ (Σ j : Fin n, Fin (b j - qb j))) :=
  (Equiv.sumCongr (spiderCut m a pa hpa)
      (Equiv.sumCongr (Equiv.refl (Fin t)) (spiderCut n b qb hqb))).trans
    (connShuffle _ _ _ _ _)

variable {pa : Fin m → ℕ} {qb : Fin n → ℕ} {hpa : ∀ i, pa i ≤ a i} {hqb : ∀ j, qb j ≤ b j}

theorem connCut_symm_core_l (p : SpiderV m pa) :
    (connCut t m a n b pa qb hpa hqb).symm (Sum.inl (Sum.inl p))
      = Sum.inl ((spiderCut m a pa hpa).symm (Sum.inl p)) := rfl

theorem connCut_symm_core_c (j : Fin t) :
    (connCut t m a n b pa qb hpa hqb).symm (Sum.inl (Sum.inr (Sum.inl j)))
      = Sum.inr (Sum.inl j) := rfl

theorem connCut_symm_core_r (q : SpiderV n qb) :
    (connCut t m a n b pa qb hpa hqb).symm (Sum.inl (Sum.inr (Sum.inr q)))
      = Sum.inr (Sum.inr ((spiderCut n b qb hqb).symm (Sum.inl q))) := rfl

theorem connCut_symm_tail_l (s : Σ i : Fin m, Fin (a i - pa i)) :
    (connCut t m a n b pa qb hpa hqb).symm (Sum.inr (Sum.inl s))
      = Sum.inl ((spiderCut m a pa hpa).symm (Sum.inr s)) := rfl

theorem connCut_symm_tail_r (s : Σ j : Fin n, Fin (b j - qb j)) :
    (connCut t m a n b pa qb hpa hqb).symm (Sum.inr (Sum.inr s))
      = Sum.inr (Sum.inr ((spiderCut n b qb hqb).symm (Sum.inr s))) := rfl

/-- **Cutting the connector tree at the active prefixes.** -/
theorem Wpoly_conn_cut (x : SpiderV m a → ℕ) (g : Fin t → ℕ) (y : SpiderV n b → ℕ)
    (hcx : ∀ (i : Fin m) (h : pa i < a i), x (some ⟨i, ⟨pa i, h⟩⟩) = 0)
    (hcy : ∀ (j : Fin n) (h : qb j < b j), y (some ⟨j, ⟨qb j, h⟩⟩) = 0) :
    Wpoly (connGraph t m a n b) (Sum.elim x (Sum.elim g y))
      = Wpoly (connGraph t m pa n qb)
          (fun z => Sum.elim x (Sum.elim g y)
            ((connCut t m a n b pa qb hpa hqb).symm (Sum.inl z)))
        * Wpoly ((ArmsG m (fun i => a i - pa i)).sum (ArmsG n (fun j => b j - qb j)))
          (fun z => Sum.elim x (Sum.elim g y)
            ((connCut t m a n b pa qb hpa hqb).symm (Sum.inr z))) := by
  refine Wpoly_split_equiv (connGraph t m a n b) (Sum.elim x (Sum.elim g y))
    (connGraph t m pa n qb) ((ArmsG m (fun i => a i - pa i)).sum (ArmsG n (fun j => b j - qb j)))
    (connCut t m a n b pa qb hpa hqb) ?_ ?_ ?_
  · rintro (p | j | q) (p' | j' | q')
    · rw [connCut_symm_core_l, connCut_symm_core_l]
      exact spiderCut_adj_core m a pa hpa p p'
    · rw [connCut_symm_core_l, connCut_symm_core_c, connGraph_adj_lc, connGraph_adj_lc,
        spiderCut_symm_inl_eq_none_iff]
    · rw [connCut_symm_core_l, connCut_symm_core_r, connGraph_adj_lr, connGraph_adj_lr,
        spiderCut_symm_inl_eq_none_iff, spiderCut_symm_inl_eq_none_iff]
    · rw [connCut_symm_core_c, connCut_symm_core_l]
      show ((spiderCut m a pa hpa).symm (Sum.inl p') = none ∧ j.val = 0) ↔ _
      rw [spiderCut_symm_inl_eq_none_iff]
      exact Iff.rfl
    · rw [connCut_symm_core_c, connCut_symm_core_c]
      exact Iff.rfl
    · rw [connCut_symm_core_c, connCut_symm_core_r, connGraph_adj_cr, connGraph_adj_cr,
        spiderCut_symm_inl_eq_none_iff]
    · rw [connCut_symm_core_r, connCut_symm_core_l]
      show ((spiderCut m a pa hpa).symm (Sum.inl p') = none
        ∧ (spiderCut n b qb hqb).symm (Sum.inl q) = none ∧ t = 0) ↔ _
      rw [spiderCut_symm_inl_eq_none_iff, spiderCut_symm_inl_eq_none_iff]
      exact Iff.rfl
    · rw [connCut_symm_core_r, connCut_symm_core_c]
      show ((spiderCut n b qb hqb).symm (Sum.inl q) = none ∧ j'.val + 1 = t) ↔ _
      rw [spiderCut_symm_inl_eq_none_iff]
      exact Iff.rfl
    · rw [connCut_symm_core_r, connCut_symm_core_r]
      exact spiderCut_adj_core n b qb hqb q q'
  · rintro (s | s) (s' | s')
    · rw [connCut_symm_tail_l, connCut_symm_tail_l, SimpleGraph.sum_adj]
      exact spiderCut_adj_tail m a pa hpa s s'
    · rw [connCut_symm_tail_l, connCut_symm_tail_r, connGraph_adj_lr, SimpleGraph.sum_adj]
      constructor
      · rintro ⟨h1, h2⟩
        exact absurd h1 (spiderCut_symm_inr_ne_none m a pa hpa s)
      · exact False.elim
    · rw [connCut_symm_tail_r, connCut_symm_tail_l, SimpleGraph.sum_adj]
      constructor
      · intro h
        exact absurd (connGraph_adj_lr.mp h.symm).1 (spiderCut_symm_inr_ne_none m a pa hpa s')
      · exact False.elim
    · rw [connCut_symm_tail_r, connCut_symm_tail_r, SimpleGraph.sum_adj]
      exact spiderCut_adj_tail n b qb hqb s s'
  · rintro (p | j | q) (s | s) hadj
    · rw [connCut_symm_core_l, connCut_symm_tail_l] at hadj ⊢
      exact spiderCut_cross m a pa hpa x hcx p s (connGraph_adj_ll.mp hadj)
    · rw [connCut_symm_core_l, connCut_symm_tail_r] at hadj
      exact absurd (connGraph_adj_lr.mp hadj).2.1 (spiderCut_symm_inr_ne_none n b qb hqb s)
    · rw [connCut_symm_core_c, connCut_symm_tail_l] at hadj
      exact absurd (connGraph_adj_lc.mp hadj.symm).1 (spiderCut_symm_inr_ne_none m a pa hpa s)
    · rw [connCut_symm_core_c, connCut_symm_tail_r] at hadj
      exact absurd (connGraph_adj_cr.mp hadj).1 (spiderCut_symm_inr_ne_none n b qb hqb s)
    · rw [connCut_symm_core_r, connCut_symm_tail_l] at hadj
      exact absurd (connGraph_adj_lr.mp hadj.symm).1 (spiderCut_symm_inr_ne_none m a pa hpa s)
    · rw [connCut_symm_core_r, connCut_symm_tail_r] at hadj ⊢
      exact spiderCut_cross n b qb hqb y hcy q s (connGraph_adj_rr.mp hadj)

/-- **Both hubs active and a unit-positive connector.**  The two hubs, the connector and
all the active prefixes form one component. -/
theorem Wpoly_conn_joint (x : SpiderV m a → ℕ) (y : SpiderV n b → ℕ)
    (hx : x none = 1) (hy : y none = 1) (hrx : AdmRun x) (hry : AdmRun y) (r : ℕ)
    (hr : Wpoly (connGraph t m (pref x) n (pref y)) (fun _ => 1) = Pw r) :
    Wpoly (connGraph t m a n b) (Sum.elim x (Sum.elim (fun _ => 1) y))
      = Pw r * ((∏ i : Fin m, tailW x i) * ∏ j : Fin n, tailW y j) := by
  classical
  have hcx : ∀ (i : Fin m) (h : pref x i < a i), x (some ⟨i, ⟨pref x i, h⟩⟩) = 0 := by
    intro i h
    have := hrx.stop i h
    rwa [armVal, dif_pos h] at this
  have hcy : ∀ (j : Fin n) (h : pref y j < b j), y (some ⟨j, ⟨pref y j, h⟩⟩) = 0 := by
    intro j h
    have := hry.stop j h
    rwa [armVal, dif_pos h] at this
  rw [Wpoly_conn_cut (pa := pref x) (qb := pref y) (hpa := pref_le x) (hqb := pref_le y)
    x (fun _ => 1) y hcx hcy]
  congr 1
  · have hones : (fun z => Sum.elim x (Sum.elim (fun _ : Fin t => 1) y)
        ((connCut t m a n b (pref x) (pref y) (pref_le x) (pref_le y)).symm (Sum.inl z)))
        = fun _ : ConnV t m (pref x) n (pref y) => 1 := by
      funext z
      match z with
      | Sum.inl none => exact hx
      | Sum.inl (some ⟨i, j⟩) =>
          rw [connCut_symm_core_l, spiderCut_symm_inl_some]
          have := pref_ones x i j.isLt
          rwa [armVal, dif_pos (show j.val < a i from by
            have := j.isLt; have := pref_le x i; omega)] at this
      | Sum.inr (Sum.inl j) => rfl
      | Sum.inr (Sum.inr none) => exact hy
      | Sum.inr (Sum.inr (some ⟨j, l⟩)) =>
          rw [connCut_symm_core_r, spiderCut_symm_inl_some]
          have := pref_ones y j l.isLt
          rwa [armVal, dif_pos (show l.val < b j from by
            have := l.isLt; have := pref_le y j; omega)] at this
    rw [hones, hr]
  · have helim : (fun z => Sum.elim x (Sum.elim (fun _ : Fin t => 1) y)
        ((connCut t m a n b (pref x) (pref y) (pref_le x) (pref_le y)).symm (Sum.inr z)))
        = Sum.elim
            (fun s : Σ i : Fin m, Fin (a i - pref x i) =>
              x ((spiderCut m a (pref x) (pref_le x)).symm (Sum.inr s)))
            (fun s : Σ j : Fin n, Fin (b j - pref y j) =>
              y ((spiderCut n b (pref y) (pref_le y)).symm (Sum.inr s))) := by
      funext z
      rcases z with s | s
      · rfl
      · rfl
    rw [helim, Wpoly_sum_elim, Wpoly_ArmsG, Wpoly_ArmsG]
    congr 1
    · refine Finset.prod_congr rfl fun i _ => ?_
      refine congrArg _ (funext fun j => ?_)
      rw [spiderCut_symm_inr, armVal,
        dif_pos (show pref x i + j.val < a i from by have := j.isLt; omega)]
    · refine Finset.prod_congr rfl fun j _ => ?_
      refine congrArg _ (funext fun l => ?_)
      rw [spiderCut_symm_inr, armVal,
        dif_pos (show pref y j + l.val < b j from by have := l.isLt; omega)]

/-- **Odd connector (`t` even).**  The joined component has imbalance `|p - q|`. -/
theorem Wpoly_conn_joint_even_t (x : SpiderV m a → ℕ) (y : SpiderV n b → ℕ)
    (hx : x none = 1) (hy : y none = 1) (hrx : AdmRun x) (hry : AdmRun y) (ht : t % 2 = 0) :
    Wpoly (connGraph t m a n b) (Sum.elim x (Sum.elim (fun _ => 1) y))
      = Pw (((pNum x : ℤ) - pNum y).natAbs)
        * ((∏ i : Fin m, tailW x i) * ∏ j : Fin n, tailW y j) :=
  Wpoly_conn_joint x y hx hy hrx hry _
    (Wpoly_connGraph_ones_ell_odd t m (pref x) n (pref y) ht (pNum x) (pNum y) _ rfl rfl rfl)

/-- **Even connector (`t` odd).**  The joined component has imbalance `|1 - p - q|`. -/
theorem Wpoly_conn_joint_odd_t (x : SpiderV m a → ℕ) (y : SpiderV n b → ℕ)
    (hx : x none = 1) (hy : y none = 1) (hrx : AdmRun x) (hry : AdmRun y) (ht : t % 2 = 1) :
    Wpoly (connGraph t m a n b) (Sum.elim x (Sum.elim (fun _ => 1) y))
      = Pw (((1 : ℤ) - pNum x - pNum y).natAbs)
        * ((∏ i : Fin m, tailW x i) * ∏ j : Fin n, tailW y j) :=
  Wpoly_conn_joint x y hx hy hrx hry _
    (Wpoly_connGraph_ones_ell_even t m (pref x) n (pref y) ht (pNum x) (pNum y) _ rfl rfl rfl)

end ClanAudit
