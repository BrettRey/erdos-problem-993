import RequestProject.ClanAudit.ImageWeight

/-!
# The normalized weight of a clan map of the adjacent two-hub tree

Two exact weight formulas, covering all clan maps whose two arm families have their
initial positive runs made of ones:

* `Wpoly_dgraph_split` — if at least one hub is switched off, the two sides are
  independent and the weight is the product of the two side weights;
* `Wpoly_dgraph_joint` — if both hubs are active, the two hubs and all the active
  prefixes form one component of imbalance `|p - q|`, and the weight is
  `(z^|p-q| + z^(-|p-q|)) * ∏ tails`.
-/

namespace ClanAudit

open Finset LaurentPolynomial SimpleGraph

variable {m n : ℕ} {a : Fin m → ℕ} {b : Fin n → ℕ}

/-! ### Splitting when a hub is switched off -/

/-- **If one of the two hubs is inactive the tree splits into its two spiders.** -/
theorem Wpoly_dgraph_split (su : SpiderV m a → ℕ) (sv : SpiderV n b → ℕ)
    (h : su none = 0 ∨ sv none = 0) :
    Wpoly (dgraph m a n b) (Sum.elim su sv)
      = Wpoly (spider m a) su * Wpoly (spider n b) sv := by
  refine Wpoly_split_sum (dgraph m a n b) (Sum.elim su sv) (spider m a) (spider n b)
    (fun _ _ => Iff.rfl) (fun _ _ => Iff.rfl) ?_
  intro x y hadj
  obtain ⟨hx, hy⟩ := dgraph_adj_lr.mp hadj
  subst hx; subst hy
  exact h

/-! ### Cutting the tree at the active prefixes -/

/-- Regrouping a sum of two sums. -/
def sumShuffle (A B C D : Type*) : (A ⊕ B) ⊕ (C ⊕ D) ≃ (A ⊕ C) ⊕ (B ⊕ D) where
  toFun x :=
    match x with
    | Sum.inl (Sum.inl u) => Sum.inl (Sum.inl u)
    | Sum.inl (Sum.inr u) => Sum.inr (Sum.inl u)
    | Sum.inr (Sum.inl u) => Sum.inl (Sum.inr u)
    | Sum.inr (Sum.inr u) => Sum.inr (Sum.inr u)
  invFun x :=
    match x with
    | Sum.inl (Sum.inl u) => Sum.inl (Sum.inl u)
    | Sum.inl (Sum.inr u) => Sum.inr (Sum.inl u)
    | Sum.inr (Sum.inl u) => Sum.inl (Sum.inr u)
    | Sum.inr (Sum.inr u) => Sum.inr (Sum.inr u)
  left_inv := by rintro ((u | u) | (u | u)) <;> rfl
  right_inv := by rintro ((u | u) | (u | u)) <;> rfl

/-- The vertex splitting of the adjacent two-hub tree into the sub-tree spanned by the
two hubs and the active prefixes of all arms, and the family of arm tails. -/
def dCut (m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) (pa : Fin m → ℕ) (qb : Fin n → ℕ)
    (hpa : ∀ i, pa i ≤ a i) (hqb : ∀ j, qb j ≤ b j) :
    DV m a n b ≃
      DV m pa n qb ⊕ ((Σ i : Fin m, Fin (a i - pa i)) ⊕ (Σ j : Fin n, Fin (b j - qb j))) :=
  (Equiv.sumCongr (spiderCut m a pa hpa) (spiderCut n b qb hqb)).trans (sumShuffle _ _ _ _)

variable {pa : Fin m → ℕ} {qb : Fin n → ℕ} {hpa : ∀ i, pa i ≤ a i} {hqb : ∀ j, qb j ≤ b j}

theorem dCut_symm_ll (x : SpiderV m pa) :
    (dCut m a n b pa qb hpa hqb).symm (Sum.inl (Sum.inl x))
      = Sum.inl ((spiderCut m a pa hpa).symm (Sum.inl x)) := rfl

theorem dCut_symm_lr (x : SpiderV n qb) :
    (dCut m a n b pa qb hpa hqb).symm (Sum.inl (Sum.inr x))
      = Sum.inr ((spiderCut n b qb hqb).symm (Sum.inl x)) := rfl

theorem dCut_symm_rl (x : Σ i : Fin m, Fin (a i - pa i)) :
    (dCut m a n b pa qb hpa hqb).symm (Sum.inr (Sum.inl x))
      = Sum.inl ((spiderCut m a pa hpa).symm (Sum.inr x)) := rfl

theorem dCut_symm_rr (x : Σ j : Fin n, Fin (b j - qb j)) :
    (dCut m a n b pa qb hpa hqb).symm (Sum.inr (Sum.inr x))
      = Sum.inr ((spiderCut n b qb hqb).symm (Sum.inr x)) := rfl

theorem spiderCut_symm_inl_eq_none_iff (k : ℕ) (len pref : Fin k → ℕ)
    (hpref : ∀ i, pref i ≤ len i) (x : SpiderV k pref) :
    (spiderCut k len pref hpref).symm (Sum.inl x) = none ↔ x = none := by
  match x with
  | none => exact ⟨fun _ => rfl, fun _ => rfl⟩
  | some ⟨i, j⟩ => exact ⟨fun h => absurd h (by rw [spiderCut_symm_inl_some]; simp), fun h => by
      exact absurd h (by simp)⟩

theorem spiderCut_symm_inr_ne_none (k : ℕ) (len pref : Fin k → ℕ)
    (hpref : ∀ i, pref i ≤ len i) (x : Σ i : Fin k, Fin (len i - pref i)) :
    (spiderCut k len pref hpref).symm (Sum.inr x) ≠ none := by
  obtain ⟨i, j⟩ := x
  rw [spiderCut_symm_inr]
  simp

/-- **Cutting the adjacent two-hub tree at the active prefixes.** -/
theorem Wpoly_dgraph_cut (su : SpiderV m a → ℕ) (sv : SpiderV n b → ℕ)
    (hcu : ∀ (i : Fin m) (h : pa i < a i), su (some ⟨i, ⟨pa i, h⟩⟩) = 0)
    (hcv : ∀ (j : Fin n) (h : qb j < b j), sv (some ⟨j, ⟨qb j, h⟩⟩) = 0) :
    Wpoly (dgraph m a n b) (Sum.elim su sv)
      = Wpoly (dgraph m pa n qb)
          (fun x => Sum.elim su sv ((dCut m a n b pa qb hpa hqb).symm (Sum.inl x)))
        * Wpoly ((ArmsG m (fun i => a i - pa i)).sum (ArmsG n (fun j => b j - qb j)))
          (fun y => Sum.elim su sv ((dCut m a n b pa qb hpa hqb).symm (Sum.inr y))) := by
  refine Wpoly_split_equiv (dgraph m a n b) (Sum.elim su sv) (dgraph m pa n qb)
    ((ArmsG m (fun i => a i - pa i)).sum (ArmsG n (fun j => b j - qb j)))
    (dCut m a n b pa qb hpa hqb) ?_ ?_ ?_
  · rintro (x | x) (y | y)
    · rw [dCut_symm_ll, dCut_symm_ll]
      exact spiderCut_adj_core m a pa hpa x y
    · rw [dCut_symm_ll, dCut_symm_lr, dgraph_adj_lr, dgraph_adj_lr]
      rw [spiderCut_symm_inl_eq_none_iff, spiderCut_symm_inl_eq_none_iff]
    · rw [dCut_symm_lr, dCut_symm_ll, dgraph_adj_rl, dgraph_adj_rl]
      rw [spiderCut_symm_inl_eq_none_iff, spiderCut_symm_inl_eq_none_iff]
    · rw [dCut_symm_lr, dCut_symm_lr]
      exact spiderCut_adj_core n b qb hqb x y
  · rintro (x | x) (y | y)
    · rw [dCut_symm_rl, dCut_symm_rl, SimpleGraph.sum_adj]
      exact spiderCut_adj_tail m a pa hpa x y
    · rw [dCut_symm_rl, dCut_symm_rr, dgraph_adj_lr, SimpleGraph.sum_adj]
      constructor
      · rintro ⟨h1, h2⟩
        exact absurd h1 (spiderCut_symm_inr_ne_none m a pa hpa x)
      · exact False.elim
    · rw [dCut_symm_rr, dCut_symm_rl, dgraph_adj_rl, SimpleGraph.sum_adj]
      constructor
      · rintro ⟨h1, h2⟩
        exact absurd h1 (spiderCut_symm_inr_ne_none n b qb hqb x)
      · exact False.elim
    · rw [dCut_symm_rr, dCut_symm_rr, SimpleGraph.sum_adj]
      exact spiderCut_adj_tail n b qb hqb x y
  · rintro (x | x) (y | y) hadj
    · rw [dCut_symm_ll, dCut_symm_rl] at hadj ⊢
      exact spiderCut_cross m a pa hpa su hcu x y (dgraph_adj_ll.mp hadj)
    · rw [dCut_symm_ll, dCut_symm_rr] at hadj
      exact absurd (dgraph_adj_lr.mp hadj).2 (spiderCut_symm_inr_ne_none n b qb hqb y)
    · rw [dCut_symm_lr, dCut_symm_rl] at hadj
      exact absurd (dgraph_adj_rl.mp hadj).2 (spiderCut_symm_inr_ne_none m a pa hpa y)
    · rw [dCut_symm_lr, dCut_symm_rr] at hadj ⊢
      exact spiderCut_cross n b qb hqb sv hcv x y (dgraph_adj_rr.mp hadj)

/-- **Both hubs active.**  The two hubs together with all the active prefixes form a
single component, of imbalance `|p - q|`. -/
theorem Wpoly_dgraph_joint (su : SpiderV m a → ℕ) (sv : SpiderV n b → ℕ)
    (hu : su none = 1) (hv : sv none = 1)
    (hsu : ∀ i : Fin m, pref su i < a i → armVal su i (pref su i) = 0)
    (hsv : ∀ j : Fin n, pref sv j < b j → armVal sv j (pref sv j) = 0) :
    Wpoly (dgraph m a n b) (Sum.elim su sv)
      = Pw (((pNum su : ℤ) - pNum sv).natAbs)
        * ((∏ i : Fin m, tailW su i) * ∏ j : Fin n, tailW sv j) := by
  classical
  have hcu : ∀ (i : Fin m) (h : pref su i < a i), su (some ⟨i, ⟨pref su i, h⟩⟩) = 0 := by
    intro i h
    have := hsu i h
    rwa [armVal, dif_pos h] at this
  have hcv : ∀ (j : Fin n) (h : pref sv j < b j), sv (some ⟨j, ⟨pref sv j, h⟩⟩) = 0 := by
    intro j h
    have := hsv j h
    rwa [armVal, dif_pos h] at this
  rw [Wpoly_dgraph_cut (pa := pref su) (qb := pref sv) (hpa := pref_le su) (hqb := pref_le sv)
    su sv hcu hcv]
  congr 1
  · -- the core carries all multiplicities one
    have hones : (fun x => Sum.elim su sv
        ((dCut m a n b (pref su) (pref sv) (pref_le su) (pref_le sv)).symm (Sum.inl x)))
        = fun _ : DV m (pref su) n (pref sv) => 1 := by
      funext x
      match x with
      | Sum.inl none => exact hu
      | Sum.inl (some ⟨i, j⟩) =>
          rw [dCut_symm_ll, spiderCut_symm_inl_some]
          have := pref_ones su i j.isLt
          rwa [armVal, dif_pos (show j.val < a i from by
            have := j.isLt; have := pref_le su i; omega)] at this
      | Sum.inr none => exact hv
      | Sum.inr (some ⟨j, l⟩) =>
          rw [dCut_symm_lr, spiderCut_symm_inl_some]
          have := pref_ones sv j l.isLt
          rwa [armVal, dif_pos (show l.val < b j from by
            have := l.isLt; have := pref_le sv j; omega)] at this
    rw [hones]
    exact Wpoly_dgraph_ones m (pref su) n (pref sv) (pNum su) (pNum sv) _ rfl rfl rfl
  · -- the tails
    have helim : (fun y => Sum.elim su sv
        ((dCut m a n b (pref su) (pref sv) (pref_le su) (pref_le sv)).symm (Sum.inr y)))
        = Sum.elim
            (fun t : Σ i : Fin m, Fin (a i - pref su i) =>
              su ((spiderCut m a (pref su) (pref_le su)).symm (Sum.inr t)))
            (fun t : Σ j : Fin n, Fin (b j - pref sv j) =>
              sv ((spiderCut n b (pref sv) (pref_le sv)).symm (Sum.inr t))) := by
      funext y
      rcases y with t | t
      · rfl
      · rfl
    rw [helim, Wpoly_sum_elim, Wpoly_ArmsG, Wpoly_ArmsG]
    congr 1
    · refine Finset.prod_congr rfl fun i _ => ?_
      refine congrArg _ (funext fun j => ?_)
      rw [spiderCut_symm_inr, armVal,
        dif_pos (show pref su i + j.val < a i from by have := j.isLt; omega)]
    · refine Finset.prod_congr rfl fun j _ => ?_
      refine congrArg _ (funext fun l => ?_)
      rw [spiderCut_symm_inr, armVal,
        dif_pos (show pref sv j + l.val < b j from by have := l.isLt; omega)]

end ClanAudit
