import RequestProject.ClanAudit.ArmSpiderR
import RequestProject.ClanAudit.ConnectorImbalance

/-!
# Cutting the connector tree

The connector tree `connGraph t m a n b` falls apart, once the vertices of multiplicity
zero are deleted, into pieces that are all of the shape already analysed:

* if the left hub is switched off, into `spider m a` and the mirrored arm spider
  `armSpiderR n b t` (`Wpoly_conn_hubL_zero`);
* if the right hub is switched off, into `armSpider m a t` and `spider n b`
  (`Wpoly_conn_hubR_zero`);
* if an interior connector vertex `i` is switched off, into `armSpider m a (i+1)` and
  `armSpiderR n b (t - i - 1)` (`Wpoly_conn_interior_zero`) — the connector prefix on
  each side becomes an ordinary ordered arm of the hub it is attached to.

The file also records the vanishing lemmas: a multiplicity at least two next to a
positive connector vertex produces a triangle, hence weight zero.
-/

namespace ClanAudit

open Finset LaurentPolynomial SimpleGraph

variable {t m n : ℕ} {a : Fin m → ℕ} {b : Fin n → ℕ}

/-! ### Reversing a path -/

/-- Reversing a path is an automorphism of the path graph. -/
def pathRevIso (L : ℕ) : pathGraph L ≃g pathGraph L where
  toEquiv := finRev L
  map_rel_iff' := by
    intro i j
    rw [pathGraph_adj, pathGraph_adj]
    have hi := i.isLt
    have hj := j.isLt
    show ((L - 1 - i.val) + 1 = L - 1 - j.val ∨ (L - 1 - j.val) + 1 = L - 1 - i.val)
      ↔ (i.val + 1 = j.val ∨ j.val + 1 = i.val)
    omega

theorem finRev_symm (L : ℕ) : (finRev L).symm = finRev L := rfl

theorem Wpoly_pathGraph_rev (L : ℕ) (g : Fin L → ℕ) :
    Wpoly (pathGraph L) (fun i => g ((finRev L).symm i)) = Wpoly (pathGraph L) g := by
  rw [← Wpoly_of_iso (pathRevIso L) g]
  rfl

/-! ### Switching off the left hub -/

/-- **The left hub is switched off.**  The connector and the whole right half form a
mirrored arm spider. -/
theorem Wpoly_conn_hubL_zero (x : SpiderV m a → ℕ) (g : Fin t → ℕ) (y : SpiderV n b → ℕ)
    (hx : x none = 0) :
    Wpoly (connGraph t m a n b) (Sum.elim x (Sum.elim g y))
      = Wpoly (spider m a) x * Wpoly (armSpiderR n b t) (Sum.elim g y) := by
  refine Wpoly_split_sum (connGraph t m a n b) (Sum.elim x (Sum.elim g y)) (spider m a)
    (armSpiderR n b t) (fun _ _ => Iff.rfl) ?_ ?_
  · rintro (i | p) (j | q) <;> exact Iff.rfl
  · rintro p (j | q) hadj
    · refine Or.inl ?_
      have h := (connGraph_adj_lc.mp hadj).1
      subst h
      exact hx
    · refine Or.inl ?_
      have h := (connGraph_adj_lr.mp hadj).1
      subst h
      exact hx

/-! ### Switching off the right hub -/

/-- Regrouping the vertex type of the connector tree. -/
def connAssoc (t m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) :
    ConnV t m a n b ≃ (SpiderV m a ⊕ Fin t) ⊕ SpiderV n b where
  toFun x :=
    match x with
    | Sum.inl p => Sum.inl (Sum.inl p)
    | Sum.inr (Sum.inl j) => Sum.inl (Sum.inr j)
    | Sum.inr (Sum.inr q) => Sum.inr q
  invFun x :=
    match x with
    | Sum.inl (Sum.inl p) => Sum.inl p
    | Sum.inl (Sum.inr j) => Sum.inr (Sum.inl j)
    | Sum.inr q => Sum.inr (Sum.inr q)
  left_inv := by rintro (p | j | q) <;> rfl
  right_inv := by rintro ((p | j) | q) <;> rfl

theorem connAssoc_symm_ll (p : SpiderV m a) :
    (connAssoc t m a n b).symm (Sum.inl (Sum.inl p)) = Sum.inl p := rfl

theorem connAssoc_symm_lr (j : Fin t) :
    (connAssoc t m a n b).symm (Sum.inl (Sum.inr j)) = Sum.inr (Sum.inl j) := rfl

theorem connAssoc_symm_r (q : SpiderV n b) :
    (connAssoc t m a n b).symm (Sum.inr q) = Sum.inr (Sum.inr q) := rfl

/-- **The right hub is switched off.**  The whole left half and the connector form an
arm spider. -/
theorem Wpoly_conn_hubR_zero (x : SpiderV m a → ℕ) (g : Fin t → ℕ) (y : SpiderV n b → ℕ)
    (hy : y none = 0) :
    Wpoly (connGraph t m a n b) (Sum.elim x (Sum.elim g y))
      = Wpoly (armSpider m a t) (Sum.elim x g) * Wpoly (spider n b) y := by
  have key := Wpoly_split_equiv (connGraph t m a n b) (Sum.elim x (Sum.elim g y))
    (armSpider m a t) (spider n b) (connAssoc t m a n b) ?_ ?_ ?_
  · rw [key]
    have h1 : (fun z => Sum.elim x (Sum.elim g y) ((connAssoc t m a n b).symm (Sum.inl z)))
        = Sum.elim x g := by
      funext z; rcases z with p | j <;> rfl
    have h2 : (fun q => Sum.elim x (Sum.elim g y) ((connAssoc t m a n b).symm (Sum.inr q)))
        = y := by
      funext q; rfl
    rw [h1, h2]
  · rintro (p | j) (q | k)
    · rw [connAssoc_symm_ll, connAssoc_symm_ll]; exact Iff.rfl
    · rw [connAssoc_symm_ll, connAssoc_symm_lr]; exact Iff.rfl
    · rw [connAssoc_symm_lr, connAssoc_symm_ll]; exact Iff.rfl
    · rw [connAssoc_symm_lr, connAssoc_symm_lr]; exact Iff.rfl
  · intro q q'
    rw [connAssoc_symm_r, connAssoc_symm_r]
    exact Iff.rfl
  · rintro (p | j) q hadj
    · rw [connAssoc_symm_ll, connAssoc_symm_r] at hadj ⊢
      refine Or.inr ?_
      have h := (connGraph_adj_lr.mp hadj).2.1
      subst h
      exact hy
    · rw [connAssoc_symm_lr, connAssoc_symm_r] at hadj ⊢
      refine Or.inr ?_
      have h := (connGraph_adj_cr.mp hadj).1
      subst h
      exact hy

/-- **Both hubs are switched off.**  The connector is an isolated path. -/
theorem Wpoly_conn_bothHubs_zero (x : SpiderV m a → ℕ) (g : Fin t → ℕ) (y : SpiderV n b → ℕ)
    (hx : x none = 0) (hy : y none = 0) :
    Wpoly (connGraph t m a n b) (Sum.elim x (Sum.elim g y))
      = Wpoly (spider m a) x * (Wpoly (pathGraph t) g * Wpoly (spider n b) y) := by
  rw [Wpoly_conn_hubL_zero x g y hx, Wpoly_armSpiderR,
    Wpoly_armSpider_hub_zero y (fun i => g ((finRev t).symm i)) hy, Wpoly_pathGraph_rev,
    mul_comm (Wpoly (spider n b) y)]

/-! ### Switching off an interior connector vertex -/

/-- The vertex splitting of the connector tree at an interior connector vertex `i`: the
left half, ending at `i`, and the right half. -/
def connSplitE (t m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) (i : ℕ) (hi : i < t) :
    ConnV t m a n b ≃ (SpiderV m a ⊕ Fin (i + 1)) ⊕ (Fin (t - i - 1) ⊕ SpiderV n b) where
  toFun x :=
    match x with
    | Sum.inl p => Sum.inl (Sum.inl p)
    | Sum.inr (Sum.inl j) =>
        if h : j.val ≤ i then Sum.inl (Sum.inr ⟨j.val, by omega⟩)
        else Sum.inr (Sum.inl ⟨j.val - i - 1, by have := j.isLt; omega⟩)
    | Sum.inr (Sum.inr q) => Sum.inr (Sum.inr q)
  invFun x :=
    match x with
    | Sum.inl (Sum.inl p) => Sum.inl p
    | Sum.inl (Sum.inr k) => Sum.inr (Sum.inl ⟨k.val, by have := k.isLt; omega⟩)
    | Sum.inr (Sum.inl k) => Sum.inr (Sum.inl ⟨k.val + i + 1, by have := k.isLt; omega⟩)
    | Sum.inr (Sum.inr q) => Sum.inr (Sum.inr q)
  left_inv := by
    rintro (p | j | q)
    · rfl
    · by_cases h : j.val ≤ i
      · simp only [dif_pos h]
      · simp only [dif_neg h]
        exact congrArg (fun z => Sum.inr (Sum.inl z)) (Fin.ext (by simp only; omega))
    · rfl
  right_inv := by
    rintro ((p | k) | (k | q))
    · rfl
    · simp only [dif_pos (show (⟨k.val, by have := k.isLt; omega⟩ : Fin t).val ≤ i from by
        have := k.isLt; simp only; omega)]
    · simp only [dif_neg (show ¬ ((⟨k.val + i + 1, by have := k.isLt; omega⟩ : Fin t).val ≤ i)
        from by simp only; omega)]
      exact congrArg (fun z => Sum.inr (Sum.inl z)) (Fin.ext (by simp only; omega))
    · rfl

section SplitLemmas
variable {i : ℕ} {hi : i < t}

theorem connSplitE_symm_ll (p : SpiderV m a) :
    (connSplitE t m a n b i hi).symm (Sum.inl (Sum.inl p)) = Sum.inl p := rfl

theorem connSplitE_symm_lr (k : Fin (i + 1)) :
    (connSplitE t m a n b i hi).symm (Sum.inl (Sum.inr k))
      = Sum.inr (Sum.inl ⟨k.val, by have := k.isLt; omega⟩) := rfl

theorem connSplitE_symm_rl (k : Fin (t - i - 1)) :
    (connSplitE t m a n b i hi).symm (Sum.inr (Sum.inl k))
      = Sum.inr (Sum.inl ⟨k.val + i + 1, by have := k.isLt; omega⟩) := rfl

theorem connSplitE_symm_rr (q : SpiderV n b) :
    (connSplitE t m a n b i hi).symm (Sum.inr (Sum.inr q)) = Sum.inr (Sum.inr q) := rfl

end SplitLemmas

/-- **An interior connector vertex is switched off.**  The tree splits into an arm
spider on the left and a mirrored arm spider on the right. -/
theorem Wpoly_conn_interior_zero (x : SpiderV m a → ℕ) (g : Fin t → ℕ) (y : SpiderV n b → ℕ)
    (i : ℕ) (hi : i < t) (hg : g ⟨i, hi⟩ = 0) :
    Wpoly (connGraph t m a n b) (Sum.elim x (Sum.elim g y))
      = Wpoly (armSpider m a (i + 1))
          (Sum.elim x (fun k : Fin (i + 1) => g ⟨k.val, by have := k.isLt; omega⟩))
        * Wpoly (armSpiderR n b (t - i - 1))
          (Sum.elim (fun k : Fin (t - i - 1) => g ⟨k.val + i + 1, by have := k.isLt; omega⟩)
            y) := by
  have key := Wpoly_split_equiv (connGraph t m a n b) (Sum.elim x (Sum.elim g y))
    (armSpider m a (i + 1)) (armSpiderR n b (t - i - 1)) (connSplitE t m a n b i hi) ?_ ?_ ?_
  · rw [key]
    congr 1
    · congr 1
      funext z
      rcases z with p | k
      · rw [connSplitE_symm_ll]; rfl
      · rw [connSplitE_symm_lr]; rfl
    · congr 1
      funext z
      rcases z with k | q
      · rw [connSplitE_symm_rl]; rfl
      · rw [connSplitE_symm_rr]; rfl
  · rintro (p | k) (q | k')
    · rw [connSplitE_symm_ll, connSplitE_symm_ll]; exact Iff.rfl
    · rw [connSplitE_symm_ll, connSplitE_symm_lr]; exact Iff.rfl
    · rw [connSplitE_symm_lr, connSplitE_symm_ll]; exact Iff.rfl
    · rw [connSplitE_symm_lr, connSplitE_symm_lr]; exact Iff.rfl
  · rintro (k | q) (k' | q')
    · rw [connSplitE_symm_rl, connSplitE_symm_rl, connGraph_adj_cc]
      show (k.val + i + 1 + 1 = k'.val + i + 1 ∨ k'.val + i + 1 + 1 = k.val + i + 1)
        ↔ (k.val + 1 = k'.val ∨ k'.val + 1 = k.val)
      omega
    · rw [connSplitE_symm_rl, connSplitE_symm_rr, connGraph_adj_cr]
      have hk := k.isLt
      show (q' = none ∧ k.val + i + 1 + 1 = t) ↔ (q' = none ∧ k.val + 1 = t - i - 1)
      constructor
      · rintro ⟨h1, h2⟩; exact ⟨h1, by omega⟩
      · rintro ⟨h1, h2⟩; exact ⟨h1, by omega⟩
    · rw [connSplitE_symm_rr, connSplitE_symm_rl]
      have hk' := k'.isLt
      show (q = none ∧ k'.val + i + 1 + 1 = t) ↔ (q = none ∧ k'.val + 1 = t - i - 1)
      constructor
      · rintro ⟨h1, h2⟩; exact ⟨h1, by omega⟩
      · rintro ⟨h1, h2⟩; exact ⟨h1, by omega⟩
    · rw [connSplitE_symm_rr, connSplitE_symm_rr]; exact Iff.rfl
  · rintro (p | k) (k' | q) hadj
    · rw [connSplitE_symm_ll, connSplitE_symm_rl] at hadj
      exact absurd (connGraph_adj_lc.mp hadj).2 (by simp only; omega)
    · rw [connSplitE_symm_ll, connSplitE_symm_rr] at hadj
      exact absurd (connGraph_adj_lr.mp hadj).2.2 (by omega)
    · rw [connSplitE_symm_lr, connSplitE_symm_rl] at hadj ⊢
      refine Or.inl ?_
      have hk := k.isLt
      have hk' := k'.isLt
      have hval : k.val = i := by
        rcases connGraph_adj_cc.mp hadj with h | h <;> simp only at h <;> omega
      show g ⟨k.val, _⟩ = 0
      rw [show (⟨k.val, show k.val < t from by omega⟩ : Fin t) = ⟨i, hi⟩ from Fin.ext hval]
      exact hg
    · rw [connSplitE_symm_lr, connSplitE_symm_rr] at hadj ⊢
      refine Or.inl ?_
      have hk := k.isLt
      have hval : k.val = i := by
        have := (connGraph_adj_cr.mp hadj).2
        simp only at this
        omega
      show g ⟨k.val, _⟩ = 0
      rw [show (⟨k.val, show k.val < t from by omega⟩ : Fin t) = ⟨i, hi⟩ from Fin.ext hval]
      exact hg

/-! ### Vanishing -/

/-- A hub of multiplicity at least two next to a positive connector vertex. -/
theorem Wpoly_conn_eq_zero_hubL_two (x : SpiderV m a → ℕ) (g : Fin t → ℕ)
    (y : SpiderV n b → ℕ) (h0 : 0 < t) (hx : 2 ≤ x none) (hg : 1 ≤ g ⟨0, h0⟩) :
    Wpoly (connGraph t m a n b) (Sum.elim x (Sum.elim g y)) = 0 :=
  Wpoly_eq_zero_of_mult_two (G := connGraph t m a n b) (u := Sum.inl none)
    (v := Sum.inr (Sum.inl ⟨0, h0⟩)) (connGraph_adj_lc.mpr ⟨rfl, rfl⟩) hx hg

theorem Wpoly_conn_eq_zero_hubR_two (x : SpiderV m a → ℕ) (g : Fin t → ℕ)
    (y : SpiderV n b → ℕ) (h0 : 0 < t) (hy : 2 ≤ y none)
    (hg : 1 ≤ g ⟨t - 1, by omega⟩) :
    Wpoly (connGraph t m a n b) (Sum.elim x (Sum.elim g y)) = 0 :=
  Wpoly_eq_zero_of_mult_two (G := connGraph t m a n b) (u := Sum.inr (Sum.inr none))
    (v := Sum.inr (Sum.inl ⟨t - 1, by omega⟩))
    (connGraph_adj_cr.mpr ⟨rfl, by simp only; omega⟩).symm hy hg

/-- A left side with a positive hub whose initial run is not made of ones. -/
theorem Wpoly_conn_eq_zero_of_not_admRun_left (x : SpiderV m a → ℕ) (g : Fin t → ℕ)
    (y : SpiderV n b → ℕ) (hhub : 0 < x none) (h : ¬ AdmRun x) :
    Wpoly (connGraph t m a n b) (Sum.elim x (Sum.elim g y)) = 0 := by
  obtain ⟨p, q, hadj, hp, hq⟩ := exists_spider_triangle_of_not_admRun hhub h
  exact Wpoly_eq_zero_of_mult_two (G := connGraph t m a n b) (u := Sum.inl p) (v := Sum.inl q)
    (connGraph_adj_ll.mpr hadj) hp hq

theorem Wpoly_conn_eq_zero_of_not_admRun_right (x : SpiderV m a → ℕ) (g : Fin t → ℕ)
    (y : SpiderV n b → ℕ) (hhub : 0 < y none) (h : ¬ AdmRun y) :
    Wpoly (connGraph t m a n b) (Sum.elim x (Sum.elim g y)) = 0 := by
  obtain ⟨p, q, hadj, hp, hq⟩ := exists_spider_triangle_of_not_admRun hhub h
  exact Wpoly_eq_zero_of_mult_two (G := connGraph t m a n b) (u := Sum.inr (Sum.inr p))
    (v := Sum.inr (Sum.inr q)) (connGraph_adj_rr.mpr hadj) hp hq

/-- A connector vertex of multiplicity at least two, in a connector all of whose
vertices are positive, next to a positive left hub. -/
theorem Wpoly_conn_eq_zero_of_conn_two_left (x : SpiderV m a → ℕ) (g : Fin t → ℕ)
    (y : SpiderV n b → ℕ) (j : Fin t) (hj : 2 ≤ g j) (hpos : ∀ k, 1 ≤ g k)
    (hx : 1 ≤ x none) :
    Wpoly (connGraph t m a n b) (Sum.elim x (Sum.elim g y)) = 0 := by
  rcases Nat.eq_zero_or_pos j.val with h0 | h0
  · refine Wpoly_eq_zero_of_mult_two (G := connGraph t m a n b) (u := Sum.inr (Sum.inl j))
      (v := Sum.inl none) (connGraph_adj_lc.mpr ⟨rfl, h0⟩).symm hj hx
  · have hlt : j.val - 1 < t := by have := j.isLt; omega
    refine Wpoly_eq_zero_of_mult_two (G := connGraph t m a n b) (u := Sum.inr (Sum.inl j))
      (v := Sum.inr (Sum.inl ⟨j.val - 1, hlt⟩)) (connGraph_adj_cc.mpr (Or.inr ?_)) hj
      (hpos _)
    simp only
    omega

theorem Wpoly_conn_eq_zero_of_conn_two_right (x : SpiderV m a → ℕ) (g : Fin t → ℕ)
    (y : SpiderV n b → ℕ) (j : Fin t) (hj : 2 ≤ g j) (hpos : ∀ k, 1 ≤ g k)
    (hy : 1 ≤ y none) :
    Wpoly (connGraph t m a n b) (Sum.elim x (Sum.elim g y)) = 0 := by
  have hjt := j.isLt
  rcases Nat.lt_or_ge (j.val + 1) t with h0 | h0
  · refine Wpoly_eq_zero_of_mult_two (G := connGraph t m a n b) (u := Sum.inr (Sum.inl j))
      (v := Sum.inr (Sum.inl ⟨j.val + 1, h0⟩)) (connGraph_adj_cc.mpr (Or.inl ?_)) hj (hpos _)
    simp only
  · refine Wpoly_eq_zero_of_mult_two (G := connGraph t m a n b) (u := Sum.inr (Sum.inl j))
      (v := Sum.inr (Sum.inr none)) (connGraph_adj_cr.mpr ⟨rfl, by omega⟩) hj hy

end ClanAudit
