import RequestProject.ClanAudit.PathClan

/-!
# Arms of a spider

`ArmsG k len` is the disjoint union of the `k` pendant paths of a spider, i.e. the spider
with its hub removed.  This file proves

* `Wpoly_ArmsG` — the normalized weight of an arm family is the product of the weights of
  the individual paths;
* `Wpoly_spider_hub_zero` — if the hub has multiplicity zero, the weight of a spider clan
  is the product of the weights of its arms;
* `Wpoly_spider_cut` — a spider clan may be cut into the sub-spider spanned by an initial
  prefix of each arm and the family of remaining arm tails, provided the multiplicity at
  the first vertex beyond each prefix is zero.
-/

namespace ClanAudit

open Finset LaurentPolynomial SimpleGraph

/-- The disjoint union of the `k` arms of a spider: the `i`-th arm is a path with
`len i` vertices. -/
def ArmsG (k : ℕ) (len : Fin k → ℕ) : SimpleGraph (Σ i : Fin k, Fin (len i)) where
  Adj p q := p.1 = q.1 ∧ (p.2.val + 1 = q.2.val ∨ q.2.val + 1 = p.2.val)
  symm := by
    rintro ⟨i, j⟩ ⟨i', j'⟩ h
    exact ⟨h.1.symm, h.2.symm⟩
  loopless := ⟨by rintro ⟨i, j⟩ h; rcases h.2 with h' | h' <;> omega⟩

theorem ArmsG_adj {k : ℕ} {len : Fin k → ℕ} {p q : Σ i : Fin k, Fin (len i)} :
    (ArmsG k len).Adj p q ↔ p.1 = q.1 ∧ (p.2.val + 1 = q.2.val ∨ q.2.val + 1 = p.2.val) :=
  Iff.rfl

/-! ### Peeling off the first arm -/

/-- `Σ i : Fin (k+1), Fin (len i)` splits off the arm `0`. -/
def sigmaFinSucc (k : ℕ) (len : Fin (k + 1) → ℕ) :
    (Σ i : Fin (k + 1), Fin (len i)) ≃ Fin (len 0) ⊕ (Σ i : Fin k, Fin (len i.succ)) where
  toFun p :=
    Fin.cases (motive := fun i => Fin (len i) → Fin (len 0) ⊕ (Σ i : Fin k, Fin (len i.succ)))
      (fun j => Sum.inl j) (fun i j => Sum.inr ⟨i, j⟩) p.1 p.2
  invFun y :=
    match y with
    | Sum.inl j => ⟨0, j⟩
    | Sum.inr ⟨i, j⟩ => ⟨i.succ, j⟩
  left_inv := by
    rintro ⟨i, j⟩
    induction i using Fin.cases with
    | zero => rfl
    | succ i => rfl
  right_inv := by
    rintro (j | ⟨i, j⟩) <;> rfl

theorem sigmaFinSucc_symm_inl (k : ℕ) (len : Fin (k + 1) → ℕ) (j : Fin (len 0)) :
    (sigmaFinSucc k len).symm (Sum.inl j) = ⟨0, j⟩ := rfl

theorem sigmaFinSucc_symm_inr (k : ℕ) (len : Fin (k + 1) → ℕ)
    (p : Σ i : Fin k, Fin (len i.succ)) :
    (sigmaFinSucc k len).symm (Sum.inr p) = ⟨p.1.succ, p.2⟩ := rfl

theorem Wpoly_ArmsG : ∀ (k : ℕ) (len : Fin k → ℕ) (alpha : (Σ i : Fin k, Fin (len i)) → ℕ),
    Wpoly (ArmsG k len) alpha
      = ∏ i : Fin k, Wpoly (pathGraph (len i)) (fun j => alpha ⟨i, j⟩)
  | 0, len, alpha => by
      haveI : IsEmpty (Σ i : Fin 0, Fin (len i)) := ⟨fun p => p.1.elim0⟩
      rw [Wpoly_isEmpty]
      simp
  | (k + 1), len, alpha => by
      have key := Wpoly_split_equiv (ArmsG (k + 1) len) alpha (pathGraph (len 0))
        (ArmsG k (fun i => len i.succ)) (sigmaFinSucc k len) ?_ ?_ ?_
      · rw [key, Wpoly_ArmsG k (fun i => len i.succ) _, Fin.prod_univ_succ]
        rfl
      · intro x y
        rw [sigmaFinSucc_symm_inl, sigmaFinSucc_symm_inl, ArmsG_adj, pathGraph_adj]
        simp
      · rintro ⟨i, j⟩ ⟨i', j'⟩
        rw [sigmaFinSucc_symm_inr, sigmaFinSucc_symm_inr, ArmsG_adj, ArmsG_adj]
        simp [Fin.succ_inj]
      · rintro x ⟨i, j⟩ hadj
        rw [sigmaFinSucc_symm_inl, sigmaFinSucc_symm_inr, ArmsG_adj] at hadj
        exact absurd hadj.1.symm (Fin.succ_ne_zero i)

/-! ### Spiders with an inactive hub -/

theorem Wpoly_of_all_zero {V : Type*} [Fintype V] [DecidableEq V] (G : SimpleGraph V)
    (alpha : V → ℕ) (h : ∀ v, alpha v = 0) : Wpoly G alpha = 1 := by
  haveI : IsEmpty (Σ v : V, Fin (alpha v)) :=
    ⟨fun x => by have hx := x.2.isLt; have hz := h x.1; omega⟩
  rw [Wpoly, imbalanceGF_isEmpty, clanFactorial]
  simp [h]

/-- `Option X` splits off its base point. -/
def optionSumUnit (X : Type*) : Option X ≃ X ⊕ Unit where
  toFun o := o.elim (Sum.inr ()) Sum.inl
  invFun y := Sum.elim some (fun _ => none) y
  left_inv := by rintro (_ | x) <;> rfl
  right_inv := by rintro (x | ⟨⟩) <;> rfl

/-- **A spider whose hub is inactive** is the disjoint union of its arms. -/
theorem Wpoly_spider_hub_zero (k : ℕ) (len : Fin k → ℕ) (s : SpiderV k len → ℕ)
    (hhub : s none = 0) :
    Wpoly (spider k len) s
      = ∏ i : Fin k, Wpoly (pathGraph (len i)) (fun j => s (some ⟨i, j⟩)) := by
  have key := Wpoly_split_equiv (spider k len) s (ArmsG k len) (⊥ : SimpleGraph Unit)
    (optionSumUnit _) ?_ ?_ ?_
  · rw [key, Wpoly_ArmsG]
    have h1 : Wpoly (⊥ : SimpleGraph Unit)
        (fun y => s ((optionSumUnit _).symm (Sum.inr y))) = 1 :=
      Wpoly_of_all_zero _ _ fun _ => hhub
    rw [h1, mul_one]
    rfl
  · intro x y
    exact Iff.rfl
  · intro x y
    simp [spider]
  · intro x y _
    right
    exact hhub

/-! ### Cutting a spider at a prefix of each arm -/

/-- The vertex splitting of a spider into the sub-spider spanned by the first `pref i`
vertices of each arm and the family of remaining arm tails. -/
def spiderCut (k : ℕ) (len pref : Fin k → ℕ) (hpref : ∀ i, pref i ≤ len i) :
    SpiderV k len ≃ SpiderV k pref ⊕ (Σ i : Fin k, Fin (len i - pref i)) where
  toFun x :=
    match x with
    | none => Sum.inl none
    | some ⟨i, j⟩ =>
        if h : j.val < pref i then Sum.inl (some ⟨i, ⟨j.val, h⟩⟩)
        else Sum.inr ⟨i, ⟨j.val - pref i, by have := j.isLt; omega⟩⟩
  invFun y :=
    match y with
    | Sum.inl none => none
    | Sum.inl (some ⟨i, j⟩) => some ⟨i, ⟨j.val, by have := j.isLt; have := hpref i; omega⟩⟩
    | Sum.inr ⟨i, j⟩ => some ⟨i, ⟨pref i + j.val, by have := j.isLt; have := hpref i; omega⟩⟩
  left_inv := by
    rintro (_ | ⟨i, j⟩)
    · rfl
    · dsimp only
      by_cases h : j.val < pref i
      · rw [dif_pos h]
      · rw [dif_neg h]
        dsimp only
        have hj : (⟨pref i + (j.val - pref i), by have := j.isLt; omega⟩ : Fin (len i)) = j :=
          Fin.ext (by simp only; omega)
        rw [hj]
  right_inv := by
    rintro ((_ | ⟨i, j⟩) | ⟨i, j⟩)
    · rfl
    · dsimp only
      rw [dif_pos (show j.val < pref i from j.isLt)]
    · dsimp only
      rw [dif_neg (show ¬ (pref i + j.val < pref i) by omega)]
      have hj : (⟨pref i + j.val - pref i, by have := j.isLt; omega⟩ :
          Fin (len i - pref i)) = j := Fin.ext (by simp only; omega)
      rw [hj]

theorem spiderCut_symm_inl_none (k : ℕ) (len pref : Fin k → ℕ) (hpref : ∀ i, pref i ≤ len i) :
    (spiderCut k len pref hpref).symm (Sum.inl none) = none := rfl

theorem spiderCut_symm_inl_some (k : ℕ) (len pref : Fin k → ℕ) (hpref : ∀ i, pref i ≤ len i)
    (i : Fin k) (j : Fin (pref i)) :
    (spiderCut k len pref hpref).symm (Sum.inl (some ⟨i, j⟩))
      = some ⟨i, ⟨j.val, by have := j.isLt; have := hpref i; omega⟩⟩ := rfl

theorem spiderCut_symm_inr (k : ℕ) (len pref : Fin k → ℕ) (hpref : ∀ i, pref i ≤ len i)
    (i : Fin k) (j : Fin (len i - pref i)) :
    (spiderCut k len pref hpref).symm (Sum.inr ⟨i, j⟩)
      = some ⟨i, ⟨pref i + j.val, by have := j.isLt; have := hpref i; omega⟩⟩ := rfl

/-- The core of the cut spans the same edges as the spider on the prefixes. -/
theorem spiderCut_adj_core (k : ℕ) (len pref : Fin k → ℕ) (hpref : ∀ i, pref i ≤ len i)
    (x y : SpiderV k pref) :
    (spider k len).Adj ((spiderCut k len pref hpref).symm (Sum.inl x))
        ((spiderCut k len pref hpref).symm (Sum.inl y))
      ↔ (spider k pref).Adj x y := by
  revert x y
  rintro (_ | ⟨i, j⟩) (_ | ⟨i', j'⟩)
  · exact Iff.rfl
  · rw [spiderCut_symm_inl_none, spiderCut_symm_inl_some]
    exact Iff.rfl
  · rw [spiderCut_symm_inl_none, spiderCut_symm_inl_some]
    exact Iff.rfl
  · rw [spiderCut_symm_inl_some, spiderCut_symm_inl_some]
    exact Iff.rfl

/-- The tails of the cut span the same edges as the family of arm tails. -/
theorem spiderCut_adj_tail (k : ℕ) (len pref : Fin k → ℕ) (hpref : ∀ i, pref i ≤ len i)
    (x y : Σ i : Fin k, Fin (len i - pref i)) :
    (spider k len).Adj ((spiderCut k len pref hpref).symm (Sum.inr x))
        ((spiderCut k len pref hpref).symm (Sum.inr y))
      ↔ (ArmsG k (fun i => len i - pref i)).Adj x y := by
  revert x y
  rintro ⟨i, j⟩ ⟨i', j'⟩
  rw [spiderCut_symm_inr, spiderCut_symm_inr, ArmsG_adj]
  show (i = i' ∧ _) ↔ (i = i' ∧ _)
  constructor
  · rintro ⟨h1, h2⟩
    subst h1
    exact ⟨rfl, by simp only at h2 ⊢; omega⟩
  · rintro ⟨h1, h2⟩
    subst h1
    exact ⟨rfl, by simp only at h2 ⊢; omega⟩

/-- Every edge of the spider joining the core of the cut to a tail has an endpoint of
multiplicity zero. -/
theorem spiderCut_cross (k : ℕ) (len pref : Fin k → ℕ) (hpref : ∀ i, pref i ≤ len i)
    (s : SpiderV k len → ℕ)
    (hcut : ∀ (i : Fin k) (h : pref i < len i), s (some ⟨i, ⟨pref i, h⟩⟩) = 0)
    (x : SpiderV k pref) (y : Σ i : Fin k, Fin (len i - pref i)) :
    (spider k len).Adj ((spiderCut k len pref hpref).symm (Sum.inl x))
        ((spiderCut k len pref hpref).symm (Sum.inr y)) →
      s ((spiderCut k len pref hpref).symm (Sum.inl x)) = 0 ∨
        s ((spiderCut k len pref hpref).symm (Sum.inr y)) = 0 := by
  revert x y
  rintro (_ | ⟨i, j⟩) ⟨i', j'⟩ hadj
  · rw [spiderCut_symm_inl_none, spiderCut_symm_inr] at hadj ⊢
    have hz : pref i' + j'.val = 0 := hadj
    have hlt : pref i' < len i' := by have := j'.isLt; omega
    right
    have hj : (⟨pref i' + j'.val, by have := j'.isLt; have := hpref i'; omega⟩ :
        Fin (len i')) = ⟨pref i', hlt⟩ := Fin.ext (by simp only; omega)
    rw [hj]
    exact hcut i' hlt
  · rw [spiderCut_symm_inl_some, spiderCut_symm_inr] at hadj ⊢
    obtain ⟨h1, h2⟩ := hadj
    simp only at h1
    subst h1
    simp only at h2
    have hj0 := j.isLt
    have hlt : pref i < len i := by have := j'.isLt; omega
    right
    have hj : (⟨pref i + j'.val, by have := j'.isLt; have := hpref i; omega⟩ :
        Fin (len i)) = ⟨pref i, hlt⟩ := Fin.ext (by simp only; omega)
    rw [hj]
    exact hcut i hlt

/-- **Cutting a spider clan at the prefixes.**  If the multiplicity of the first vertex
beyond the prefix of each arm is zero, the weight factors into the weight of the
sub-spider on the prefixes and the weight of the family of arm tails. -/
theorem Wpoly_spider_cut (k : ℕ) (len pref : Fin k → ℕ) (hpref : ∀ i, pref i ≤ len i)
    (s : SpiderV k len → ℕ)
    (hcut : ∀ (i : Fin k) (h : pref i < len i), s (some ⟨i, ⟨pref i, h⟩⟩) = 0) :
    Wpoly (spider k len) s
      = Wpoly (spider k pref)
          (fun x => s ((spiderCut k len pref hpref).symm (Sum.inl x)))
        * Wpoly (ArmsG k (fun i => len i - pref i))
          (fun y => s ((spiderCut k len pref hpref).symm (Sum.inr y))) :=
  Wpoly_split_equiv (spider k len) s (spider k pref)
    (ArmsG k (fun i => len i - pref i)) (spiderCut k len pref hpref)
    (spiderCut_adj_core k len pref hpref) (spiderCut_adj_tail k len pref hpref)
    (spiderCut_cross k len pref hpref s hcut)

end ClanAudit
