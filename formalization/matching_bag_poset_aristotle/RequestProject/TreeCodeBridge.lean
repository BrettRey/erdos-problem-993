import RequestProject.BagPoset
import RequestProject.CodeInvariance
import RequestProject.PascalBridge

/-!
# The code of a tree with a maximum matching is a forest-poset code

This file proves equation (2) of `matching_bag_poset_reduction_2026-08-20.md` and derives
equation (6) and the log-concavity statements of §3 for an arbitrary finite forest with a
maximum matching.

Let `D` be a `TreeMatching`: a forest with a proper 2-colouring and a maximum matching
oriented by it.  Write `P = D.Free` for the poset of free variables (whose cover graph is a
forest, `coverGraph_isAcyclic`) and let `D.Const` be the remaining, constant, coordinates:
the forced matched bags and the singleton (unmatched) bags.

* `treeCode_eq_codeRelabel` is equation (2): the code of minimum covers on the bags is,
  after the coordinate bijection `bagEquiv` and the coordinate flips `bagFlip`, exactly the
  order-ideal code `Ω(P)` together with one constant word on the `c = |D.Const|` remaining
  coordinates.
* `codeP_treeCode` is equation (6) coefficientwise: the profile polynomial of the tree code
  is `(1+t)^c I(B(P);t)`.
* `erasure_eq_erasureProfile` identifies the erasure profile `e_d = p_{M-d}` of the tree code
  with the poset erasure profile `E_d` of `PascalBridge.lean`, where `M = |bags| = α(T)`.
* `erasure_log_concave_depth_le_eight`, `erasure_log_concave_of_le_33`,
  `erasure_depth_three`, `erasure_depth_three_reserve` are the resulting inequalities.
-/

open scoped BigOperators
open Finset Polynomial

namespace MatchingBag

attribute [local instance 1] Classical.propDecidable

/-! ### Coefficients of the independence polynomial -/

section Coeff

variable {P : Type*} [Fintype P] [DecidableEq P] [PartialOrder P]
  [DecidableRel ((· ≤ ·) : P → P → Prop)]

/-- The coefficients of `I(B(P);t)` are the retained-coordinate profile of `Ω(P)`. -/
lemma coeff_indepPoly (j : ℕ) : (indepPoly P).coeff j = (codeP (idealCode P) j : ℤ) := by
  classical
  rw [indepPoly_eq_sum_codeP, Polynomial.finset_sum_coeff]
  by_cases hj : j < 2 * Fintype.card P + 1
  · rw [Finset.sum_eq_single j]
    · simp
    · intro k _ hk
      simp [Polynomial.coeff_X_pow, Ne.symm hk]
    · intro h
      exact absurd (Finset.mem_range.2 hj) h
  · have hzero : ∀ k ∈ range (2 * Fintype.card P + 1),
        (((codeP (idealCode P) k : ℤ) : Polynomial ℤ) * X ^ k).coeff j = 0 := by
      intro k hk
      rw [Finset.mem_range] at hk
      have hkj : j ≠ k := by omega
      simp [Polynomial.coeff_X_pow, hkj]
    rw [Finset.sum_eq_zero hzero, codeP_eq_zero_of_lt (by omega)]
    simp

/-- The coefficients of `(1+t)^c I(B(P);t)`. -/
lemma coeff_pow_mul_indepPoly (c k : ℕ) :
    (((1 : Polynomial ℤ) + X) ^ c * indepPoly P).coeff k
      = ∑ j ∈ range (k + 1), (codeP (idealCode P) j : ℤ) * (c.choose (k - j) : ℤ) := by
  rw [Polynomial.coeff_mul, Finset.Nat.sum_antidiagonal_eq_sum_range_succ_mk]
  conv_rhs => rw [← Finset.sum_range_reflect]
  refine Finset.sum_congr rfl fun i hi => ?_
  rw [Finset.mem_range] at hi
  have h1 : k + 1 - 1 - i = k - i := by omega
  have h2 : k - (k - i) = i := by omega
  rw [h1, h2, Polynomial.coeff_one_add_X_pow, coeff_indepPoly]
  ring

/-- The profile polynomial `∑_k p_k(Ω) t^k` of a code. -/
noncomputable def codePoly {ι : Type*} [Fintype ι] [DecidableEq ι] (Ω : Finset (ι → Bool)) :
    Polynomial ℤ :=
  ∑ k ∈ range (Fintype.card ι + 1), ((codeP Ω k : ℤ) : Polynomial ℤ) * (X : Polynomial ℤ) ^ k

end Coeff

variable {V : Type*} [Fintype V] [DecidableEq V]

namespace TreeMatching

variable (D : TreeMatching V)

/-! ### Splitting the bags into free and constant coordinates -/

/-- The constant coordinates: the forced matched bags together with the singleton bags. -/
def Const : Type _ := {i : D.Idx // D.Forced i} ⊕ D.Unm

instance : Finite D.Const := by unfold Const; infer_instance
noncomputable instance : Fintype D.Const := Fintype.ofFinite _
noncomputable instance : DecidableEq D.Const := Classical.decEq _

/-- The number `c` of constant coordinates. -/
noncomputable def numConst : ℕ := Fintype.card D.Const

/-- The bags are the free variables together with the constant coordinates. -/
noncomputable def bagEquiv : D.Bag ≃ D.Free ⊕ D.Const where
  toFun := Sum.elim
    (fun i => if h : D.Forced i then Sum.inr (Sum.inl ⟨i, h⟩) else Sum.inl ⟨i, h⟩)
    (fun v => Sum.inr (Sum.inr v))
  invFun := Sum.elim (fun a => Sum.inl a.1) (Sum.elim (fun a => Sum.inl a.1) (fun v => Sum.inr v))
  left_inv := by
    rintro (i | v)
    · by_cases h : D.Forced i <;> simp [h]
    · simp
  right_inv := by
    rintro (a | (a | v))
    · simp [a.2]
    · simp [a.2]
    · simp

@[simp] lemma bagEquiv_inl_forced {i : D.Idx} (h : D.Forced i) :
    D.bagEquiv (Sum.inl i) = Sum.inr (Sum.inl ⟨i, h⟩) := by
  simp [bagEquiv, h]

@[simp] lemma bagEquiv_inl_free {i : D.Idx} (h : ¬ D.Forced i) :
    D.bagEquiv (Sum.inl i) = Sum.inl ⟨i, h⟩ := by
  simp [bagEquiv, h]

@[simp] lemma bagEquiv_inr (v : D.Unm) : D.bagEquiv (Sum.inr v) = Sum.inr (Sum.inr v) := rfl

/-- The coordinate flips relating the cover code to the order-ideal code: a forced bag whose
forced value is `0` and a singleton bag are flipped. -/
noncomputable def bagFlip : D.Bag → Bool :=
  Sum.elim (fun i => if D.Forced i then !(D.forcedVal i) else false) (fun _ => true)

theorem card_Bag_eq_card_Free_add_numConst :
    Fintype.card D.Bag = Fintype.card D.Free + D.numConst := by
  rw [Fintype.card_congr D.bagEquiv, Fintype.card_sum, numConst]

/-! ### Equation (2) -/

variable {D}

/-- The word of the cover attached to a solution, written through the coordinate bijection
and the coordinate flips. -/
lemma coverWord_coverOf_eq {x : D.Idx → Bool} (hx : x ∈ D.Sol) :
    D.coverWord (D.coverOf x)
      = fun k => xor (D.bagFlip k) (Sum.elim (restrictFree x) (fun _ => true) (D.bagEquiv k)) := by
  funext k
  cases k with
  | inl i =>
      have hL : D.coverWord (D.coverOf x) (Sum.inl i) = x i := by
        simp [coverWord, Lv_mem_coverOf]
      rw [hL]
      by_cases h : D.Forced i
      · rw [bagEquiv_inl_forced D h]
        simp only [bagFlip, Sum.elim_inl, Sum.elim_inr, if_pos h]
        rw [forcedVal_spec h hx]
        cases D.forcedVal i <;> simp
      · rw [bagEquiv_inl_free D h]
        simp [bagFlip, h, restrictFree]
  | inr v =>
      have hR : D.coverWord (D.coverOf x) (Sum.inr v) = false := by
        simp [coverWord, unmatched_not_mem_coverOf v.2]
      rw [hR, bagEquiv_inr D]
      simp [bagFlip]

variable (D)

/-- **Equation (2)**: the code of the tree is the order-ideal code of the forest poset `P`
together with one constant word on the remaining `c` coordinates, read through the
coordinate bijection `bagEquiv` and the coordinate flips `bagFlip`. -/
theorem treeCode_eq_codeRelabel :
    D.treeCode = codeRelabel D.bagEquiv D.bagFlip (appendConsts (idealCode D.Free) D.Const) := by
  ext w
  constructor
  · intro hw
    obtain ⟨C, hC, rfl⟩ := Finset.mem_image.1 hw
    have hx : D.solOf C ∈ D.Sol := solOf_mem_Sol hC
    refine Finset.mem_image.2 ⟨Sum.elim (restrictFree (D.solOf C)) (fun _ => true), ?_, ?_⟩
    · exact Finset.mem_image.2 ⟨restrictFree (D.solOf C),
        by simpa [idealCode] using isIdealIndicator_restrictFree hx, rfl⟩
    · refine Eq.symm ?_
      conv_lhs => rw [← coverOf_solOf hC]
      exact coverWord_coverOf_eq hx
  · intro hw
    obtain ⟨u, hu, rfl⟩ := Finset.mem_image.1 hw
    obtain ⟨y, hy, rfl⟩ := Finset.mem_image.1 hu
    have hy' : IsIdealIndicator y := by simpa [idealCode] using hy
    have hx : D.extendFree y ∈ D.Sol := extendFree_mem_Sol hy'
    refine Finset.mem_image.2
      ⟨D.coverOf (D.extendFree y), coverOf_mem_minCovers (mem_Sol.1 hx), ?_⟩
    rw [coverWord_coverOf_eq hx, restrictFree_extendFree]

/-! ### Equation (6) -/

/-- **Equation (6)**, coefficientwise: the profile polynomial of the tree code is
`(1+t)^c I(B(P);t)`. -/
theorem codeP_treeCode (k : ℕ) :
    codeP D.treeCode k
      = ∑ j ∈ range (k + 1), codeP (idealCode D.Free) j * D.numConst.choose (k - j) := by
  rw [treeCode_eq_codeRelabel, codeP_codeRelabel, codeP_appendConsts, numConst]

/-- **The code of maximum independent sets** is the code of minimum covers with every
coordinate flipped; in particular the two have the same retained-coordinate profile. -/
theorem maxIndepCode_eq_flip_treeCode :
    D.maxIndepSets.image D.coverWord = codeRelabel (Equiv.refl D.Bag) (fun _ => true) D.treeCode :=
  by
  rw [maxIndepSets_eq_image_compl, treeCode, codeRelabel, Finset.image_image, Finset.image_image]
  refine Finset.image_congr ?_
  intro C _
  funext k
  simp [coverWord_compl]

theorem codeP_maxIndepCode (k : ℕ) :
    codeP (D.maxIndepSets.image D.coverWord) k = codeP D.treeCode k := by
  rw [maxIndepCode_eq_flip_treeCode, codeP_codeRelabel]

/-- **Equation (6)** in polynomial form: the profile polynomial of the tree code is
`(1+t)^c I(B(P);t)`. -/
theorem codePoly_treeCode :
    codePoly D.treeCode = ((1 : Polynomial ℤ) + X) ^ D.numConst * indepPoly D.Free := by
  have hM : Fintype.card D.Free + D.numConst = Fintype.card D.Bag :=
    (card_Bag_eq_card_Free_add_numConst D).symm
  refine Polynomial.ext fun m => ?_
  rw [coeff_pow_mul_indepPoly, codePoly, Polynomial.finset_sum_coeff]
  by_cases hm : m < Fintype.card D.Bag + 1
  · rw [Finset.sum_eq_single m]
    · have hco : (((codeP D.treeCode m : ℤ) : Polynomial ℤ) * X ^ m).coeff m
          = (codeP D.treeCode m : ℤ) := by simp
      rw [hco, codeP_treeCode]
      push_cast
      rfl
    · intro j _ hj
      simp [Polynomial.coeff_X_pow, Ne.symm hj]
    · intro h
      exact absurd (Finset.mem_range.2 hm) h
  · have hzero : ∀ j ∈ range (Fintype.card D.Bag + 1),
        (((codeP D.treeCode j : ℤ) : Polynomial ℤ) * X ^ j).coeff m = 0 := by
      intro j hj
      rw [Finset.mem_range] at hj
      have : m ≠ j := by omega
      simp [Polynomial.coeff_X_pow, this]
    rw [Finset.sum_eq_zero hzero]
    refine (Finset.sum_eq_zero ?_).symm
    intro j hj
    rw [Finset.mem_range] at hj
    by_cases hjr : Fintype.card D.Free < j
    · rw [codeP_eq_zero_of_lt hjr]
      simp
    · have : D.numConst < m - j := by omega
      rw [Nat.choose_eq_zero_of_lt this]
      simp

/-- The erasure profile of the tree code: `e_d = p_{M-d}` with `M = α(T)` the number of
bags. -/
noncomputable def erasure (d : ℕ) : ℕ := codeP D.treeCode (Fintype.card D.Bag - d)

/-- The erasure profile of the tree code is the poset erasure profile of `PascalBridge`. -/
theorem erasure_eq_erasureProfile (d : ℕ) (hd : d ≤ Fintype.card D.Bag) :
    D.erasure d = erasureProfile D.Free D.numConst d := by
  have hM : Fintype.card D.Free + D.numConst = Fintype.card D.Bag :=
    (card_Bag_eq_card_Free_add_numConst D).symm
  have hcoeff := erasureProfile_eq_coeff (P := D.Free) D.numConst (d := d) (by omega)
  rw [coeff_pow_mul_indepPoly, hM] at hcoeff
  have hsum : ((D.erasure d : ℕ) : ℤ)
      = ∑ j ∈ range (Fintype.card D.Bag - d + 1),
          (codeP (idealCode D.Free) j : ℤ) *
            (D.numConst.choose (Fintype.card D.Bag - d - j) : ℤ) := by
    rw [erasure, codeP_treeCode]
    push_cast
    rfl
  exact_mod_cast hsum.trans hcoeff

/-! ### The inequalities of §3 for the tree code -/

/-- Log-concavity of the erasure profile of a tree through defect depth eight. -/
theorem erasure_log_concave_depth_le_eight {d : ℕ} (hd1 : 1 ≤ d)
    (hd2 : d ≤ Fintype.card D.Bag - 1) (hd8 : d ≤ 8) :
    D.erasure (d - 1) * D.erasure (d + 1) ≤ D.erasure d ^ 2 := by
  have hM : Fintype.card D.Free + D.numConst = Fintype.card D.Bag :=
    (card_Bag_eq_card_Free_add_numConst D).symm
  have hd1' : d + 1 ≤ Fintype.card D.Bag := by omega
  rw [erasure_eq_erasureProfile D (d - 1) (by omega), erasure_eq_erasureProfile D d (by omega),
    erasure_eq_erasureProfile D (d + 1) hd1']
  exact erasureProfile_log_concave_depth_le_eight (P := D.Free) D.numConst hd1 (by omega) hd8

/-- Log-concavity of the erasure profile of a tree at every interior defect, when the
independence number is at most 33. -/
theorem erasure_log_concave_of_le_33 (hM33 : Fintype.card D.Bag ≤ 33) {d : ℕ} (hd1 : 1 ≤ d)
    (hd2 : d ≤ Fintype.card D.Bag - 1) :
    D.erasure (d - 1) * D.erasure (d + 1) ≤ D.erasure d ^ 2 := by
  have hM : Fintype.card D.Free + D.numConst = Fintype.card D.Bag :=
    (card_Bag_eq_card_Free_add_numConst D).symm
  have hd1' : d + 1 ≤ Fintype.card D.Bag := by omega
  rw [erasure_eq_erasureProfile D (d - 1) (by omega), erasure_eq_erasureProfile D d (by omega),
    erasure_eq_erasureProfile D (d + 1) hd1']
  exact erasureProfile_log_concave_of_M_le_33 (P := D.Free) D.numConst (by omega) hd1 (by omega)

/-- **Equation (7)** for a tree: `e_3^2 ≥ e_2 e_4`. -/
theorem erasure_depth_three (hM4 : 4 ≤ Fintype.card D.Bag) :
    D.erasure 2 * D.erasure 4 ≤ D.erasure 3 ^ 2 := by
  have hM : Fintype.card D.Free + D.numConst = Fintype.card D.Bag :=
    (card_Bag_eq_card_Free_add_numConst D).symm
  rw [erasure_eq_erasureProfile D 2 (by omega), erasure_eq_erasureProfile D 3 (by omega),
    erasure_eq_erasureProfile D 4 (by omega)]
  exact erasureProfile_depth_three_log_concave (P := D.Free) D.numConst (by omega)

/-- **Equation (7a)** for a tree: the quantitative depth-three reserve. -/
theorem erasure_depth_three_reserve (hM4 : 4 ≤ Fintype.card D.Bag) :
    32 * (Fintype.card D.Bag - 2) * (D.erasure 2 * D.erasure 4)
      ≤ 27 * (Fintype.card D.Bag - 3) * D.erasure 3 ^ 2 := by
  have hM : Fintype.card D.Free + D.numConst = Fintype.card D.Bag :=
    (card_Bag_eq_card_Free_add_numConst D).symm
  rw [erasure_eq_erasureProfile D 2 (by omega), erasure_eq_erasureProfile D 3 (by omega),
    erasure_eq_erasureProfile D 4 (by omega), ← hM]
  exact erasureProfile_depth_three (P := D.Free) D.numConst (by omega)

end TreeMatching

end MatchingBag
