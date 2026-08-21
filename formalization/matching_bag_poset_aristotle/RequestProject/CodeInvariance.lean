import RequestProject.Codes

/-!
# Relabelling coordinates, flipping coordinates, and adjoining constant coordinates

The code produced by a tree with a maximum matching lives on the bags, and matches the
order-ideal code of the forest poset only after renaming the coordinates and flipping some of
them ("harmlessly flipping coordinates" in §1 of
`matching_bag_poset_reduction_2026-08-20.md`).  This file records that the
retained-coordinate profile `p_k` is invariant under such an isomorphism of codes, and
computes the profile of a code with several constant coordinates adjoined.

* `MatchingBag.codeRelabel σ s Ω`: the code `Ω` read through the coordinate bijection `σ`,
  with the coordinates in `s` flipped.
* `MatchingBag.codeP_codeRelabel`: `p_k` is invariant.
* `MatchingBag.appendConsts Ω κ`: `Ω` with one constant coordinate for each element of `κ`.
* `MatchingBag.codeP_appendConsts`: `p_k(Ω × const) = ∑_j p_j(Ω) C(|κ|, k - j)`, the
  coefficientwise form of multiplication of the profile polynomial by `(1+t)^{|κ|}`.
-/

open scoped BigOperators
open Finset

namespace MatchingBag

variable {ι κ : Type*} [Fintype ι] [DecidableEq ι] [Fintype κ] [DecidableEq κ]

/-! ### Relabelling and flipping coordinates -/

/-- The code `Ω` read through the coordinate bijection `σ : κ ≃ ι`, with the coordinates
where `s` is `true` flipped. -/
def codeRelabel (σ : κ ≃ ι) (s : κ → Bool) (Ω : Finset (ι → Bool)) : Finset (κ → Bool) :=
  Ω.image (fun f k => xor (s k) (f (σ k)))

/-- The retained-coordinate profile is invariant under renaming and flipping coordinates. -/
theorem codeP_codeRelabel (σ : κ ≃ ι) (s : κ → Bool) (Ω : Finset (ι → Bool)) (k : ℕ) :
    codeP (codeRelabel σ s Ω) k = codeP Ω k := by
  classical
  have key : ∀ K' : Finset κ, (codeProj K' (codeRelabel σ s Ω)).card
      = (codeProj (K'.image σ) Ω).card := by
    intro K'
    set ψ : PartialWord ι → PartialWord κ :=
      fun w x => (w (σ x)).map (fun b => xor (s x) b) with hψ
    have hψinj : Function.Injective ψ := by
      intro w1 w2 h
      funext i
      have hi := congrFun h (σ.symm i)
      simp only [hψ, Equiv.apply_symm_apply] at hi
      have hb : Function.Injective (fun b => xor (s (σ.symm i)) b) := by
        intro a b hab; cases s (σ.symm i) <;> simpa using hab
      exact Option.map_injective hb hi
    have himg : codeProj K' (codeRelabel σ s Ω) = (codeProj (K'.image σ) Ω).image ψ := by
      rw [codeProj, codeRelabel, codeProj, Finset.image_image, Finset.image_image]
      refine Finset.image_congr ?_
      intro f _
      funext x
      by_cases hx : x ∈ K'
      · simp [restrictTo, hψ, hx, Finset.mem_image.2 ⟨x, hx, rfl⟩]
      · have hnot : σ x ∉ K'.image σ := by
          simp only [Finset.mem_image, not_exists]
          rintro y ⟨hy, hxy⟩
          exact hx (by rwa [σ.injective hxy] at hy)
        simp [restrictTo, hψ, hx, hnot]
    rw [himg, Finset.card_image_of_injective _ hψinj]
  rw [codeP, codeP]
  refine Finset.sum_nbij' (i := fun K' => K'.image σ) (j := fun K => K.image σ.symm) ?_ ?_ ?_ ?_ ?_
  · intro K' hK'
    rw [Finset.mem_powersetCard] at hK' ⊢
    exact ⟨Finset.subset_univ _, by rw [Finset.card_image_of_injective _ σ.injective, hK'.2]⟩
  · intro K hK
    rw [Finset.mem_powersetCard] at hK ⊢
    exact ⟨Finset.subset_univ _, by rw [Finset.card_image_of_injective _ σ.symm.injective, hK.2]⟩
  · intro K' _
    simp [Finset.image_image]
  · intro K _
    simp [Finset.image_image]
  · intro K' _
    exact key K'

/-! ### Constant coordinates -/

/-- The code `Ω` with one extra coordinate for each element of `κ`, on which every word of
the code is `1`. -/
def appendConsts (Ω : Finset (ι → Bool)) (κ : Type*) [Fintype κ] [DecidableEq κ] :
    Finset ((ι ⊕ κ) → Bool) :=
  Ω.image (fun f => Sum.elim f (fun _ => true))

/-- The profile vanishes above the number of coordinates. -/
theorem codeP_eq_zero_of_lt {Ω : Finset (ι → Bool)} {k : ℕ} (hk : Fintype.card ι < k) :
    codeP Ω k = 0 := by
  rw [codeP, Finset.powersetCard_eq_empty.2 (by simpa using hk), Finset.sum_empty]

/-- Projecting a code with constant coordinates adjoined only sees the original
coordinates. -/
theorem card_codeProj_appendConsts (Ω : Finset (ι → Bool)) (K : Finset (ι ⊕ κ)) :
    (codeProj K (appendConsts Ω κ)).card = (codeProj K.toLeft Ω).card := by
  classical
  set Φ : PartialWord ι → PartialWord (ι ⊕ κ) :=
    fun w x => Sum.elim w (fun y => if Sum.inr y ∈ K then some true else none) x with hΦ
  have hinj : Function.Injective Φ := by
    intro w1 w2 h
    funext i
    exact congrFun h (Sum.inl i)
  have himg : codeProj K (appendConsts Ω κ) = (codeProj K.toLeft Ω).image Φ := by
    rw [codeProj, appendConsts, codeProj, Finset.image_image, Finset.image_image]
    refine Finset.image_congr ?_
    intro f _
    funext x
    cases x with
    | inl i => by_cases h : Sum.inl i ∈ K <;> simp [restrictTo, hΦ, h, Finset.mem_toLeft]
    | inr j => by_cases h : Sum.inr j ∈ K <;> simp [restrictTo, hΦ, h]
  rw [himg, Finset.card_image_of_injective _ hinj]

omit [DecidableEq ι] [DecidableEq κ] in
/-- Splitting a sum over the `k`-subsets of `ι ⊕ κ` according to the trace on `ι`. -/
theorem sum_powersetCard_sum (f : Finset ι → ℕ) (k : ℕ) :
    ∑ K ∈ (univ : Finset (ι ⊕ κ)).powersetCard k, f K.toLeft
      = ∑ j ∈ range (k + 1), ∑ A ∈ (univ : Finset ι).powersetCard j,
          ∑ _B ∈ (univ : Finset κ).powersetCard (k - j), f A := by
  classical
  have h1 : ∑ j ∈ range (k + 1), ∑ A ∈ (univ : Finset ι).powersetCard j,
        ∑ _B ∈ (univ : Finset κ).powersetCard (k - j), f A
      = ∑ x ∈ (range (k + 1)).sigma (fun j => ((univ : Finset ι).powersetCard j) ×ˢ
          ((univ : Finset κ).powersetCard (k - j))), f x.2.1 := by
    rw [Finset.sum_sigma]
    refine Finset.sum_congr rfl fun j _ => ?_
    rw [Finset.sum_product]
  rw [h1]
  refine Finset.sum_nbij'
    (i := fun K => (⟨K.toLeft.card, (K.toLeft, K.toRight)⟩ : (_ : ℕ) × (Finset ι × Finset κ)))
    (j := fun x => x.2.1.disjSum x.2.2) ?_ ?_ ?_ ?_ ?_
  · intro K hK
    rw [Finset.mem_powersetCard] at hK
    have hsum := Finset.card_toLeft_add_card_toRight (u := K)
    simp only [Finset.mem_sigma, Finset.mem_product, Finset.mem_powersetCard, Finset.mem_range]
    refine ⟨by omega, ⟨Finset.subset_univ _, ?_⟩, Finset.subset_univ _, by omega⟩
    trivial
  · rintro ⟨j, A, B⟩ hx
    simp only [Finset.mem_sigma, Finset.mem_product, Finset.mem_powersetCard,
      Finset.mem_range] at hx
    rw [Finset.mem_powersetCard]
    refine ⟨Finset.subset_univ _, ?_⟩
    rw [Finset.card_disjSum, hx.2.1.2, hx.2.2.2]
    omega
  · intro K _
    simp [Finset.toLeft_disjSum_toRight]
  · rintro ⟨j, A, B⟩ hx
    simp only [Finset.mem_sigma, Finset.mem_product, Finset.mem_powersetCard,
      Finset.mem_range] at hx
    simp [hx.2.1.2]
  · intro K _
    rfl

/-- Adjoining `|κ|` constant coordinates multiplies the profile polynomial by
`(1+t)^{|κ|}`. -/
theorem codeP_appendConsts (Ω : Finset (ι → Bool)) (k : ℕ) :
    codeP (appendConsts Ω κ) k
      = ∑ j ∈ range (k + 1), codeP Ω j * (Fintype.card κ).choose (k - j) := by
  classical
  rw [codeP, Finset.sum_congr rfl (fun K _ => card_codeProj_appendConsts Ω K),
    sum_powersetCard_sum (fun A => (codeProj A Ω).card) k]
  refine Finset.sum_congr rfl fun j _ => ?_
  rw [codeP, Finset.sum_mul]
  refine Finset.sum_congr rfl fun A _ => ?_
  rw [Finset.sum_const, Finset.card_powersetCard, Finset.card_univ, smul_eq_mul, mul_comm]

end MatchingBag
