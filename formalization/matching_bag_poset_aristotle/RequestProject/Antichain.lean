import RequestProject.PosetCode

/-!
# The antichain form of the independence polynomial of `B(P)`

This file formalizes Section 3 of `matching_bag_poset_reduction_2026-08-20.md`:
formulas (8), (9) and (10).

* `MatchingBag.indepPoly P = I(B(P); t) = ∑_{S independent} t^{|S|}`.
* `indepPoly_eq_sum_codeP` is formula (5) in polynomial form.
* `indepPoly_eq_sum_downClosure` is (8):
  `I(B(P);t) = ∑_{O ⊆ P} t^{|O|} (1+t)^{r - |↓O|}`.
* `indepPoly_eq_sum_antichains` is (9):
  `I(B(P);t) = ∑_{A antichain} t^{|A|} (1+t)^{r - |A|} = (1+t)^r A_P(t/(1+t))`.
* `coeff_erasure` is (10): `E_d = ∑_j a_j (M - j).choose d`.
-/

open scoped BigOperators
open Finset Polynomial

namespace MatchingBag

variable {P : Type*} [Fintype P] [DecidableEq P] [PartialOrder P]
  [DecidableRel ((· ≤ ·) : P → P → Prop)]

/-- The down-closure `↓O = {x | ∃ i ∈ O, x ≤ i}` of a subset of the poset. -/
def downClosure (O : Finset P) : Finset P := Finset.univ.filter (fun x => ∃ i ∈ O, x ≤ i)

/-- The antichains of `P`. -/
def antichainsOf (P : Type*) [Fintype P] [DecidableEq P] [PartialOrder P]
    [DecidableRel ((· ≤ ·) : P → P → Prop)] : Finset (Finset P) :=
  (Finset.univ : Finset (Finset P)).filter (fun A => ∀ i ∈ A, ∀ j ∈ A, i ≤ j → i = j)

/-- `a_j`, the number of antichains of `P` with `j` elements; the coefficients of the
antichain polynomial `A_P(u) = ∑_j a_j u^j`. -/
def antichainCount (P : Type*) [Fintype P] [DecidableEq P] [PartialOrder P]
    [DecidableRel ((· ≤ ·) : P → P → Prop)] (j : ℕ) : ℕ :=
  ((antichainsOf P).filter (fun A => A.card = j)).card

/-- The independence polynomial `I(B(P); t) = ∑_{S independent} t^{|S|}`. -/
noncomputable def indepPoly (P : Type*) [Fintype P] [DecidableEq P] [PartialOrder P]
    [DecidableRel ((· ≤ ·) : P → P → Prop)] : Polynomial ℤ :=
  ∑ S ∈ indepSets P, X ^ S.card

omit [DecidableEq P] in
lemma mem_downClosure {O : Finset P} {x : P} : x ∈ downClosure O ↔ ∃ i ∈ O, x ≤ i := by
  simp [downClosure]

omit [DecidableEq P] in
lemma subset_downClosure (O : Finset P) : O ⊆ downClosure O :=
  fun x hx => mem_downClosure.2 ⟨x, hx, le_refl x⟩

omit [Fintype P] [PartialOrder P] [DecidableRel ((· ≤ ·) : P → P → Prop)] in
lemma sum_powerset_pow_card (s : Finset P) :
    ∑ Z ∈ s.powerset, (X : Polynomial ℤ) ^ Z.card = (1 + X) ^ s.card := by
  have h := Finset.prod_add (fun _ : P => (X : Polynomial ℤ)) (fun _ : P => (1 : Polynomial ℤ)) s
  simpa [Finset.prod_const, add_comm] using h.symm

/-- **Formula (5)**, polynomial form: `∑_k p_k(Ω(P)) t^k = I(B(P); t)`. -/
theorem indepPoly_eq_sum_codeP :
    indepPoly P = ∑ k ∈ range (2 * Fintype.card P + 1),
      (codeP (idealCode P) k : ℤ) * (X : Polynomial ℤ) ^ k := by
  have hmaps : ∀ S ∈ indepSets P, S.card ∈ range (2 * Fintype.card P + 1) := by
    intro S _
    rw [Finset.mem_range]
    have h1 : S.card ≤ (Finset.univ : Finset (P × Bool)).card := Finset.card_le_univ S
    rw [Finset.card_univ, Fintype.card_prod, Fintype.card_bool] at h1
    omega
  rw [indepPoly, ← Finset.sum_fiberwise_of_maps_to hmaps
    (f := fun S : Finset (P × Bool) => (X : Polynomial ℤ) ^ S.card)]
  refine Finset.sum_congr rfl fun k _ => ?_
  have hpow : ∀ S ∈ (indepSets P).filter (fun S => S.card = k),
      (X : Polynomial ℤ) ^ S.card = X ^ k := by
    intro S hS
    rw [(Finset.mem_filter.1 hS).2]
  rw [Finset.sum_congr rfl hpow, Finset.sum_const, nsmul_eq_mul, codeP_eq_indepSets_card,
    indepSetsCard]
  norm_cast

/-- The set of coordinates a partial assignment sends to `1`. -/
def onesOf (S : Finset (P × Bool)) : Finset P := (S.filter (fun a => a.2 = true)).image Prod.fst

/-- The set of coordinates a partial assignment sends to `0`. -/
def zerosOf (S : Finset (P × Bool)) : Finset P := (S.filter (fun a => a.2 = false)).image Prod.fst

/-- The partial assignment with `1`s on `O` and `0`s on `Z`, as a vertex set of `B(P)`. -/
def pairSet (O Z : Finset P) : Finset (P × Bool) :=
  O.image (fun i => (i, true)) ∪ Z.image (fun i => (i, false))

section
omit [Fintype P] [PartialOrder P] [DecidableRel ((· ≤ ·) : P → P → Prop)]

lemma mem_onesOf {S : Finset (P × Bool)} {i : P} : i ∈ onesOf S ↔ (i, true) ∈ S := by
  simp only [onesOf, Finset.mem_image, Finset.mem_filter]
  constructor
  · rintro ⟨a, ⟨ha, h2⟩, rfl⟩
    rwa [← h2]
  · intro h
    exact ⟨(i, true), ⟨h, rfl⟩, rfl⟩

lemma mem_zerosOf {S : Finset (P × Bool)} {i : P} : i ∈ zerosOf S ↔ (i, false) ∈ S := by
  simp only [zerosOf, Finset.mem_image, Finset.mem_filter]
  constructor
  · rintro ⟨a, ⟨ha, h2⟩, rfl⟩
    rwa [← h2]
  · intro h
    exact ⟨(i, false), ⟨h, rfl⟩, rfl⟩

lemma mem_pairSet {O Z : Finset P} {a : P × Bool} :
    a ∈ pairSet O Z ↔ (a.2 = true ∧ a.1 ∈ O) ∨ (a.2 = false ∧ a.1 ∈ Z) := by
  simp only [pairSet, Finset.mem_union, Finset.mem_image]
  constructor
  · rintro (⟨i, hi, rfl⟩ | ⟨i, hi, rfl⟩)
    · exact Or.inl ⟨rfl, hi⟩
    · exact Or.inr ⟨rfl, hi⟩
  · rintro (⟨h2, h1⟩ | ⟨h2, h1⟩)
    · exact Or.inl ⟨a.1, h1, Prod.ext_iff.2 ⟨rfl, h2.symm⟩⟩
    · exact Or.inr ⟨a.1, h1, Prod.ext_iff.2 ⟨rfl, h2.symm⟩⟩

lemma onesOf_pairSet (O Z : Finset P) : onesOf (pairSet O Z) = O := by
  ext i
  rw [mem_onesOf, mem_pairSet]
  simp

lemma zerosOf_pairSet (O Z : Finset P) : zerosOf (pairSet O Z) = Z := by
  ext i
  rw [mem_zerosOf, mem_pairSet]
  simp

lemma pairSet_onesOf_zerosOf (S : Finset (P × Bool)) : pairSet (onesOf S) (zerosOf S) = S := by
  ext a
  rw [mem_pairSet, mem_onesOf, mem_zerosOf]
  constructor
  · rintro (⟨h2, h1⟩ | ⟨h2, h1⟩) <;> rwa [← h2] at h1
  · intro ha
    cases hA : a.2 with
    | true => exact Or.inl ⟨rfl, by rwa [← hA]⟩
    | false => exact Or.inr ⟨rfl, by rwa [← hA]⟩

lemma card_pairSet (O Z : Finset P) : (pairSet O Z).card = O.card + Z.card := by
  have hdisj : Disjoint (O.image (fun i => (i, true))) (Z.image (fun i : P => (i, false))) := by
    rw [Finset.disjoint_left]
    rintro a ha hb
    obtain ⟨i, -, rfl⟩ := Finset.mem_image.1 ha
    obtain ⟨j, -, hj⟩ := Finset.mem_image.1 hb
    exact Bool.noConfusion (congrArg Prod.snd hj)
  rw [pairSet, Finset.card_union_of_disjoint hdisj,
    Finset.card_image_of_injective _ (fun i j hij => congrArg Prod.fst hij),
    Finset.card_image_of_injective _ (fun i j hij => congrArg Prod.fst hij)]

end

lemma pairSet_mem_indepSets_iff (O Z : Finset P) :
    pairSet O Z ∈ indepSets P ↔ Z ⊆ Finset.univ \ downClosure O := by
  rw [indepSets, Finset.mem_filter]
  simp only [Finset.mem_univ, true_and]
  constructor
  · intro h j hj
    rw [Finset.mem_sdiff]
    refine ⟨Finset.mem_univ _, ?_⟩
    intro hjd
    obtain ⟨i, hi, hij⟩ := mem_downClosure.1 hjd
    exact h (i, true) (mem_pairSet.2 (Or.inl ⟨rfl, hi⟩)) (j, false)
      (mem_pairSet.2 (Or.inr ⟨rfl, hj⟩)) (Or.inl ⟨rfl, rfl, hij⟩)
  · rintro h a ha b hb hab
    rcases hab with ⟨h1, h2, h3⟩ | ⟨h1, h2, h3⟩
    · rcases mem_pairSet.1 ha with ⟨-, ha1⟩ | ⟨ha2, -⟩
      · rcases mem_pairSet.1 hb with ⟨hb2, -⟩ | ⟨-, hb1⟩
        · rw [h2] at hb2; exact Bool.noConfusion hb2
        · exact (Finset.mem_sdiff.1 (h hb1)).2 (mem_downClosure.2 ⟨a.1, ha1, h3⟩)
      · rw [h1] at ha2; exact Bool.noConfusion ha2
    · rcases mem_pairSet.1 hb with ⟨-, hb1⟩ | ⟨hb2, -⟩
      · rcases mem_pairSet.1 ha with ⟨ha2, -⟩ | ⟨-, ha1⟩
        · rw [h1] at ha2; exact Bool.noConfusion ha2
        · exact (Finset.mem_sdiff.1 (h ha1)).2 (mem_downClosure.2 ⟨b.1, hb1, h3⟩)
      · rw [h2] at hb2; exact Bool.noConfusion hb2

/-- The independent sets of `B(P)` with a prescribed set `O` of `1`-coordinates are exactly
the subsets of `P \ ↓O` used as `0`-coordinates. -/
lemma sum_indepSets_fiber_onesOf (O : Finset P) :
    ∑ S ∈ (indepSets P).filter (fun S => onesOf S = O), (X : Polynomial ℤ) ^ S.card
      = ∑ Z ∈ (Finset.univ \ downClosure O).powerset, (X : Polynomial ℤ) ^ (O.card + Z.card) := by
  refine Finset.sum_nbij' (i := zerosOf) (j := pairSet O) ?_ ?_ ?_ ?_ ?_
  · intro S hS
    rw [Finset.mem_filter] at hS
    rw [Finset.mem_powerset, ← hS.2, ← pairSet_onesOf_zerosOf S]
    rw [onesOf_pairSet, zerosOf_pairSet]
    exact (pairSet_mem_indepSets_iff (onesOf S) (zerosOf S)).1
      (by rw [pairSet_onesOf_zerosOf]; exact hS.1)
  · intro Z hZ
    rw [Finset.mem_powerset] at hZ
    rw [Finset.mem_filter]
    exact ⟨(pairSet_mem_indepSets_iff O Z).2 hZ, onesOf_pairSet O Z⟩
  · intro S hS
    rw [Finset.mem_filter] at hS
    rw [← hS.2, pairSet_onesOf_zerosOf]
  · intro Z _
    exact zerosOf_pairSet O Z
  · intro S hS
    rw [Finset.mem_filter] at hS
    conv_lhs => rw [← pairSet_onesOf_zerosOf S]
    rw [card_pairSet, hS.2]

/-- **Formula (8)**: grouping independent sets of `B(P)` by their set `O` of coordinates
assigned `1`, the coordinates assigned `0` may be any subset of `P \ ↓O`. -/
theorem indepPoly_eq_sum_downClosure :
    indepPoly P = ∑ O ∈ (Finset.univ : Finset P).powerset,
      (X : Polynomial ℤ) ^ O.card * (1 + X) ^ (Fintype.card P - (downClosure O).card) := by
  have hmaps : ∀ S ∈ indepSets P, onesOf S ∈ (Finset.univ : Finset P).powerset := by
    intro S _
    simp
  rw [indepPoly, ← Finset.sum_fiberwise_of_maps_to hmaps
    (f := fun S : Finset (P × Bool) => (X : Polynomial ℤ) ^ S.card)]
  refine Finset.sum_congr rfl fun O _ => ?_
  rw [sum_indepSets_fiber_onesOf O,
    Finset.sum_congr rfl (fun Z _ => pow_add (X : Polynomial ℤ) O.card Z.card), ← Finset.mul_sum,
    sum_powerset_pow_card, Finset.card_sdiff_of_subset (Finset.subset_univ _), Finset.card_univ]

/-- The set of maximal elements of `O`. -/
def maxOf (O : Finset P) : Finset P := O.filter (fun i => ∀ j ∈ O, i ≤ j → j = i)

lemma maxOf_subset (O : Finset P) : maxOf O ⊆ O := Finset.filter_subset _ _

lemma exists_maximal_ge {O : Finset P} {x : P} (hx : x ∈ O) : ∃ m ∈ maxOf O, x ≤ m := by
  have : WellFoundedGT P := Finite.to_wellFoundedGT
  obtain ⟨b, hab, hb⟩ := exists_maximal_ge_of_wellFoundedGT (fun y => y ∈ O) x hx
  exact ⟨b, Finset.mem_filter.2 ⟨hb.1, fun j hj hbj => le_antisymm (hb.2 hj hbj) hbj⟩, hab⟩

lemma maxOf_mem_antichains (O : Finset P) : maxOf O ∈ antichainsOf P := by
  rw [antichainsOf, Finset.mem_filter]
  refine ⟨Finset.mem_univ _, ?_⟩
  intro i hi j hj hij
  exact ((Finset.mem_filter.1 hi).2 j (maxOf_subset O hj) hij).symm

lemma downClosure_maxOf (O : Finset P) : downClosure (maxOf O) = downClosure O := by
  ext x
  rw [mem_downClosure, mem_downClosure]
  constructor
  · rintro ⟨i, hi, hxi⟩
    exact ⟨i, maxOf_subset O hi, hxi⟩
  · rintro ⟨i, hi, hxi⟩
    obtain ⟨m, hm, him⟩ := exists_maximal_ge hi
    exact ⟨m, hm, le_trans hxi him⟩

lemma maxOf_eq_iff {A O : Finset P} (hA : A ∈ antichainsOf P) :
    maxOf O = A ↔ (A ⊆ O ∧ O ⊆ downClosure A) := by
  rw [antichainsOf, Finset.mem_filter] at hA
  constructor
  · rintro rfl
    refine ⟨maxOf_subset O, ?_⟩
    intro x hx
    obtain ⟨m, hm, hxm⟩ := exists_maximal_ge hx
    exact mem_downClosure.2 ⟨m, hm, hxm⟩
  · rintro ⟨hAO, hOA⟩
    apply Finset.Subset.antisymm
    · intro m hm
      have hmO : m ∈ O := maxOf_subset O hm
      obtain ⟨a, haA, hma⟩ := mem_downClosure.1 (hOA hmO)
      have hmam : a = m := (Finset.mem_filter.1 hm).2 a (hAO haA) hma
      exact hmam ▸ haA
    · intro a haA
      refine Finset.mem_filter.2 ⟨hAO haA, ?_⟩
      intro j hjO haj
      obtain ⟨a', ha'A, hja'⟩ := mem_downClosure.1 (hOA hjO)
      have haa' : a = a' := hA.2 a haA a' ha'A (le_trans haj hja')
      exact le_antisymm (by rw [haa']; exact hja') haj

/-- **Formula (9)**: `I(B(P);t) = (1+t)^r A_P(t/(1+t))`, written without denominators. -/
theorem indepPoly_eq_sum_antichains :
    indepPoly P = ∑ A ∈ antichainsOf P,
      (X : Polynomial ℤ) ^ A.card * (1 + X) ^ (Fintype.card P - A.card) := by
  rw [indepPoly_eq_sum_downClosure]
  have hmaps : ∀ O ∈ (Finset.univ : Finset P).powerset, maxOf O ∈ antichainsOf P :=
    fun O _ => maxOf_mem_antichains O
  rw [← Finset.sum_fiberwise_of_maps_to hmaps
    (f := fun O : Finset P => (X : Polynomial ℤ) ^ O.card *
      (1 + X) ^ (Fintype.card P - (downClosure O).card))]
  refine Finset.sum_congr rfl fun A hA => ?_
  have hAd : A ⊆ downClosure A := subset_downClosure A
  have hcard1 : A.card ≤ (downClosure A).card := Finset.card_le_card hAd
  have hcard2 : (downClosure A).card ≤ Fintype.card P := by
    rw [← Finset.card_univ]
    exact Finset.card_le_card (Finset.subset_univ _)
  have step : ∑ O ∈ ((Finset.univ : Finset P).powerset).filter (fun O => maxOf O = A),
      (X : Polynomial ℤ) ^ O.card * (1 + X) ^ (Fintype.card P - (downClosure O).card)
      = ∑ Z ∈ (downClosure A \ A).powerset, (X : Polynomial ℤ) ^ (A.card + Z.card) *
        (1 + X) ^ (Fintype.card P - (downClosure A).card) := by
    refine Finset.sum_nbij' (i := fun O => O \ A) (j := fun Z => A ∪ Z) ?_ ?_ ?_ ?_ ?_
    · intro O hO
      rw [Finset.mem_filter] at hO
      obtain ⟨hAO, hOA⟩ := (maxOf_eq_iff hA).1 hO.2
      rw [Finset.mem_powerset]
      intro x hx
      rw [Finset.mem_sdiff] at hx ⊢
      exact ⟨hOA hx.1, hx.2⟩
    · intro Z hZ
      rw [Finset.mem_powerset] at hZ
      rw [Finset.mem_filter, Finset.mem_powerset]
      refine ⟨Finset.subset_univ _, (maxOf_eq_iff hA).2 ⟨Finset.subset_union_left, ?_⟩⟩
      intro x hx
      rcases Finset.mem_union.1 hx with h | h
      · exact hAd h
      · exact (Finset.mem_sdiff.1 (hZ h)).1
    · intro O hO
      rw [Finset.mem_filter] at hO
      obtain ⟨hAO, -⟩ := (maxOf_eq_iff hA).1 hO.2
      exact Finset.union_sdiff_of_subset hAO
    · intro Z hZ
      rw [Finset.mem_powerset] at hZ
      ext x
      simp only [Finset.mem_sdiff, Finset.mem_union]
      constructor
      · rintro ⟨h | h, hx⟩
        · exact absurd h hx
        · exact h
      · intro hx
        exact ⟨Or.inr hx, (Finset.mem_sdiff.1 (hZ hx)).2⟩
    · intro O hO
      rw [Finset.mem_filter] at hO
      obtain ⟨hAO, -⟩ := (maxOf_eq_iff hA).1 hO.2
      have hcardO : O.card = A.card + (O \ A).card := by
        rw [← Finset.card_union_of_disjoint (Finset.disjoint_sdiff),
          Finset.union_sdiff_of_subset hAO]
      have hdO : downClosure O = downClosure A := by
        rw [← hO.2, downClosure_maxOf]
      rw [hcardO, hdO]
  have hexp : (downClosure A \ A).card + (Fintype.card P - (downClosure A).card)
      = Fintype.card P - A.card := by
    rw [Finset.card_sdiff_of_subset hAd]
    omega
  rw [step, Finset.sum_congr rfl
      (fun Z _ => by rw [pow_add (X : Polynomial ℤ) A.card Z.card]),
    ← Finset.sum_mul, ← Finset.mul_sum, sum_powerset_pow_card, mul_assoc, ← pow_add, hexp]

/-- **Formula (9)**, collected by the size of the antichain:
`I(B(P);t) = ∑_j a_j t^j (1+t)^{r-j}`. -/
theorem indepPoly_eq_antichainCount :
    indepPoly P = ∑ j ∈ range (Fintype.card P + 1),
      Polynomial.C (antichainCount P j : ℤ) * X ^ j * (1 + X) ^ (Fintype.card P - j) := by
  rw [indepPoly_eq_sum_antichains]
  have hmaps : ∀ A ∈ antichainsOf P, A.card ∈ range (Fintype.card P + 1) := by
    intro A _
    rw [Finset.mem_range]
    have : A.card ≤ (Finset.univ : Finset P).card := Finset.card_le_univ A
    rw [Finset.card_univ] at this
    omega
  rw [← Finset.sum_fiberwise_of_maps_to hmaps
    (f := fun A : Finset P => (X : Polynomial ℤ) ^ A.card * (1 + X) ^ (Fintype.card P - A.card))]
  refine Finset.sum_congr rfl fun j _ => ?_
  have hterm : ∀ A ∈ (antichainsOf P).filter (fun A => A.card = j),
      (X : Polynomial ℤ) ^ A.card * (1 + X) ^ (Fintype.card P - A.card)
        = X ^ j * (1 + X) ^ (Fintype.card P - j) := by
    intro A hA
    rw [(Finset.mem_filter.1 hA).2]
  rw [Finset.sum_congr rfl hterm, Finset.sum_const, nsmul_eq_mul, antichainCount, mul_assoc,
    Polynomial.C_eq_natCast]

/-- **Formula (10)**: with `M = r + c` and `E_d = [t^{M-d}] (1+t)^c I(B(P);t)`, one has
`E_d = ∑_j a_j (M - j).choose d`. -/
theorem coeff_erasure (c d : ℕ) (hd : d ≤ Fintype.card P + c) :
    (((1 : Polynomial ℤ) + X) ^ c * indepPoly P).coeff (Fintype.card P + c - d) =
      ∑ j ∈ range (Fintype.card P + 1),
        (antichainCount P j : ℤ) * ((Fintype.card P + c - j).choose d : ℤ) := by
  rw [indepPoly_eq_antichainCount, Finset.mul_sum, Polynomial.finset_sum_coeff]
  refine Finset.sum_congr rfl fun j hj => ?_
  rw [Finset.mem_range] at hj
  have hjr : j ≤ Fintype.card P := by omega
  have hrw : ((1 : Polynomial ℤ) + X) ^ c *
      (Polynomial.C (antichainCount P j : ℤ) * X ^ j * (1 + X) ^ (Fintype.card P - j))
      = Polynomial.C (antichainCount P j : ℤ) *
        ((1 + X) ^ (Fintype.card P + c - j) * X ^ j) := by
    have hc : Fintype.card P + c - j = c + (Fintype.card P - j) := by omega
    rw [hc, pow_add]
    ring
  rw [hrw, Polynomial.coeff_C_mul, Polynomial.coeff_mul_X_pow']
  by_cases hcase : j ≤ Fintype.card P + c - d
  · rw [if_pos hcase, Polynomial.coeff_one_add_X_pow]
    have h1 : Fintype.card P + c - d - j = (Fintype.card P + c - j) - d := by omega
    have h2 : d ≤ Fintype.card P + c - j := by omega
    rw [h1, Nat.choose_symm h2]
  · rw [if_neg hcase]
    have : Fintype.card P + c - j < d := by omega
    rw [Nat.choose_eq_zero_of_lt this]
    simp

end MatchingBag
