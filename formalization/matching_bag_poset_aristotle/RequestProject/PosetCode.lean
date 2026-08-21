import RequestProject.Codes

/-!
# The order-ideal code of a finite poset and its bipartite graph

This file formalizes Section 2 of `matching_bag_poset_reduction_2026-08-20.md`.

* `MatchingBag.idealCode P` is the binary code `Ω(P)` whose words are the indicator
  functions of the order ideals (down-closed subsets) of a finite poset `P`.
* `MatchingBag.numInducedIdeals K` is `|J(P[K])|`, the number of order ideals of the
  induced subposet on `K`.
* Theorem `codeP_idealCode` is formula (3):
  `p_k(Ω(P)) = ∑_{|K| = k} |J(P[K])|`.
* `MatchingBag.posetBip P` is the graph `B(P)` of formula (4): two vertices `i^1`, `i^0`
  for each `i ∈ P`, with an edge `i^1 j^0` exactly when the partial assignment
  `x_i = 1, x_j = 0` violates the order (here, when `j ≤ i`, since a `1` propagates
  downwards in the ideal convention).
* Theorem `codeP_eq_indepSets_card` is formula (5) in coefficientwise form: the number
  of independent sets of `B(P)` of size `k` equals `p_k(Ω(P))`.
-/

open scoped BigOperators
open Finset

namespace MatchingBag

variable {P : Type*} [Fintype P] [DecidableEq P] [PartialOrder P]
  [DecidableRel ((· ≤ ·) : P → P → Prop)]

/-- `f : P → Bool` is the indicator of an order ideal (down-closed subset) of `P`. -/
def IsIdealIndicator (f : P → Bool) : Prop := ∀ i j : P, i ≤ j → f j = true → f i = true

instance (f : P → Bool) : Decidable (IsIdealIndicator f) := by
  unfold IsIdealIndicator; infer_instance

/-- The order-ideal code `Ω(P)`: all indicators of order ideals of `P`. -/
def idealCode (P : Type*) [Fintype P] [DecidableEq P] [PartialOrder P]
    [DecidableRel ((· ≤ ·) : P → P → Prop)] : Finset (P → Bool) :=
  Finset.univ.filter IsIdealIndicator

/-- The order ideals of the induced subposet `P[K]`, as indicator functions on `K`. -/
def inducedIdeals (K : Finset P) : Finset ({x // x ∈ K} → Bool) :=
  Finset.univ.filter (fun g => ∀ i j : {x // x ∈ K}, (i : P) ≤ (j : P) → g j = true → g i = true)

/-- `|J(P[K])|`, the number of order ideals of the induced subposet on `K`. -/
def numInducedIdeals (K : Finset P) : ℕ := (inducedIdeals K).card

/-- The partial word on `K` determined by an assignment `g : K → Bool`. -/
def extendPartial (K : Finset P) (g : {x // x ∈ K} → Bool) : PartialWord P :=
  fun i => if h : i ∈ K then some (g ⟨i, h⟩) else none

omit [Fintype P] [PartialOrder P] [DecidableRel ((· ≤ ·) : P → P → Prop)] in
lemma extendPartial_injective (K : Finset P) :
    Function.Injective (extendPartial K) := by
  intro g₁ g₂ h
  funext i
  have hi := congrFun h (i : P)
  simpa [extendPartial, i.2] using hi

/-- The projection of the order-ideal code to `K` consists exactly of the order ideals of
the induced subposet `P[K]`. -/
theorem codeProj_idealCode (K : Finset P) :
    codeProj K (idealCode P) = (inducedIdeals K).image (extendPartial K) := by
  ext w
  simp only [codeProj, idealCode, inducedIdeals, Finset.mem_image, Finset.mem_filter,
    Finset.mem_univ, true_and, IsIdealIndicator]
  constructor
  · rintro ⟨f, hf, rfl⟩
    refine ⟨fun i => f (i : P), fun i j hij hj => hf i j hij hj, ?_⟩
    funext i
    by_cases h : i ∈ K <;> simp [restrictTo, extendPartial, h]
  · rintro ⟨g, hg, rfl⟩
    set F : P → Bool := fun x => decide (∃ j : {y // y ∈ K}, x ≤ (j : P) ∧ g j = true) with hF
    have hFiff : ∀ x : P, F x = true ↔ ∃ j : {y // y ∈ K}, x ≤ (j : P) ∧ g j = true := by
      intro x; simp only [hF, decide_eq_true_eq]
    refine ⟨F, ?_, ?_⟩
    · intro i j hij hj
      rw [hFiff] at hj ⊢
      obtain ⟨k, hk1, hk2⟩ := hj
      exact ⟨k, le_trans hij hk1, hk2⟩
    · funext i
      by_cases h : i ∈ K
      · have key : F i = g ⟨i, h⟩ := by
          cases hgi : g ⟨i, h⟩ with
          | false =>
            have hne : ¬ (F i = true) := by
              rw [hFiff]
              rintro ⟨k, hk1, hk2⟩
              have := hg ⟨i, h⟩ k hk1 hk2
              rw [hgi] at this
              exact Bool.noConfusion this
            simpa using hne
          | true => exact (hFiff i).mpr ⟨⟨i, h⟩, le_refl (i : P), hgi⟩
        simp [restrictTo, extendPartial, h, key]
      · simp [restrictTo, extendPartial, h]

/-- Formula (3), pointwise: `|π_K(Ω(P))| = |J(P[K])|`. -/
theorem card_codeProj_idealCode (K : Finset P) :
    (codeProj K (idealCode P)).card = numInducedIdeals K := by
  rw [codeProj_idealCode, Finset.card_image_of_injective _ (extendPartial_injective K)]
  rfl

/-- **Formula (3)**: the retained-coordinate profile of the order-ideal code counts order
ideals of induced subposets, `p_k(Ω(P)) = ∑_{|K| = k} |J(P[K])|`. -/
theorem codeP_idealCode (k : ℕ) :
    codeP (idealCode P) k = ∑ K ∈ (Finset.univ : Finset P).powersetCard k, numInducedIdeals K := by
  unfold codeP
  exact Finset.sum_congr rfl fun K _ => card_codeProj_idealCode K

/-- The graph `B(P)` of formula (4). Vertices are `(i, b)` with `b = true` standing for
`i^1` and `b = false` for `i^0`; there is an edge `i^1 j^0` exactly when `j ≤ i`, i.e.
exactly when assigning `1` to `i` and `0` to `j` violates down-closedness. -/
def posetBip (P : Type*) [Fintype P] [DecidableEq P] [PartialOrder P]
    [DecidableRel ((· ≤ ·) : P → P → Prop)] : SimpleGraph (P × Bool) where
  Adj a b := (a.2 = true ∧ b.2 = false ∧ b.1 ≤ a.1) ∨ (a.2 = false ∧ b.2 = true ∧ a.1 ≤ b.1)
  symm := by
    rintro a b (⟨h1, h2, h3⟩ | ⟨h1, h2, h3⟩)
    · exact Or.inr ⟨h2, h1, h3⟩
    · exact Or.inl ⟨h2, h1, h3⟩
  loopless := ⟨by rintro a (⟨h1, h2, -⟩ | ⟨h1, h2, -⟩) <;> simp [h1] at h2⟩

instance : DecidableRel (posetBip P).Adj := fun a b => by
  unfold posetBip; infer_instance

/-- The independent sets of `B(P)`. -/
def indepSets (P : Type*) [Fintype P] [DecidableEq P] [PartialOrder P]
    [DecidableRel ((· ≤ ·) : P → P → Prop)] : Finset (Finset (P × Bool)) :=
  (Finset.univ : Finset (Finset (P × Bool))).filter
    (fun S => ∀ a ∈ S, ∀ b ∈ S, ¬ (posetBip P).Adj a b)

/-- The independent sets of `B(P)` of size `k`. -/
def indepSetsCard (P : Type*) [Fintype P] [DecidableEq P] [PartialOrder P]
    [DecidableRel ((· ≤ ·) : P → P → Prop)] (k : ℕ) : Finset (Finset (P × Bool)) :=
  (indepSets P).filter (fun S => S.card = k)

/-- The vertex set of `B(P)` encoding the partial assignment `g` on `K`. -/
def graphOf (K : Finset P) (g : {x // x ∈ K} → Bool) : Finset (P × Bool) :=
  K.attach.image (fun i => (i.1, g i))

section
omit [Fintype P] [PartialOrder P] [DecidableRel ((· ≤ ·) : P → P → Prop)]

lemma mem_graphOf {K : Finset P} {g : {x // x ∈ K} → Bool} {a : P × Bool} :
    a ∈ graphOf K g ↔ ∃ h : a.1 ∈ K, g ⟨a.1, h⟩ = a.2 := by
  constructor
  · intro ha
    rw [graphOf, Finset.mem_image] at ha
    obtain ⟨i, -, hi⟩ := ha
    subst hi
    exact ⟨i.2, by simp⟩
  · rintro ⟨h, hg⟩
    rw [graphOf, Finset.mem_image]
    exact ⟨⟨a.1, h⟩, Finset.mem_attach _ _, by simp [hg]⟩

lemma image_fst_graphOf (K : Finset P) (g : {x // x ∈ K} → Bool) :
    (graphOf K g).image Prod.fst = K := by
  ext x
  simp only [Finset.mem_image]
  constructor
  · rintro ⟨a, ha, rfl⟩
    exact (mem_graphOf.1 ha).1
  · intro hx
    exact ⟨(x, g ⟨x, hx⟩), mem_graphOf.2 ⟨hx, rfl⟩, rfl⟩

lemma card_graphOf (K : Finset P) (g : {x // x ∈ K} → Bool) :
    (graphOf K g).card = K.card := by
  have hinj : Function.Injective (fun i : {x // x ∈ K} => (i.1, g i)) := by
    intro i j hij
    exact Subtype.ext (congrArg Prod.fst hij)
  rw [graphOf, Finset.card_image_of_injective _ hinj, Finset.card_attach]

end

lemma indepSets_injOn_fst {S : Finset (P × Bool)} (hS : S ∈ indepSets P) :
    Set.InjOn Prod.fst (S : Set (P × Bool)) := by
  rw [indepSets, Finset.mem_filter] at hS
  intro a ha b hb hab
  by_contra hne
  have h2 : a.2 ≠ b.2 := fun h => hne (Prod.ext_iff.2 ⟨hab, h⟩)
  have ha' : a ∈ S := ha
  have hb' : b ∈ S := hb
  cases hA : a.2 with
  | true =>
    cases hB : b.2 with
    | true => exact h2 (by rw [hA, hB])
    | false => exact hS.2 a ha' b hb' (Or.inl ⟨hA, hB, le_of_eq hab.symm⟩)
  | false =>
    cases hB : b.2 with
    | true => exact hS.2 a ha' b hb' (Or.inr ⟨hA, hB, le_of_eq hab⟩)
    | false => exact h2 (by rw [hA, hB])

lemma graphOf_mem_indepSets_iff (K : Finset P) (g : {x // x ∈ K} → Bool) :
    graphOf K g ∈ indepSets P ↔ g ∈ inducedIdeals K := by
  rw [indepSets, Finset.mem_filter, inducedIdeals, Finset.mem_filter]
  simp only [Finset.mem_univ, true_and]
  constructor
  · intro h i j hij hj
    by_contra hi
    rw [Bool.not_eq_true] at hi
    exact h ((j : P), true) (mem_graphOf.2 ⟨j.2, by simpa using hj⟩)
      ((i : P), false) (mem_graphOf.2 ⟨i.2, by simpa using hi⟩)
      (Or.inl ⟨rfl, rfl, hij⟩)
  · rintro hg a ha b hb hab
    obtain ⟨ha1, ha2⟩ := mem_graphOf.1 ha
    obtain ⟨hb1, hb2⟩ := mem_graphOf.1 hb
    rcases hab with ⟨h1, h2, h3⟩ | ⟨h1, h2, h3⟩
    · have := hg ⟨b.1, hb1⟩ ⟨a.1, ha1⟩ h3 (by rw [ha2, h1])
      rw [hb2, h2] at this
      exact Bool.noConfusion this
    · have := hg ⟨a.1, ha1⟩ ⟨b.1, hb1⟩ h3 (by rw [hb2, h2])
      rw [ha2, h1] at this
      exact Bool.noConfusion this

/-- An independent set of `B(P)` is the graph of the partial assignment it defines on its
set of used coordinates. -/
lemma eq_graphOf_of_indep {S : Finset (P × Bool)} {K : Finset P} (hS1 : S ∈ indepSets P)
    (hS3 : S.image Prod.fst = K) :
    S = graphOf K (fun i => decide (((i.1, true) : P × Bool) ∈ S)) := by
  apply Finset.ext
  intro a
  constructor
  · intro ha
    have ha1 : a.1 ∈ K := by rw [← hS3]; exact Finset.mem_image_of_mem _ ha
    refine mem_graphOf.2 ⟨ha1, ?_⟩
    cases hA : a.2 with
    | true =>
      have hmem : ((a.1, true) : P × Bool) ∈ S := by rwa [← hA]
      simp [hmem]
    | false =>
      have hnot : ((a.1, true) : P × Bool) ∉ S := by
        intro hmem
        have heq := indepSets_injOn_fst hS1 (x₁ := ((a.1, true) : P × Bool)) (x₂ := a)
          (by exact_mod_cast hmem) (by exact_mod_cast ha) rfl
        rw [Prod.ext_iff] at heq
        rw [hA] at heq
        exact Bool.noConfusion heq.2
      simp [hnot]
  · intro ha
    obtain ⟨ha1, ha2⟩ := mem_graphOf.1 ha
    cases hA : a.2 with
    | true =>
      have hd : decide (((a.1 : P), true) ∈ S) = true := by rw [ha2, hA]
      have hmem : ((a.1, true) : P × Bool) ∈ S := by simpa using hd
      rwa [← hA] at hmem
    | false =>
      have hd : decide (((a.1 : P), true) ∈ S) = false := by rw [ha2, hA]
      have hnot : ((a.1, true) : P × Bool) ∉ S := by simpa using hd
      have ha1' : a.1 ∈ S.image Prod.fst := by rw [hS3]; exact ha1
      obtain ⟨b, hb, hb1⟩ := Finset.mem_image.1 ha1'
      cases hB : b.2 with
      | true =>
        exact absurd (by rwa [← hb1, ← hB] : ((a.1, true) : P × Bool) ∈ S) hnot
      | false =>
        have hba : b = a := Prod.ext_iff.2 ⟨hb1, by rw [hB, hA]⟩
        rwa [hba] at hb

/-- The independent sets of `B(P)` whose set of used coordinates is exactly `K` are in
bijection with the order ideals of the induced subposet `P[K]`. -/
lemma card_indepSets_fiber (K : Finset P) :
    ((indepSetsCard P K.card).filter (fun S => S.image Prod.fst = K)).card
      = numInducedIdeals K := by
  rw [numInducedIdeals]
  refine Finset.card_bij' (fun S _ => fun i => decide (((i.1, true) : P × Bool) ∈ S))
    (fun g _ => graphOf K g) ?_ ?_ ?_ ?_
  · intro S hS
    rw [Finset.mem_filter, indepSetsCard, Finset.mem_filter] at hS
    obtain ⟨⟨hS1, -⟩, hS3⟩ := hS
    have hSg := eq_graphOf_of_indep hS1 hS3
    have := (graphOf_mem_indepSets_iff K (fun i => decide (((i.1, true) : P × Bool) ∈ S))).1
      (by rw [← hSg]; exact hS1)
    exact this
  · intro g hg
    rw [Finset.mem_filter, indepSetsCard, Finset.mem_filter]
    exact ⟨⟨(graphOf_mem_indepSets_iff K g).2 hg, card_graphOf K g⟩, image_fst_graphOf K g⟩
  · intro S hS
    rw [Finset.mem_filter, indepSetsCard, Finset.mem_filter] at hS
    obtain ⟨⟨hS1, -⟩, hS3⟩ := hS
    exact (eq_graphOf_of_indep hS1 hS3).symm
  · intro g _
    funext i
    by_cases h : g i = true
    · have : ((i : P), true) ∈ graphOf K g := mem_graphOf.2 ⟨i.2, by simpa using h⟩
      simp [this, h]
    · rw [Bool.not_eq_true] at h
      have hnot : ((i : P), true) ∉ graphOf K g := by
        rw [mem_graphOf]
        rintro ⟨hi, hgi⟩
        rw [show (⟨(i : P), hi⟩ : {x // x ∈ K}) = i from Subtype.ext rfl, h] at hgi
        exact Bool.noConfusion hgi
      simp [hnot, h]

/-- **Formula (5)**, coefficientwise: the independent sets of `B(P)` of size `k` are exactly
the extendable partial assignments on `k` retained coordinates, so their number is
`p_k(Ω(P))`. -/
theorem codeP_eq_indepSets_card (k : ℕ) :
    codeP (idealCode P) k = (indepSetsCard P k).card := by
  rw [codeP_idealCode,
    Finset.card_eq_sum_card_fiberwise (f := fun S => S.image Prod.fst)
      (t := (Finset.univ : Finset P).powersetCard k) ?maps]
  case maps =>
    intro S hS
    rw [Finset.mem_coe, indepSetsCard, Finset.mem_filter] at hS
    rw [Finset.mem_coe, Finset.mem_powersetCard]
    refine ⟨Finset.subset_univ _, ?_⟩
    rw [Finset.card_image_of_injOn (indepSets_injOn_fst hS.1), hS.2]
  refine Finset.sum_congr rfl fun K hK => ?_
  rw [Finset.mem_powersetCard] at hK
  rw [← hK.2]
  exact (card_indepSets_fiber K).symm

end MatchingBag
