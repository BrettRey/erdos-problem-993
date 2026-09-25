import RequestProject.Codes

/-!
# Adding a constant coordinate multiplies the profile polynomial by `1 + t`

Section 2 of `matching_bag_poset_reduction_2026-08-20.md` uses the fact that adjoining a
coordinate on which the code is constant multiplies `∑_k p_k(Ω) t^k` by `1 + t`.
Coefficientwise this is `p_{k+1}(Ω') = p_{k+1}(Ω) + p_k(Ω)` and `p_0(Ω') = p_0(Ω)`.
-/

open scoped BigOperators
open Finset

namespace MatchingBag

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- The code `Ω` with one extra coordinate (indexed by `none`) on which every word is `1`. -/
def appendConst (Ω : Finset (ι → Bool)) : Finset (Option ι → Bool) :=
  Ω.image (fun f o => o.elim true f)

/-- Projecting the extended code onto a set of retained coordinates whose non-`none` part is
`K'` gives as many partial words as projecting `Ω` onto `K'`: the extra coordinate carries no
information. -/
lemma card_codeProj_appendConst (Ω : Finset (ι → Bool)) (K : Finset (Option ι)) (K' : Finset ι)
    (hK : ∀ i : ι, (some i ∈ K ↔ i ∈ K')) :
    (codeProj K (appendConst Ω)).card = (codeProj K' Ω).card := by
  classical
  set b : Option Bool := if none ∈ K then some true else none with hb
  set Phi : PartialWord ι → PartialWord (Option ι) := fun w o => o.elim b w with hPhi
  have hPhiinj : Function.Injective Phi := by
    intro w1 w2 h
    funext i
    exact congrFun h (some i)
  have hcomp : ∀ f : ι → Bool,
      restrictTo K (fun o => o.elim true f) = Phi (restrictTo K' f) := by
    intro f
    funext o
    cases o with
    | none => simp [restrictTo, hPhi, hb]
    | some i => simp [restrictTo, hPhi, hK i]
  rw [codeProj, appendConst, Finset.image_image,
    show ((restrictTo K) ∘ (fun (f : ι → Bool) (o : Option ι) => o.elim true f))
      = Phi ∘ (restrictTo K') from funext hcomp,
    ← Finset.image_image, Finset.card_image_of_injective _ hPhiinj, codeProj]

omit [Fintype ι] in
lemma card_image_some (K' : Finset ι) : (K'.image some).card = K'.card :=
  Finset.card_image_of_injective _ (Option.some_injective ι)

omit [Fintype ι] in
lemma none_notMem_image_some (K' : Finset ι) : none ∉ K'.image some := by
  intro h
  obtain ⟨i, -, hi⟩ := Finset.mem_image.1 h
  exact Option.some_ne_none i hi

/-- The profile of a code with a constant coordinate adjoined, in degree `0`. -/
theorem codeP_appendConst_zero (Ω : Finset (ι → Bool)) :
    codeP (appendConst Ω) 0 = codeP Ω 0 := by
  rw [codeP, codeP, Finset.powersetCard_zero, Finset.powersetCard_zero, Finset.sum_singleton,
    Finset.sum_singleton]
  exact card_codeProj_appendConst Ω ∅ ∅ (by simp)

/-- Adjoining a constant coordinate multiplies `∑_k p_k t^k` by `1 + t`. -/
theorem codeP_appendConst_succ (Ω : Finset (ι → Bool)) (k : ℕ) :
    codeP (appendConst Ω) (k + 1) = codeP Ω (k + 1) + codeP Ω k := by
  classical
  rw [codeP, ← Finset.sum_filter_add_sum_filter_not
    ((Finset.univ : Finset (Option ι)).powersetCard (k + 1)) (fun K => none ∈ K)]
  rw [add_comm]
  congr 1
  · -- the retained sets not containing the constant coordinate
    refine (Finset.sum_nbij (i := fun K' : Finset ι => K'.image some) ?_ ?_ ?_ ?_).symm
    · intro K' hK'
      rw [Finset.mem_powersetCard] at hK'
      rw [Finset.mem_filter, Finset.mem_powersetCard]
      exact ⟨⟨Finset.subset_univ _, by rw [card_image_some, hK'.2]⟩,
        none_notMem_image_some K'⟩
    · intro K1 _ K2 _ h
      exact Finset.image_injective (Option.some_injective ι) h
    · intro K hK
      simp only [Finset.coe_filter, Set.mem_setOf_eq, Finset.mem_powersetCard] at hK
      obtain ⟨⟨-, hcard⟩, hnone⟩ := hK
      refine ⟨(Finset.univ : Finset ι).filter (fun i => some i ∈ K), ?_, ?_⟩
      · rw [Finset.mem_coe, Finset.mem_powersetCard]
        refine ⟨Finset.subset_univ _, ?_⟩
        have himg : ((Finset.univ : Finset ι).filter (fun i => some i ∈ K)).image some = K := by
          ext o
          cases o with
          | none => simp [hnone]
          | some i => simp
        rw [← card_image_some, himg, hcard]
      · ext o
        cases o with
        | none => simp [hnone]
        | some i => simp
    · intro K' _
      exact (card_codeProj_appendConst Ω (K'.image some) K' (by simp)).symm
  · -- the retained sets containing the constant coordinate
    refine (Finset.sum_nbij (i := fun K' : Finset ι => insert none (K'.image some)) ?_ ?_ ?_ ?_).symm
    · intro K' hK'
      rw [Finset.mem_powersetCard] at hK'
      rw [Finset.mem_filter, Finset.mem_powersetCard]
      refine ⟨⟨Finset.subset_univ _, ?_⟩, Finset.mem_insert_self _ _⟩
      rw [Finset.card_insert_of_notMem (none_notMem_image_some K'), card_image_some, hK'.2]
    · intro K1 _ K2 _ h
      have h' : K1.image some = K2.image some := by
        have h2 := congrArg (fun s : Finset (Option ι) => s.erase none) h
        simpa [Finset.erase_insert, none_notMem_image_some] using h2
      exact Finset.image_injective (Option.some_injective ι) h'
    · intro K hK
      simp only [Finset.coe_filter, Set.mem_setOf_eq, Finset.mem_powersetCard] at hK
      obtain ⟨⟨-, hcard⟩, hnone⟩ := hK
      have himg : insert none (((Finset.univ : Finset ι).filter
          (fun i => some i ∈ K)).image some) = K := by
        ext o
        cases o with
        | none => simp [hnone]
        | some i => simp
      refine ⟨(Finset.univ : Finset ι).filter (fun i => some i ∈ K), ?_, himg⟩
      rw [Finset.mem_coe, Finset.mem_powersetCard]
      refine ⟨Finset.subset_univ _, ?_⟩
      have := congrArg Finset.card himg
      rw [Finset.card_insert_of_notMem (none_notMem_image_some _), card_image_some, hcard] at this
      omega
    · intro K' _
      exact (card_codeProj_appendConst Ω (insert none (K'.image some)) K' (by simp)).symm

end MatchingBag
