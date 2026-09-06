import Mathlib

/-!
# Minimum covers meet a maximum matching in exactly one endpoint per edge

Section 1 of `matching_bag_poset_reduction_2026-08-20.md` uses the following elementary
pigeonhole fact, which is what makes the variables `x_i` of the note well defined:

> a vertex cover whose size equals the size of a matching `M` contains exactly one endpoint
> of every edge of `M`, and no vertex outside the edges of `M`.

The statement below is graph-free: `M` is a finset of pairs with distinct endpoints which are
pairwise disjoint, and `C` is a finset of vertices meeting every edge of `M`.
-/

open Finset

namespace MatchingBag

variable {V : Type*} [DecidableEq V]

/-- If a set `C` meeting every edge of a matching `M` has exactly `|M|` elements, then it
contains exactly one endpoint of each edge of `M`, and every element of `C` is an endpoint of
an edge of `M`. -/
theorem cover_card_eq_matching_card
    (M : Finset (V × V)) (C : Finset V)
    (hne : ∀ e ∈ M, e.1 ≠ e.2)
    (hdisj : ∀ e ∈ M, ∀ e' ∈ M, e ≠ e' →
      e.1 ≠ e'.1 ∧ e.1 ≠ e'.2 ∧ e.2 ≠ e'.1 ∧ e.2 ≠ e'.2)
    (hcov : ∀ e ∈ M, e.1 ∈ C ∨ e.2 ∈ C)
    (hcard : C.card = M.card) :
    (∀ e ∈ M, ¬(e.1 ∈ C ∧ e.2 ∈ C)) ∧ (∀ v ∈ C, ∃ e ∈ M, v = e.1 ∨ v = e.2) := by
  classical
  set f : V × V → V := fun e => if e.1 ∈ C then e.1 else e.2 with hf
  have hfmem : ∀ e ∈ M, f e = e.1 ∨ f e = e.2 := by
    intro e _
    by_cases h : e.1 ∈ C
    · exact Or.inl (by simp [hf, h])
    · exact Or.inr (by simp [hf, h])
  have hfC : ∀ e ∈ M, f e ∈ C := by
    intro e he
    by_cases h : e.1 ∈ C
    · simp [hf, h]
    · rcases hcov e he with h1 | h2
      · exact absurd h1 h
      · simpa [hf, h] using h2
  have hinj : Set.InjOn f (M : Set (V × V)) := by
    intro e he e' he' hee'
    by_contra hne'
    obtain ⟨h1, h2, h3, h4⟩ := hdisj e he e' he' hne'
    rcases hfmem e he with hA | hA <;> rcases hfmem e' he' with hB | hB <;>
      rw [hA, hB] at hee'
    · exact h1 hee'
    · exact h2 hee'
    · exact h3 hee'
    · exact h4 hee'
  have himg : M.image f = C := by
    refine Finset.eq_of_subset_of_card_le ?_ ?_
    · intro v hv
      obtain ⟨e, he, rfl⟩ := Finset.mem_image.1 hv
      exact hfC e he
    · rw [Finset.card_image_of_injOn hinj, hcard]
  refine ⟨?_, ?_⟩
  · rintro e he ⟨h1, h2⟩
    have h2' : e.2 ∈ M.image f := by rw [himg]; exact h2
    obtain ⟨e', he', hfe'⟩ := Finset.mem_image.1 h2'
    by_cases hee' : e' = e
    · subst hee'
      have : f e' = e'.1 := by simp [hf, h1]
      rw [this] at hfe'
      exact hne e' he' hfe'
    · obtain ⟨-, -, h3, h4⟩ := hdisj e he e' he' (Ne.symm hee')
      rcases hfmem e' he' with hA | hA <;> rw [hA] at hfe'
      · exact h3 hfe'.symm
      · exact h4 hfe'.symm
  · intro v hv
    rw [← himg] at hv
    obtain ⟨e, he, rfl⟩ := Finset.mem_image.1 hv
    exact ⟨e, he, hfmem e he⟩

end MatchingBag
