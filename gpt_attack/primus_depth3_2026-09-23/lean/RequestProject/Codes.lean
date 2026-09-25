import Mathlib

/-!
# Binary codes, projections and the retained-coordinate profile

This file sets up the basic combinatorial notions used in the note
`matching_bag_poset_reduction_2026-08-20.md`:

* a *binary code* on a finite index type `ι` is a `Finset (ι → Bool)`;
* the projection `codeProj K Ω` of a code onto a set `K` of retained coordinates;
* the *retained-coordinate profile* `codeP Ω k = ∑_{|K| = k} |π_K(Ω)|`.
-/

open scoped BigOperators

namespace MatchingBag

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- A partial word: `none` at a coordinate means that the coordinate has been erased. -/
abbrev PartialWord (ι : Type*) := ι → Option Bool

/-- Restriction of a word to the set `K` of retained coordinates. -/
def restrictTo (K : Finset ι) (f : ι → Bool) : PartialWord ι :=
  fun i => if i ∈ K then some (f i) else none

/-- The projection `π_K(Ω)` of a code `Ω` onto the retained coordinates `K`. -/
def codeProj (K : Finset ι) (Ω : Finset (ι → Bool)) : Finset (PartialWord ι) :=
  Ω.image (restrictTo K)

/-- The retained-coordinate profile `p_k(Ω) = ∑_{|K| = k} |π_K(Ω)|`. -/
def codeP (Ω : Finset (ι → Bool)) (k : ℕ) : ℕ :=
  ∑ K ∈ (Finset.univ : Finset ι).powersetCard k, (codeProj K Ω).card

end MatchingBag
