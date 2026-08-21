import RequestProject.TreeMatching
import RequestProject.PosetCode

/-!
# The forest poset of a tree with a maximum matching

This file carries out the second half of §1 of `matching_bag_poset_reduction_2026-08-20.md`:
propagating and deleting the forced variables of the constraint system, and identifying the
solutions of the remaining system with the order ideals of a poset whose cover graph is a
forest.

* `TreeMatching.Forced i` says that the variable `x_i` takes the same value in *every*
  solution; `TreeMatching.forcedVal i` is that value.  (Since the solution set is nonempty by
  König's theorem, this is exactly the propagation of the unary boundary values, and its
  consistency is automatic.)
* `TreeMatching.Free` is the type of remaining (free) variables.  It carries a partial order:
  `i ≤ j` iff there is a chain of comparisons `j → ⋯ → i`.  Antisymmetry is *proved*, from
  acyclicity of the forest (`no_closed_matching_chain`), not assumed.
* `TreeMatching.Sol_image_restrict`: restricting a solution to the free variables is a
  bijection from the solution set onto the order-ideal code `Ω(P)` of `P = Free`.
* `TreeMatching.compGraph_isAcyclic`: the undirected comparison graph on all matched bags is
  a forest (this is the contraction statement `bagGraph_isAcyclic`).
* `TreeMatching.coverGraph_isAcyclic`: the cover graph of `P` is a forest.
-/

open Finset SimpleGraph

namespace MatchingBag

attribute [local instance 1] Classical.propDecidable

variable {V : Type*} [Fintype V] [DecidableEq V]

namespace TreeMatching

variable (D : TreeMatching V)

/-! ### Forced and free variables -/

/-- A variable is *forced* when it takes the same value in every solution. -/
def Forced (i : D.Idx) : Prop := (∀ x ∈ D.Sol, x i = true) ∨ (∀ x ∈ D.Sol, x i = false)

/-- The value of a forced variable. -/
noncomputable def forcedVal (i : D.Idx) : Bool := decide (∀ x ∈ D.Sol, x i = true)

/-- The free variables: the carrier of the poset `P` of the note. -/
def Free : Type _ := {i : D.Idx // ¬ D.Forced i}

noncomputable instance : DecidableEq D.Free := Classical.decEq _
instance : Finite D.Free := by unfold Free; infer_instance
noncomputable instance : Fintype D.Free := Fintype.ofFinite _

variable {D}

lemma forcedVal_eq_true_iff {i : D.Idx} : D.forcedVal i = true ↔ ∀ x ∈ D.Sol, x i = true := by
  simp [forcedVal]

lemma forcedVal_spec {i : D.Idx} (h : D.Forced i) {x : D.Idx → Bool} (hx : x ∈ D.Sol) :
    x i = D.forcedVal i := by
  rcases h with hall | hall
  · rw [hall x hx]
    exact (forcedVal_eq_true_iff.2 hall).symm
  · rw [hall x hx]
    symm
    simp only [forcedVal, decide_eq_false_iff_not]
    push_neg
    exact ⟨x, hx, by rw [hall x hx]; simp⟩

lemma Free.coe_injective : Function.Injective (fun i : D.Free => i.1) := fun _ _ h => Subtype.ext h

/-- A closed chain of the comparison relation between free variables gives a closed chain of
the comparison relation on all matched bags. -/
lemma transGen_rel_of_free {i j : D.Free}
    (h : Relation.TransGen (fun a b : D.Free => D.rel b.1 a.1) i j) :
    Relation.TransGen D.rel j.1 i.1 := by
  induction h with
  | single hs => exact Relation.TransGen.single hs
  | tail _ hbc ih => exact Relation.TransGen.head hbc ih

/-! ### The order on the free variables -/

variable (D)

/-- The order on the free variables: `i ≤ j` when there is a chain of comparisons from `j`
down to `i`, so that `x_j = 1` forces `x_i = 1`.  With this convention the solutions restrict
to indicators of order *ideals*. -/
noncomputable instance instPartialOrderFree : PartialOrder D.Free where
  le i j := Relation.ReflTransGen (fun a b : D.Free => D.rel b.1 a.1) i j
  le_refl _ := Relation.ReflTransGen.refl
  le_trans _ _ _ h₁ h₂ := h₁.trans h₂
  le_antisymm := by
    intro i j h1 h2
    by_contra hne
    have ht1 : Relation.TransGen (fun a b : D.Free => D.rel b.1 a.1) i j := by
      rcases Relation.reflTransGen_iff_eq_or_transGen.1 h1 with hh | hh
      · exact absurd hh.symm hne
      · exact hh
    have ht : Relation.TransGen (fun a b : D.Free => D.rel b.1 a.1) i i := ht1.trans_left h2
    exact no_closed_matching_chain D.acyclic disjointEdges i.1 (transGen_rel_of_free ht)

noncomputable instance : DecidableRel ((· ≤ ·) : D.Free → D.Free → Prop) := fun _ _ =>
  Classical.dec _

variable {D}

lemma le_free_iff {i j : D.Free} :
    i ≤ j ↔ Relation.ReflTransGen (fun a b : D.Free => D.rel b.1 a.1) i j := Iff.rfl

/-- One step of the comparison relation between free variables is a relation of the order. -/
lemma le_of_rel {i j : D.Free} (h : D.rel j.1 i.1) : i ≤ j :=
  Relation.ReflTransGen.single h

/-! ### Solutions and order ideals -/

/-- The restriction of an assignment to the free variables. -/
def restrictFree (x : D.Idx → Bool) : D.Free → Bool := fun i => x i.1

/-- The assignment obtained from an assignment of the free variables by filling in the forced
values. -/
noncomputable def extendFree (D : TreeMatching V) (y : D.Free → Bool) : D.Idx → Bool :=
  fun i => if h : D.Forced i then D.forcedVal i else y ⟨i, h⟩

lemma extendFree_of_forced {y : D.Free → Bool} {i : D.Idx} (h : D.Forced i) :
    D.extendFree y i = D.forcedVal i := dif_pos h

lemma extendFree_of_not_forced {y : D.Free → Bool} {i : D.Idx} (h : ¬ D.Forced i) :
    D.extendFree y i = y ⟨i, h⟩ := dif_neg h

/-- Restricting a solution to the free variables gives the indicator of an order ideal. -/
theorem isIdealIndicator_restrictFree {x : D.Idx → Bool} (hx : x ∈ D.Sol) :
    IsIdealIndicator (restrictFree x) := by
  intro a b hab
  induction hab with
  | refl => exact id
  | tail _ hcb ih =>
      intro hb
      exact ih ((mem_Sol.1 hx).1 _ _ hcb hb)

/-- Conversely, every order ideal of `P` extends (uniquely) to a solution: the forced
variables take their forced values.  This is the step of the note which shows that no
comparison between two free variables can pass through a forced variable. -/
theorem extendFree_mem_Sol {y : D.Free → Bool} (hy : IsIdealIndicator y) :
    D.extendFree y ∈ D.Sol := by
  refine mem_Sol.2 ⟨?_, ?_, ?_⟩
  · intro i j hij hxi
    by_cases hfi : D.Forced i
    · rw [extendFree_of_forced hfi] at hxi
      have hall : ∀ z ∈ D.Sol, z i = true := forcedVal_eq_true_iff.1 hxi
      have hallj : ∀ z ∈ D.Sol, z j = true := fun z hz => (mem_Sol.1 hz).1 i j hij (hall z hz)
      by_cases hfj : D.Forced j
      · rw [extendFree_of_forced hfj]
        exact forcedVal_eq_true_iff.2 hallj
      · exact absurd (Or.inl hallj) hfj
    · rw [extendFree_of_not_forced hfi] at hxi
      by_cases hfj : D.Forced j
      · rw [extendFree_of_forced hfj]
        refine forcedVal_eq_true_iff.2 ?_
        rcases hfj with hall | hall
        · exact hall
        · refine absurd (Or.inr ?_) hfi
          intro z hz
          by_contra hzi
          simp only [Bool.not_eq_false] at hzi
          have := (mem_Sol.1 hz).1 i j hij hzi
          rw [hall z hz] at this
          exact Bool.noConfusion this
      · rw [extendFree_of_not_forced hfj]
        exact hy ⟨j, hfj⟩ ⟨i, hfi⟩ (le_of_rel hij) hxi
  · intro i hz
    have hall : ∀ z ∈ D.Sol, z i = false := fun z hzs => (mem_Sol.1 hzs).2.1 i hz
    rw [extendFree_of_forced (Or.inr hall : D.Forced i)]
    simp only [forcedVal, decide_eq_false_iff_not]
    push_neg
    obtain ⟨z, hzS⟩ := Sol_nonempty D
    exact ⟨z, hzS, by rw [hall z hzS]; simp⟩
  · intro i hone
    have hall : ∀ z ∈ D.Sol, z i = true := fun z hzs => (mem_Sol.1 hzs).2.2 i hone
    rw [extendFree_of_forced (Or.inl hall : D.Forced i)]
    exact forcedVal_eq_true_iff.2 hall

theorem restrictFree_extendFree (y : D.Free → Bool) : restrictFree (D.extendFree y) = y := by
  funext i
  obtain ⟨i, hi⟩ := i
  show D.extendFree y i = y ⟨i, hi⟩
  rw [extendFree_of_not_forced hi]

theorem extendFree_restrictFree {x : D.Idx → Bool} (hx : x ∈ D.Sol) :
    D.extendFree (restrictFree x) = x := by
  funext i
  by_cases h : D.Forced i
  · rw [extendFree_of_forced h, ← forcedVal_spec h hx]
  · rw [extendFree_of_not_forced h]
    rfl

/-- **The solutions of the constraint system are exactly the order ideals of `P`.** -/
theorem Sol_image_restrict (D : TreeMatching V) :
    D.Sol.image restrictFree = idealCode D.Free := by
  ext y
  simp only [Finset.mem_image, idealCode, Finset.mem_filter, Finset.mem_univ, true_and]
  constructor
  · rintro ⟨x, hx, rfl⟩
    exact isIdealIndicator_restrictFree hx
  · intro hy
    exact ⟨D.extendFree y, extendFree_mem_Sol hy, restrictFree_extendFree y⟩

theorem restrictFree_injOn (D : TreeMatching V) :
    Set.InjOn restrictFree (D.Sol : Set (D.Idx → Bool)) := by
  intro x hx x' hx' h
  rw [← extendFree_restrictFree (Finset.mem_coe.1 hx),
    ← extendFree_restrictFree (Finset.mem_coe.1 hx'), h]

/-! ### The comparison graph and the cover graph are forests -/

variable (D)

/-- The undirected comparison graph on the matched bags. -/
def compGraph : SimpleGraph D.Idx := SimpleGraph.fromRel D.rel

/-- The cover graph of the poset `P`. -/
noncomputable def coverGraph : SimpleGraph D.Free := SimpleGraph.fromRel (fun i j : D.Free => i ⋖ j)

/-- **The undirected comparison graph is a forest**: this is exactly the statement that
contracting the matching edges of the tree gives a tree. -/
theorem compGraph_isAcyclic : (D.compGraph).IsAcyclic := by
  have h := bagGraph_isAcyclic D.acyclic (disjointEdges (D := D))
  have hg : D.compGraph = SimpleGraph.fromRel (fun i j : D.Idx => D.G.Adj (D.Rv i) (D.Lv j)) := by
    ext a b
    simp only [compGraph, SimpleGraph.fromRel_adj, rel]
    constructor
    · rintro ⟨hne, h1 | h1⟩
      · exact ⟨hne, Or.inl h1.2⟩
      · exact ⟨hne, Or.inr h1.2⟩
    · rintro ⟨hne, h1 | h1⟩
      · exact ⟨hne, Or.inl ⟨hne, h1⟩⟩
      · exact ⟨hne, Or.inr ⟨Ne.symm hne, h1⟩⟩
  rw [hg]
  exact h

variable {D}

/-- A cover relation of `P` is a single comparison, not a longer chain. -/
theorem rel_of_covBy {i j : D.Free} (h : i ⋖ j) : D.rel j.1 i.1 := by
  have hlt := h.1
  have hne : i ≠ j := ne_of_lt hlt
  have htg : Relation.TransGen (fun a b : D.Free => D.rel b.1 a.1) i j := by
    rcases Relation.reflTransGen_iff_eq_or_transGen.1 (le_free_iff.1 hlt.le) with hh | hh
    · exact absurd hh.symm hne
    · exact hh
  obtain ⟨k, hik, hkj⟩ := Relation.TransGen.head'_iff.1 htg
  rcases eq_or_ne i k with rfl | hik2
  · exact absurd rfl hik.1
  · have hkeq : k = j := by
      by_contra hkj2
      exact h.2 (lt_of_le_of_ne (le_free_iff.2 (Relation.ReflTransGen.single hik)) hik2)
        (lt_of_le_of_ne (le_free_iff.2 hkj) hkj2)
    subst hkeq
    exact hik

variable (D)

/-- **The cover graph of the forest poset `P` is a forest.** -/
theorem coverGraph_isAcyclic : (D.coverGraph).IsAcyclic := by
  refine SimpleGraph.IsAcyclic.comap (G' := D.compGraph) ⟨fun i => i.1, ?_⟩
    Free.coe_injective (compGraph_isAcyclic D)
  · intro a b hab
    simp only [coverGraph, SimpleGraph.fromRel_adj] at hab
    obtain ⟨hne, hcov⟩ := hab
    simp only [compGraph, SimpleGraph.fromRel_adj]
    refine ⟨fun hh => hne (Subtype.ext hh), ?_⟩
    rcases hcov with hc | hc
    · exact Or.inr (rel_of_covBy hc)
    · exact Or.inl (rel_of_covBy hc)

end TreeMatching

end MatchingBag
