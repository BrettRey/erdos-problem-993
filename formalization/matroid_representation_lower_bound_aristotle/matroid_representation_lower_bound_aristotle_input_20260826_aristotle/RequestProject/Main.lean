import Mathlib

open scoped BigOperators
open scoped Real
open scoped Nat
open scoped Classical
open scoped Pointwise

set_option maxHeartbeats 8000000
set_option maxRecDepth 4000
set_option synthInstance.maxHeartbeats 20000
set_option synthInstance.maxSize 128

set_option relaxedAutoImplicit false
set_option autoImplicit false

set_option grind.warning false

/-!
# The matroid-representation lower bound for graph independence systems

Main result (`MatroidRep.le_card_of_represents`): if the independence system of a
simple graph `G` on `α` is the intersection of the independence systems of `k`
matroids `M : Fin k → Matroid α`, then `deg v ≤ k` for every vertex `v` whose
neighbourhood is an independent set.
-/

namespace MatroidRep

variable {α : Type*} {k : ℕ}

/-- A finite family of matroids on `α` *represents* the graph `G` when the independent
sets of `G` are exactly the sets that are independent in every matroid of the family. -/
def Represents (G : SimpleGraph α) (M : Fin k → Matroid α) : Prop :=
  ∀ S : Set α, G.IsIndepSet S ↔ ∀ i, (M i).Indep S

/-- Step 1: every singleton is independent in each matroid of a representing family
(so no matroid of the family has a loop). -/
theorem indep_singleton_of_represents {G : SimpleGraph α} {M : Fin k → Matroid α}
    (hrep : Represents G M) (x : α) (i : Fin k) : (M i).Indep {x} :=
  (hrep {x}).1 (Set.pairwise_singleton _ _) i

/-- Step 1': every element of `α` is a nonloop of each matroid of a representing family. -/
theorem isNonloop_of_represents {G : SimpleGraph α} {M : Fin k → Matroid α}
    (hrep : Represents G M) (x : α) (i : Fin k) : (M i).IsNonloop x :=
  Matroid.indep_singleton.1 (indep_singleton_of_represents hrep x i)

/-- A dependent pair in a loopless matroid is a circuit. -/
theorem isCircuit_pair_of_represents {G : SimpleGraph α} {M : Fin k → Matroid α}
    (hrep : Represents G M) {x y : α} {i : Fin k} (hdep : ¬ (M i).Indep {x, y}) :
    (M i).IsCircuit {x, y} := by
  have hx : (M i).Indep {x} := indep_singleton_of_represents hrep x i
  have hy : (M i).Indep {y} := indep_singleton_of_represents hrep y i
  refine Matroid.isCircuit_iff_forall_ssubset.2 ⟨?_, ?_⟩
  · refine Matroid.dep_iff.2 ⟨hdep, ?_⟩
    rintro z (rfl | rfl)
    · exact hx.subset_ground rfl
    · exact hy.subset_ground rfl
  · intro I hI
    rcases Set.subset_pair_iff_eq.mp hI.subset with rfl | rfl | rfl | rfl
    · exact (M i).empty_indep
    · exact hx
    · exact hy
    · exact absurd rfl hI.ne

/-- Step 3, key step: if two distinct neighbours `u`, `w` of `v` are non-adjacent and
both pairs `{v, u}` and `{v, w}` are dependent in the *same* matroid `M i`, we get a
contradiction, via circuit elimination. -/
theorem not_dep_pair_of_dep_pair {G : SimpleGraph α} {M : Fin k → Matroid α}
    (hrep : Represents G M) {v u w : α} {i : Fin k}
    (huv : u ≠ v) (hwv : w ≠ v) (huw : u ≠ w) (hG : ¬ G.Adj u w)
    (hu : ¬ (M i).Indep {v, u}) (hw : ¬ (M i).Indep {v, w}) : False := by
  have hC₁ : (M i).IsCircuit {v, u} := isCircuit_pair_of_represents hrep hu
  have hC₂ : (M i).IsCircuit {v, w} := isCircuit_pair_of_represents hrep hw
  have hvC₁ : v ∈ ({v, u} : Set α) := by simp
  have hvC₂ : v ∈ ({v, w} : Set α) := by simp
  have huC₁ : u ∈ ({v, u} : Set α) := by simp
  have huC₂ : u ∉ ({v, w} : Set α) := by
    simp only [Set.mem_insert_iff, Set.mem_singleton_iff]
    exact fun h => h.elim huv huw
  obtain ⟨C, hCsub, hC, huC⟩ := hC₁.strong_elimination hC₂ hvC₁ hvC₂ huC₁ huC₂
  -- `C ⊆ {u, w}` and `u ∈ C`
  have hCsub' : C ⊆ {u, w} := by
    intro z hz
    have hz' := hCsub hz
    simp only [Set.mem_diff, Set.mem_union, Set.mem_insert_iff, Set.mem_singleton_iff] at hz'
    simp only [Set.mem_insert_iff, Set.mem_singleton_iff]
    tauto
  -- the pair `{u, w}` is independent in `M i`, since `u` and `w` are non-adjacent
  have hindep : (M i).Indep {u, w} := by
    refine (hrep {u, w}).1 ?_ i
    rw [SimpleGraph.isIndepSet_iff]
    intro a ha b hb hab
    simp only [Set.mem_insert_iff, Set.mem_singleton_iff] at ha hb
    rcases ha with rfl | rfl <;> rcases hb with rfl | rfl <;> simp_all [G.adj_comm]
  exact hC.dep.not_indep (hindep.subset hCsub')

/-- **G0.** If a family of `k` matroids on `α` represents the independence system of `G`,
then the degree of any vertex `v` whose neighbourhood is independent is at most `k`. -/
theorem le_card_of_represents {G : SimpleGraph α} {M : Fin k → Matroid α}
    (hrep : Represents G M) (v : α) (hv : G.IsIndepSet (G.neighborSet v))
    (hfin : (G.neighborSet v).Finite) : hfin.toFinset.card ≤ k := by
  classical
  -- Step 2: each edge at `v` is killed by some matroid.
  have key : ∀ u ∈ G.neighborSet v, ∃ i : Fin k, ¬ (M i).Indep {v, u} := by
    intro u hu
    have hadj : G.Adj v u := hu
    have : ¬ G.IsIndepSet ({v, u} : Set α) := by
      rw [SimpleGraph.isIndepSet_iff]
      intro h
      exact h (by simp) (by simp) hadj.ne hadj
    rw [hrep] at this
    push_neg at this
    exact this
  rcases Set.eq_empty_or_nonempty (G.neighborSet v) with hemp | ⟨u₀, hu₀⟩
  · simp [hemp]
  · obtain ⟨i₀, -⟩ := key u₀ hu₀
    have : Nonempty (Fin k) := ⟨i₀⟩
    choose! f hf using key
    have hinj : Set.InjOn f (G.neighborSet v) := by
      intro u hu w hw hfe
      by_contra hne
      exact not_dep_pair_of_dep_pair hrep (G.ne_of_adj hu).symm (G.ne_of_adj hw).symm hne
        (hv hu hw hne) (hf u hu) (hfe ▸ hf w hw)
    calc hfin.toFinset.card
        ≤ (Finset.univ : Finset (Fin k)).card := by
          refine Finset.card_le_card_of_injOn f (fun a _ => by simp) ?_
          intro a ha b hb hab
          exact hinj (by simpa using ha) (by simpa using hb) hab
      _ = k := by simp

/-- **G0'**, `Fintype` form: with `α` a fintype and `G` locally finite, the degree of a
vertex with independent neighbourhood is at most the number of representing matroids. -/
theorem degree_le_of_represents [Fintype α] {G : SimpleGraph α}
    [DecidableRel G.Adj] {M : Fin k → Matroid α} (hrep : Represents G M) (v : α)
    (hv : G.IsIndepSet (G.neighborSet v)) : G.degree v ≤ k := by
  classical
  have hfin : (G.neighborSet v).Finite := Set.toFinite _
  have := le_card_of_represents hrep v hv hfin
  have hcard : hfin.toFinset.card = G.degree v := by
    rw [← SimpleGraph.card_neighborFinset_eq_degree]
    congr 1
    ext x
    simp
  omega

/-! ## G1: the triangle-free case -/

/-- In a triangle-free graph every neighbourhood is an independent set. -/
theorem isIndepSet_neighborSet_of_cliqueFree {G : SimpleGraph α} (h : G.CliqueFree 3) (v : α) :
    G.IsIndepSet (G.neighborSet v) := by
  rw [SimpleGraph.isIndepSet_iff]
  intro u hu w hw hne hadj
  refine h {v, u, w} ?_
  have hvu : G.Adj v u := hu
  have hvw : G.Adj v w := hw
  rw [SimpleGraph.isNClique_iff]
  refine ⟨?_, ?_⟩
  · rw [SimpleGraph.isClique_iff]
    intro a ha b hb hab
    simp only [Finset.coe_insert, Set.mem_insert_iff, Finset.coe_singleton,
      Set.mem_singleton_iff] at ha hb
    rcases ha with rfl | rfl | rfl <;> rcases hb with rfl | rfl | rfl <;>
      simp_all [G.adj_comm]
  · have h1 : v ≠ u := hvu.ne
    have h2 : v ≠ w := hvw.ne
    rw [Finset.card_insert_of_notMem (by simp [h1, h2]),
      Finset.card_insert_of_notMem (by simp [hne]), Finset.card_singleton]

/-- **G1.** In a triangle-free graph, every vertex degree is a lower bound for the number
of matroids in any representing family. -/
theorem le_card_of_represents_cliqueFree {G : SimpleGraph α} {M : Fin k → Matroid α}
    (hrep : Represents G M) (h3 : G.CliqueFree 3) (v : α)
    (hfin : (G.neighborSet v).Finite) : hfin.toFinset.card ≤ k :=
  le_card_of_represents hrep v (isIndepSet_neighborSet_of_cliqueFree h3 v) hfin

/-- An acyclic graph is triangle-free. -/
theorem cliqueFree_three_of_isAcyclic {G : SimpleGraph α} (h : G.IsAcyclic) : G.CliqueFree 3 := by
  classical
  intro s hs
  rw [SimpleGraph.isNClique_iff] at hs
  obtain ⟨hclique, hcard⟩ := hs
  obtain ⟨v, u, w, hvu, hvw, huw, hse⟩ := Finset.card_eq_three.mp hcard
  subst hse
  have a1 : G.Adj v u := hclique (by simp) (by simp) hvu
  have a2 : G.Adj u w := hclique (by simp) (by simp) huw
  have a3 : G.Adj w v := (hclique (by simp) (by simp) hvw).symm
  exact h (SimpleGraph.Walk.cons a1 (SimpleGraph.Walk.cons a2
      (SimpleGraph.Walk.cons a3 SimpleGraph.Walk.nil)))
    (by simp [SimpleGraph.Walk.isCycle_def, SimpleGraph.Walk.isTrail_def, hvu, hvw, huw,
      hvu.symm, hvw.symm])

/-- **G1, tree form.** A tree is triangle-free, so the degree bound applies. -/
theorem le_card_of_represents_isTree {G : SimpleGraph α} {M : Fin k → Matroid α}
    (hrep : Represents G M) (ht : G.IsTree) (v : α)
    (hfin : (G.neighborSet v).Finite) : hfin.toFinset.card ≤ k :=
  le_card_of_represents_cliqueFree hrep (cliqueFree_three_of_isAcyclic ht.IsAcyclic) v hfin

/-! ## G2: robustness under ground-set enlargement (restrictions and minors) -/

/-- **G2, reduction lemma.** A representation of `G` by matroids living on a *larger*
ground set (the vertices `α` embedded into `β` via `e`) yields, by comparison along `e`
(`Matroid.comap`), a representation by matroids on `α` with the same number of matroids. -/
theorem represents_comap {β : Type*} {G : SimpleGraph α} {N : Fin k → Matroid β} (e : α ↪ β)
    (hrep : ∀ S : Set α, G.IsIndepSet S ↔ ∀ i, (N i).Indep (e '' S)) :
    Represents G (fun i => (N i).comap e) := by
  intro S
  rw [hrep S]
  constructor
  · intro h i
    exact Matroid.comap_indep_iff.2 ⟨h i, e.injective.injOn⟩
  · intro h i
    exact (Matroid.comap_indep_iff.1 (h i)).1

/-- **G2.** The lower bound survives ground-set enlargement: if the independence system of
`G` is recovered from matroids on a bigger ground set, the degree bound still holds. -/
theorem le_card_of_represents_of_embedding {β : Type*} {G : SimpleGraph α}
    {N : Fin k → Matroid β} (e : α ↪ β)
    (hrep : ∀ S : Set α, G.IsIndepSet S ↔ ∀ i, (N i).Indep (e '' S))
    (v : α) (hv : G.IsIndepSet (G.neighborSet v)) (hfin : (G.neighborSet v).Finite) :
    hfin.toFinset.card ≤ k :=
  le_card_of_represents (represents_comap e hrep) v hv hfin

/-- **G2, minor form.** Since every minor of a matroid is again a matroid, the bound is
robust to representations by arbitrary minors `(P i / C i) ↾ R i` of matroids on a larger
ground set. -/
theorem le_card_of_represents_minor {β : Type*} {G : SimpleGraph α}
    {P : Fin k → Matroid β} {C R : Fin k → Set β} (e : α ↪ β)
    (hrep : ∀ S : Set α, G.IsIndepSet S ↔ ∀ i, (((P i).contract (C i)).restrict (R i)).Indep
      (e '' S))
    (v : α) (hv : G.IsIndepSet (G.neighborSet v)) (hfin : (G.neighborSet v).Finite) :
    hfin.toFinset.card ≤ k :=
  le_card_of_represents_of_embedding e hrep v hv hfin

/-! ## Axiom audit -/

section AxiomAudit

#print axioms MatroidRep.indep_singleton_of_represents
#print axioms MatroidRep.isNonloop_of_represents
#print axioms MatroidRep.isCircuit_pair_of_represents
#print axioms MatroidRep.not_dep_pair_of_dep_pair
#print axioms MatroidRep.le_card_of_represents
#print axioms MatroidRep.degree_le_of_represents
#print axioms MatroidRep.isIndepSet_neighborSet_of_cliqueFree
#print axioms MatroidRep.le_card_of_represents_cliqueFree
#print axioms MatroidRep.cliqueFree_three_of_isAcyclic
#print axioms MatroidRep.le_card_of_represents_isTree
#print axioms MatroidRep.represents_comap
#print axioms MatroidRep.le_card_of_represents_of_embedding
#print axioms MatroidRep.le_card_of_represents_minor

end AxiomAudit

end MatroidRep
