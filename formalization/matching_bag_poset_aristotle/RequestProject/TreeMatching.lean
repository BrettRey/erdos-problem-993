import RequestProject.KonigHall
import RequestProject.ForestLemmas
import RequestProject.Codes

/-!
# Trees with a maximum matching: bags, minimum covers and the constraint system

This file sets up §1 of `matching_bag_poset_reduction_2026-08-20.md`.

A `MatchingBag.TreeMatching V` is a forest `G` on the finite vertex type `V`, a proper
2-colouring `col` (which exists for every forest, see `exists_treeMatching`), and a maximum
matching `M` whose edges are *oriented* by the colouring: `e.1` is the `true`-endpoint
("`L`") and `e.2` the `false`-endpoint ("`R`").

* `TreeMatching.Idx` are the *matched bags* (the edges of `M`), `TreeMatching.Unm` the
  unmatched vertices (the *singleton bags*) and `TreeMatching.Bag` their sum.  The number of
  bags is `|V| - |M| = α(T)`, the independence number (`card_Bag`).
* `TreeMatching.rel i j` is the comparison `x_i ≤ x_j` imposed by a non-matching edge
  `R_i — L_j` (equation (1) of the note), and `ForcedZero`/`ForcedOne` are the unary
  constraints coming from edges at unmatched vertices.
* `TreeMatching.IsSol` is the resulting Boolean constraint system and `TreeMatching.Sol` its
  solution set; `coverOf`/`solOf` are mutually inverse bijections between `Sol` and the
  minimum vertex covers, and `Sol` is nonempty by König's theorem.
* `TreeMatching.treeCode` is the binary code of minimum covers on the bags; equivalently,
  by `maxIndepSets_eq_image_compl`, the code of maximum independent sets.
-/

open Finset SimpleGraph

namespace MatchingBag

attribute [local instance] Classical.propDecidable

variable {V : Type*} [Fintype V] [DecidableEq V]

/-- A finite forest with a proper 2-colouring and a maximum matching whose edges are oriented
from the `true` side to the `false` side. -/
structure TreeMatching (V : Type*) [Fintype V] [DecidableEq V] where
  /-- The graph, which is assumed to be a forest (a tree is a connected forest). -/
  G : SimpleGraph V
  /-- The graph has no cycles. -/
  acyclic : G.IsAcyclic
  /-- A proper 2-colouring, i.e. the bipartition. -/
  col : V → Bool
  /-- The colouring is proper. -/
  colProper : ∀ ⦃u v⦄, G.Adj u v → col u ≠ col v
  /-- The matching, as a set of ordered pairs. -/
  M : Finset (V × V)
  /-- `M` is a matching. -/
  isMatching : IsMatchingSet G M
  /-- Every matched edge is oriented from its `true`-endpoint to its `false`-endpoint. -/
  colLeft : ∀ e ∈ M, col e.1 = true
  /-- `M` is a maximum matching. -/
  maximum : ∀ S : Finset (V × V), IsMatchingSet G S → S.card ≤ M.card

/-- Every forest with a maximum matching gives rise to a `TreeMatching`: a forest is
bipartite, and the matched edges may be oriented by the bipartition without changing their
number. -/
theorem exists_treeMatching (G : SimpleGraph V) (hG : G.IsAcyclic) (M₀ : Finset (V × V))
    (hM₀ : IsMatchingSet G M₀) (hmax : ∀ S : Finset (V × V), IsMatchingSet G S → S.card ≤ M₀.card) :
    ∃ D : TreeMatching V, D.G = G ∧ D.M.card = M₀.card := by
  classical
  set c := hG.coloringTwo with hc
  set col : V → Bool := fun v => decide (c v = 0) with hcol
  have hcolProper : ∀ ⦃u v⦄, G.Adj u v → col u ≠ col v := by
    intro u v huv hcc
    have h : c u ≠ c v := c.valid huv
    have h2 : ∀ a b : Fin 2, ((a = 0) ↔ (b = 0)) → a = b := by decide
    simp only [hcol, decide_eq_decide] at hcc
    exact h (h2 _ _ hcc)
  set o : V × V → V × V := fun e => if col e.1 = true then e else (e.2, e.1) with ho
  have hcomp : ∀ e : V × V, ((o e).1 = e.1 ∧ (o e).2 = e.2) ∨ ((o e).1 = e.2 ∧ (o e).2 = e.1) := by
    intro e
    by_cases h : col e.1 = true
    · exact Or.inl ⟨by simp [ho, h], by simp [ho, h]⟩
    · exact Or.inr ⟨by simp [ho, h], by simp [ho, h]⟩
  have hadj : ∀ e ∈ M₀, G.Adj (o e).1 (o e).2 := by
    intro e he
    rcases hcomp e with ⟨h1, h2⟩ | ⟨h1, h2⟩
    · rw [h1, h2]; exact hM₀.1 e he
    · rw [h1, h2]; exact (hM₀.1 e he).symm
  have hleft : ∀ e ∈ M₀, col (o e).1 = true := by
    intro e he
    by_cases h : col e.1 = true
    · simp [ho, h]
    · simp only [Bool.not_eq_true] at h
      have hne := hcolProper (hM₀.1 e he)
      have h2t : col e.2 = true := by
        cases hh : col e.2
        · exact absurd (by rw [h, hh]) hne
        · rfl
      simp [ho, h, h2t]
  have hinj : Set.InjOn o (M₀ : Set (V × V)) := by
    intro e he e' he' hoe
    by_contra hne
    obtain ⟨d1, d2, d3, d4⟩ := hM₀.2 e he e' he' hne
    have h1 := congrArg Prod.fst hoe
    have h2 := congrArg Prod.snd hoe
    rcases hcomp e with ⟨a1, a2⟩ | ⟨a1, a2⟩ <;> rcases hcomp e' with ⟨b1, b2⟩ | ⟨b1, b2⟩ <;>
      rw [a1] at h1 <;> rw [b1] at h1 <;> rw [a2] at h2 <;> rw [b2] at h2
    · exact d1 h1
    · exact d2 h1
    · exact d3 h1
    · exact d4 h1
  have hmatch : IsMatchingSet G (M₀.image o) := by
    constructor
    · intro f hf
      obtain ⟨e, he, rfl⟩ := Finset.mem_image.1 hf
      exact hadj e he
    · intro f hf f' hf' hne
      obtain ⟨e, he, rfl⟩ := Finset.mem_image.1 hf
      obtain ⟨e', he', rfl⟩ := Finset.mem_image.1 hf'
      have hee : e ≠ e' := by rintro rfl; exact hne rfl
      obtain ⟨d1, d2, d3, d4⟩ := hM₀.2 e he e' he' hee
      rcases hcomp e with ⟨a1, a2⟩ | ⟨a1, a2⟩ <;> rcases hcomp e' with ⟨b1, b2⟩ | ⟨b1, b2⟩ <;>
        rw [a1, a2, b1, b2] <;>
        exact ⟨by assumption, by assumption, by assumption, by assumption⟩
  have hcard : (M₀.image o).card = M₀.card := Finset.card_image_of_injOn hinj
  have hcolLeft : ∀ f ∈ M₀.image o, col f.1 = true := by
    intro f hf
    obtain ⟨e, he, rfl⟩ := Finset.mem_image.1 hf
    exact hleft e he
  have hmaximum : ∀ S : Finset (V × V), IsMatchingSet G S → S.card ≤ (M₀.image o).card := by
    intro S hS
    rw [hcard]
    exact hmax S hS
  exact ⟨⟨G, hG, col, hcolProper, M₀.image o, hmatch, hcolLeft, hmaximum⟩, rfl, hcard⟩

namespace TreeMatching

variable (D : TreeMatching V)

/-! ### Bags -/

/-- The matched bags: the edges of the matching. -/
abbrev Idx : Type _ := {e : V × V // e ∈ D.M}

/-- The `L`-endpoint (the `true`-coloured one) of a matched bag. -/
def Lv (i : D.Idx) : V := i.1.1

/-- The `R`-endpoint (the `false`-coloured one) of a matched bag. -/
def Rv (i : D.Idx) : V := i.1.2

/-- The vertices covered by the matching. -/
def matchedVerts : Finset V := D.M.image Prod.fst ∪ D.M.image Prod.snd

/-- The unmatched vertices, i.e. the singleton bags. -/
def unmatchedVerts : Finset V := Finset.univ \ D.matchedVerts

/-- The singleton bags. -/
abbrev Unm : Type _ := {v : V // v ∈ D.unmatchedVerts}

/-- All bags: the matched edges together with the unmatched vertices. -/
abbrev Bag : Type _ := D.Idx ⊕ D.Unm

variable {D}

@[simp] lemma adj_Lv_Rv (i : D.Idx) : D.G.Adj (D.Lv i) (D.Rv i) := D.isMatching.1 i.1 i.2

lemma Lv_ne_Rv (i : D.Idx) : D.Lv i ≠ D.Rv i := (adj_Lv_Rv i).ne

@[simp] lemma col_Lv (i : D.Idx) : D.col (D.Lv i) = true := D.colLeft i.1 i.2

@[simp] lemma col_Rv (i : D.Idx) : D.col (D.Rv i) = false := by
  have := D.colProper (adj_Lv_Rv i)
  rw [col_Lv] at this
  simpa using (Bool.eq_false_iff.2 (Ne.symm this))

/-- Distinct bags are disjoint. -/
lemma bags_disjoint {i j : D.Idx} (h : i ≠ j) :
    D.Lv i ≠ D.Lv j ∧ D.Lv i ≠ D.Rv j ∧ D.Rv i ≠ D.Lv j ∧ D.Rv i ≠ D.Rv j :=
  D.isMatching.2 i.1 i.2 j.1 j.2 (fun hh => h (Subtype.ext hh))

/-- The bags of a `TreeMatching` are disjoint edges of the forest. -/
lemma disjointEdges : DisjointEdges D.G D.Lv D.Rv where
  adj := adj_Lv_Rv
  disj := fun _ _ h => bags_disjoint h

lemma mem_matchedVerts_iff {v : V} :
    v ∈ D.matchedVerts ↔ ∃ i : D.Idx, v = D.Lv i ∨ v = D.Rv i := by
  simp only [matchedVerts, Finset.mem_union, Finset.mem_image]
  constructor
  · rintro (⟨e, he, rfl⟩ | ⟨e, he, rfl⟩)
    · exact ⟨⟨e, he⟩, Or.inl rfl⟩
    · exact ⟨⟨e, he⟩, Or.inr rfl⟩
  · rintro ⟨i, rfl | rfl⟩
    · exact Or.inl ⟨i.1, i.2, rfl⟩
    · exact Or.inr ⟨i.1, i.2, rfl⟩

lemma mem_unmatchedVerts_iff {v : V} :
    v ∈ D.unmatchedVerts ↔ ∀ i : D.Idx, v ≠ D.Lv i ∧ v ≠ D.Rv i := by
  simp only [unmatchedVerts, Finset.mem_sdiff, Finset.mem_univ, true_and,
    mem_matchedVerts_iff, not_exists, not_or]

/-- Every vertex is either the `L`-endpoint of a bag, or the `R`-endpoint of a bag, or
unmatched. -/
lemma vertex_cases (D : TreeMatching V) (v : V) :
    (∃ i : D.Idx, v = D.Lv i) ∨ (∃ i : D.Idx, v = D.Rv i) ∨ v ∈ D.unmatchedVerts := by
  by_cases h : v ∈ D.matchedVerts
  · obtain ⟨i, hi | hi⟩ := mem_matchedVerts_iff.1 h
    · exact Or.inl ⟨i, hi⟩
    · exact Or.inr (Or.inl ⟨i, hi⟩)
  · refine Or.inr (Or.inr ?_)
    simp only [unmatchedVerts, Finset.mem_sdiff, Finset.mem_univ, true_and]
    exact h

/-- Maximality of the matching forbids edges between two unmatched vertices. -/
lemma not_adj_of_unmatched {u v : V} (hu : u ∈ D.unmatchedVerts) (hv : v ∈ D.unmatchedVerts) :
    ¬ D.G.Adj u v := by
  intro h
  have hu' := mem_unmatchedVerts_iff.1 hu
  have hv' := mem_unmatchedVerts_iff.1 hv
  have hnotmem : (u, v) ∉ D.M := fun hmem => (hu' ⟨(u, v), hmem⟩).1 rfl
  have hS : IsMatchingSet D.G (insert (u, v) D.M) := by
    constructor
    · intro e he
      rcases Finset.mem_insert.1 he with rfl | he
      · exact h
      · exact D.isMatching.1 e he
    · intro e he e' he' hne
      rcases Finset.mem_insert.1 he with rfl | he <;>
        rcases Finset.mem_insert.1 he' with rfl | he'
      · exact absurd rfl hne
      · exact ⟨(hu' ⟨e', he'⟩).1, (hu' ⟨e', he'⟩).2, (hv' ⟨e', he'⟩).1, (hv' ⟨e', he'⟩).2⟩
      · exact ⟨fun hh => (hu' ⟨e, he⟩).1 hh.symm, fun hh => (hv' ⟨e, he⟩).1 hh.symm,
          fun hh => (hu' ⟨e, he⟩).2 hh.symm, fun hh => (hv' ⟨e, he⟩).2 hh.symm⟩
      · exact D.isMatching.2 e he e' he' hne
  have hle := D.maximum _ hS
  rw [Finset.card_insert_of_notMem hnotmem] at hle
  omega

variable (D)

lemma card_matchedVerts : D.matchedVerts.card = 2 * D.M.card := by
  classical
  have hne : ∀ e ∈ D.M, e.1 ≠ e.2 := fun e he => (D.isMatching.1 e he).ne
  have hfst : (D.M.image Prod.fst).card = D.M.card := by
    refine Finset.card_image_of_injOn ?_
    intro e he e' he' hh
    by_contra hne'
    exact (D.isMatching.2 e he e' he' hne').1 hh
  have hsnd : (D.M.image Prod.snd).card = D.M.card := by
    refine Finset.card_image_of_injOn ?_
    intro e he e' he' hh
    by_contra hne'
    exact (D.isMatching.2 e he e' he' hne').2.2.2 hh
  have hdisj : Disjoint (D.M.image Prod.fst) (D.M.image Prod.snd) := by
    rw [Finset.disjoint_left]
    rintro v hv hv'
    obtain ⟨e, he, rfl⟩ := Finset.mem_image.1 hv
    obtain ⟨e', he', hh⟩ := Finset.mem_image.1 hv'
    by_cases hee : e = e'
    · subst hee; exact hne e he hh.symm
    · exact (D.isMatching.2 e he e' he' hee).2.1 hh.symm
  rw [matchedVerts, Finset.card_union_of_disjoint hdisj, hfst, hsnd]
  ring

lemma two_mul_card_M_le : 2 * D.M.card ≤ Fintype.card V := by
  have hle : D.matchedVerts.card ≤ Fintype.card V := by
    simpa [Finset.card_univ] using Finset.card_le_card (Finset.subset_univ D.matchedVerts)
  rw [card_matchedVerts] at hle
  exact hle

lemma card_unmatchedVerts : D.unmatchedVerts.card = Fintype.card V - 2 * D.M.card := by
  rw [unmatchedVerts, Finset.card_sdiff_of_subset (Finset.subset_univ _), Finset.card_univ,
    card_matchedVerts]

lemma card_Idx : Fintype.card D.Idx = D.M.card := Fintype.card_coe _

/-- The number of bags is the independence number `α(T) = |V| - |M|`. -/
theorem card_Bag : Fintype.card D.Bag = Fintype.card V - D.M.card := by
  have hsum : Fintype.card D.Bag = Fintype.card D.Idx + Fintype.card D.Unm :=
    Fintype.card_sum
  have h1 : Fintype.card D.Idx = D.M.card := card_Idx D
  have h2 : Fintype.card D.Unm = D.unmatchedVerts.card := Fintype.card_coe _
  have h3 := card_unmatchedVerts D
  have h4 := two_mul_card_M_le D
  omega

/-! ### The constraint system -/

variable {D}

/-- The comparison `x_i ≤ x_j` forced by a non-matching edge `R_i — L_j`; this is
equation (1) of the note. -/
def rel (D : TreeMatching V) (i j : D.Idx) : Prop := i ≠ j ∧ D.G.Adj (D.Rv i) (D.Lv j)

/-- `x_i = 0` is forced: some unmatched vertex is adjacent to `R_i`. -/
def ForcedZero (D : TreeMatching V) (i : D.Idx) : Prop := ∃ v ∈ D.unmatchedVerts, D.G.Adj v (D.Rv i)

/-- `x_i = 1` is forced: some unmatched vertex is adjacent to `L_i`. -/
def ForcedOne (D : TreeMatching V) (i : D.Idx) : Prop := ∃ v ∈ D.unmatchedVerts, D.G.Adj v (D.Lv i)

/-- The Boolean constraint system attached to the tree and the matching. -/
def IsSol (D : TreeMatching V) (x : D.Idx → Bool) : Prop :=
  (∀ i j, D.rel i j → x i = true → x j = true) ∧
    (∀ i, D.ForcedZero i → x i = false) ∧ (∀ i, D.ForcedOne i → x i = true)

variable (D)

/-- The solution set of the constraint system. -/
noncomputable def Sol : Finset (D.Idx → Bool) := Finset.univ.filter D.IsSol

variable {D}

@[simp] lemma mem_Sol {x : D.Idx → Bool} : x ∈ D.Sol ↔ D.IsSol x := by simp [Sol]

/-! ### Minimum covers -/

variable (D)

/-- The minimum vertex covers: the covers of size `|M|`. -/
noncomputable def minCovers : Finset (Finset V) :=
  Finset.univ.filter (fun C => IsVertexCover D.G C ∧ C.card = D.M.card)

/-- The maximum independent sets: the independent sets of size `|V| - |M|`. -/
noncomputable def maxIndepSets : Finset (Finset V) :=
  Finset.univ.filter
    (fun I => (∀ u ∈ I, ∀ v ∈ I, ¬ D.G.Adj u v) ∧ I.card = Fintype.card V - D.M.card)

variable {D}

@[simp] lemma mem_minCovers {C : Finset V} :
    C ∈ D.minCovers ↔ IsVertexCover D.G C ∧ C.card = D.M.card := by simp [minCovers]

/-- **König's theorem** for the tree: a cover of size `|M|` exists. -/
theorem exists_minCover (D : TreeMatching V) : ∃ C : Finset V, C ∈ D.minCovers := by
  obtain ⟨C, hC, hcard⟩ := konig D.colProper D.isMatching D.maximum
  exact ⟨C, mem_minCovers.2 ⟨hC, hcard⟩⟩

/-- The maximum independent sets are exactly the complements of the minimum covers. -/
theorem maxIndepSets_eq_image_compl (D : TreeMatching V) :
    D.maxIndepSets = D.minCovers.image (fun C => Cᶜ) := by
  classical
  have h2 := two_mul_card_M_le D
  ext I
  simp only [maxIndepSets, minCovers, Finset.mem_filter, Finset.mem_univ, true_and,
    Finset.mem_image]
  constructor
  · rintro ⟨hind, hcard⟩
    refine ⟨Iᶜ, ⟨?_, ?_⟩, by simp⟩
    · intro u v huv
      by_contra hcon
      push_neg at hcon
      obtain ⟨hu, hv⟩ := hcon
      simp only [Finset.mem_compl, not_not] at hu hv
      exact hind u hu v hv huv
    · rw [Finset.card_compl, hcard]
      omega
  · rintro ⟨C, ⟨hcov, hcard⟩, rfl⟩
    refine ⟨?_, ?_⟩
    · intro u hu v hv huv
      simp only [Finset.mem_compl] at hu hv
      rcases hcov huv with h | h
      · exact hu h
      · exact hv h
    · rw [Finset.card_compl, hcard]

/-- A minimum cover contains exactly one endpoint of each matched bag, and no unmatched
vertex.  (This is the pigeonhole step `cover_card_eq_matching_card` of §1.) -/
theorem minCover_structure {C : Finset V} (hC : C ∈ D.minCovers) :
    (∀ i : D.Idx, ¬ (D.Lv i ∈ C ∧ D.Rv i ∈ C)) ∧
      (∀ i : D.Idx, D.Lv i ∈ C ∨ D.Rv i ∈ C) ∧
      (∀ v ∈ D.unmatchedVerts, v ∉ C) := by
  obtain ⟨hcov, hcard⟩ := mem_minCovers.1 hC
  have h := cover_card_eq_matching_card D.M C (fun e he => (D.isMatching.1 e he).ne)
    (fun e he e' he' hne => D.isMatching.2 e he e' he' hne)
    (fun e he => hcov (D.isMatching.1 e he)) hcard
  refine ⟨fun i => h.1 i.1 i.2, fun i => hcov (adj_Lv_Rv i), ?_⟩
  intro v hv hvC
  obtain ⟨e, he, hve⟩ := h.2 v hvC
  have hvi := mem_unmatchedVerts_iff.1 hv ⟨e, he⟩
  rcases hve with hh | hh
  · exact hvi.1 hh
  · exact hvi.2 hh

/-! ### Covers and solutions -/

/-- The cover determined by a solution of the constraint system: take `L_i` when `x_i = 1`
and `R_i` when `x_i = 0`. -/
noncomputable def coverOf (D : TreeMatching V) (x : D.Idx → Bool) : Finset V :=
  Finset.univ.filter
    (fun v => ∃ i : D.Idx, (x i = true ∧ v = D.Lv i) ∨ (x i = false ∧ v = D.Rv i))

/-- The Boolean assignment determined by a cover: `x_i = 1` exactly when `L_i` is in the
cover. -/
noncomputable def solOf (D : TreeMatching V) (C : Finset V) : D.Idx → Bool :=
  fun i => decide (D.Lv i ∈ C)

lemma mem_coverOf {x : D.Idx → Bool} {v : V} :
    v ∈ D.coverOf x ↔ ∃ i : D.Idx, (x i = true ∧ v = D.Lv i) ∨ (x i = false ∧ v = D.Rv i) := by
  simp [coverOf]

lemma Lv_mem_coverOf {x : D.Idx → Bool} {i : D.Idx} : D.Lv i ∈ D.coverOf x ↔ x i = true := by
  rw [mem_coverOf]
  constructor
  · rintro ⟨j, ⟨hj, hij⟩ | ⟨hj, hij⟩⟩
    · rcases eq_or_ne i j with rfl | hne
      · exact hj
      · exact absurd hij (bags_disjoint hne).1
    · rcases eq_or_ne i j with rfl | hne
      · exact absurd hij (Lv_ne_Rv i)
      · exact absurd hij (bags_disjoint hne).2.1
  · intro h
    exact ⟨i, Or.inl ⟨h, rfl⟩⟩

lemma Rv_mem_coverOf {x : D.Idx → Bool} {i : D.Idx} : D.Rv i ∈ D.coverOf x ↔ x i = false := by
  rw [mem_coverOf]
  constructor
  · rintro ⟨j, ⟨hj, hij⟩ | ⟨hj, hij⟩⟩
    · rcases eq_or_ne i j with rfl | hne
      · exact absurd hij.symm (Lv_ne_Rv i)
      · exact absurd hij (bags_disjoint hne).2.2.1
    · rcases eq_or_ne i j with rfl | hne
      · exact hj
      · exact absurd hij (bags_disjoint hne).2.2.2
  · intro h
    exact ⟨i, Or.inr ⟨h, rfl⟩⟩

lemma unmatched_not_mem_coverOf {x : D.Idx → Bool} {v : V} (hv : v ∈ D.unmatchedVerts) :
    v ∉ D.coverOf x := by
  rw [mem_coverOf]
  rintro ⟨i, ⟨_, hh⟩ | ⟨_, hh⟩⟩
  · exact (mem_unmatchedVerts_iff.1 hv i).1 hh
  · exact (mem_unmatchedVerts_iff.1 hv i).2 hh

lemma coverOf_eq_image (x : D.Idx → Bool) :
    D.coverOf x = Finset.univ.image (fun i : D.Idx => if x i = true then D.Lv i else D.Rv i) := by
  ext v
  simp only [mem_coverOf, Finset.mem_image, Finset.mem_univ, true_and]
  constructor
  · rintro ⟨i, ⟨hx, rfl⟩ | ⟨hx, rfl⟩⟩
    · exact ⟨i, by simp [hx]⟩
    · exact ⟨i, by simp [hx]⟩
  · rintro ⟨i, rfl⟩
    by_cases hx : x i = true
    · exact ⟨i, Or.inl ⟨hx, by simp [hx]⟩⟩
    · simp only [Bool.not_eq_true] at hx
      exact ⟨i, Or.inr ⟨hx, by simp [hx]⟩⟩

lemma card_coverOf (x : D.Idx → Bool) : (D.coverOf x).card = D.M.card := by
  classical
  rw [coverOf_eq_image, Finset.card_image_of_injective _ ?inj, Finset.card_univ, card_Idx]
  case inj =>
    intro i j hij
    by_contra hne
    obtain ⟨d1, d2, d3, d4⟩ := bags_disjoint hne
    simp only [] at hij
    rcases Bool.eq_false_or_eq_true (x i) with hi | hi <;>
      rcases Bool.eq_false_or_eq_true (x j) with hj | hj
    · rw [if_pos hi, if_pos hj] at hij
      exact d1 hij
    · rw [if_pos hi, if_neg (by simp [hj])] at hij
      exact d2 hij
    · rw [if_neg (by simp [hi]), if_pos hj] at hij
      exact d3 hij
    · rw [if_neg (by simp [hi]), if_neg (by simp [hj])] at hij
      exact d4 hij

theorem coverOf_mem_minCovers {x : D.Idx → Bool} (hx : D.IsSol x) : D.coverOf x ∈ D.minCovers := by
  obtain ⟨hrel, hz, ho⟩ := hx
  refine mem_minCovers.2 ⟨?_, card_coverOf x⟩
  intro u v huv
  rcases vertex_cases D u with ⟨i, rfl⟩ | ⟨i, rfl⟩ | hu <;>
    rcases vertex_cases D v with ⟨j, rfl⟩ | ⟨j, rfl⟩ | hv
  -- Lv i — Lv j : impossible, both colours true
  · exact absurd (by rw [col_Lv, col_Lv]) (D.colProper huv)
  -- Lv i — Rv j
  · rcases eq_or_ne i j with rfl | hne
    · by_cases hxi : x i = true
      · exact Or.inl (Lv_mem_coverOf.2 hxi)
      · simp only [Bool.not_eq_true] at hxi
        exact Or.inr (Rv_mem_coverOf.2 hxi)
    · by_cases hxj : x j = true
      · have : x i = true := hrel j i ⟨Ne.symm hne, huv.symm⟩ hxj
        exact Or.inl (Lv_mem_coverOf.2 this)
      · simp only [Bool.not_eq_true] at hxj
        exact Or.inr (Rv_mem_coverOf.2 hxj)
  -- Lv i — unmatched
  · exact Or.inl (Lv_mem_coverOf.2 (ho i ⟨v, hv, huv.symm⟩))
  -- Rv i — Lv j
  · rcases eq_or_ne i j with rfl | hne
    · by_cases hxi : x i = true
      · exact Or.inr (Lv_mem_coverOf.2 hxi)
      · simp only [Bool.not_eq_true] at hxi
        exact Or.inl (Rv_mem_coverOf.2 hxi)
    · by_cases hxi : x i = true
      · have : x j = true := hrel i j ⟨hne, huv⟩ hxi
        exact Or.inr (Lv_mem_coverOf.2 this)
      · simp only [Bool.not_eq_true] at hxi
        exact Or.inl (Rv_mem_coverOf.2 hxi)
  -- Rv i — Rv j : impossible
  · exact absurd (by rw [col_Rv, col_Rv]) (D.colProper huv)
  -- Rv i — unmatched
  · exact Or.inl (Rv_mem_coverOf.2 (hz i ⟨v, hv, huv.symm⟩))
  -- unmatched — Lv j
  · exact Or.inr (Lv_mem_coverOf.2 (ho j ⟨u, hu, huv⟩))
  -- unmatched — Rv j
  · exact Or.inr (Rv_mem_coverOf.2 (hz j ⟨u, hu, huv⟩))
  -- unmatched — unmatched : impossible
  · exact absurd huv (not_adj_of_unmatched hu hv)

theorem solOf_mem_Sol {C : Finset V} (hC : C ∈ D.minCovers) : D.solOf C ∈ D.Sol := by
  obtain ⟨hcov, hcard⟩ := mem_minCovers.1 hC
  obtain ⟨hone, hcovbag, hunm⟩ := minCover_structure hC
  refine mem_Sol.2 ⟨?_, ?_, ?_⟩
  · intro i j hij hxi
    simp only [solOf, decide_eq_true_eq] at hxi ⊢
    have hRi : D.Rv i ∉ C := fun hh => hone i ⟨hxi, hh⟩
    rcases hcov hij.2 with hh | hh
    · exact absurd hh hRi
    · exact hh
  · intro i ⟨v, hv, hadj⟩
    simp only [solOf, decide_eq_false_iff_not]
    intro hLi
    have hRi : D.Rv i ∈ C := by
      rcases hcov hadj with hh | hh
      · exact absurd hh (hunm v hv)
      · exact hh
    exact hone i ⟨hLi, hRi⟩
  · intro i ⟨v, hv, hadj⟩
    simp only [solOf, decide_eq_true_eq]
    rcases hcov hadj with hh | hh
    · exact absurd hh (hunm v hv)
    · exact hh

/-- **The inequality system, equation (1)**, read directly on covers: for a non-matching
edge `R_i — L_j`, a minimum cover containing `L_i` also contains `L_j`, i.e. `x_i ≤ x_j`. -/
theorem minCover_rel {C : Finset V} (hC : C ∈ D.minCovers) {i j : D.Idx} (hij : D.rel i j)
    (hi : D.Lv i ∈ C) : D.Lv j ∈ C := by
  have h := (mem_Sol.1 (solOf_mem_Sol hC)).1 i j hij (by simpa [solOf] using hi)
  simpa [solOf] using h

/-- The unary constraint attached to an unmatched vertex adjacent to `R_i`: `x_i = 0`. -/
theorem minCover_forcedZero {C : Finset V} (hC : C ∈ D.minCovers) {i : D.Idx}
    (h : D.ForcedZero i) : D.Lv i ∉ C := by
  have h' := (mem_Sol.1 (solOf_mem_Sol hC)).2.1 i h
  simpa [solOf] using h'

/-- The unary constraint attached to an unmatched vertex adjacent to `L_i`: `x_i = 1`. -/
theorem minCover_forcedOne {C : Finset V} (hC : C ∈ D.minCovers) {i : D.Idx}
    (h : D.ForcedOne i) : D.Lv i ∈ C := by
  have h' := (mem_Sol.1 (solOf_mem_Sol hC)).2.2 i h
  simpa [solOf] using h'

theorem solOf_coverOf {x : D.Idx → Bool} : D.solOf (D.coverOf x) = x := by
  funext i
  simp only [solOf, Lv_mem_coverOf]
  cases x i <;> simp

theorem coverOf_solOf {C : Finset V} (hC : C ∈ D.minCovers) : D.coverOf (D.solOf C) = C := by
  obtain ⟨hone, hcovbag, hunm⟩ := minCover_structure hC
  ext v
  constructor
  · intro hv
    rcases mem_coverOf.1 hv with ⟨i, ⟨hx, rfl⟩ | ⟨hx, rfl⟩⟩
    · simpa [solOf] using hx
    · simp only [solOf, decide_eq_false_iff_not] at hx
      rcases hcovbag i with hh | hh
      · exact absurd hh hx
      · exact hh
  · intro hv
    rcases vertex_cases D v with ⟨i, rfl⟩ | ⟨i, rfl⟩ | hvu
    · exact Lv_mem_coverOf.2 (by simpa [solOf] using hv)
    · refine Rv_mem_coverOf.2 ?_
      simp only [solOf, decide_eq_false_iff_not]
      exact fun hh => hone i ⟨hh, hv⟩
    · exact absurd hv (hunm v hvu)

/-- The minimum covers are in explicit bijection with the solutions of the constraint
system, via the mutually inverse maps `solOf` and `coverOf`. -/
theorem minCovers_image_solOf (D : TreeMatching V) : D.minCovers.image D.solOf = D.Sol := by
  ext x
  constructor
  · rintro hx
    obtain ⟨C, hC, rfl⟩ := Finset.mem_image.1 hx
    exact solOf_mem_Sol hC
  · intro hx
    exact Finset.mem_image.2 ⟨D.coverOf x, coverOf_mem_minCovers (mem_Sol.1 hx), solOf_coverOf⟩

/-- The solution set of the constraint system is nonempty: this is the *consistency* of the
unary boundary values, and it comes from maximality of the matching (through König). -/
theorem Sol_nonempty (D : TreeMatching V) : D.Sol.Nonempty := by
  obtain ⟨C, hC⟩ := exists_minCover D
  exact ⟨D.solOf C, solOf_mem_Sol hC⟩

/-! ### The code of the tree -/

/-- The word on the bags attached to a cover: on a matched bag the bit says whether the
`L`-endpoint is in the cover, on a singleton bag whether the vertex is in the cover. -/
noncomputable def coverWord (D : TreeMatching V) (C : Finset V) : D.Bag → Bool :=
  Sum.elim (fun i => decide (D.Lv i ∈ C)) (fun v => decide (v.1 ∈ C))

/-- The binary code of the tree: the minimum covers (equivalently, by
`maxIndepSets_eq_image_compl`, the maximum independent sets), read on the bags. -/
noncomputable def treeCode (D : TreeMatching V) : Finset (D.Bag → Bool) :=
  D.minCovers.image D.coverWord

/-- An independent set contains at most one endpoint of each matched bag. -/
lemma indep_not_both {I : Finset V} (hI : ∀ u ∈ I, ∀ v ∈ I, ¬ D.G.Adj u v) (i : D.Idx) :
    ¬ (D.Lv i ∈ I ∧ D.Rv i ∈ I) := fun h => hI _ h.1 _ h.2 (adj_Lv_Rv i)

/-- Complementing a cover flips every bit of its word. -/
lemma coverWord_compl (D : TreeMatching V) (C : Finset V) (k : D.Bag) :
    D.coverWord Cᶜ k = !(D.coverWord C k) := by
  cases k <;> simp [coverWord]

/-- The code of the tree is nonempty (König). -/
theorem treeCode_nonempty (D : TreeMatching V) : D.treeCode.Nonempty := by
  obtain ⟨C, hC⟩ := exists_minCover D
  exact ⟨D.coverWord C, Finset.mem_image_of_mem _ hC⟩

end TreeMatching

end MatchingBag
