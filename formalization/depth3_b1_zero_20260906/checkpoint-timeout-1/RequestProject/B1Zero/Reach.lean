import RequestProject.B1Zero.Convolution

/-!
# Connectivity inside a finite vertex set

For the corona decomposition we peel one connected component of an induced subgraph at a
time.  This file develops the necessary bookkeeping in `Finset` form:

* `ReachIn G X a b` — `b` is reachable from `a` by a walk all of whose vertices lie in `X`;
* `compIn G X v` — the connected component of `v` in the subgraph induced on `X`;
* `IsTransversal G X Y` — `Y` meets every component of `X` in exactly one vertex;
* `transversal_erase` and `transversal_children`: the two transversals produced when a
  component and then its root are peeled off;
* `exists_transversal_subset`: a set meeting every component contains a transversal.

The only place where acyclicity enters is `not_reachIn_of_adj`: two distinct neighbours of a
vertex `v` of a forest are not joined by a walk avoiding `v`.
-/

open Finset SimpleGraph

namespace MatchingBag

namespace B1Zero

attribute [local instance 1] Classical.propDecidable

variable {V : Type*} [Fintype V] [DecidableEq V] {G : SimpleGraph V} {X : Finset V}

/-- One step of a walk confined to `X`. -/
def StepIn (G : SimpleGraph V) (X : Finset V) (a b : V) : Prop := a ∈ X ∧ b ∈ X ∧ G.Adj a b

/-- `b` is reachable from `a` by a walk inside `X`. -/
def ReachIn (G : SimpleGraph V) (X : Finset V) (a b : V) : Prop :=
  Relation.ReflTransGen (StepIn G X) a b

namespace ReachIn

@[refl] lemma refl {a : V} : ReachIn G X a a := Relation.ReflTransGen.refl

lemma single {a b : V} (h : StepIn G X a b) : ReachIn G X a b := Relation.ReflTransGen.single h

lemma trans {a b c : V} (h : ReachIn G X a b) (h' : ReachIn G X b c) : ReachIn G X a c :=
  Relation.ReflTransGen.trans h h'

lemma symm {a b : V} (h : ReachIn G X a b) : ReachIn G X b a := by
  refine Relation.ReflTransGen.symmetric (fun x y hxy => ?_) h
  exact ⟨hxy.2.1, hxy.1, hxy.2.2.symm⟩

lemma mem_right {a b : V} (h : ReachIn G X a b) (hne : b ≠ a) : b ∈ X := by
  rcases Relation.ReflTransGen.cases_tail h with rfl | ⟨c, -, hstep⟩
  · exact absurd rfl hne
  · exact hstep.2.1

lemma mono {X' : Finset V} (hXX : X ⊆ X') {a b : V} (h : ReachIn G X a b) : ReachIn G X' a b := by
  induction h with
  | refl => exact ReachIn.refl
  | tail _ hstep ih => exact ih.trans (ReachIn.single ⟨hXX hstep.1, hXX hstep.2.1, hstep.2.2⟩)

/-- Reachability inside a closed subset is reachability inside the ambient set. -/
lemma restrict {Cl : Finset V} (hclosed : ∀ u ∈ Cl, ∀ v ∈ X, G.Adj u v → v ∈ Cl)
    {a b : V} (ha : a ∈ Cl) (h : ReachIn G X a b) : ReachIn G Cl a b := by
  induction h with
  | refl => exact ReachIn.refl
  | @tail c d _ hstep ih =>
      have hc : c ∈ Cl := by
        by_cases hca : c = a
        · exact hca ▸ ha
        · exact ih.mem_right hca
      exact ih.trans (ReachIn.single ⟨hc, hclosed c hc d hstep.2.1 hstep.2.2, hstep.2.2⟩)

end ReachIn

/-- A walk realising reachability inside `X`. -/
lemma exists_walk_of_reachIn {a b : V} (ha : a ∈ X) (h : ReachIn G X a b) :
    ∃ w : G.Walk a b, ∀ x ∈ w.support, x ∈ X := by
  induction h with
  | refl => exact ⟨SimpleGraph.Walk.nil, by simpa using ha⟩
  | tail _ hstep ih =>
      obtain ⟨w, hw⟩ := ih
      refine ⟨w.append (SimpleGraph.Walk.cons hstep.2.2 SimpleGraph.Walk.nil), ?_⟩
      intro x hx
      rw [SimpleGraph.Walk.support_append] at hx
      simp only [List.mem_append, SimpleGraph.Walk.support_cons, SimpleGraph.Walk.support_nil,
        List.tail_cons, List.mem_singleton] at hx
      rcases hx with hx | rfl
      · exact hw x hx
      · exact hstep.2.1

/-- **Acyclicity.**  Two distinct neighbours of `v` are not joined inside a set avoiding `v`. -/
theorem not_reachIn_of_adj (hac : G.IsAcyclic) {v c c' : V} (hvX : v ∉ X)
    (hvc : G.Adj v c) (hvc' : G.Adj v c') (hcX : c ∈ X) (hne : c ≠ c') :
    ¬ ReachIn G X c c' := by
  intro h
  obtain ⟨w, hw⟩ := exists_walk_of_reachIn hcX h
  have hbridge : s(v, c) ∈ (SimpleGraph.Walk.cons hvc' w.reverse).edges :=
    (SimpleGraph.isBridge_iff_adj_and_forall_walk_mem_edges.1
      ((SimpleGraph.isAcyclic_iff_forall_adj_isBridge.1 hac) hvc)).2 _
  rw [SimpleGraph.Walk.edges_cons, List.mem_cons] at hbridge
  rcases hbridge with he | he
  · rw [Sym2.eq_iff] at he
    rcases he with ⟨-, h2⟩ | ⟨-, h2⟩
    · exact hne h2
    · exact hvc.ne h2.symm
  · have : v ∈ w.reverse.support := SimpleGraph.Walk.fst_mem_support_of_mem_edges _ he
    rw [SimpleGraph.Walk.support_reverse, List.mem_reverse] at this
    exact hvX (hw v this)

/-- Every vertex reachable from `v` other than `v` is reachable from a neighbour of `v`
without using `v`. -/
theorem exists_adj_reachIn_erase {v x : V} (h : ReachIn G X v x) (hne : x ≠ v) :
    ∃ c, G.Adj v c ∧ c ∈ X ∧ ReachIn G (X.erase v) c x := by
  have key : ∀ {b : V}, ReachIn G X v b →
      (b = v ∨ ∃ c, G.Adj v c ∧ c ∈ X ∧ ReachIn G (X.erase v) c b) := by
    intro b hb
    induction hb with
    | refl => exact Or.inl rfl
    | @tail y d _ hstep ih =>
        by_cases hdv : d = v
        · exact Or.inl hdv
        · have hd : d ∈ X.erase v := Finset.mem_erase.2 ⟨hdv, hstep.2.1⟩
          by_cases hyv : y = v
          · subst hyv
            exact Or.inr ⟨d, hstep.2.2, hstep.2.1, ReachIn.refl⟩
          · rcases ih with rfl | ⟨c, hadj, hcX, hreach⟩
            · exact absurd rfl hyv
            · exact Or.inr ⟨c, hadj, hcX,
                hreach.trans (ReachIn.single
                  ⟨Finset.mem_erase.2 ⟨hyv, hstep.1⟩, hd, hstep.2.2⟩)⟩
  rcases key h with hx | hx
  · exact absurd hx hne
  · exact hx

/-! ### Components -/

/-- The connected component of `v` inside `X`. -/
noncomputable def compIn (G : SimpleGraph V) (X : Finset V) (v : V) : Finset V :=
  X.filter (fun x => ReachIn G X v x)

lemma mem_compIn {v x : V} : x ∈ compIn G X v ↔ x ∈ X ∧ ReachIn G X v x := by
  rw [compIn, Finset.mem_filter]

lemma compIn_subset {v : V} : compIn G X v ⊆ X := Finset.filter_subset _ _

lemma self_mem_compIn {v : V} (hv : v ∈ X) : v ∈ compIn G X v :=
  mem_compIn.2 ⟨hv, ReachIn.refl⟩

lemma compIn_closed {v x y : V} (hx : x ∈ compIn G X v) (hy : y ∈ X) (hadj : G.Adj x y) :
    y ∈ compIn G X v := by
  rw [mem_compIn] at hx ⊢
  exact ⟨hy, hx.2.trans (ReachIn.single ⟨hx.1, hy, hadj⟩)⟩

/-- The complement of a component is closed. -/
lemma sdiff_compIn_closed {v : V} :
    ∀ u ∈ X \ compIn G X v, ∀ w ∈ X, G.Adj u w → w ∈ X \ compIn G X v := by
  intro u hu w hw hadj
  rw [Finset.mem_sdiff] at hu ⊢
  refine ⟨hw, fun hwc => hu.2 ?_⟩
  exact compIn_closed hwc hu.1 hadj.symm

/-! ### Transversals -/

/-- `Y` meets every connected component of `X` in exactly one vertex. -/
def IsTransversal (G : SimpleGraph V) (X Y : Finset V) : Prop :=
  Y ⊆ X ∧ (∀ x ∈ X, ∃ y ∈ Y, ReachIn G X x y) ∧
    (∀ y ∈ Y, ∀ y' ∈ Y, ReachIn G X y y' → y = y')

lemma IsTransversal.nonempty (hY : IsTransversal G X Y) (hX : X.Nonempty) : Y.Nonempty := by
  obtain ⟨x, hx⟩ := hX
  obtain ⟨y, hy, -⟩ := hY.2.1 x hx
  exact ⟨y, hy⟩

lemma IsTransversal.eq_empty (hY : IsTransversal G X Y) (hX : X = ∅) : Y = ∅ := by
  rw [Finset.eq_empty_iff_forall_notMem]
  intro y hy
  have := hY.1 hy
  rw [hX] at this
  exact absurd this (Finset.notMem_empty y)

/-- Peeling off the component of a transversal vertex leaves a transversal. -/
theorem transversal_erase {Y : Finset V} (hY : IsTransversal G X Y) {v : V} (hv : v ∈ Y) :
    IsTransversal G (X \ compIn G X v) (Y.erase v) := by
  obtain ⟨hYX, hex, huniq⟩ := hY
  have hvX : v ∈ X := hYX hv
  refine ⟨?_, ?_, ?_⟩
  · intro y hy
    rw [Finset.mem_erase] at hy
    rw [Finset.mem_sdiff]
    refine ⟨hYX hy.2, fun hyc => hy.1 ?_⟩
    exact (huniq y hy.2 v hv ((mem_compIn.1 hyc).2.symm)).symm ▸ rfl
  · intro x hx
    rw [Finset.mem_sdiff] at hx
    obtain ⟨y, hy, hreach⟩ := hex x hx.1
    have hyv : y ≠ v := by
      rintro rfl
      exact hx.2 (mem_compIn.2 ⟨hx.1, hreach.symm⟩)
    refine ⟨y, Finset.mem_erase.2 ⟨hyv, hy⟩, ?_⟩
    exact ReachIn.restrict sdiff_compIn_closed (Finset.mem_sdiff.2 hx) hreach
  · intro y hy y' hy' hreach
    rw [Finset.mem_erase] at hy hy'
    exact huniq y hy.2 y' hy'.2 (hreach.mono Finset.sdiff_subset)

/-- After removing the root, the neighbours of the root form a transversal of what is left of
its component. -/
theorem transversal_children (hac : G.IsAcyclic) {Y : Finset V} (hY : IsTransversal G X Y)
    {v : V} (hv : v ∈ Y) :
    IsTransversal G ((compIn G X v).erase v)
      (((compIn G X v).erase v).filter (fun c => G.Adj v c)) := by
  obtain ⟨hYX, hex, huniq⟩ := hY
  have hvX : v ∈ X := hYX hv
  set T := compIn G X v with hT
  set Ch := T.erase v with hCh
  have hChX : Ch ⊆ X := (Finset.erase_subset _ _).trans compIn_subset
  have hvCh : v ∉ Ch := fun h => (Finset.mem_erase.1 h).1 rfl
  have hclosed : ∀ u ∈ Ch, ∀ w ∈ X.erase v, G.Adj u w → w ∈ Ch := by
    intro u hu w hw hadj
    rw [Finset.mem_erase] at hw ⊢
    exact ⟨hw.1, compIn_closed (Finset.mem_of_mem_erase hu) hw.2 hadj⟩
  refine ⟨Finset.filter_subset _ _, ?_, ?_⟩
  · intro x hx
    have hxT : x ∈ T := Finset.mem_of_mem_erase hx
    have hxv : x ≠ v := (Finset.mem_erase.1 hx).1
    obtain ⟨c, hadj, hcX, hreach⟩ := exists_adj_reachIn_erase (mem_compIn.1 hxT).2 hxv
    have hcCh : c ∈ Ch := by
      rw [hCh, Finset.mem_erase]
      exact ⟨hadj.ne', compIn_closed (self_mem_compIn hvX) hcX hadj⟩
    refine ⟨c, Finset.mem_filter.2 ⟨hcCh, hadj⟩, ?_⟩
    exact (ReachIn.restrict hclosed hcCh hreach).symm
  · intro y hy y' hy' hreach
    rw [Finset.mem_filter] at hy hy'
    by_contra hne
    exact not_reachIn_of_adj hac hvCh hy.2 hy'.2 hy.1 hne hreach

/-- A set meeting every component contains a transversal. -/
theorem exists_transversal_subset {R : Finset V} (hRX : R ⊆ X)
    (hR : ∀ x ∈ X, ∃ y ∈ R, ReachIn G X x y) :
    ∃ Y ⊆ R, IsTransversal G X Y := by
  classical
  set idx : V → ℕ := fun v => ((Fintype.equivFin V) v : ℕ) with hidx
  have hinj : Function.Injective idx := fun a b hab => by
    have : (Fintype.equivFin V) a = (Fintype.equivFin V) b := Fin.ext hab
    exact (Fintype.equivFin V).injective this
  refine ⟨R.filter (fun y => ∀ z ∈ R, ReachIn G X y z → idx y ≤ idx z),
    Finset.filter_subset _ _, (Finset.filter_subset _ _).trans hRX, ?_, ?_⟩
  · intro x hx
    have hne : (R.filter (fun y => ReachIn G X x y)).Nonempty := by
      obtain ⟨y, hy, hreach⟩ := hR x hx
      exact ⟨y, Finset.mem_filter.2 ⟨hy, hreach⟩⟩
    obtain ⟨y₀, hy₀, hmin⟩ := Finset.exists_min_image _ idx hne
    rw [Finset.mem_filter] at hy₀
    refine ⟨y₀, Finset.mem_filter.2 ⟨hy₀.1, ?_⟩, hy₀.2⟩
    intro z hz hreach
    exact hmin z (Finset.mem_filter.2 ⟨hz, hy₀.2.trans hreach⟩)
  · intro y hy y' hy' hreach
    rw [Finset.mem_filter] at hy hy'
    have h1 := hy.2 y' hy'.1 hreach
    have h2 := hy'.2 y hy.1 hreach.symm
    exact hinj (le_antisymm h1 h2)

end B1Zero

end MatchingBag
