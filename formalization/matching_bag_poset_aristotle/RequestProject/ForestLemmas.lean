import Mathlib

/-!
# Forests: closed non-backtracking walks, and contraction of a matching

This file contains the graph-theoretic input for the tree-to-forest-poset bridge of
`matching_bag_poset_reduction_2026-08-20.md` (§1).  Everything here is about an arbitrary
acyclic simple graph (a forest); nothing refers to codes or posets.

Main results.

* `MatchingBag.unique_closer_neighbour`: in a forest, a vertex has at most one neighbour
  which is strictly closer to a fixed base point.
* `MatchingBag.no_periodic_two_block`: a forest carries no closed non-backtracking walk.
  The walk is presented in the *two-block* form which the contraction argument produces: a
  periodic pair of sequences `a k` ("the vertex at which block `k` is entered") and `b k`
  ("the vertex at which block `k` is left").
* `MatchingBag.bagGraph_isAcyclic`: contracting a set of pairwise disjoint edges of a forest
  yields a forest.  This is the statement that *contracting the matching edges of a tree
  gives a tree*, i.e. the port-labelled quotient tree of the note.
* `MatchingBag.no_closed_matching_chain`: the directed comparison relation
  `i → j  ⟺  i ≠ j ∧ T.Adj (R i) (L j)` obtained from the non-matching edges has no
  closed chain; this is the acyclicity used to get antisymmetry of the induced order.
-/

open SimpleGraph

namespace MatchingBag

variable {V : Type*}

/-! ### Arithmetic of periodic indexing -/

lemma mod_succ_eq {m : ℕ} (k : ℕ) : (k + 1) % m = (k % m + 1) % m := by
  simp [Nat.add_mod]

lemma mod_add_two_eq {m : ℕ} (k : ℕ) : (k + 2) % m = (k % m + 2) % m := by
  simp [Nat.add_mod]

/-! ### Distances in a forest -/

/-- In a forest, a vertex `a` has at most one neighbour at distance `dist u a - 1` from a base
point `u`: the neighbour on the unique path from `u` to `a`. -/
theorem unique_closer_neighbour {G : SimpleGraph V} (hG : G.IsAcyclic) {u a b c : V}
    (hab : G.Adj a b) (hac : G.Adj a c) (hub : G.Reachable u b) (huc : G.Reachable u c)
    (hb : G.dist u b + 1 = G.dist u a) (hc : G.dist u c + 1 = G.dist u a) : b = c := by
  obtain ⟨p, hp, hpl⟩ := hub.exists_path_of_dist
  obtain ⟨q, hq, hql⟩ := huc.exists_path_of_dist
  have hP : (p.concat hab.symm).IsPath := by
    refine Walk.isPath_of_length_eq_dist _ ?_
    rw [Walk.length_concat, hpl, hb]
  have hQ : (q.concat hac.symm).IsPath := by
    refine Walk.isPath_of_length_eq_dist _ ?_
    rw [Walk.length_concat, hql, hc]
  have hEq : p.concat hab.symm = q.concat hac.symm :=
    congrArg Subtype.val (hG.path_unique ⟨_, hP⟩ ⟨_, hQ⟩)
  calc b = (p.concat hab.symm).penultimate := (Walk.penultimate_concat p hab.symm).symm
    _ = (q.concat hac.symm).penultimate := by rw [hEq]
    _ = c := Walk.penultimate_concat q hac.symm

/-! ### No closed non-backtracking walk in a forest -/

/-- A forest carries no closed non-backtracking walk.

The walk is given in *two-block* form: the `k`-th block of the walk consists of the vertex
`a k` followed (if different) by the vertex `b k`; consecutive blocks are joined by an edge
`b k — a (k+1)`.  Both sequences are periodic with period `n > 0`, so the walk is closed.
The last two hypotheses express that the walk never backtracks. -/
theorem no_periodic_two_block {G : SimpleGraph V} (hG : G.IsAcyclic) {n : ℕ} (hn : 0 < n)
    (a b : ℕ → V) (hpa : ∀ k, a (k + n) = a k) (hpb : ∀ k, b (k + n) = b k)
    (hin : ∀ k, a k = b k ∨ G.Adj (a k) (b k))
    (hstep : ∀ k, G.Adj (b k) (a (k + 1)))
    (hnb_stay : ∀ k, a (k + 1) = b (k + 1) → b k ≠ a (k + 2))
    (hnb_move : ∀ k, a (k + 1) ≠ b (k + 1) → b k ≠ b (k + 1) ∧ a (k + 1) ≠ a (k + 2)) :
    False := by
  classical
  set u := a 1 with hu
  -- every vertex of the walk is reachable from the base point
  have hra1 : ∀ j, G.Reachable u (a (j + 1)) ∧ G.Reachable u (b (j + 1)) := by
    intro j
    induction j with
    | zero =>
        refine ⟨Reachable.refl _, ?_⟩
        rcases hin 1 with h | h
        · exact h ▸ Reachable.refl _
        · exact h.reachable
    | succ j ih =>
        have h1 : G.Reachable u (a (j + 1 + 1)) := ih.2.trans (hstep (j + 1)).reachable
        refine ⟨h1, ?_⟩
        rcases hin (j + 1 + 1) with h | h
        · exact h ▸ h1
        · exact h1.trans h.reachable
  have hra : ∀ j, G.Reachable u (a j) := by
    intro j
    match j with
    | 0 =>
        have h0 : a 0 = a (n - 1 + 1) := by
          have hp := hpa 0
          have hn1 : n - 1 + 1 = 0 + n := by omega
          rw [hn1, hp]
        rw [h0]
        exact (hra1 (n - 1)).1
    | (j + 1) => exact (hra1 j).1
  have hrb : ∀ j, G.Reachable u (b j) := by
    intro j
    match j with
    | 0 =>
        have h0 : b 0 = b (n - 1 + 1) := by
          have hp := hpb 0
          have hn1 : n - 1 + 1 = 0 + n := by omega
          rw [hn1, hp]
        rw [h0]
        exact (hra1 (n - 1)).2
    | (j + 1) => exact (hra1 j).2
  -- adjacent vertices at no greater distance are strictly closer
  have key : ∀ w w' : V, G.Adj w w' → G.Reachable u w → G.dist u w' ≤ G.dist u w →
      G.dist u w' + 1 = G.dist u w := by
    intro w w' hadj hre hle
    rcases hG.dist_eq_dist_add_one_of_adj_of_reachable u hadj hre with h | h
    · omega
    · omega
  -- a vertex of the walk at maximal distance from the base point
  set F : ℕ → ℕ := fun j => max (G.dist u (a j)) (G.dist u (b j)) with hF
  have hFper : ∀ j, F (j + n) = F j := by
    intro j
    simp only [hF, hpa, hpb]
  have hFmod : ∀ j, F j = F (j % n) := by
    intro j
    induction j using Nat.strong_induction_on with
    | _ j ih =>
      by_cases hj : j < n
      · rw [Nat.mod_eq_of_lt hj]
      · have h1 : j - n + n = j := by omega
        have h2 : F j = F (j - n) := by rw [← hFper (j - n), h1]
        rw [h2, ih (j - n) (by omega)]
        conv_rhs => rw [Nat.mod_eq_sub_mod (show n ≤ j by omega)]
  obtain ⟨k, hk_mem, hk_max⟩ :=
    Finset.exists_max_image (Finset.Icc 1 n) F ⟨n, by simp [Finset.mem_Icc]; omega⟩
  have hall : ∀ j, F j ≤ F k := by
    intro j
    by_cases hj0 : j % n = 0
    · have hFn : F j = F n := by
        rw [hFmod j, hj0]
        have h := hFper 0
        rw [Nat.zero_add] at h
        exact h.symm
      rw [hFn]
      exact hk_max n (by simp [Finset.mem_Icc]; omega)
    · have hlt : j % n < n := Nat.mod_lt _ hn
      rw [hFmod j]
      refine hk_max (j % n) ?_
      simp only [Finset.mem_Icc]
      omega
  obtain ⟨k', rfl⟩ : ∃ k', k = k' + 1 := ⟨k - 1, by have := (Finset.mem_Icc.1 hk_mem).1; omega⟩
  have hda : ∀ j, G.dist u (a j) ≤ F (k' + 1) := by
    intro j
    refine le_trans ?_ (hall j)
    simp only [hF]
    exact le_max_left _ _
  have hdb : ∀ j, G.dist u (b j) ≤ F (k' + 1) := by
    intro j
    refine le_trans ?_ (hall j)
    simp only [hF]
    exact le_max_right _ _
  by_cases hcase : G.dist u (b (k' + 1)) ≤ G.dist u (a (k' + 1))
  · -- the maximum is attained at `a (k'+1)`
    have hMa : F (k' + 1) = G.dist u (a (k' + 1)) := by
      simp only [hF]
      exact max_eq_left hcase
    have h1 : G.dist u (b k') + 1 = G.dist u (a (k' + 1)) := by
      refine key _ _ (hstep k').symm (hra _) ?_
      rw [← hMa]
      exact hdb k'
    by_cases heq : a (k' + 1) = b (k' + 1)
    · have hadj2 : G.Adj (a (k' + 1)) (a (k' + 1 + 1)) := by
        rw [heq]
        exact hstep (k' + 1)
      have h2 : G.dist u (a (k' + 1 + 1)) + 1 = G.dist u (a (k' + 1)) := by
        refine key _ _ hadj2 (hra _) ?_
        rw [← hMa]
        exact hda _
      exact hnb_stay k' heq
        (unique_closer_neighbour hG (hstep k').symm hadj2 (hrb k') (hra _) h1 h2)
    · have hadj2 : G.Adj (a (k' + 1)) (b (k' + 1)) := (hin (k' + 1)).resolve_left heq
      have h2 : G.dist u (b (k' + 1)) + 1 = G.dist u (a (k' + 1)) := by
        refine key _ _ hadj2 (hra _) ?_
        rw [← hMa]
        exact hdb _
      exact (hnb_move k' heq).1
        (unique_closer_neighbour hG (hstep k').symm hadj2 (hrb k') (hrb _) h1 h2)
  · -- the maximum is attained at `b (k'+1)`
    push_neg at hcase
    have hMb : F (k' + 1) = G.dist u (b (k' + 1)) := by
      simp only [hF]
      exact max_eq_right (le_of_lt hcase)
    have hne : a (k' + 1) ≠ b (k' + 1) := by
      intro h
      rw [h] at hcase
      omega
    have hadj2 : G.Adj (a (k' + 1)) (b (k' + 1)) := (hin (k' + 1)).resolve_left hne
    have h1 : G.dist u (a (k' + 1)) + 1 = G.dist u (b (k' + 1)) := by
      refine key _ _ hadj2.symm (hrb _) ?_
      rw [← hMb]
      exact hda _
    have h2 : G.dist u (a (k' + 1 + 1)) + 1 = G.dist u (b (k' + 1)) := by
      refine key _ _ (hstep (k' + 1)) (hrb _) ?_
      rw [← hMb]
      exact hda _
    exact (hnb_move k' hne).2
      (unique_closer_neighbour hG hadj2.symm (hstep (k' + 1)) (hra _) (hra _) h1 h2)

/-! ### Closed chains of a relation -/

/-- A closed chain of a relation can be presented as a periodic sequence. -/
theorem exists_periodic_of_transGen {ι : Type*} {r : ι → ι → Prop} {i : ι}
    (h : Relation.TransGen r i i) :
    ∃ n : ℕ, 0 < n ∧ ∃ f : ℕ → ι, (∀ k, f (k + n) = f k) ∧ ∀ k, r (f k) (f (k + 1)) := by
  classical
  have key : ∀ {x y : ι}, Relation.TransGen r x y →
      ∃ m : ℕ, 0 < m ∧ ∃ g : ℕ → ι, g 0 = x ∧ g m = y ∧ ∀ j < m, r (g j) (g (j + 1)) := by
    intro x y hxy
    induction hxy with
    | @single y hs =>
        refine ⟨1, one_pos, fun j => if j = 0 then x else y, by simp, by simp, ?_⟩
        intro j hj
        have hj0 : j = 0 := by omega
        subst hj0
        simpa using hs
    | @tail y z hxy hyz ih =>
        obtain ⟨m, hm, g, h0, hmy, hstep⟩ := ih
        refine ⟨m + 1, by omega, fun j => if j ≤ m then g j else z, ?_, by simp, ?_⟩
        · simp only [if_pos (Nat.zero_le m)]
          exact h0
        · intro j hj
          by_cases hjm : j < m
          · simp only [if_pos (show j ≤ m by omega), if_pos (show j + 1 ≤ m by omega)]
            exact hstep j hjm
          · have hjm' : j = m := by omega
            subst hjm'
            simp only [if_pos (le_refl j), if_neg (show ¬ (j + 1 ≤ j) by omega), hmy]
            exact hyz
  obtain ⟨m, hm, g, h0, hmi, hstep⟩ := key h
  refine ⟨m, hm, fun k => g (k % m), ?_, ?_⟩
  · intro k
    simp [Nat.add_mod_right]
  · intro k
    have hs : k % m < m := Nat.mod_lt _ hm
    by_cases hcase : k % m + 1 = m
    · have hz : (k + 1) % m = 0 := by
        rw [mod_succ_eq, hcase, Nat.mod_self]
      simp only [hz, h0, ← hmi]
      have := hstep (k % m) hs
      rw [hcase] at this
      simpa [h0, hmi] using this
    · have hlt : k % m + 1 < m := by omega
      have hz : (k + 1) % m = k % m + 1 := by
        rw [mod_succ_eq, Nat.mod_eq_of_lt hlt]
      simp only [hz]
      exact hstep (k % m) hs

/-- A cycle of a graph gives a periodic sequence of vertices with distinct consecutive and
next-to-consecutive entries. -/
theorem exists_periodic_of_isCycle {G : SimpleGraph V} {u : V} {c : G.Walk u u}
    (hc : c.IsCycle) :
    ∃ n : ℕ, 3 ≤ n ∧ ∃ B : ℕ → V, (∀ k, B (k + n) = B k) ∧ (∀ k, G.Adj (B k) (B (k + 1)))
      ∧ (∀ k, B k ≠ B (k + 2)) := by
  classical
  set n := c.length with hn
  have hn3 : 3 ≤ n := hc.three_le_length
  refine ⟨n, hn3, fun k => c.getVert (k % n), ?_, ?_, ?_⟩
  · intro k
    simp [Nat.add_mod_right]
  · intro k
    show G.Adj (c.getVert (k % n)) (c.getVert ((k + 1) % n))
    have hs : k % n < n := Nat.mod_lt _ (by omega)
    by_cases hcase : k % n + 1 = n
    · have hz : (k + 1) % n = 0 := by rw [mod_succ_eq, hcase, Nat.mod_self]
      rw [hz, c.getVert_zero]
      have h := c.adj_getVert_succ (i := k % n) (by omega)
      rw [hcase, hn, c.getVert_length] at h
      exact h
    · have hlt : k % n + 1 < n := by omega
      have hz : (k + 1) % n = k % n + 1 := by rw [mod_succ_eq, Nat.mod_eq_of_lt hlt]
      rw [hz]
      exact c.adj_getVert_succ (by omega)
  · intro k hcon
    have hs : k % n < n := Nat.mod_lt _ (by omega)
    have h2 : (k + 2) % n < n := Nat.mod_lt _ (by omega)
    have hne : k % n ≠ (k + 2) % n := by
      rcases lt_trichotomy (k % n + 2) n with hlt | heq | hgt
      · rw [mod_add_two_eq, Nat.mod_eq_of_lt hlt]
        omega
      · rw [mod_add_two_eq, heq, Nat.mod_self]
        omega
      · have hk : k % n + 2 = n + 1 := by omega
        have : (k % n + 2) % n = 1 := by
          rw [hk, Nat.add_mod, Nat.mod_self, Nat.zero_add,
            Nat.mod_eq_of_lt (show (1 : ℕ) < n by omega),
            Nat.mod_eq_of_lt (show (1 : ℕ) < n by omega)]
        rw [mod_add_two_eq, this]
        omega
    exact hne (hc.getVert_injOn' (by simp only [Set.mem_setOf_eq]; omega)
      (by simp only [Set.mem_setOf_eq]; omega) hcon)

/-! ### Contracting a matching in a forest -/

section Contraction

variable {ι : Type*} {T : SimpleGraph V} {Lv Rv : ι → V}

/-- The hypotheses of the contraction lemmas: `Lv i — Rv i` are pairwise disjoint edges of the
forest `T`. -/
structure DisjointEdges (T : SimpleGraph V) (Lv Rv : ι → V) : Prop where
  adj : ∀ i, T.Adj (Lv i) (Rv i)
  disj : ∀ i j : ι, i ≠ j → Lv i ≠ Lv j ∧ Lv i ≠ Rv j ∧ Rv i ≠ Lv j ∧ Rv i ≠ Rv j

/-- **Contracting a matching of a forest gives a forest.**  The vertices of the contracted
graph are the matched edges `i`, and two of them are adjacent when the forest has an edge
between them. -/
theorem bagGraph_isAcyclic (hT : T.IsAcyclic) (hD : DisjointEdges T Lv Rv) :
    (SimpleGraph.fromRel (fun i j : ι => T.Adj (Rv i) (Lv j))).IsAcyclic := by
  classical
  intro v c hc
  obtain ⟨n, hn3, B, hBper, hBadj, hBne2⟩ := exists_periodic_of_isCycle hc
  have hBne : ∀ k, B k ≠ B (k + 1) := fun k => (hBadj k).ne
  have hBor : ∀ k, T.Adj (Rv (B k)) (Lv (B (k + 1))) ∨ T.Adj (Rv (B (k + 1))) (Lv (B k)) :=
    fun k => (SimpleGraph.fromRel_adj _ _ _ |>.1 (hBadj k)).2
  set dir : ℕ → Bool := fun k => decide (T.Adj (Rv (B k)) (Lv (B (k + 1)))) with hdir
  have hdirT : ∀ k, dir k = true → T.Adj (Rv (B k)) (Lv (B (k + 1))) := by
    intro k h
    simpa [hdir] using h
  have hdirF : ∀ k, dir k = false → T.Adj (Rv (B (k + 1))) (Lv (B k)) := by
    intro k h
    have hnot : ¬ T.Adj (Rv (B k)) (Lv (B (k + 1))) := by simpa [hdir] using h
    exact (hBor k).resolve_left hnot
  have hdirper : ∀ k, dir (k + n) = dir k := by
    intro k
    have h1 : B (k + n) = B k := hBper k
    have h2 : B (k + n + 1) = B (k + 1) := by
      have := hBper (k + 1)
      rw [show k + 1 + n = k + n + 1 by omega] at this
      exact this
    simp only [hdir, h1, h2]
  set A : ℕ → V := fun k => if dir k = true then Lv (B (k + 1)) else Rv (B (k + 1)) with hA
  set Bl : ℕ → V := fun k => if dir (k + 1) = true then Rv (B (k + 1)) else Lv (B (k + 1))
    with hBl
  refine no_periodic_two_block hT (n := n) (by omega) A Bl ?_ ?_ ?_ ?_ ?_ ?_
  · intro k
    have h1 : B (k + n + 1) = B (k + 1) := by
      have := hBper (k + 1)
      rw [show k + 1 + n = k + n + 1 by omega] at this
      exact this
    simp only [hA, hdirper, h1]
  · intro k
    have h1 : B (k + n + 1) = B (k + 1) := by
      have := hBper (k + 1)
      rw [show k + 1 + n = k + n + 1 by omega] at this
      exact this
    have h2 : dir (k + n + 1) = dir (k + 1) := by
      have := hdirper (k + 1)
      rwa [show k + 1 + n = k + n + 1 by omega] at this
    simp only [hBl, h1, h2]
  · intro k
    cases hd1 : dir k <;> cases hd2 : dir (k + 1) <;>
      simp only [hA, hBl, hd1, hd2, if_true, if_false, Bool.false_eq_true]
    · exact Or.inr ((hD.adj (B (k + 1))).symm)
    · exact Or.inl trivial
    · exact Or.inl trivial
    · exact Or.inr (hD.adj (B (k + 1)))
  · intro k
    cases hd : dir (k + 1) <;> simp only [hA, hBl, hd, if_true, if_false, Bool.false_eq_true]
    · exact (hdirF (k + 1) hd).symm
    · exact hdirT (k + 1) hd
  · intro k hstay
    cases hd1 : dir (k + 1) <;> cases hd2 : dir (k + 2) <;>
      simp only [hA, hBl, hd1, hd2, if_true, if_false, Bool.false_eq_true,
        show k + 1 + 1 = k + 2 from rfl, show k + 2 + 1 = k + 3 from rfl] at hstay ⊢
    · exact absurd hstay (hD.adj (B (k + 2))).ne'
    · intro hcon
      exact (hD.disj (B (k + 1)) (B (k + 3)) (by
        intro hh
        exact hBne2 (k + 1) (by rw [hh])) ).1 hcon
    · intro hcon
      exact (hD.disj (B (k + 1)) (B (k + 3)) (by
        intro hh
        exact hBne2 (k + 1) (by rw [hh])) ).2.2.2 hcon
    · exact absurd hstay (hD.adj (B (k + 2))).ne
  · intro k hmove
    cases hd1 : dir (k + 1) <;> cases hd2 : dir (k + 2) <;>
      simp only [hA, hBl, hd1, hd2, if_true, if_false, Bool.false_eq_true,
        show k + 1 + 1 = k + 2 from rfl, show k + 2 + 1 = k + 3 from rfl] at hmove ⊢
    · refine ⟨?_, ?_⟩
      · exact (hD.disj (B (k + 1)) (B (k + 2)) (hBne (k + 1))).1
      · exact (hD.disj (B (k + 2)) (B (k + 3)) (hBne (k + 2))).2.2.2
    · exact absurd rfl hmove
    · exact absurd rfl hmove
    · refine ⟨?_, ?_⟩
      · exact (hD.disj (B (k + 1)) (B (k + 2)) (hBne (k + 1))).2.2.2
      · exact (hD.disj (B (k + 2)) (B (k + 3)) (hBne (k + 2))).1

/-- **The comparison relation has no closed chain.**  Here `i → j` means that the forest has a
non-matching edge from the right endpoint of `i` to the left endpoint of `j`. -/
theorem no_closed_matching_chain (hT : T.IsAcyclic) (hD : DisjointEdges T Lv Rv) (i : ι) :
    ¬ Relation.TransGen (fun a b : ι => a ≠ b ∧ T.Adj (Rv a) (Lv b)) i i := by
  intro h
  obtain ⟨n, hn, f, hper, hstep⟩ := exists_periodic_of_transGen h
  refine no_periodic_two_block hT hn (fun k => Lv (f k)) (fun k => Rv (f k)) ?_ ?_ ?_ ?_ ?_ ?_
  · intro k; simp only [hper]
  · intro k; simp only [hper]
  · intro k; exact Or.inr (hD.adj (f k))
  · intro k; exact (hstep k).2
  · intro k hk
    exact absurd hk (hD.adj (f (k + 1))).ne
  · intro k _
    exact ⟨(hD.disj (f k) (f (k + 1)) (hstep k).1).2.2.2,
      (hD.disj (f (k + 1)) (f (k + 2)) (hstep (k + 1)).1).1⟩

end Contraction

end MatchingBag
