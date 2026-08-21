import RequestProject.ClanAudit.DoubleSpider

/-!
# The arbitrary-connector two-hub tree

`connGraph t m a n b` is the tree consisting of

* a *left hub* carrying `m` pendant paths of lengths `a i`;
* a *right hub* carrying `n` pendant paths of lengths `b j`;
* a connector path of `t` interior vertices joining the two hubs, i.e. `ell = t + 1`
  connector edges.

`connectorGraph ell m a n b` is the same graph written in terms of the number `ell ≥ 1`
of connector edges.  For `ell = 1` (`t = 0`) it is the adjacent two-hub tree.

Proved here:

* `connGraph_adj_*` — the graph has exactly the advertised edges;
* `connGraph_degree_le_two` — every vertex other than the two hubs has degree at most
  two, so the only possible branch vertices are the two hubs;
* `connGraph_connected` — the graph is connected.
-/

namespace ClanAudit

open Finset SimpleGraph

open scoped Classical

/-- Vertices of the connector tree: the left spider, the `t` interior connector vertices,
and the right spider. -/
abbrev ConnV (t m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) :=
  SpiderV m a ⊕ Fin t ⊕ SpiderV n b

/-- The connector tree with `t` interior connector vertices (that is, with `t + 1`
connector edges). -/
def connGraph (t m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) :
    SimpleGraph (ConnV t m a n b) where
  Adj x y :=
    match x, y with
    | Sum.inl p, Sum.inl q => (spider m a).Adj p q
    | Sum.inr (Sum.inr p), Sum.inr (Sum.inr q) => (spider n b).Adj p q
    | Sum.inl p, Sum.inr (Sum.inl j) => p = none ∧ j.val = 0
    | Sum.inr (Sum.inl j), Sum.inl p => p = none ∧ j.val = 0
    | Sum.inl p, Sum.inr (Sum.inr q) => p = none ∧ q = none ∧ t = 0
    | Sum.inr (Sum.inr q), Sum.inl p => p = none ∧ q = none ∧ t = 0
    | Sum.inr (Sum.inl i), Sum.inr (Sum.inl j) => i.val + 1 = j.val ∨ j.val + 1 = i.val
    | Sum.inr (Sum.inl i), Sum.inr (Sum.inr q) => q = none ∧ i.val + 1 = t
    | Sum.inr (Sum.inr q), Sum.inr (Sum.inl i) => q = none ∧ i.val + 1 = t
  symm := by
    rintro (p | i | p) (q | j | q) h <;> first | exact h | exact h.symm
  loopless := ⟨by
    rintro (p | i | p) h
    · exact (spider m a).irrefl h
    · rcases h with h | h <;> omega
    · exact (spider n b).irrefl h⟩

/-- The connector tree written in terms of the number `ell ≥ 1` of connector edges. -/
abbrev connectorGraph (ell m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) :
    SimpleGraph (ConnV (ell - 1) m a n b) := connGraph (ell - 1) m a n b

variable {t m n : ℕ} {a : Fin m → ℕ} {b : Fin n → ℕ}

/-! ### The advertised edges -/

theorem connGraph_adj_ll {p q : SpiderV m a} :
    (connGraph t m a n b).Adj (Sum.inl p) (Sum.inl q) ↔ (spider m a).Adj p q := Iff.rfl

theorem connGraph_adj_rr {p q : SpiderV n b} :
    (connGraph t m a n b).Adj (Sum.inr (Sum.inr p)) (Sum.inr (Sum.inr q))
      ↔ (spider n b).Adj p q := Iff.rfl

theorem connGraph_adj_lc {p : SpiderV m a} {j : Fin t} :
    (connGraph t m a n b).Adj (Sum.inl p) (Sum.inr (Sum.inl j))
      ↔ (p = none ∧ j.val = 0) := Iff.rfl

theorem connGraph_adj_lr {p : SpiderV m a} {q : SpiderV n b} :
    (connGraph t m a n b).Adj (Sum.inl p) (Sum.inr (Sum.inr q))
      ↔ (p = none ∧ q = none ∧ t = 0) := Iff.rfl

theorem connGraph_adj_cc {i j : Fin t} :
    (connGraph t m a n b).Adj (Sum.inr (Sum.inl i)) (Sum.inr (Sum.inl j))
      ↔ (i.val + 1 = j.val ∨ j.val + 1 = i.val) := Iff.rfl

theorem connGraph_adj_cr {i : Fin t} {q : SpiderV n b} :
    (connGraph t m a n b).Adj (Sum.inr (Sum.inl i)) (Sum.inr (Sum.inr q))
      ↔ (q = none ∧ i.val + 1 = t) := Iff.rfl

/-- For `t = 0` the connector tree is the adjacent two-hub tree. -/
theorem connGraph_zero_adj_lr {p : SpiderV m a} {q : SpiderV n b} :
    (connGraph 0 m a n b).Adj (Sum.inl p) (Sum.inr (Sum.inr q))
      ↔ (p = none ∧ q = none) := by
  rw [connGraph_adj_lr]
  simp

/-! ### The level function -/

/-- The distance from the left hub. -/
def clevel : ConnV t m a n b → ℕ
  | Sum.inl p => spiderLevel p
  | Sum.inr (Sum.inl i) => i.val + 1
  | Sum.inr (Sum.inr q) => spiderLevel q + t + 1

theorem clevel_adj {x y : ConnV t m a n b} (h : (connGraph t m a n b).Adj x y) :
    clevel x + 1 = clevel y ∨ clevel y + 1 = clevel x := by
  match x, y with
  | Sum.inl p, Sum.inl q => exact spiderLevel_adj (connGraph_adj_ll.mp h)
  | Sum.inr (Sum.inr p), Sum.inr (Sum.inr q) =>
      have := spiderLevel_adj (connGraph_adj_rr.mp h)
      simp only [clevel]
      omega
  | Sum.inl p, Sum.inr (Sum.inl j) =>
      obtain ⟨hp, hj⟩ := connGraph_adj_lc.mp h
      subst hp
      simp only [clevel, spiderLevel]
      omega
  | Sum.inr (Sum.inl j), Sum.inl p =>
      obtain ⟨hp, hj⟩ := h
      subst hp
      simp only [clevel, spiderLevel]
      omega
  | Sum.inl p, Sum.inr (Sum.inr q) =>
      obtain ⟨hp, hq, ht⟩ := connGraph_adj_lr.mp h
      subst hp; subst hq; subst ht
      simp [clevel, spiderLevel]
  | Sum.inr (Sum.inr q), Sum.inl p =>
      obtain ⟨hp, hq, ht⟩ := h
      subst hp; subst hq; subst ht
      simp [clevel, spiderLevel]
  | Sum.inr (Sum.inl i), Sum.inr (Sum.inl j) =>
      have := connGraph_adj_cc.mp h
      simp only [clevel]
      omega
  | Sum.inr (Sum.inl i), Sum.inr (Sum.inr q) =>
      obtain ⟨hq, hi⟩ := connGraph_adj_cr.mp h
      subst hq
      simp only [clevel, spiderLevel]
      omega
  | Sum.inr (Sum.inr q), Sum.inr (Sum.inl i) =>
      obtain ⟨hq, hi⟩ := h
      subst hq
      simp only [clevel, spiderLevel]
      omega

/-! ### Only the two hubs can be branch vertices -/

/-- The left hub. -/
abbrev chubL (t m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) : ConnV t m a n b :=
  Sum.inl none

/-- The right hub. -/
abbrev chubR (t m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) : ConnV t m a n b :=
  Sum.inr (Sum.inr none)

/-- Two neighbours of a non-hub vertex at the same distance from the left hub coincide. -/
theorem eq_of_adj_of_clevel_eq {x y z : ConnV t m a n b}
    (hxL : x ≠ chubL t m a n b) (hxR : x ≠ chubR t m a n b)
    (hy : (connGraph t m a n b).Adj x y) (hz : (connGraph t m a n b).Adj x z)
    (hlev : clevel y = clevel z) : y = z := by
  match x with
  | Sum.inl none => exact absurd rfl hxL
  | Sum.inr (Sum.inr none) => exact absurd rfl hxR
  | Sum.inl (some ⟨i, j⟩) =>
      -- neighbours lie in the left spider
      have hcase : ∀ w : ConnV t m a n b,
          (connGraph t m a n b).Adj (Sum.inl (some ⟨i, j⟩)) w →
            ∃ p : SpiderV m a, w = Sum.inl p ∧ (spider m a).Adj (some ⟨i, j⟩) p := by
        rintro (p | k | p) hw
        · exact ⟨p, rfl, hw⟩
        · exact absurd hw.1 (by simp)
        · exact absurd hw.1 (by simp)
      obtain ⟨p, rfl, hp⟩ := hcase y hy
      obtain ⟨q, rfl, hq⟩ := hcase z hz
      have hlev' : spiderLevel p = spiderLevel q := hlev
      match p, q with
      | none, none => rfl
      | none, some ⟨i', j'⟩ => exact absurd hlev' (by simp [spiderLevel])
      | some ⟨i', j'⟩, none => exact absurd hlev' (by simp [spiderLevel])
      | some ⟨i', j'⟩, some ⟨i'', j''⟩ =>
          obtain ⟨h1, -⟩ := hp
          obtain ⟨h2, -⟩ := hq
          simp only at h1 h2
          subst h1; subst h2
          have : j'.val = j''.val := by simpa [spiderLevel] using hlev'
          simp only [Sum.inl.injEq, Option.some.injEq, Sigma.mk.injEq, heq_eq_eq, true_and]
          exact Fin.ext this
  | Sum.inr (Sum.inr (some ⟨i, j⟩)) =>
      have hcase : ∀ w : ConnV t m a n b,
          (connGraph t m a n b).Adj (Sum.inr (Sum.inr (some ⟨i, j⟩))) w →
            ∃ p : SpiderV n b, w = Sum.inr (Sum.inr p) ∧ (spider n b).Adj (some ⟨i, j⟩) p := by
        rintro (p | k | p) hw
        · exact absurd hw.2.1 (by simp)
        · exact absurd hw.1 (by simp)
        · exact ⟨p, rfl, hw⟩
      obtain ⟨p, rfl, hp⟩ := hcase y hy
      obtain ⟨q, rfl, hq⟩ := hcase z hz
      have hlev' : spiderLevel p = spiderLevel q := by
        simpa [clevel] using hlev
      match p, q with
      | none, none => rfl
      | none, some ⟨i', j'⟩ => exact absurd hlev' (by simp [spiderLevel])
      | some ⟨i', j'⟩, none => exact absurd hlev' (by simp [spiderLevel])
      | some ⟨i', j'⟩, some ⟨i'', j''⟩ =>
          obtain ⟨h1, -⟩ := hp
          obtain ⟨h2, -⟩ := hq
          simp only at h1 h2
          subst h1; subst h2
          have : j'.val = j''.val := by simpa [spiderLevel] using hlev'
          simp only [Sum.inr.injEq, Option.some.injEq, Sigma.mk.injEq, heq_eq_eq, true_and]
          exact Fin.ext this
  | Sum.inr (Sum.inl i) =>
      -- an interior connector vertex: its neighbours are determined by their level
      have hlevy := clevel_adj hy
      have hlevz := clevel_adj hz
      match y, z with
      | Sum.inl p, Sum.inl q =>
          obtain ⟨hp, -⟩ := hy
          obtain ⟨hq, -⟩ := hz
          rw [hp, hq]
      | Sum.inl p, Sum.inr (Sum.inl k) =>
          obtain ⟨hp, hi⟩ := hy
          subst hp
          exact absurd hlev (by simp [clevel, spiderLevel])
      | Sum.inl p, Sum.inr (Sum.inr q) =>
          obtain ⟨hp, hi⟩ := hy
          obtain ⟨hq, ht⟩ := hz
          subst hp; subst hq
          simp only [clevel, spiderLevel] at hlev
          omega
      | Sum.inr (Sum.inl k), Sum.inl q =>
          obtain ⟨hq, -⟩ := hz
          subst hq
          exact absurd hlev (by simp [clevel, spiderLevel])
      | Sum.inr (Sum.inl k), Sum.inr (Sum.inl k') =>
          have : k.val = k'.val := by simpa [clevel] using hlev
          simp only [Sum.inr.injEq, Sum.inl.injEq]
          exact Fin.ext this
      | Sum.inr (Sum.inl k), Sum.inr (Sum.inr q) =>
          obtain ⟨hq, ht⟩ := hz
          subst hq
          simp only [clevel, spiderLevel] at hlev
          have hk : k.val < t := k.isLt
          omega
      | Sum.inr (Sum.inr p), Sum.inl q =>
          obtain ⟨hp, ht⟩ := hy
          obtain ⟨hq, -⟩ := hz
          subst hp; subst hq
          simp only [clevel, spiderLevel] at hlev
          omega
      | Sum.inr (Sum.inr p), Sum.inr (Sum.inl k) =>
          obtain ⟨hp, ht⟩ := hy
          subst hp
          simp only [clevel, spiderLevel] at hlev
          have hk : k.val < t := k.isLt
          omega
      | Sum.inr (Sum.inr p), Sum.inr (Sum.inr q) =>
          obtain ⟨hp, ht⟩ := hy
          obtain ⟨hq, -⟩ := hz
          rw [hp, hq]

/-- **Only the two hubs can be branch vertices.**  Every other vertex of the connector
tree has degree at most two. -/
theorem connGraph_degree_le_two (x : ConnV t m a n b) (hxL : x ≠ chubL t m a n b)
    (hxR : x ≠ chubR t m a n b) : (connGraph t m a n b).degree x ≤ 2 := by
  classical
  have hmap : Set.MapsTo (fun y => decide (clevel y = clevel x + 1))
      (((connGraph t m a n b).neighborFinset x : Finset (ConnV t m a n b)) :
        Set (ConnV t m a n b)) ((Finset.univ : Finset Bool) : Set Bool) := by
    intro y _
    exact Finset.mem_coe.mpr (Finset.mem_univ _)
  have hinj : Set.InjOn (fun y => decide (clevel y = clevel x + 1))
      (((connGraph t m a n b).neighborFinset x : Finset (ConnV t m a n b)) :
        Set (ConnV t m a n b)) := by
    intro y hy z hz hyz
    rw [Finset.mem_coe, SimpleGraph.mem_neighborFinset] at hy hz
    have hly := clevel_adj hy
    have hlz := clevel_adj hz
    refine eq_of_adj_of_clevel_eq hxL hxR hy hz ?_
    by_cases hy1 : clevel y = clevel x + 1
    · have hz1 : clevel z = clevel x + 1 := by
        by_contra hc
        simp [hy1, hc] at hyz
      omega
    · have hz1 : ¬ clevel z = clevel x + 1 := by
        by_contra hc
        simp [hy1, hc] at hyz
      omega
  have hcard := Finset.card_le_card_of_injOn _ hmap hinj
  simpa using hcard

/-- The only possible vertices of degree at least three are the two hubs. -/
theorem connGraph_branch_vertex (x : ConnV t m a n b)
    (hx : 3 ≤ (connGraph t m a n b).degree x) :
    x = chubL t m a n b ∨ x = chubR t m a n b := by
  by_contra hc
  push_neg at hc
  have := connGraph_degree_le_two x hc.1 hc.2
  omega

/-! ### Connectivity -/

/-- The inclusion of the left spider. -/
def cinl (t m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) :
    spider m a →g connGraph t m a n b where
  toFun := Sum.inl
  map_rel' := id

/-- The inclusion of the right spider. -/
def cinr (t m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) :
    spider n b →g connGraph t m a n b where
  toFun := fun q => Sum.inr (Sum.inr q)
  map_rel' := id

theorem creachable_interior (t m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) :
    ∀ (k : ℕ) (i : Fin t), i.val = k →
      (connGraph t m a n b).Reachable (chubL t m a n b) (Sum.inr (Sum.inl i)) := by
  intro k
  induction k with
  | zero =>
      intro i hi
      exact SimpleGraph.Adj.reachable (connGraph_adj_lc.mpr ⟨rfl, hi⟩)
  | succ k ih =>
      intro i hi
      have hk : k < t := by omega
      have h1 := ih ⟨k, hk⟩ rfl
      exact h1.trans (SimpleGraph.Adj.reachable (connGraph_adj_cc.mpr (Or.inl (by simpa using hi.symm))))

theorem creachable_hubR (t m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) :
    (connGraph t m a n b).Reachable (chubL t m a n b) (chubR t m a n b) := by
  match t with
  | 0 => exact SimpleGraph.Adj.reachable (connGraph_adj_lr.mpr ⟨rfl, rfl, rfl⟩)
  | (s + 1) =>
      have h1 := creachable_interior (s + 1) m a n b s ⟨s, by omega⟩ rfl
      exact h1.trans (SimpleGraph.Adj.reachable (connGraph_adj_cr.mpr ⟨rfl, rfl⟩))

theorem connGraph_connected (t m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) :
    (connGraph t m a n b).Connected := by
  have hall : ∀ x : ConnV t m a n b,
      (connGraph t m a n b).Reachable (chubL t m a n b) x := by
    rintro (p | i | p)
    · exact ((spider_connected m a).preconnected none p).map (cinl t m a n b)
    · exact creachable_interior t m a n b i.val i rfl
    · exact (creachable_hubR t m a n b).trans
        (((spider_connected n b).preconnected none p).map (cinr t m a n b))
  rw [SimpleGraph.connected_iff]
  exact ⟨fun x y => (hall x).symm.trans (hall y), ⟨chubL t m a n b⟩⟩

end ClanAudit
