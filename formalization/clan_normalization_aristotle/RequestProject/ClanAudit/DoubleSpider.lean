import RequestProject.ClanAudit.Arms

/-!
# The adjacent two-hub tree

`dgraph m a n b` is the tree `D(a,b)` of the target: two adjacent hubs `u`, `v`, the
first carrying `m` pendant paths of lengths `a i`, the second carrying `n` pendant paths
of lengths `b j`.

The new clan component of the adjacent two-hub problem is the one containing *both*
hubs.  Here we prove the exact imbalance law asserted in the proof notes:

* `dgraph_connected` — with all multiplicities one the graph is connected;
* `dgraph_imbalance` — its colour classes differ by `|p - q|`, where `p` and `q` are the
  numbers of odd arms at the two hubs;
* `Wpoly_dgraph_ones` — hence its normalized weight is `z^|p-q| + z^(-|p-q|)`.
-/

namespace ClanAudit

open Finset LaurentPolynomial SimpleGraph

/-- Vertices of the adjacent two-hub tree. -/
abbrev DV (m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) := SpiderV m a ⊕ SpiderV n b

/-- The adjacent two-hub tree: two spiders joined by an edge between their hubs. -/
def dgraph (m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) : SimpleGraph (DV m a n b) where
  Adj x y :=
    match x, y with
    | Sum.inl p, Sum.inl q => (spider m a).Adj p q
    | Sum.inr p, Sum.inr q => (spider n b).Adj p q
    | Sum.inl p, Sum.inr q => p = none ∧ q = none
    | Sum.inr p, Sum.inl q => p = none ∧ q = none
  symm := by
    rintro (p | p) (q | q) h
    · exact h.symm
    · exact ⟨h.2, h.1⟩
    · exact ⟨h.2, h.1⟩
    · exact h.symm
  loopless := ⟨by
    rintro (p | p) h
    · exact (spider m a).irrefl h
    · exact (spider n b).irrefl h⟩

variable {m n : ℕ} {a : Fin m → ℕ} {b : Fin n → ℕ}

theorem dgraph_adj_ll {p q : SpiderV m a} :
    (dgraph m a n b).Adj (Sum.inl p) (Sum.inl q) ↔ (spider m a).Adj p q := Iff.rfl

theorem dgraph_adj_rr {p q : SpiderV n b} :
    (dgraph m a n b).Adj (Sum.inr p) (Sum.inr q) ↔ (spider n b).Adj p q := Iff.rfl

theorem dgraph_adj_lr {p : SpiderV m a} {q : SpiderV n b} :
    (dgraph m a n b).Adj (Sum.inl p) (Sum.inr q) ↔ (p = none ∧ q = none) := Iff.rfl

theorem dgraph_adj_rl {p : SpiderV n b} {q : SpiderV m a} :
    (dgraph m a n b).Adj (Sum.inr p) (Sum.inl q) ↔ (p = none ∧ q = none) := Iff.rfl

/-- The distance from the first hub. -/
def dlevel : DV m a n b → ℕ
  | Sum.inl x => spiderLevel x
  | Sum.inr y => spiderLevel y + 1

theorem dlevel_adj {x y : DV m a n b} (h : (dgraph m a n b).Adj x y) :
    dlevel x + 1 = dlevel y ∨ dlevel y + 1 = dlevel x := by
  match x, y with
  | Sum.inl p, Sum.inl q => exact spiderLevel_adj (dgraph_adj_ll.mp h)
  | Sum.inr p, Sum.inr q =>
      have := spiderLevel_adj (dgraph_adj_rr.mp h)
      simp only [dlevel]
      omega
  | Sum.inl p, Sum.inr q =>
      obtain ⟨hp, hq⟩ := dgraph_adj_lr.mp h
      subst hp; subst hq
      exact Or.inl rfl
  | Sum.inr p, Sum.inl q =>
      obtain ⟨hp, hq⟩ := dgraph_adj_rl.mp h
      subst hp; subst hq
      exact Or.inr rfl

/-! ### Connectivity -/

/-- The inclusion of the first spider. -/
def dinl (m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) : spider m a →g dgraph m a n b where
  toFun := Sum.inl
  map_rel' := id

/-- The inclusion of the second spider. -/
def dinr (m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) : spider n b →g dgraph m a n b where
  toFun := Sum.inr
  map_rel' := id

theorem dgraph_connected (m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) :
    (dgraph m a n b).Connected := by
  have hhub : (dgraph m a n b).Reachable (Sum.inl none) (Sum.inr none) :=
    SimpleGraph.Adj.reachable (dgraph_adj_lr.mpr ⟨rfl, rfl⟩)
  have hall : ∀ x : DV m a n b, (dgraph m a n b).Reachable (Sum.inl none) x := by
    rintro (p | p)
    · exact ((spider_connected m a).preconnected none p).map (dinl m a n b)
    · exact hhub.trans (((spider_connected n b).preconnected none p).map (dinr m a n b))
  rw [SimpleGraph.connected_iff]
  exact ⟨fun x y => (hall x).symm.trans (hall y), ⟨Sum.inl none⟩⟩

/-! ### The parity bipartition -/

/-- The colour class of vertices at even distance from the first hub. -/
noncomputable def dEven (m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) :
    Finset (DV m a n b) := Finset.univ.filter (fun x => dlevel x % 2 = 0)

/-- The colour class of vertices at odd distance from the first hub. -/
noncomputable def dOdd (m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) :
    Finset (DV m a n b) := Finset.univ.filter (fun x => dlevel x % 2 = 1)

theorem dgraph_isStablePair (m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) :
    IsStablePair (dgraph m a n b) (dEven m a n b) (dOdd m a n b) := by
  classical
  refine ⟨?_, ?_, ?_, ?_⟩
  · intro x hx y hy hadj
    simp only [dEven, Finset.mem_filter] at hx hy
    have := dlevel_adj hadj
    omega
  · intro x hx y hy hadj
    simp only [dOdd, Finset.mem_filter] at hx hy
    have := dlevel_adj hadj
    omega
  · rw [Finset.disjoint_left]
    intro x hx hx'
    simp only [dEven, dOdd, Finset.mem_filter] at hx hx'
    omega
  · ext x
    simp only [dEven, dOdd, Finset.mem_union, Finset.mem_filter, Finset.mem_univ, true_and,
      iff_true]
    omega

theorem card_dEven (m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) :
    (dEven m a n b).card = (spiderEven m a).card + (spiderOdd n b).card := by
  classical
  rw [dEven, Finset.card_filter, Fintype.sum_sum_type, spiderEven, spiderOdd,
    Finset.card_filter, Finset.card_filter]
  congr 1
  refine Finset.sum_congr rfl fun y _ => ?_
  have : (dlevel (Sum.inr y : DV m a n b)) = spiderLevel y + 1 := rfl
  rw [this]
  by_cases h : spiderLevel y % 2 = 1
  · rw [if_pos (by omega), if_pos h]
  · rw [if_neg (by omega), if_neg h]

theorem card_dOdd (m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) :
    (dOdd m a n b).card = (spiderOdd m a).card + (spiderEven n b).card := by
  classical
  rw [dOdd, Finset.card_filter, Fintype.sum_sum_type, spiderEven, spiderOdd,
    Finset.card_filter, Finset.card_filter]
  congr 1
  refine Finset.sum_congr rfl fun y _ => ?_
  have : (dlevel (Sum.inr y : DV m a n b)) = spiderLevel y + 1 := rfl
  rw [this]
  by_cases h : spiderLevel y % 2 = 0
  · rw [if_pos (by omega), if_pos h]
  · rw [if_neg (by omega), if_neg h]

/-- **The exact imbalance law for a two-hub component.**  Colouring the first hub
positive, the difference of the two colour class sizes is `q - p`, where `p` and `q` are
the numbers of odd arms at the first and the second hub. -/
theorem dgraph_imbalance (m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) :
    ((dEven m a n b).card : ℤ) - (dOdd m a n b).card
      = ((Finset.univ.filter (fun j : Fin n => b j % 2 = 1)).card : ℤ)
        - ((Finset.univ.filter (fun i : Fin m => a i % 2 = 1)).card : ℤ) := by
  classical
  have h1 := spider_imbalance m a
  have h2 := spider_imbalance n b
  rw [card_dEven, card_dOdd]
  push_cast
  omega

/-- **The normalized weight of a two-hub component** is `z^|p-q| + z^(-|p-q|)`. -/
theorem Wpoly_dgraph_ones (m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) (p q r : ℕ)
    (hp : (Finset.univ.filter (fun i : Fin m => a i % 2 = 1)).card = p)
    (hq : (Finset.univ.filter (fun j : Fin n => b j % 2 = 1)).card = q)
    (hr : ((p : ℤ) - q).natAbs = r) :
    Wpoly (dgraph m a n b) (fun _ => 1) = Pw r := by
  classical
  have himb := dgraph_imbalance m a n b
  rw [hp, hq] at himb
  have hgf : imbalanceGF (dgraph m a n b) = Pw r := by
    rw [imbalanceGF_connected (dgraph_connected m a n b) (dgraph_isStablePair m a n b), Pw]
    rcases le_total (p : ℤ) q with h | h
    · have e1 : ((dEven m a n b).card : ℤ) - (dOdd m a n b).card = (r : ℤ) := by omega
      have e2 : ((dOdd m a n b).card : ℤ) - (dEven m a n b).card = -(r : ℤ) := by omega
      rw [e1, e2]
    · have e1 : ((dEven m a n b).card : ℤ) - (dOdd m a n b).card = -(r : ℤ) := by omega
      have e2 : ((dOdd m a n b).card : ℤ) - (dEven m a n b).card = (r : ℤ) := by omega
      rw [e1, e2]
      ring
  rw [Wpoly, clanFactorial_one, imbalanceGF_iso (clanOneIso _), hgf]
  simp

end ClanAudit
