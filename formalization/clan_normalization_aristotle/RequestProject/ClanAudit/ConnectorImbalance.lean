import RequestProject.ClanAudit.Connector

/-!
# The exact imbalance law of the connector component

For the connector tree `connGraph t m a n b` (two hubs joined by a path with `t` interior
vertices, i.e. `ell = t + 1` connector edges) with all multiplicities one, the clan graph
is the tree itself; it is connected and bipartite, and the difference of its two colour
classes is computed here exactly:

* `connGraph_imbalance_ell_odd` (`t` even, `ell = t+1` odd): the imbalance is `q - p`;
* `connGraph_imbalance_ell_even` (`t` odd, `ell` even): the imbalance is `1 - p - q`;

where `p` and `q` are the numbers of odd arms at the two hubs.  Consequently

* `Wpoly_connGraph_ones_ell_odd` — the normalized weight is `z^|p-q| + z^(-|p-q|)`;
* `Wpoly_connGraph_ones_ell_even` — it is `z^|1-p-q| + z^(-|1-p-q|)`.

This is item 3 of the connector target: the joined component imbalance is `|p-q|` for odd
`ell` and `|1-p-q|` for even `ell`, with no restriction on `p`, `q` (the endpoint,
zero-active and one-active cases are all included, since `p` and `q` are arbitrary).
-/

namespace ClanAudit

open Finset LaurentPolynomial SimpleGraph

open scoped Classical

/-! ### A general weight formula for a connected bipartite clan graph of multiplicity one -/

/-- If a connected graph has a stable two-block partition with colour classes differing by
`r` in absolute value, then its normalized weight at multiplicity one is `z^r + z^(-r)`. -/
theorem Wpoly_ones_of_imbalance {V : Type*} [Fintype V] [DecidableEq V] {G : SimpleGraph V}
    (hc : G.Connected) {A B : Finset V} (h : IsStablePair G A B) (r : ℕ)
    (hr : ((A.card : ℤ) - B.card).natAbs = r) :
    Wpoly G (fun _ => 1) = Pw r := by
  have hgf : imbalanceGF G = Pw r := by
    rw [imbalanceGF_connected hc h, Pw]
    rcases le_total (A.card : ℤ) B.card with hle | hle
    · have e1 : ((A.card : ℤ) - B.card) = -(r : ℤ) := by omega
      have e2 : ((B.card : ℤ) - A.card) = (r : ℤ) := by omega
      rw [e1, e2, add_comm]
    · have e1 : ((A.card : ℤ) - B.card) = (r : ℤ) := by omega
      have e2 : ((B.card : ℤ) - A.card) = -(r : ℤ) := by omega
      rw [e1, e2]
  rw [Wpoly, clanFactorial_one, imbalanceGF_iso (clanOneIso _), hgf]
  simp

/-! ### The parity bipartition of the connector tree -/

variable {t m n : ℕ} {a : Fin m → ℕ} {b : Fin n → ℕ}

/-- The colour class of vertices at even distance from the left hub. -/
noncomputable def cEven (t m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) :
    Finset (ConnV t m a n b) := Finset.univ.filter (fun x => clevel x % 2 = 0)

/-- The colour class of vertices at odd distance from the left hub. -/
noncomputable def cOdd (t m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) :
    Finset (ConnV t m a n b) := Finset.univ.filter (fun x => clevel x % 2 = 1)

theorem connGraph_isStablePair (t m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) :
    IsStablePair (connGraph t m a n b) (cEven t m a n b) (cOdd t m a n b) := by
  refine ⟨?_, ?_, ?_, ?_⟩
  · intro x hx y hy hadj
    simp only [cEven, Finset.mem_filter] at hx hy
    have := clevel_adj hadj
    omega
  · intro x hx y hy hadj
    simp only [cOdd, Finset.mem_filter] at hx hy
    have := clevel_adj hadj
    omega
  · rw [Finset.disjoint_left]
    intro x hx hx'
    simp only [cEven, cOdd, Finset.mem_filter] at hx hx'
    omega
  · ext x
    simp only [cEven, cOdd, Finset.mem_union, Finset.mem_filter, Finset.mem_univ, true_and,
      iff_true]
    omega

/-! ### Counting the two colour classes -/

theorem card_interior_even (t : ℕ) :
    (Finset.univ.filter (fun i : Fin t => (i.val + 1) % 2 = 0)).card = t / 2 := by
  rw [Finset.card_filter,
    Fin.sum_univ_eq_sum_range (fun i => if (i + 1) % 2 = 0 then 1 else 0) t]
  induction t with
  | zero => simp
  | succ s ih =>
      rw [Finset.sum_range_succ, ih]
      split_ifs with h <;> omega

theorem card_interior_odd (t : ℕ) :
    (Finset.univ.filter (fun i : Fin t => (i.val + 1) % 2 = 1)).card = (t + 1) / 2 := by
  rw [Finset.card_filter,
    Fin.sum_univ_eq_sum_range (fun i => if (i + 1) % 2 = 1 then 1 else 0) t]
  induction t with
  | zero => simp
  | succ s ih =>
      rw [Finset.sum_range_succ, ih]
      split_ifs with h <;> omega

theorem card_cEven (t m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) :
    (cEven t m a n b).card
      = (spiderEven m a).card + t / 2
        + (if t % 2 = 1 then (spiderEven n b).card else (spiderOdd n b).card) := by
  rw [cEven, Finset.card_filter, Fintype.sum_sum_type, Fintype.sum_sum_type]
  have h1 : ∑ x : SpiderV m a, (if clevel (Sum.inl x : ConnV t m a n b) % 2 = 0 then 1 else 0)
      = (spiderEven m a).card := by
    rw [spiderEven, Finset.card_filter]
    rfl
  have h2 : ∑ i : Fin t, (if clevel (Sum.inr (Sum.inl i) : ConnV t m a n b) % 2 = 0 then 1 else 0)
      = t / 2 := by
    rw [← card_interior_even t, Finset.card_filter]
    rfl
  have h3 : ∑ y : SpiderV n b,
        (if clevel (Sum.inr (Sum.inr y) : ConnV t m a n b) % 2 = 0 then 1 else 0)
      = (if t % 2 = 1 then (spiderEven n b).card else (spiderOdd n b).card) := by
    rw [spiderEven, spiderOdd, Finset.card_filter, Finset.card_filter]
    by_cases ht : t % 2 = 1
    · rw [if_pos ht]
      refine Finset.sum_congr rfl fun y _ => ?_
      have hl : clevel (Sum.inr (Sum.inr y) : ConnV t m a n b) = spiderLevel y + t + 1 := rfl
      rw [hl]
      by_cases h : spiderLevel y % 2 = 0
      · rw [if_pos (by omega), if_pos h]
      · rw [if_neg (by omega), if_neg h]
    · rw [if_neg ht]
      refine Finset.sum_congr rfl fun y _ => ?_
      have hl : clevel (Sum.inr (Sum.inr y) : ConnV t m a n b) = spiderLevel y + t + 1 := rfl
      rw [hl]
      by_cases h : spiderLevel y % 2 = 1
      · rw [if_pos (by omega), if_pos h]
      · rw [if_neg (by omega), if_neg h]
  rw [h1, h2, h3]
  exact (add_assoc _ _ _).symm

theorem card_cOdd (t m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ) :
    (cOdd t m a n b).card
      = (spiderOdd m a).card + (t + 1) / 2
        + (if t % 2 = 1 then (spiderOdd n b).card else (spiderEven n b).card) := by
  rw [cOdd, Finset.card_filter, Fintype.sum_sum_type, Fintype.sum_sum_type]
  have h1 : ∑ x : SpiderV m a, (if clevel (Sum.inl x : ConnV t m a n b) % 2 = 1 then 1 else 0)
      = (spiderOdd m a).card := by
    rw [spiderOdd, Finset.card_filter]
    rfl
  have h2 : ∑ i : Fin t, (if clevel (Sum.inr (Sum.inl i) : ConnV t m a n b) % 2 = 1 then 1 else 0)
      = (t + 1) / 2 := by
    rw [← card_interior_odd t, Finset.card_filter]
    rfl
  have h3 : ∑ y : SpiderV n b,
        (if clevel (Sum.inr (Sum.inr y) : ConnV t m a n b) % 2 = 1 then 1 else 0)
      = (if t % 2 = 1 then (spiderOdd n b).card else (spiderEven n b).card) := by
    rw [spiderEven, spiderOdd, Finset.card_filter, Finset.card_filter]
    by_cases ht : t % 2 = 1
    · rw [if_pos ht]
      refine Finset.sum_congr rfl fun y _ => ?_
      have hl : clevel (Sum.inr (Sum.inr y) : ConnV t m a n b) = spiderLevel y + t + 1 := rfl
      rw [hl]
      by_cases h : spiderLevel y % 2 = 1
      · rw [if_pos (by omega), if_pos h]
      · rw [if_neg (by omega), if_neg h]
    · rw [if_neg ht]
      refine Finset.sum_congr rfl fun y _ => ?_
      have hl : clevel (Sum.inr (Sum.inr y) : ConnV t m a n b) = spiderLevel y + t + 1 := rfl
      rw [hl]
      by_cases h : spiderLevel y % 2 = 0
      · rw [if_pos (by omega), if_pos h]
      · rw [if_neg (by omega), if_neg h]
  rw [h1, h2, h3]
  exact (add_assoc _ _ _).symm

/-! ### The imbalance law -/

/-- **Odd connector (`ell = t + 1` odd, i.e. `t` even): the imbalance is `q - p`.** -/
theorem connGraph_imbalance_ell_odd (t m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ)
    (ht : t % 2 = 0) :
    ((cEven t m a n b).card : ℤ) - (cOdd t m a n b).card
      = ((Finset.univ.filter (fun j : Fin n => b j % 2 = 1)).card : ℤ)
        - ((Finset.univ.filter (fun i : Fin m => a i % 2 = 1)).card : ℤ) := by
  have h1 := spider_imbalance m a
  have h2 := spider_imbalance n b
  rw [card_cEven, card_cOdd, if_neg (by omega), if_neg (by omega)]
  have h3 : t / 2 = (t + 1) / 2 := by omega
  push_cast
  omega

/-- **Even connector (`ell = t + 1` even, i.e. `t` odd): the imbalance is `1 - p - q`.** -/
theorem connGraph_imbalance_ell_even (t m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ)
    (ht : t % 2 = 1) :
    ((cEven t m a n b).card : ℤ) - (cOdd t m a n b).card
      = 1 - ((Finset.univ.filter (fun i : Fin m => a i % 2 = 1)).card : ℤ)
        - ((Finset.univ.filter (fun j : Fin n => b j % 2 = 1)).card : ℤ) := by
  have h1 := spider_imbalance m a
  have h2 := spider_imbalance n b
  rw [card_cEven, card_cOdd, if_pos ht, if_pos ht]
  have h3 : (t + 1) / 2 = t / 2 + 1 := by omega
  push_cast
  omega

/-! ### The normalized weight of the joined component -/

/-- **Odd connector: the normalized weight of the joined component is
`z^|p-q| + z^(-|p-q|)`.** -/
theorem Wpoly_connGraph_ones_ell_odd (t m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ)
    (ht : t % 2 = 0) (p q r : ℕ)
    (hp : (Finset.univ.filter (fun i : Fin m => a i % 2 = 1)).card = p)
    (hq : (Finset.univ.filter (fun j : Fin n => b j % 2 = 1)).card = q)
    (hr : ((p : ℤ) - q).natAbs = r) :
    Wpoly (connGraph t m a n b) (fun _ => 1) = Pw r := by
  refine Wpoly_ones_of_imbalance (connGraph_connected t m a n b)
    (connGraph_isStablePair t m a n b) r ?_
  rw [connGraph_imbalance_ell_odd t m a n b ht, hp, hq]
  omega

/-- **Even connector: the normalized weight of the joined component is
`z^|1-p-q| + z^(-|1-p-q|)`.** -/
theorem Wpoly_connGraph_ones_ell_even (t m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ)
    (ht : t % 2 = 1) (p q r : ℕ)
    (hp : (Finset.univ.filter (fun i : Fin m => a i % 2 = 1)).card = p)
    (hq : (Finset.univ.filter (fun j : Fin n => b j % 2 = 1)).card = q)
    (hr : ((1 : ℤ) - p - q).natAbs = r) :
    Wpoly (connGraph t m a n b) (fun _ => 1) = Pw r := by
  refine Wpoly_ones_of_imbalance (connGraph_connected t m a n b)
    (connGraph_isStablePair t m a n b) r ?_
  rw [connGraph_imbalance_ell_even t m a n b ht, hp, hq]
  omega

end ClanAudit
