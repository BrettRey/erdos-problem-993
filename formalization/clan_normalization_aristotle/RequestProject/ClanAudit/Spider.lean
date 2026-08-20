import RequestProject.ClanAudit.Weight
import RequestProject.ClanAudit.Blocks

/-!
# The hub component of an active spider state and its imbalance `r = p - 1`

For a spider with hub `w`, a clan map `alpha` with `alpha w = 1` and no neighbour of
multiplicity `≥ 2`, the hub component of the clan graph is the hub together with a
positive prefix of each active arm, all of multiplicity one (see the request,
"Published spider transformation").

This file constructs that component as a graph `spider n len` (a hub with `n` arms, the
`i`-th of length `len i`), and proves:

* `spider_connected` — it is connected;
* `imbalanceGF_spider` — its imbalance polynomial is `z^r + z^(-r)` where `r + 1` is the
  number `p` of arms of odd length;
* `Wpoly_spider` — the same after normalization (the denominator is `1`).

This *derives* the term `z^r + z^(-r)` of the proposed block
`A(r,c;z) = c (z+z⁻¹)^r + z^r + z^(-r)` together with the identification `r = p - 1`
claimed in the request.
-/

namespace ClanAudit

open Finset LaurentPolynomial

/-- Vertices of a spider: the hub (`none`) and position `j` of arm `i`. -/
abbrev SpiderV (n : ℕ) (len : Fin n → ℕ) := Option (Σ i : Fin n, Fin (len i))

/-- A spider: a hub joined to the first vertex of each of `n` pendant paths, the `i`-th
of length `len i`. -/
def spider (n : ℕ) (len : Fin n → ℕ) : SimpleGraph (SpiderV n len) where
  Adj x y :=
    match x, y with
    | none, none => False
    | none, some q => q.2.val = 0
    | some p, none => p.2.val = 0
    | some p, some q => p.1 = q.1 ∧ (p.2.val + 1 = q.2.val ∨ q.2.val + 1 = p.2.val)
  symm := by
    rintro (_ | ⟨i, j⟩) (_ | ⟨i', j'⟩) h
    · exact h
    · exact h
    · exact h
    · exact ⟨h.1.symm, h.2.symm⟩
  loopless := ⟨by
    rintro (_ | ⟨i, j⟩) h
    · exact h
    · rcases h.2 with h' | h' <;> omega⟩

theorem spider_adj_hub {n : ℕ} {len : Fin n → ℕ} {q : Σ i : Fin n, Fin (len i)}
    (hq : q.2.val = 0) : (spider n len).Adj none (some q) := hq

theorem spider_adj_arm {n : ℕ} {len : Fin n → ℕ} {i : Fin n} {j j' : Fin (len i)}
    (h : j.val + 1 = j'.val) : (spider n len).Adj (some ⟨i, j⟩) (some ⟨i, j'⟩) :=
  ⟨rfl, Or.inl h⟩

/-- The distance of a vertex from the hub. -/
def spiderLevel {n : ℕ} {len : Fin n → ℕ} : SpiderV n len → ℕ
  | none => 0
  | some p => p.2.val + 1

theorem spiderLevel_adj {n : ℕ} {len : Fin n → ℕ} {x y : SpiderV n len}
    (h : (spider n len).Adj x y) :
    spiderLevel x + 1 = spiderLevel y ∨ spiderLevel y + 1 = spiderLevel x := by
  match x, y with
  | none, none => exact absurd h (by simp [spider])
  | none, some q => exact Or.inl (by simp [spiderLevel, show q.2.val = 0 from h])
  | some p, none => exact Or.inr (by simp [spiderLevel, show p.2.val = 0 from h])
  | some p, some q =>
      simp only [spiderLevel]
      rcases h.2 with h' | h' <;> omega

/-! ### Connectivity -/

theorem spider_connected (n : ℕ) (len : Fin n → ℕ) : (spider n len).Connected := by
  have key : ∀ (k : ℕ) (i : Fin n) (j : Fin (len i)), j.val = k →
      (spider n len).Reachable none (some ⟨i, j⟩) := by
    intro k
    induction k with
    | zero => intro i j hj; exact (spider_adj_hub (q := ⟨i, j⟩) hj).reachable
    | succ k ih =>
        intro i j hj
        have hk : k < len i := by omega
        refine (ih i ⟨k, hk⟩ rfl).trans (SimpleGraph.Adj.reachable ?_)
        exact spider_adj_arm (j := ⟨k, hk⟩) (j' := j) (by simpa using hj.symm)
  have hall : ∀ x : SpiderV n len, (spider n len).Reachable none x := by
    rintro (_ | ⟨i, j⟩)
    · rfl
    · exact key j.val i j rfl
  rw [SimpleGraph.connected_iff]
  exact ⟨fun x y => ((hall x).symm).trans (hall y), ⟨none⟩⟩

/-! ### The parity bipartition -/

/-- The colour class of vertices at even distance from the hub (including the hub). -/
noncomputable def spiderEven (n : ℕ) (len : Fin n → ℕ) : Finset (SpiderV n len) :=
  Finset.univ.filter (fun x => spiderLevel x % 2 = 0)

/-- The colour class of vertices at odd distance from the hub. -/
noncomputable def spiderOdd (n : ℕ) (len : Fin n → ℕ) : Finset (SpiderV n len) :=
  Finset.univ.filter (fun x => spiderLevel x % 2 = 1)

theorem spider_isStablePair (n : ℕ) (len : Fin n → ℕ) :
    IsStablePair (spider n len) (spiderEven n len) (spiderOdd n len) := by
  classical
  refine ⟨?_, ?_, ?_, ?_⟩
  · intro x hx y hy hadj
    simp only [spiderEven, Finset.mem_filter] at hx hy
    have := spiderLevel_adj hadj
    omega
  · intro x hx y hy hadj
    simp only [spiderOdd, Finset.mem_filter] at hx hy
    have := spiderLevel_adj hadj
    omega
  · rw [Finset.disjoint_left]
    intro x hx hx'
    simp only [spiderEven, spiderOdd, Finset.mem_filter] at hx hx'
    omega
  · ext x
    simp only [spiderEven, spiderOdd, Finset.mem_union, Finset.mem_filter, Finset.mem_univ,
      true_and, iff_true]
    omega

/-! ### Counting the two classes -/

private theorem card_filter_even_range (m : ℕ) :
    ((Finset.range m).filter (fun j => j % 2 = 0)).card = (m + 1) / 2 := by
  induction m with
  | zero => simp
  | succ m ih =>
      rw [Finset.range_add_one, Finset.filter_insert]
      by_cases h : m % 2 = 0
      · rw [if_pos h, Finset.card_insert_of_notMem (by simp), ih]
        omega
      · rw [if_neg h, ih]
        omega

private theorem card_filter_odd_range (m : ℕ) :
    ((Finset.range m).filter (fun j => j % 2 = 1)).card = m / 2 := by
  induction m with
  | zero => simp
  | succ m ih =>
      rw [Finset.range_add_one, Finset.filter_insert]
      by_cases h : m % 2 = 1
      · rw [if_pos h, Finset.card_insert_of_notMem (by simp), ih]
        omega
      · rw [if_neg h, ih]
        omega

private theorem card_filter_fin (m : ℕ) (b : ℕ) :
    ((Finset.univ : Finset (Fin m)).filter (fun j => j.val % 2 = b)).card
      = ((Finset.range m).filter (fun j => j % 2 = b)).card := by
  classical
  rw [Finset.card_filter, Finset.card_filter, Fin.sum_univ_eq_sum_range (fun j => if j % 2 = b then 1 else 0)]

private theorem card_filter_sigma (n : ℕ) (len : Fin n → ℕ) (b : ℕ) :
    ((Finset.univ : Finset (Σ i : Fin n, Fin (len i))).filter (fun p => p.2.val % 2 = b)).card
      = ∑ i : Fin n, ((Finset.univ : Finset (Fin (len i))).filter (fun j => j.val % 2 = b)).card := by
  classical
  rw [Finset.card_filter]
  rw [← Finset.univ_sigma_univ, Finset.sum_sigma]
  exact Finset.sum_congr rfl fun i _ => (Finset.card_filter _ _).symm

theorem card_spiderEven (n : ℕ) (len : Fin n → ℕ) :
    (spiderEven n len).card = 1 + ∑ i : Fin n, (len i) / 2 := by
  classical
  have h : (spiderEven n len).card
      = 1 + ((Finset.univ : Finset (Σ i : Fin n, Fin (len i))).filter
          (fun p => p.2.val % 2 = 1)).card := by
    rw [spiderEven]
    rw [show (Finset.univ.filter (fun x : SpiderV n len => spiderLevel x % 2 = 0))
        = insert none (((Finset.univ : Finset (Σ i : Fin n, Fin (len i))).filter
            (fun p => p.2.val % 2 = 1)).image some) from ?_]
    · rw [Finset.card_insert_of_notMem (by simp), Finset.card_image_of_injective _
        (Option.some_injective _)]
      omega
    · ext x
      match x with
      | none =>
          constructor
          · intro _
            exact Finset.mem_insert_self _ _
          · intro _
            exact Finset.mem_filter.mpr ⟨Finset.mem_univ _, rfl⟩
      | some p =>
          constructor
          · intro hx
            have h0 : (p.2.val + 1) % 2 = 0 := (Finset.mem_filter.mp hx).2
            exact Finset.mem_insert_of_mem (Finset.mem_image_of_mem _
              (Finset.mem_filter.mpr ⟨Finset.mem_univ _, by omega⟩))
          · intro hx
            rcases Finset.mem_insert.mp hx with h | h
            · exact absurd h (by simp)
            · obtain ⟨q, hq, hqp⟩ := Finset.mem_image.mp h
              have hqp' : q = p := Option.some_injective _ hqp
              subst hqp'
              have h1 : q.2.val % 2 = 1 := (Finset.mem_filter.mp hq).2
              exact Finset.mem_filter.mpr ⟨Finset.mem_univ _,
                show (q.2.val + 1) % 2 = 0 by omega⟩
  rw [h, card_filter_sigma]
  congr 1
  exact Finset.sum_congr rfl fun i _ => by
    rw [card_filter_fin, card_filter_odd_range]

theorem card_spiderOdd (n : ℕ) (len : Fin n → ℕ) :
    (spiderOdd n len).card = ∑ i : Fin n, (len i + 1) / 2 := by
  classical
  have h : (spiderOdd n len).card
      = ((Finset.univ : Finset (Σ i : Fin n, Fin (len i))).filter
          (fun p => p.2.val % 2 = 0)).card := by
    rw [spiderOdd]
    rw [show (Finset.univ.filter (fun x : SpiderV n len => spiderLevel x % 2 = 1))
        = (((Finset.univ : Finset (Σ i : Fin n, Fin (len i))).filter
            (fun p => p.2.val % 2 = 0)).image some) from ?_]
    · rw [Finset.card_image_of_injective _ (Option.some_injective _)]
    · ext x
      match x with
      | none =>
          constructor
          · intro hx
            have h0 : (0 : ℕ) % 2 = 1 := (Finset.mem_filter.mp hx).2
            omega
          · intro hx
            obtain ⟨q, -, hq⟩ := Finset.mem_image.mp hx
            exact absurd hq (by simp)
      | some p =>
          constructor
          · intro hx
            have h0 : (p.2.val + 1) % 2 = 1 := (Finset.mem_filter.mp hx).2
            exact Finset.mem_image_of_mem _
              (Finset.mem_filter.mpr ⟨Finset.mem_univ _, by omega⟩)
          · intro hx
            obtain ⟨q, hq, hqp⟩ := Finset.mem_image.mp hx
            have hqp' : q = p := Option.some_injective _ hqp
            subst hqp'
            have h1 : q.2.val % 2 = 0 := (Finset.mem_filter.mp hq).2
            exact Finset.mem_filter.mpr ⟨Finset.mem_univ _,
              show (q.2.val + 1) % 2 = 1 by omega⟩
  rw [h, card_filter_sigma]
  exact Finset.sum_congr rfl fun i _ => by
    rw [card_filter_fin, card_filter_even_range]

/-- The imbalance of the hub component is `p - 1`, where `p` is the number of arms of
odd length. -/
theorem spider_imbalance (n : ℕ) (len : Fin n → ℕ) :
    ((spiderOdd n len).card : ℤ) - (spiderEven n len).card
      = ((Finset.univ.filter (fun i : Fin n => len i % 2 = 1)).card : ℤ) - 1 := by
  classical
  rw [card_spiderEven, card_spiderOdd, Finset.card_filter]
  have key : ∀ i : Fin n, (len i + 1) / 2 = len i / 2 + (if len i % 2 = 1 then 1 else 0) := by
    intro i
    by_cases h : len i % 2 = 1
    · rw [if_pos h]; omega
    · rw [if_neg h]; omega
  rw [Finset.sum_congr rfl fun i _ => key i, Finset.sum_add_distrib]
  push_cast
  ring

/-- **The hub component of an active spider state has imbalance polynomial
`z^r + z^(-r)`, with `r = p - 1`.** -/
theorem imbalanceGF_spider (n : ℕ) (len : Fin n → ℕ) (r : ℕ)
    (hp : (Finset.univ.filter (fun i : Fin n => len i % 2 = 1)).card = r + 1) :
    imbalanceGF (spider n len) = Pw r := by
  classical
  have himb := spider_imbalance n len
  rw [hp] at himb
  have h1 : ((spiderOdd n len).card : ℤ) - (spiderEven n len).card = (r : ℤ) := by
    rw [himb]; push_cast; ring
  have h2 : ((spiderEven n len).card : ℤ) - (spiderOdd n len).card = -(r : ℤ) := by omega
  rw [imbalanceGF_connected (spider_connected n len) (spider_isStablePair n len), h1, h2, Pw]
  ring

/-! ### The normalized weight of the active hub component -/

variable {V : Type*} [Fintype V] [DecidableEq V]

/-- With all multiplicities one, the clan graph is the graph itself. -/
def clanOneIso (G : SimpleGraph V) : clan G (fun _ => 1) ≃g G where
  toFun x := x.1
  invFun v := ⟨v, ⟨0, by norm_num⟩⟩
  left_inv := by
    rintro ⟨v, k⟩
    have : k = ⟨0, by norm_num⟩ := Subsingleton.elim _ _
    rw [this]
  right_inv := by intro v; rfl
  map_rel_iff' := by
    rintro ⟨u, i⟩ ⟨v, j⟩
    constructor
    · intro h
      exact Or.inr h
    · rintro (⟨h1, h2⟩ | h)
      · exact absurd (by subst h1; rw [Subsingleton.elim i j]) h2
      · exact h

omit [DecidableEq V] in
theorem clanFactorial_one : clanFactorial (fun _ : V => 1) = 1 := by
  simp [clanFactorial]

/-- **The normalized weight of the active hub component is `z^r + z^(-r)`.** -/
theorem Wpoly_spider (n : ℕ) (len : Fin n → ℕ) (r : ℕ)
    (hp : (Finset.univ.filter (fun i : Fin n => len i % 2 = 1)).card = r + 1) :
    Wpoly (spider n len) (fun _ => 1) = Pw r := by
  rw [Wpoly, clanFactorial_one, imbalanceGF_iso (clanOneIso (spider n len)),
    imbalanceGF_spider n len r hp]
  simp

end ClanAudit
