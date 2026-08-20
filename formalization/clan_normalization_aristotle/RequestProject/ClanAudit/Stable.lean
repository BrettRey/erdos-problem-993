import Mathlib
import RequestProject.ClanAudit.Laurent

/-!
# Stable two-block partitions and the imbalance generating function

For a finite simple graph `H` we count *ordered* pairs `(A,B)` of independent sets
partitioning `V(H)`.  This is exactly the "semi-ordered" convention of the request:
for `k ≠ l` an ordered pair is determined by the unordered partition, while for
`k = l` the two equal-sized blocks are ordered, i.e. counted twice.

The imbalance generating function `∑ z^(|A|-|B|)` packages all of these counts at
once; `stableCount H k l` is one of its coefficients, and the two-row coefficient of
Li–Li–Yang–Zhang Proposition 2.4 is the difference of two consecutive coefficients.
-/

namespace ClanAudit

open Finset LaurentPolynomial

variable {W W₁ W₂ : Type*} [Fintype W] [DecidableEq W] [Fintype W₁] [DecidableEq W₁]
  [Fintype W₂] [DecidableEq W₂]

/-- An independent (stable) set of vertices. -/
def IsIndep (H : SimpleGraph W) (A : Finset W) : Prop := ∀ x ∈ A, ∀ y ∈ A, ¬ H.Adj x y

/-- `(A, B)` is an ordered stable two-block partition of `H`.  Empty blocks are allowed. -/
def IsStablePair (H : SimpleGraph W) (A B : Finset W) : Prop :=
  IsIndep H A ∧ IsIndep H B ∧ Disjoint A B ∧ A ∪ B = Finset.univ

open Classical in
/-- All ordered stable two-block partitions of `H`. -/
noncomputable def stablePairs (H : SimpleGraph W) : Finset (Finset W × Finset W) :=
  Finset.univ.filter (fun p => IsStablePair H p.1 p.2)

theorem mem_stablePairs {H : SimpleGraph W} {p : Finset W × Finset W} :
    p ∈ stablePairs H ↔ IsStablePair H p.1 p.2 := by
  classical
  simp [stablePairs]

open Classical in
/-- `stableCount H k l` is the number of ordered stable two-block partitions with
`|A| = k` and `|B| = l`. -/
noncomputable def stableCount (H : SimpleGraph W) (k l : ℕ) : ℕ :=
  ((stablePairs H).filter (fun p => p.1.card = k ∧ p.2.card = l)).card

/-- The imbalance generating function `∑ z^(|A|-|B|)` over ordered stable two-block
partitions. -/
noncomputable def imbalanceGF (H : SimpleGraph W) : LaurentPolynomial ℚ :=
  ∑ p ∈ stablePairs H, T ((p.1.card : ℤ) - p.2.card)

theorem card_add_card_of_stablePair {H : SimpleGraph W} {A B : Finset W}
    (h : IsStablePair H A B) : A.card + B.card = Fintype.card W := by
  rw [← Finset.card_union_of_disjoint h.2.2.1, h.2.2.2, Finset.card_univ]

/-- In a stable two-block partition the second block is the complement of the first. -/
theorem mem_snd_iff {H : SimpleGraph W} {A B : Finset W} (h : IsStablePair H A B) (w : W) :
    w ∈ B ↔ w ∉ A := by
  obtain ⟨-, -, hd, hu⟩ := h
  constructor
  · intro hw hwA; exact (Finset.disjoint_left.mp hd hwA) hw
  · intro hw
    have : w ∈ A ∪ B := by rw [hu]; exact Finset.mem_univ w
    rcases Finset.mem_union.mp this with h1 | h1
    · exact absurd h1 hw
    · exact h1

/-- Along an edge, the two blocks alternate. -/
theorem mem_fst_iff_adj {H : SimpleGraph W} {A B : Finset W} (h : IsStablePair H A B)
    {u v : W} (huv : H.Adj u v) : (u ∈ A ↔ v ∉ A) := by
  have hA := h.1
  have hB := h.2.1
  constructor
  · intro huA hvA; exact hA u huA v hvA huv
  · intro hvA
    by_contra huA
    exact hB u ((mem_snd_iff h u).mpr huA) v ((mem_snd_iff h v).mpr hvA) huv

/-- The coefficients of the imbalance generating function are exactly the counts
`stableCount`. -/
theorem coeffL_imbalanceGF (H : SimpleGraph W) (k l : ℕ) (h : k + l = Fintype.card W) :
    coeffL (imbalanceGF H) ((k : ℤ) - l) = (stableCount H k l : ℚ) := by
  classical
  rw [imbalanceGF, coeffL_sum, stableCount, Finset.card_filter]
  push_cast
  refine Finset.sum_congr rfl fun p hp => ?_
  have hp' := card_add_card_of_stablePair (mem_stablePairs.mp hp)
  rw [coeffL_T]
  by_cases hc : p.1.card = k ∧ p.2.card = l
  · rw [if_pos (by omega), if_pos hc]
  · rw [if_neg (by omega), if_neg hc]

/-- The two-row coefficient of Li–Li–Yang–Zhang Proposition 2.4, taken as a definition
(as the request permits). -/
noncomputable def twoRowCoeff (H : SimpleGraph W) (k l : ℕ) : ℤ :=
  (stableCount H k l : ℤ) - (stableCount H (k + 1) (l - 1) : ℤ)

/-- The two-row coefficient is a difference of two coefficients of the imbalance
generating function, two apart. -/
theorem twoRowCoeff_eq_coeff (H : SimpleGraph W) (k l : ℕ) (hl : 1 ≤ l)
    (h : k + l = Fintype.card W) :
    (twoRowCoeff H k l : ℚ)
      = coeffL (imbalanceGF H) ((k : ℤ) - l) - coeffL (imbalanceGF H) ((k : ℤ) - l + 2) := by
  rw [twoRowCoeff, coeffL_imbalanceGF H k l h]
  have hcast : ((k : ℤ) - l + 2) = ((k + 1 : ℕ) : ℤ) - ((l - 1 : ℕ) : ℤ) := by
    push_cast [Nat.cast_sub hl]; ring
  rw [hcast, coeffL_imbalanceGF H (k + 1) (l - 1) (by omega)]
  push_cast
  ring

/-! ### Bipartiteness -/

/-- A graph admits a stable two-block partition iff it is 2-colorable. -/
theorem stablePairs_nonempty_iff (H : SimpleGraph W) :
    (stablePairs H).Nonempty ↔ H.Colorable 2 := by
  classical
  constructor
  · rintro ⟨p, hp⟩
    obtain ⟨hA, hB, hd, hu⟩ := mem_stablePairs.mp hp
    refine ⟨SimpleGraph.Coloring.mk (fun w => if w ∈ p.1 then (0 : Fin 2) else 1) ?_⟩
    intro u v huv
    have hstable : IsStablePair H p.1 p.2 := ⟨hA, hB, hd, hu⟩
    have halt := mem_fst_iff_adj hstable huv
    by_cases hu1 : u ∈ p.1
    · have hv1 : v ∉ p.1 := halt.mp hu1
      simp [hu1, hv1]
    · have hv1 : v ∈ p.1 := by
        by_contra hv1
        exact hu1 (halt.mpr hv1)
      simp [hu1, hv1]
  · rintro ⟨c⟩
    refine ⟨(Finset.univ.filter (fun w => c w = 0), Finset.univ.filter (fun w => c w ≠ 0)), ?_⟩
    rw [mem_stablePairs]
    refine ⟨?_, ?_, ?_, ?_⟩
    · intro x hx y hy hxy
      simp only [Finset.mem_filter] at hx hy
      exact c.valid hxy (by rw [hx.2, hy.2])
    · intro x hx y hy hxy
      simp only [Finset.mem_filter] at hx hy
      refine c.valid hxy ?_
      have h1 : c x = 1 := by omega
      have h2 : c y = 1 := by omega
      rw [h1, h2]
    · rw [Finset.disjoint_left]
      intro a ha ha'
      simp only [Finset.mem_filter] at ha ha'
      exact ha'.2 ha.2
    · ext a
      simp only [Finset.mem_union, Finset.mem_filter, Finset.mem_univ, true_and, iff_true]
      exact em _

/-- A non-2-colorable (i.e. non-bipartite) graph has vanishing imbalance generating
function, hence all its two-row coefficients vanish. -/
theorem imbalanceGF_eq_zero_of_not_colorable {H : SimpleGraph W} (h : ¬ H.Colorable 2) :
    imbalanceGF H = 0 := by
  have : stablePairs H = ∅ := by
    by_contra hne
    exact h ((stablePairs_nonempty_iff H).mp (Finset.nonempty_iff_ne_empty.mpr hne))
  rw [imbalanceGF, this, Finset.sum_empty]

/-- A graph containing a triangle has vanishing imbalance generating function.
This is the "a neighbour of the hub has multiplicity ≥ 2 forces a triangle, so the
two-row contribution is zero" step of the spider transformation. -/
theorem imbalanceGF_eq_zero_of_triangle {H : SimpleGraph W} {x y z : W}
    (hxy : H.Adj x y) (hyz : H.Adj y z) (hxz : H.Adj x z) : imbalanceGF H = 0 := by
  have hempty : stablePairs H = ∅ := by
    rw [Finset.eq_empty_iff_forall_notMem]
    rintro ⟨A, B⟩ hp
    have hst : IsStablePair H A B := mem_stablePairs.mp hp
    have hA := hst.1
    have hB := hst.2.1
    by_cases hx : x ∈ A <;> by_cases hy : y ∈ A <;> by_cases hz : z ∈ A
    · exact hA x hx y hy hxy
    · exact hA x hx y hy hxy
    · exact hA x hx z hz hxz
    · exact hB y ((mem_snd_iff hst y).mpr hy) z ((mem_snd_iff hst z).mpr hz) hyz
    · exact hA y hy z hz hyz
    · exact hB x ((mem_snd_iff hst x).mpr hx) z ((mem_snd_iff hst z).mpr hz) hxz
    · exact hB x ((mem_snd_iff hst x).mpr hx) y ((mem_snd_iff hst y).mpr hy) hxy
    · exact hB x ((mem_snd_iff hst x).mpr hx) y ((mem_snd_iff hst y).mpr hy) hxy
  rw [imbalanceGF, hempty, Finset.sum_empty]

/-! ### Connected bipartite graphs -/

/-- For a connected graph, a stable two-block partition is unique up to swapping the
two blocks. -/
theorem stablePair_unique {H : SimpleGraph W} (hc : H.Connected) {A B A' B' : Finset W}
    (h : IsStablePair H A B) (h' : IsStablePair H A' B') :
    (A' = A ∧ B' = B) ∨ (A' = B ∧ B' = A) := by
  have step : ∀ u v : W, H.Adj u v → ((u ∈ A ↔ u ∈ A') ↔ (v ∈ A ↔ v ∈ A')) := by
    intro u v huv
    have e1 := mem_fst_iff_adj h huv
    have e2 := mem_fst_iff_adj h' huv
    constructor
    · intro hP
      constructor
      · intro hv
        by_contra hv'
        exact absurd (e2.mpr hv') (by tauto)
      · intro hv'
        by_contra hv
        exact absurd (e1.mpr hv) (by tauto)
    · intro hP
      constructor
      · intro hu
        by_contra hu'
        have hv : v ∉ A := e1.mp hu
        have hv' : v ∈ A' := by
          by_contra hv'
          exact hu' (e2.mpr hv')
        exact hv (hP.mpr hv')
      · intro hu'
        by_contra hu
        have hv : v ∈ A := by
          by_contra hv
          exact hu (e1.mpr hv)
        have hv' : v ∉ A' := e2.mp hu'
        exact hv' (hP.mp hv)
  have const : ∀ u v : W, ((u ∈ A ↔ u ∈ A') ↔ (v ∈ A ↔ v ∈ A')) := by
    intro u v
    obtain ⟨w⟩ := hc.preconnected u v
    induction w with
    | nil => exact Iff.rfl
    | cons hadj p ih => exact (step _ _ hadj).trans ih
  by_cases hex : ∃ w : W, (w ∈ A ↔ w ∈ A')
  · obtain ⟨w, hw⟩ := hex
    left
    have hAA : A' = A := by
      ext a
      exact ((const w a).mp hw).symm
    refine ⟨hAA, ?_⟩
    ext a
    rw [mem_snd_iff h' a, mem_snd_iff h a, hAA]
  · push_neg at hex
    right
    have hAB : A' = B := by
      ext a
      rw [mem_snd_iff h a]
      have := hex a
      tauto
    refine ⟨hAB, ?_⟩
    ext a
    rw [mem_snd_iff h' a, hAB, mem_snd_iff h a]
    tauto

/-- **Weighted bipartition formula for a connected component.**  A connected bipartite
graph with colour classes of sizes `a` and `b` has imbalance generating function
`z^(a-b) + z^(b-a)`; in particular a balanced component contributes `2 = z^0 + z^0`,
and an isolated vertex contributes `z + z⁻¹`. -/
theorem imbalanceGF_connected {H : SimpleGraph W} (hc : H.Connected) {A B : Finset W}
    (h : IsStablePair H A B) :
    imbalanceGF H = T ((A.card : ℤ) - B.card) + T ((B.card : ℤ) - A.card) := by
  classical
  obtain ⟨hA, hB, hd, hu⟩ := h
  have hst : IsStablePair H A B := ⟨hA, hB, hd, hu⟩
  have hst' : IsStablePair H B A := ⟨hB, hA, hd.symm, by rw [Finset.union_comm]; exact hu⟩
  have hne : A ≠ B := by
    intro heq
    obtain ⟨w⟩ := hc.nonempty
    have hw : w ∈ A ∨ w ∈ B := by
      have : w ∈ A ∪ B := by rw [hu]; exact Finset.mem_univ w
      exact Finset.mem_union.mp this
    rcases hw with hw | hw
    · exact (Finset.disjoint_left.mp hd hw) (heq ▸ hw)
    · exact (Finset.disjoint_left.mp hd (heq ▸ hw)) hw
  have hpairs : stablePairs H = {(A, B), (B, A)} := by
    ext p
    simp only [Finset.mem_insert, Finset.mem_singleton, mem_stablePairs]
    constructor
    · intro hp
      rcases stablePair_unique hc hst hp with ⟨h1, h2⟩ | ⟨h1, h2⟩
      · left; exact Prod.ext h1 h2
      · right; exact Prod.ext h1 h2
    · rintro (rfl | rfl)
      · exact hst
      · exact hst'
  rw [imbalanceGF, hpairs, Finset.sum_pair (by simp [Prod.ext_iff, hne])]

/-! ### Products over disjoint unions -/

private theorem isStablePair_toLeft {H₁ : SimpleGraph W₁} {H₂ : SimpleGraph W₂}
    {A B : Finset (W₁ ⊕ W₂)} (h : IsStablePair (H₁.sum H₂) A B) :
    IsStablePair H₁ A.toLeft B.toLeft := by
  obtain ⟨hA, hB, hd, hu⟩ := h
  refine ⟨?_, ?_, ?_, ?_⟩
  · intro x hx y hy hxy
    exact hA _ (Finset.mem_toLeft.mp hx) _ (Finset.mem_toLeft.mp hy) (by simpa using hxy)
  · intro x hx y hy hxy
    exact hB _ (Finset.mem_toLeft.mp hx) _ (Finset.mem_toLeft.mp hy) (by simpa using hxy)
  · rw [Finset.disjoint_left]
    intro a ha ha'
    exact Finset.disjoint_left.mp hd (Finset.mem_toLeft.mp ha) (Finset.mem_toLeft.mp ha')
  · ext a
    simp only [Finset.mem_union, Finset.mem_toLeft, Finset.mem_univ, iff_true]
    have : (Sum.inl a : W₁ ⊕ W₂) ∈ A ∪ B := by rw [hu]; exact Finset.mem_univ _
    exact Finset.mem_union.mp this

private theorem isStablePair_toRight {H₁ : SimpleGraph W₁} {H₂ : SimpleGraph W₂}
    {A B : Finset (W₁ ⊕ W₂)} (h : IsStablePair (H₁.sum H₂) A B) :
    IsStablePair H₂ A.toRight B.toRight := by
  obtain ⟨hA, hB, hd, hu⟩ := h
  refine ⟨?_, ?_, ?_, ?_⟩
  · intro x hx y hy hxy
    exact hA _ (Finset.mem_toRight.mp hx) _ (Finset.mem_toRight.mp hy) (by simpa using hxy)
  · intro x hx y hy hxy
    exact hB _ (Finset.mem_toRight.mp hx) _ (Finset.mem_toRight.mp hy) (by simpa using hxy)
  · rw [Finset.disjoint_left]
    intro a ha ha'
    exact Finset.disjoint_left.mp hd (Finset.mem_toRight.mp ha) (Finset.mem_toRight.mp ha')
  · ext a
    simp only [Finset.mem_union, Finset.mem_toRight, Finset.mem_univ, iff_true]
    have : (Sum.inr a : W₁ ⊕ W₂) ∈ A ∪ B := by rw [hu]; exact Finset.mem_univ _
    exact Finset.mem_union.mp this

private theorem isStablePair_disjSum {H₁ : SimpleGraph W₁} {H₂ : SimpleGraph W₂}
    {A₁ B₁ : Finset W₁} {A₂ B₂ : Finset W₂} (h1 : IsStablePair H₁ A₁ B₁)
    (h2 : IsStablePair H₂ A₂ B₂) :
    IsStablePair (H₁.sum H₂) (A₁.disjSum A₂) (B₁.disjSum B₂) := by
  obtain ⟨hA1, hB1, hd1, hu1⟩ := h1
  obtain ⟨hA2, hB2, hd2, hu2⟩ := h2
  refine ⟨?_, ?_, ?_, ?_⟩
  · rintro (x | x) hx (y | y) hy hxy <;>
      simp only [Finset.mem_disjSum] at hx hy <;> simp only [SimpleGraph.sum_adj] at hxy
    · exact hA1 x (by simpa using hx) y (by simpa using hy) hxy
    · exact hA2 x (by simpa using hx) y (by simpa using hy) hxy
  · rintro (x | x) hx (y | y) hy hxy <;>
      simp only [Finset.mem_disjSum] at hx hy <;> simp only [SimpleGraph.sum_adj] at hxy
    · exact hB1 x (by simpa using hx) y (by simpa using hy) hxy
    · exact hB2 x (by simpa using hx) y (by simpa using hy) hxy
  · rw [Finset.disjoint_left]
    rintro (a | a) ha ha' <;> simp only [Finset.mem_disjSum] at ha ha'
    · exact Finset.disjoint_left.mp hd1 (by simpa using ha) (by simpa using ha')
    · exact Finset.disjoint_left.mp hd2 (by simpa using ha) (by simpa using ha')
  · ext x
    simp only [Finset.mem_union, Finset.mem_univ, iff_true]
    rcases x with a | a
    · have hmem : a ∈ A₁ ∪ B₁ := by rw [hu1]; exact Finset.mem_univ _
      rcases Finset.mem_union.mp hmem with hb | hb
      · exact Or.inl (Finset.inl_mem_disjSum.mpr hb)
      · exact Or.inr (Finset.inl_mem_disjSum.mpr hb)
    · have hmem : a ∈ A₂ ∪ B₂ := by rw [hu2]; exact Finset.mem_univ _
      rcases Finset.mem_union.mp hmem with hb | hb
      · exact Or.inl (Finset.inr_mem_disjSum.mpr hb)
      · exact Or.inr (Finset.inr_mem_disjSum.mpr hb)

/-- **Multiplicativity.**  The imbalance generating function of a disjoint union is the
product of the imbalance generating functions.  This is what lets outside components be
factored out of a local cancellation block. -/
theorem imbalanceGF_sum (H₁ : SimpleGraph W₁) (H₂ : SimpleGraph W₂) :
    imbalanceGF (H₁.sum H₂) = imbalanceGF H₁ * imbalanceGF H₂ := by
  classical
  rw [imbalanceGF, imbalanceGF, imbalanceGF, Finset.sum_mul_sum, ← Finset.sum_product']
  refine Finset.sum_nbij' (i := fun p => ((p.1.toLeft, p.2.toLeft), (p.1.toRight, p.2.toRight)))
    (j := fun q => (q.1.1.disjSum q.2.1, q.1.2.disjSum q.2.2)) ?_ ?_ ?_ ?_ ?_
  · intro p hp
    rw [mem_stablePairs] at hp
    simp only [Finset.mem_product, mem_stablePairs]
    exact ⟨isStablePair_toLeft hp, isStablePair_toRight hp⟩
  · intro q hq
    simp only [Finset.mem_product, mem_stablePairs] at hq
    rw [mem_stablePairs]
    exact isStablePair_disjSum hq.1 hq.2
  · intro p _
    simp [Finset.toLeft_disjSum_toRight]
  · intro q _
    simp
  · intro p hp
    rw [← T_add]
    congr 1
    have h1 := Finset.card_toLeft_add_card_toRight (u := p.1)
    have h2 := Finset.card_toLeft_add_card_toRight (u := p.2)
    push_cast [← h1, ← h2]
    ring

end ClanAudit
