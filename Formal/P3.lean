import Mathlib

/-!
# P3: Support-Root Domination

P3 (tail domination) states that at a support vertex r of a graph with pendant (ℓ, r),
the number of size-k independent sets excluding r is at least the number including r.
Equivalently, e_k ≥ j_{k-1} for all k, where E = dp[r][0] and J = dp[r][1]/x.

This file provides:
1. **Algebraic P3**: if A ≥ B coefficientwise with nonneg coefficients,
   then (1+x)A - xB ≥ 0 coefficientwise.
2. **Combinatorial P3**: the leaf-swap injection φ(S) = (S \ {r}) ∪ {ℓ}
   is a size-preserving injection from {IS ∋ r} to {IS ∌ r}.
3. **Tree corollary**: every nontrivial tree has a pendant pair.
-/

namespace Formal.P3

/-! ## Section 1: Algebraic P3 -/

section Algebraic

/-- If `a ≥ b` coefficientwise (over ℕ), then `b k ≤ a k + a (k+1)`.
    This captures: the k-th coefficient of `(1+x)A - xB` is nonneg. -/
theorem p3_algebraic (a b : ℕ → ℕ)
    (hab : ∀ k, b k ≤ a k)
    (k : ℕ) : b k ≤ a k + a (k + 1) :=
  Nat.le_add_right_of_le (hab k)

/-- Wrapper with DP names: E_k = a_k + a_{k-1}, (xJ)_k = b_{k-1}.
    Given a ≥ b coefficientwise, E_{k+1} ≥ (xJ)_{k+1} = b_k. -/
theorem p3_support_vertex (a b : ℕ → ℕ)
    (hab : ∀ k, b k ≤ a k)
    (k : ℕ) : b k ≤ a (k + 1) + a k := by
  have := hab k; omega

end Algebraic

/-! ## Section 2: Combinatorial P3 -/

section Combinatorial

set_option linter.unusedSectionVars false
set_option linter.unusedFintypeInType false
set_option linter.unusedDecidableInType false

variable {V : Type*} [Fintype V] [DecidableEq V]
variable (G : SimpleGraph V) [DecidableRel G.Adj]

/-- A set S of vertices is independent if no two distinct members are adjacent. -/
def IsIndepSet (S : Finset V) : Prop :=
  ∀ u ∈ S, ∀ v ∈ S, u ≠ v → ¬G.Adj u v

instance IsIndepSet.decidable (S : Finset V) : Decidable (IsIndepSet G S) :=
  inferInstanceAs (Decidable (∀ u ∈ S, ∀ v ∈ S, u ≠ v → ¬G.Adj u v))

/-- Vertex ℓ is a pendant (leaf) with unique neighbor r. -/
def IsPendant (ℓ r : V) : Prop :=
  G.Adj ℓ r ∧ ∀ w, G.Adj ℓ w → w = r

/-- The leaf-swap map: replace r with ℓ in S. -/
def leafSwap (r ℓ : V) (S : Finset V) : Finset V :=
  insert ℓ (S.erase r)

/-- Membership in leafSwap: `x ∈ leafSwap r ℓ S ↔ x = ℓ ∨ (x ≠ r ∧ x ∈ S)`. -/
theorem mem_leafSwap {r ℓ x : V} {S : Finset V} :
    x ∈ leafSwap r ℓ S ↔ x = ℓ ∨ (x ≠ r ∧ x ∈ S) := by
  simp [leafSwap, Finset.mem_insert, Finset.mem_erase]

/-- ℓ is not in an independent set containing r, when (ℓ, r) is pendant. -/
theorem leaf_not_mem_of_root_mem {ℓ r : V} {S : Finset V}
    (hpend : IsPendant G ℓ r) (hind : IsIndepSet G S) (hr : r ∈ S) :
    ℓ ∉ S := by
  intro hℓ
  have hne : ℓ ≠ r := fun h => by subst h; exact G.irrefl hpend.1
  exact hind ℓ hℓ r hr hne hpend.1

/-- Claim 1: leafSwap preserves independence when (ℓ, r) is a pendant pair. -/
theorem leafSwap_isIndep {ℓ r : V} {S : Finset V}
    (hpend : IsPendant G ℓ r)
    (hind : IsIndepSet G S) :
    IsIndepSet G (leafSwap r ℓ S) := by
  intro u hu v hv huv
  rw [mem_leafSwap] at hu hv
  rcases hu with rfl | ⟨hu_ne, hu_mem⟩ <;> rcases hv with rfl | ⟨hv_ne, hv_mem⟩
  · exact absurd rfl huv
  · intro hadj
    exact hv_ne (hpend.2 v hadj)
  · intro hadj
    exact hu_ne (hpend.2 u (G.symm hadj))
  · exact hind u hu_mem v hv_mem huv

/-- Claim 2: leafSwap preserves cardinality when r ∈ S and ℓ ∉ S. -/
theorem leafSwap_card {r ℓ : V} {S : Finset V}
    (hr : r ∈ S) (hℓ : ℓ ∉ S) :
    (leafSwap r ℓ S).card = S.card := by
  have hℓ_erase : ℓ ∉ S.erase r := by
    rw [Finset.mem_erase]; exact fun ⟨_, hm⟩ => hℓ hm
  simp only [leafSwap]
  rw [Finset.card_insert_of_notMem hℓ_erase, Finset.card_erase_of_mem hr]
  exact Nat.succ_pred_eq_of_pos (Finset.card_pos.mpr ⟨r, hr⟩)

/-- Claim 3: r is not in the image of leafSwap (when r ≠ ℓ). -/
theorem leafSwap_not_mem_root {r ℓ : V} {S : Finset V}
    (hne : r ≠ ℓ) : r ∉ leafSwap r ℓ S := by
  rw [mem_leafSwap]
  rintro (rfl | ⟨h, _⟩)
  · exact hne rfl
  · exact h rfl

/-- Claim 4: leafSwap is injective (on sets containing r and not ℓ). -/
theorem leafSwap_injective {r ℓ : V} (hne : r ≠ ℓ)
    {S₁ S₂ : Finset V} (hr₁ : r ∈ S₁) (hℓ₁ : ℓ ∉ S₁) (hr₂ : r ∈ S₂) (hℓ₂ : ℓ ∉ S₂)
    (heq : leafSwap r ℓ S₁ = leafSwap r ℓ S₂) : S₁ = S₂ := by
  ext x
  by_cases hxr : x = r
  · subst hxr; simp [hr₁, hr₂]
  · by_cases hxℓ : x = ℓ
    · subst hxℓ; simp [hℓ₁, hℓ₂]
    · have h1 : x ∈ leafSwap r ℓ S₁ ↔ x ∈ leafSwap r ℓ S₂ := heq ▸ Iff.rfl
      simp only [mem_leafSwap] at h1
      constructor <;> intro hx
      · exact ((h1.mp (Or.inr ⟨hxr, hx⟩)).resolve_left hxℓ).2
      · exact ((h1.mpr (Or.inr ⟨hxr, hx⟩)).resolve_left hxℓ).2

/-- **Main theorem (combinatorial P3):** For a pendant pair (ℓ, r), the number of
    size-k independent sets containing r is at most the number excluding r.

    This is: |{S : IS(G) | r ∈ S, |S| = k}| ≤ |{S : IS(G) | r ∉ S, |S| = k}|. -/
theorem support_root_domination {ℓ r : V} (hpend : IsPendant G ℓ r) (k : ℕ) :
    (Finset.univ.filter (fun S : Finset V =>
      IsIndepSet G S ∧ r ∈ S ∧ S.card = k)).card ≤
    (Finset.univ.filter (fun S : Finset V =>
      IsIndepSet G S ∧ r ∉ S ∧ S.card = k)).card := by
  have hne : r ≠ ℓ := fun h => by subst h; exact G.irrefl hpend.1
  -- The image of the source under leafSwap lands in the target
  have himage :
      (Finset.univ.filter (fun S : Finset V =>
        IsIndepSet G S ∧ r ∈ S ∧ S.card = k)).image (leafSwap r ℓ) ⊆
      Finset.univ.filter (fun S : Finset V =>
        IsIndepSet G S ∧ r ∉ S ∧ S.card = k) := by
    intro T hT
    rw [Finset.mem_image] at hT
    obtain ⟨S, hS, rfl⟩ := hT
    rw [Finset.mem_filter] at hS ⊢
    obtain ⟨_, hind, hr_mem, hcard⟩ := hS
    exact ⟨Finset.mem_univ _,
           leafSwap_isIndep G hpend hind,
           leafSwap_not_mem_root hne,
           (leafSwap_card hr_mem (leaf_not_mem_of_root_mem G hpend hind hr_mem)).trans hcard⟩
  -- leafSwap is injective on the source
  have hinj : Set.InjOn (leafSwap r ℓ)
      ((Finset.univ.filter (fun S : Finset V =>
        IsIndepSet G S ∧ r ∈ S ∧ S.card = k)) : Set (Finset V)) := by
    intro S₁ hS₁ S₂ hS₂ heq
    rw [Finset.mem_coe, Finset.mem_filter] at hS₁ hS₂
    exact leafSwap_injective hne
      hS₁.2.2.1 (leaf_not_mem_of_root_mem G hpend hS₁.2.1 hS₁.2.2.1)
      hS₂.2.2.1 (leaf_not_mem_of_root_mem G hpend hS₂.2.1 hS₂.2.2.1)
      heq
  calc (Finset.univ.filter (fun S : Finset V =>
          IsIndepSet G S ∧ r ∈ S ∧ S.card = k)).card
      = ((Finset.univ.filter (fun S : Finset V =>
          IsIndepSet G S ∧ r ∈ S ∧ S.card = k)).image (leafSwap r ℓ)).card :=
        (Finset.card_image_of_injOn hinj).symm
    _ ≤ (Finset.univ.filter (fun S : Finset V =>
          IsIndepSet G S ∧ r ∉ S ∧ S.card = k)).card :=
        Finset.card_le_card himage

end Combinatorial

/-! ## Section 3: Tree corollary -/

section TreeCorollary

set_option linter.unusedSectionVars false
set_option linter.unusedDecidableInType false

variable {V : Type*} [Fintype V] [DecidableEq V]
variable (G : SimpleGraph V) [DecidableRel G.Adj]

/-- Every tree on ≥ 2 vertices has a pendant pair.
    This is a standard fact but requires nontrivial Lean graph theory
    infrastructure (finite connected acyclic ⟹ has leaf). -/
theorem tree_has_pendant (hT : G.IsTree) (hcard : 1 < Fintype.card V) :
    ∃ ℓ r : V, IsPendant G ℓ r := by
  obtain ⟨ℓ, hℓ⟩ : ∃ ℓ : V, G.degree ℓ = 1 := by
    apply_rules [SimpleGraph.IsTree.exists_vert_degree_one_of_nontrivial]
    exact Fintype.one_lt_card_iff_nontrivial.mp hcard
  obtain ⟨r, hℓr, hr_unique⟩ := SimpleGraph.degree_eq_one_iff_existsUnique_adj.mp hℓ
  exact ⟨ℓ, r, hℓr, hr_unique⟩

/-- P3 holds at some rooting of every nontrivial tree. -/
theorem tree_p3 (hT : G.IsTree) (hcard : 1 < Fintype.card V) (k : ℕ) :
    ∃ r : V, (Finset.univ.filter (fun S : Finset V =>
      IsIndepSet G S ∧ r ∈ S ∧ S.card = k)).card ≤
    (Finset.univ.filter (fun S : Finset V =>
      IsIndepSet G S ∧ r ∉ S ∧ S.card = k)).card := by
  obtain ⟨ℓ, r, hpend⟩ := tree_has_pendant G hT hcard
  exact ⟨r, support_root_domination G hpend k⟩

end TreeCorollary

end Formal.P3
