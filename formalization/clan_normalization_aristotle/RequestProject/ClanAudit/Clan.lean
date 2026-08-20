import Mathlib

/-!
# Clan graphs and the clan normalization identity

The clan graph `clan G alpha` replaces each vertex `v` of `G` by a clique of size
`alpha v` and joins all clones of `u` to all clones of `v` whenever `uv ∈ E(G)`.

The *clan normalization* of Li–Li–Yang–Zhang says that proper colourings of the clan
graph are in `∏ v, (alpha v)!`-to-one correspondence with proper multicolourings of
type `alpha`, with equal monomials.  We prove exactly that statement here, in the
finite-counting form the request allows (no symmetric-function library is used):

* `toMulti_isMulticoloring` : the colour sets of a proper colouring of the clan graph
  form a proper multicolouring of type `alpha`;
* `clan_fiber_card` : each proper multicolouring has exactly `∏ v, (alpha v)!`
  preimages;
* `monomial_toMulti` : a colouring and its multicolouring have the same monomial.
-/

namespace ClanAudit

open Finset

variable {V : Type*} [Fintype V] [DecidableEq V]

/-- The clan graph of `G` with multiplicities `alpha`: each vertex `v` is blown up into a
clique of `alpha v` clones, and clones of adjacent vertices are joined. -/
def clan (G : SimpleGraph V) (alpha : V → ℕ) : SimpleGraph (Σ v : V, Fin (alpha v)) where
  Adj x y := (x.1 = y.1 ∧ x ≠ y) ∨ G.Adj x.1 y.1
  symm := by
    rintro x y (⟨h1, h2⟩ | h)
    · exact Or.inl ⟨h1.symm, h2.symm⟩
    · exact Or.inr h.symm
  loopless := ⟨by
    rintro x (⟨-, h⟩ | h)
    · exact h rfl
    · exact G.irrefl h⟩

/-- A proper colouring of a graph by natural-number colours. -/
def IsProperColoring {U : Type*} (H : SimpleGraph U) (f : U → ℕ) : Prop :=
  ∀ x y, H.Adj x y → f x ≠ f y

/-- A proper multicolouring of type `alpha`: each vertex `v` receives a set of exactly
`alpha v` colours, and adjacent vertices receive disjoint colour sets. -/
structure IsMulticoloring (G : SimpleGraph V) (alpha : V → ℕ) (m : V → Finset ℕ) : Prop where
  card : ∀ v, (m v).card = alpha v
  disj : ∀ u v, G.Adj u v → Disjoint (m u) (m v)

/-- The colour sets induced by a colouring of the clan graph. -/
def toMulti (alpha : V → ℕ) (f : (Σ v : V, Fin (alpha v)) → ℕ) : V → Finset ℕ :=
  fun v => Finset.image (fun i => f ⟨v, i⟩) Finset.univ

/-- The monomial of a colouring: how many clan vertices carry the colour `c`. -/
noncomputable def monomialOfColoring {alpha : V → ℕ} (f : (Σ v : V, Fin (alpha v)) → ℕ)
    (c : ℕ) : ℕ := by
  classical exact (Finset.univ.filter (fun x => f x = c)).card

/-- The monomial of a multicolouring: how many original vertices use the colour `c`. -/
noncomputable def monomialOfMulti (m : V → Finset ℕ) (c : ℕ) : ℕ := by
  classical exact (Finset.univ.filter (fun v => c ∈ m v)).card

omit [Fintype V] [DecidableEq V] in
/-- Within one clan clique a proper colouring is injective. -/
theorem clone_injective {G : SimpleGraph V} {alpha : V → ℕ}
    {f : (Σ v : V, Fin (alpha v)) → ℕ} (hf : IsProperColoring (clan G alpha) f) (v : V) :
    Function.Injective (fun i : Fin (alpha v) => f ⟨v, i⟩) := by
  intro i j hij
  by_contra hne
  refine hf ⟨v, i⟩ ⟨v, j⟩ (Or.inl ⟨rfl, ?_⟩) hij
  simp only [ne_eq, Sigma.mk.injEq, heq_eq_eq, true_and]
  exact hne

omit [Fintype V] [DecidableEq V] in
/-- Proper colourings of the clan graph induce proper multicolourings of type `alpha`. -/
theorem toMulti_isMulticoloring {G : SimpleGraph V} {alpha : V → ℕ}
    {f : (Σ v : V, Fin (alpha v)) → ℕ} (hf : IsProperColoring (clan G alpha) f) :
    IsMulticoloring G alpha (toMulti alpha f) where
  card v := by
    rw [toMulti, Finset.card_image_of_injective _ (clone_injective hf v),
      Finset.card_univ, Fintype.card_fin]
  disj u v huv := by
    rw [Finset.disjoint_left]
    rintro c hcu hcv
    rw [toMulti, Finset.mem_image] at hcu hcv
    obtain ⟨i, -, hi⟩ := hcu
    obtain ⟨j, -, hj⟩ := hcv
    exact hf ⟨u, i⟩ ⟨v, j⟩ (Or.inr huv) (by rw [hi, hj])

omit [Fintype V] [DecidableEq V] in
/-- Membership in the fibre over a multicolouring is a purely local condition. -/
theorem fiber_iff {G : SimpleGraph V} {alpha : V → ℕ} {m : V → Finset ℕ}
    (hm : IsMulticoloring G alpha m) (f : (Σ v : V, Fin (alpha v)) → ℕ) :
    (IsProperColoring (clan G alpha) f ∧ toMulti alpha f = m)
      ↔ ∀ v, Finset.image (fun i : Fin (alpha v) => f ⟨v, i⟩) Finset.univ = m v := by
  constructor
  · rintro ⟨-, h2⟩ v
    exact congrFun h2 v
  · intro h
    have hinj : ∀ v, Function.Injective (fun i : Fin (alpha v) => f ⟨v, i⟩) := by
      intro v
      have hcard : (Finset.image (fun i : Fin (alpha v) => f ⟨v, i⟩) Finset.univ).card
          = (Finset.univ : Finset (Fin (alpha v))).card := by
        rw [h v, hm.card v, Finset.card_univ, Fintype.card_fin]
      have := Finset.card_image_iff.mp hcard
      intro i j hij
      exact this (Finset.mem_coe.mpr (Finset.mem_univ i)) (Finset.mem_coe.mpr (Finset.mem_univ j))
        hij
    refine ⟨?_, funext h⟩
    rintro ⟨u, i⟩ ⟨v, j⟩ (⟨huv, hne⟩ | hadj)
    · subst huv
      intro hcol
      apply hne
      have := hinj u hcol
      simp only [Sigma.mk.injEq, heq_eq_eq, true_and]
      exact this
    · intro hcol
      have hu : f ⟨u, i⟩ ∈ m u := by rw [← h u]; exact Finset.mem_image_of_mem _ (mem_univ i)
      have hv : f ⟨v, j⟩ ∈ m v := by rw [← h v]; exact Finset.mem_image_of_mem _ (mem_univ j)
      rw [hcol] at hu
      exact Finset.disjoint_left.mp (hm.disj u v hadj) hu hv

/-- Functions with a prescribed image on `Fin k` biject with embeddings into that image. -/
noncomputable def imageFiberEquiv {k : ℕ} {S : Finset ℕ} (hS : S.card = k) :
    {g : Fin k → ℕ // Finset.image g Finset.univ = S} ≃ (Fin k ↪ {x // x ∈ S}) where
  toFun g :=
    ⟨fun i => ⟨g.1 i, by
      have hmem := Finset.mem_image_of_mem g.1 (mem_univ i)
      rwa [g.2] at hmem⟩, by
      have hcard : (Finset.image g.1 Finset.univ).card
          = (Finset.univ : Finset (Fin k)).card := by
        rw [g.2, hS, Finset.card_univ, Fintype.card_fin]
      have hinj := Finset.card_image_iff.mp hcard
      intro i j hij
      exact hinj (Finset.mem_coe.mpr (Finset.mem_univ i)) (Finset.mem_coe.mpr (Finset.mem_univ j))
        (by simpa using hij)⟩
  invFun e :=
    ⟨fun i => (e i).1, by
      have hsub : Finset.image (fun i => (e i).1) Finset.univ ⊆ S := by
        intro c hc
        rw [Finset.mem_image] at hc
        obtain ⟨i, -, rfl⟩ := hc
        exact (e i).2
      have hinj : Function.Injective (fun i => ((e i).1 : ℕ)) := by
        intro i j hij
        exact e.injective (Subtype.ext hij)
      have hcard : (Finset.image (fun i => (e i).1) Finset.univ).card = S.card := by
        rw [Finset.card_image_of_injective _ hinj, Finset.card_univ, Fintype.card_fin, hS]
      exact Finset.eq_of_subset_of_card_le hsub (le_of_eq hcard.symm)⟩
  left_inv g := by ext i; rfl
  right_inv e := by ext i; rfl

theorem card_image_fiber {k : ℕ} {S : Finset ℕ} (hS : S.card = k) :
    Nat.card {g : Fin k → ℕ // Finset.image g Finset.univ = S} = Nat.factorial k := by
  classical
  rw [Nat.card_congr (imageFiberEquiv hS), Nat.card_eq_fintype_card, Fintype.card_embedding_eq,
    Fintype.card_coe, Fintype.card_fin, hS, Nat.descFactorial_self]

/-- Currying the colourings of a clan graph. -/
def clanCurryEquiv {alpha : V → ℕ} (Q : ∀ v, (Fin (alpha v) → ℕ) → Prop) :
    {f : (Σ v : V, Fin (alpha v)) → ℕ // ∀ v, Q v (fun i => f ⟨v, i⟩)}
      ≃ ∀ v : V, {g : Fin (alpha v) → ℕ // Q v g} where
  toFun f v := ⟨fun i => f.1 ⟨v, i⟩, f.2 v⟩
  invFun F := ⟨fun x => (F x.1).1 x.2, fun v => (F v).2⟩
  left_inv f := by
    apply Subtype.ext
    funext x
    cases x
    rfl
  right_inv F := by
    funext v
    rfl

omit [DecidableEq V] in
/-- **Clan normalization (fibre count).**  Every proper multicolouring of type `alpha`
is the colour-set family of exactly `∏ v, (alpha v)!` proper colourings of the clan
graph. -/
theorem clan_fiber_card {G : SimpleGraph V} {alpha : V → ℕ} {m : V → Finset ℕ}
    (hm : IsMulticoloring G alpha m) :
    Nat.card {f : (Σ v : V, Fin (alpha v)) → ℕ //
        IsProperColoring (clan G alpha) f ∧ toMulti alpha f = m}
      = ∏ v : V, Nat.factorial (alpha v) := by
  classical
  have hequiv : {f : (Σ v : V, Fin (alpha v)) → ℕ //
      IsProperColoring (clan G alpha) f ∧ toMulti alpha f = m}
      ≃ {f : (Σ v : V, Fin (alpha v)) → ℕ //
        ∀ v, Finset.image (fun i : Fin (alpha v) => f ⟨v, i⟩) Finset.univ = m v} :=
    Equiv.subtypeEquivRight (fiber_iff hm)
  rw [Nat.card_congr hequiv,
    Nat.card_congr (clanCurryEquiv (alpha := alpha)
      (fun v g => Finset.image g Finset.univ = m v)),
    Nat.card_pi]
  exact Finset.prod_congr rfl fun v _ => card_image_fiber (hm.card v)

omit [DecidableEq V] in
/-- **Clan normalization (monomials agree).**  A proper colouring of the clan graph and
its multicolouring have the same monomial. -/
theorem monomial_toMulti {G : SimpleGraph V} {alpha : V → ℕ}
    {f : (Σ v : V, Fin (alpha v)) → ℕ} (hf : IsProperColoring (clan G alpha) f) (c : ℕ) :
    monomialOfColoring f c = monomialOfMulti (toMulti alpha f) c := by
  classical
  rw [monomialOfColoring, monomialOfMulti]
  refine Finset.card_bij (fun x _ => x.1) ?_ ?_ ?_
  · intro x hx
    simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hx ⊢
    rw [← hx, toMulti]
    exact Finset.mem_image_of_mem _ (mem_univ x.2)
  · intro x hx y hy hxy
    simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hx hy
    have hcol : f x = f y := by rw [hx, hy]
    obtain ⟨u, i⟩ := x
    obtain ⟨v, j⟩ := y
    simp only at hxy
    subst hxy
    have hij : i = j := clone_injective hf u hcol
    subst hij
    rfl
  · intro v hv
    simp only [Finset.mem_filter, Finset.mem_univ, true_and, toMulti, Finset.mem_image] at hv
    obtain ⟨i, hi⟩ := hv
    exact ⟨⟨v, i⟩, Finset.mem_filter.mpr ⟨Finset.mem_univ _, hi⟩, rfl⟩

end ClanAudit
