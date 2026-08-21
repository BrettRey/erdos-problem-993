import Mathlib

/-!
# A rank/parent criterion for being a tree

A connected finite graph is a tree as soon as it carries a *rank* function `f : V → ℕ`
minimized at a root `r` and a *parent* function `par` such that

* `par v` is adjacent to `v` and has strictly smaller rank, for every `v ≠ r`;
* adjacent vertices have different ranks;
* the parent is the *only* neighbour of smaller rank.

Indeed `v ↦ s(par v, v)` is then a bijection from the non-root vertices onto the edges, so
the graph has exactly `|V| - 1` edges.

This is the criterion used to check that the explicit spider and connector models of this
development really are trees.
-/

namespace ClanAudit

open SimpleGraph

variable {V : Type*} [Fintype V] [DecidableEq V] {G : SimpleGraph V}

/-- The edges of a graph with a rank and a parent function are in bijection with the
non-root vertices. -/
noncomputable def parentEdgeEquiv (r : V) (f : V → ℕ) (par : V → V)
    (hlt : ∀ v, v ≠ r → f (par v) < f v)
    (hadj : ∀ v, v ≠ r → G.Adj (par v) v)
    (hne : ∀ u v, G.Adj u v → f u ≠ f v)
    (hroot : ∀ v, f r ≤ f v)
    (huniq : ∀ u v, G.Adj u v → f u < f v → u = par v) :
    {v : V // v ≠ r} ≃ G.edgeSet := by
  refine Equiv.ofBijective (fun v => ⟨s(par v.1, v.1), (hadj v.1 v.2)⟩) ⟨?_, ?_⟩
  · rintro ⟨u, hu⟩ ⟨v, hv⟩ h
    simp only [Subtype.mk.injEq, Sym2.eq, Sym2.rel_iff', Prod.mk.injEq, Prod.swap_prod_mk] at h
    rcases h with ⟨_, h2⟩ | ⟨h1, h2⟩
    · exact Subtype.ext h2
    · exfalso
      have e1 := hlt u hu
      have e2 := hlt v hv
      rw [h1] at e1
      rw [← h2] at e2
      omega
  · rintro ⟨e, he⟩
    induction e with
    | _ u v =>
      have hadj' : G.Adj u v := he
      rcases lt_or_gt_of_ne (hne u v hadj') with hlt' | hlt'
      · have hvr : v ≠ r := by
          intro hc
          have := hroot u
          rw [hc] at hlt'
          omega
        refine ⟨⟨v, hvr⟩, ?_⟩
        rw [Subtype.mk.injEq, huniq u v hadj' hlt']
      · have hur : u ≠ r := by
          intro hc
          have := hroot v
          rw [hc] at hlt'
          omega
        refine ⟨⟨u, hur⟩, ?_⟩
        rw [Subtype.mk.injEq, huniq v u hadj'.symm hlt', Sym2.eq_swap]

/-- **Rank/parent criterion.**  A connected finite graph with a rank function and a
parent function as above is a tree. -/
theorem isTree_of_rank_parent (hc : G.Connected) (r : V) (f : V → ℕ) (par : V → V)
    (hlt : ∀ v, v ≠ r → f (par v) < f v)
    (hadj : ∀ v, v ≠ r → G.Adj (par v) v)
    (hne : ∀ u v, G.Adj u v → f u ≠ f v)
    (hroot : ∀ v, f r ≤ f v)
    (huniq : ∀ u v, G.Adj u v → f u < f v → u = par v) :
    G.IsTree := by
  classical
  rw [isTree_iff_connected_and_card]
  refine ⟨hc, ?_⟩
  have hbij := parentEdgeEquiv (G := G) r f par hlt hadj hne hroot huniq
  have hcard : Nat.card G.edgeSet = Nat.card {v : V // v ≠ r} := (Nat.card_eq_of_bijective _
    (Equiv.bijective hbij.symm))
  have hsub : Fintype.card {v : V // v ≠ r} = Fintype.card V - 1 := by
    simp [Fintype.card_subtype_compl (p := fun v : V => v = r)]
  have hpos : 1 ≤ Fintype.card V := Fintype.card_pos_iff.mpr ⟨r⟩
  rw [hcard, Nat.card_eq_fintype_card, Nat.card_eq_fintype_card, hsub]
  omega

end ClanAudit
