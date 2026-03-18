/-
  J ≤ E coefficientwise: At any vertex v of a graph G,
  the number of size-k independent sets avoiding the closed neighborhood N[v]
  is at most the number of size-k independent sets avoiding just v.

  This is because V \ N[v] ⊆ V \ {v}, so every IS in G[V \ N[v]]
  is also an IS in G[V \ {v}].
-/
import Formal.P3

open SimpleGraph Classical Finset

noncomputable section

set_option linter.unusedSectionVars false
set_option linter.unusedDecidableInType false

variable {V : Type*} [Fintype V] [DecidableEq V]
variable (G : SimpleGraph V) [DecidableRel G.Adj]

namespace Formal.JleE

/-- Number of independent sets of size k whose vertices all lie in a given set W. -/
def numIndSetsIn (W : Finset V) (k : ℕ) : ℕ :=
  (Finset.univ.powerset.filter
    (fun S : Finset V => S.card = k ∧ S ⊆ W ∧ Formal.P3.IsIndepSet G S)).card

/-- The closed neighborhood N[v] = {v} ∪ N(v). -/
def closedNeighborFinset (v : V) : Finset V :=
  insert v (G.neighborFinset v)

/-- J_k(v) = number of size-k IS in V \ N[v]. -/
def numJ (v : V) (k : ℕ) : ℕ :=
  numIndSetsIn G (Finset.univ \ closedNeighborFinset G v) k

/-- E_k(v) = number of size-k IS in V \ {v}. -/
def numE (v : V) (k : ℕ) : ℕ :=
  numIndSetsIn G (Finset.univ \ {v}) k

/-- Subgraph monotonicity: if W₁ ⊆ W₂, IS count in W₁ ≤ IS count in W₂. -/
theorem numIndSetsIn_mono
    (W₁ W₂ : Finset V) (k : ℕ) (h : W₁ ⊆ W₂) :
    numIndSetsIn G W₁ k ≤ numIndSetsIn G W₂ k := by
  convert Finset.card_mono ?_ using 2 ; aesop_cat;

/-- V \ N[v] ⊆ V \ {v}, because v ∈ N[v]. -/
theorem sdiff_closedNbhd_subset_sdiff_singleton (v : V) :
    Finset.univ \ closedNeighborFinset G v ⊆ Finset.univ \ {v} := by
  exact Finset.sdiff_subset_sdiff ( Finset.Subset.refl _ ) ( Finset.singleton_subset_iff.2 ( Finset.mem_insert_self _ _ ) )

/-- J ≤ E coefficientwise: numJ G v k ≤ numE G v k. -/
theorem J_le_E (v : V) (k : ℕ) :
    numJ G v k ≤ numE G v k := by
  apply numIndSetsIn_mono; exact sdiff_closedNbhd_subset_sdiff_singleton G v;

end Formal.JleE
