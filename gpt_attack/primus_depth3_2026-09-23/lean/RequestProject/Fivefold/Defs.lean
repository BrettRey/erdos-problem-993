import Mathlib.Combinatorics.SimpleGraph.Clique
import Mathlib.Combinatorics.SimpleGraph.Acyclic
import Mathlib.Tactic

/-!
# Fivefold bounds for blocked sets in finite forests: core definitions

The whole development works with a *fixed* ambient graph `G : SimpleGraph V`
and studies its induced subgraphs, which are encoded simply as finsets
`S : Finset V` of vertices.  This makes vertex deletion, disjoint unions and
induction on the number of vertices painless, and every induced subgraph of a
forest is automatically a forest.

At the very end (`RequestProject.Main`) the counts defined here, taken over
`S = Finset.univ`, are identified with the counts from the problem statement,
which are phrased with `SimpleGraph.IsMaximumIndepSet` and
`SimpleGraph.indepNum`.
-/

open Finset

namespace FivefoldForest

open scoped Classical

variable {V : Type*} (G : SimpleGraph V)

/-- The independent subsets of the vertex set `S`. -/
noncomputable def indFam (S : Finset V) : Finset (Finset V) :=
  S.powerset.filter fun A => G.IsIndepSet (A : Set V)

/-- The independence number of the subgraph induced on `S`. -/
noncomputable def alpha (S : Finset V) : ℕ := (indFam G S).sup Finset.card

/-- `A` is contained in a maximum independent set of the subgraph induced on `S`. -/
def Ext (S A : Finset V) : Prop :=
  ∃ M ∈ indFam G S, M.card = alpha G S ∧ A ⊆ M

/-- Number of independent subsets of `S` of size `k`. -/
noncomputable def icnt (S : Finset V) (k : ℕ) : ℕ :=
  ((indFam G S).filter fun A => A.card = k).card

/-- Number of *extendable* independent subsets of `S` of size `k`. -/
noncomputable def ecnt (S : Finset V) (k : ℕ) : ℕ :=
  ((indFam G S).filter fun A => A.card = k ∧ Ext G S A).card

/-- Number of *blocked* independent subsets of `S` of size `k`. -/
noncomputable def bcnt (S : Finset V) (k : ℕ) : ℕ :=
  ((indFam G S).filter fun A => A.card = k ∧ ¬ Ext G S A).card

/-- The extendable count at defect `d`: independent sets of size `α - d`
contained in a maximum independent set.  Defects larger than `α` give `0`. -/
noncomputable def eD (S : Finset V) (d : ℕ) : ℕ :=
  if d ≤ alpha G S then ecnt G S (alpha G S - d) else 0

/-- The blocked count at defect `d`: independent sets of size `α - d`
*not* contained in a maximum independent set.  Defects larger than `α` give `0`. -/
noncomputable def bD (S : Finset V) (d : ℕ) : ℕ :=
  if d ≤ alpha G S then bcnt G S (alpha G S - d) else 0

/-- The margin `e_d - 5 b_d`, as an integer. -/
noncomputable def margin (S : Finset V) (d : ℕ) : ℤ :=
  (eD G S d : ℤ) - 5 * (bD G S d : ℤ)

end FivefoldForest
