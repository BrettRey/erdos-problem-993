import Mathlib.Combinatorics.SimpleGraph.Clique
import Mathlib.Combinatorics.SimpleGraph.Acyclic

namespace FivefoldForest

variable {V : Type*} [Fintype V]

/-- Containment in a maximum independent set, not merely a maximal one. -/
def Extendable (G : SimpleGraph V) (S : Finset V) : Prop :=
  ∃ M : Finset V, G.IsMaximumIndepSet M ∧ S ⊆ M

noncomputable def extendableCount (G : SimpleGraph V) (d : ℕ) : ℕ := by
  classical
  exact if d ≤ G.indepNum then
    (Finset.univ.filter fun S : Finset V =>
      S.card = G.indepNum - d ∧ G.IsIndepSet (S : Set V) ∧ Extendable G S).card
  else 0

noncomputable def blockedCount (G : SimpleGraph V) (d : ℕ) : ℕ := by
  classical
  exact if d ≤ G.indepNum then
    (Finset.univ.filter fun S : Finset V =>
      S.card = G.indepNum - d ∧ G.IsIndepSet (S : Set V) ∧ ¬ Extendable G S).card
  else 0

/-- Specification only; this definition is not a proof of the assertion. -/
def FivefoldConclusion (G : SimpleGraph V) : Prop :=
  5 * blockedCount G 3 ≤ extendableCount G 3 ∧
  5 * blockedCount G 4 ≤ extendableCount G 4

#check @FivefoldForest.extendableCount
#check @FivefoldForest.blockedCount
#check @FivefoldForest.FivefoldConclusion

end FivefoldForest
