import Mathlib.Combinatorics.SimpleGraph.Clique
import Mathlib.Combinatorics.SimpleGraph.Acyclic

/-! Frozen specification. Definitions below state propositions, not proofs. -/
namespace DepthThree

variable {V : Type*} [Fintype V]

def Extendable (G : SimpleGraph V) (S : Finset V) : Prop :=
  ∃ M : Finset V, G.IsMaximumIndepSet M ∧ S ⊆ M

noncomputable def e (G : SimpleGraph V) (d : ℕ) : ℕ := by
  classical
  exact if d ≤ G.indepNum then
    (Finset.univ.filter fun S : Finset V =>
      S.card = G.indepNum - d ∧ G.IsIndepSet (S : Set V) ∧ Extendable G S).card
  else 0

noncomputable def b (G : SimpleGraph V) (d : ℕ) : ℕ := by
  classical
  exact if d ≤ G.indepNum then
    (Finset.univ.filter fun S : Finset V =>
      S.card = G.indepNum - d ∧ G.IsIndepSet (S : Set V) ∧ ¬ Extendable G S).card
  else 0

noncomputable def s (G : SimpleGraph V) (d : ℕ) : ℕ := e G d + b G d

def DepthThreeStrict (G : SimpleGraph V) : Prop := s G 2 * s G 4 < s G 3 ^ 2

def BlockedShadowTarget : Prop :=
  ∀ {V : Type*} [Fintype V] (G : SimpleGraph V), G.IsAcyclic →
    4 ≤ G.indepNum → (G.indepNum - 4) * b G 2 ≤ 6 * b G 3

def LowDensityTarget : Prop :=
  ∀ {V : Type*} [Fintype V] (G : SimpleGraph V), G.IsAcyclic →
    17 ≤ G.indepNum → G.indepNum ≤ 19 →
    3 * (G.indepNum - 3) * b G 4 ≤ (G.indepNum - 7) * e G 4 →
    DepthThreeStrict G

def B1ZeroWindowTarget : Prop :=
  ∀ {V : Type*} [Fintype V] (G : SimpleGraph V), G.IsTree →
    33 ≤ Fintype.card V → Fintype.card V ≤ 38 →
    17 ≤ G.indepNum → G.indepNum ≤ 19 →
    2 * G.indepNum ≤ Fintype.card V + 5 → b G 1 = 0 →
    DepthThreeStrict G

#check BlockedShadowTarget
#check LowDensityTarget
#check B1ZeroWindowTarget

end DepthThree
