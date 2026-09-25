import RequestProject.DepthThreeSpec

/-!
New open specifications for the Primus handoff of 23 September 2026.
These definitions are propositions, not proofs or axioms asserting truth.
The historical DepthThree definitions are imported unchanged.
-/

namespace Primus993

open DepthThree

variable {V : Type*} [Fintype V]

def Window (G : SimpleGraph V) : Prop :=
  G.IsTree ∧ 33 ≤ Fintype.card V ∧ Fintype.card V ≤ 38 ∧
    17 ≤ G.indepNum ∧ G.indepNum ≤ 19 ∧
    2 * G.indepNum ≤ Fintype.card V + 5

noncomputable def correction (G : SimpleGraph V) : ℤ :=
  (b G 3 : ℤ) ^ 2 - (b G 2 : ℤ) * (b G 4 : ℤ) +
    2 * (e G 3 : ℤ) * (b G 3 : ℤ) -
    (e G 2 : ℤ) * (b G 4 : ℤ) - (e G 4 : ℤ) * (b G 2 : ℤ)

def Residual (G : SimpleGraph V) : Prop :=
  0 < b G 1 ∧
    (G.indepNum - 7) * e G 4 < 3 * (G.indepNum - 3) * b G 4 ∧
    correction G < 0

def WindowTarget : Prop :=
  ∀ {V : Type*} [Fintype V] (G : SimpleGraph V), Window G →
    s G 2 * s G 4 ≤ s G 3 ^ 2

def ResidualTarget : Prop :=
  ∀ {V : Type*} [Fintype V] (G : SimpleGraph V), Window G → Residual G →
    s G 2 * s G 4 ≤ s G 3 ^ 2

def SufficientJointTarget : Prop :=
  ∀ {V : Type*} [Fintype V] (G : SimpleGraph V), Window G → Residual G →
    0 ≤ (5 * (G.indepNum : ℤ) + 17) * (e G 2 : ℤ) * (e G 4 : ℤ) +
      27 * ((G.indepNum : ℤ) - 3) * correction G

#check WindowTarget
#check ResidualTarget
#check SufficientJointTarget
#print Window
#print correction
#print Residual
#print WindowTarget
#print ResidualTarget
#print SufficientJointTarget

end Primus993
