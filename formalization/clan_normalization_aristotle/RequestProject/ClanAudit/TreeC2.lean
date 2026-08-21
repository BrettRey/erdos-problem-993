import RequestProject.ClanAudit.C2Family
import RequestProject.ClanAudit.SpiderRecognition
import RequestProject.ClanAudit.ConnectorRecognition

/-!
# Abstract trees with at most two branch vertices

A *branch vertex* of a graph is a vertex of degree at least three.  Combining the two
recognition theorems `spiderIsoOfTree` and `connIsoOfTree` we show that every finite tree
with at most two branch vertices is a `C₂` model, and therefore has a log-concave
independence-count sequence and a log-concave independence polynomial.
-/

namespace ClanAudit

open SimpleGraph Finset

variable {V : Type*} [Fintype V] [DecidableEq V] {G : SimpleGraph V} [DecidableRel G.Adj]

variable (G) in
/-- The branch vertices of `G`: the vertices of degree at least three. -/
def branchVerts : Finset V := univ.filter (fun v => 3 ≤ G.degree v)

omit [DecidableEq V] in
theorem mem_branchVerts {v : V} : v ∈ branchVerts G ↔ 3 ≤ G.degree v := by
  simp [branchVerts]

/-- **Recognition.**  A finite tree with at most two branch vertices is a `C₂` model:
it is isomorphic either to a spider (this covers the single vertex, paths and all
one-hub trees) or to a two-hub connector tree. -/
theorem isC2Model_of_isTree_of_branchVerts_card_le_two (hG : G.IsTree)
    (hb : (branchVerts G).card ≤ 2) : IsC2Model G := by
  by_cases hpair : ∃ a b : V, a ∈ branchVerts G ∧ b ∈ branchVerts G ∧ a ≠ b
  · obtain ⟨a, b, ha, hbmem, hab⟩ := hpair
    have hsub : ({a, b} : Finset V) ⊆ branchVerts G := by
      intro x hx
      rcases Finset.mem_insert.mp hx with rfl | hx
      · exact ha
      · rw [Finset.mem_singleton] at hx
        exact hx ▸ hbmem
    have hcard : (branchVerts G).card ≤ ({a, b} : Finset V).card := by
      rw [Finset.card_insert_of_notMem (by simpa using hab), Finset.card_singleton]
      exact hb
    have heq : branchVerts G = {a, b} := (Finset.eq_of_subset_of_card_le hsub hcard).symm
    have hdeg : ∀ w : V, w ≠ a → w ≠ b → G.degree w ≤ 2 := by
      intro w hw1 hw2
      by_contra hcon
      have : w ∈ branchVerts G := mem_branchVerts.mpr (by omega)
      rw [heq] at this
      rcases Finset.mem_insert.mp this with h | h
      · exact hw1 h
      · exact hw2 (Finset.mem_singleton.mp h)
    exact Or.inr ⟨_, _, _, _, _, ⟨connIsoOfTree hG hab hdeg⟩⟩
  · push_neg at hpair
    have hV : Nonempty V := hG.isConnected.nonempty
    obtain ⟨h, hh⟩ : ∃ h : V, ∀ v : V, v ≠ h → G.degree v ≤ 2 := by
      by_cases hne : (branchVerts G).Nonempty
      · obtain ⟨a, ha⟩ := hne
        refine ⟨a, fun v hv => ?_⟩
        by_contra hcon
        exact hv (hpair v a (mem_branchVerts.mpr (by omega)) ha)
      · obtain ⟨x⟩ := hV
        refine ⟨x, fun v _ => ?_⟩
        by_contra hcon
        exact hne ⟨v, mem_branchVerts.mpr (by omega)⟩
    exact Or.inl ⟨_, _, ⟨spiderIsoOfTree hG h hh⟩⟩

/-- **Log-concavity of independence counts** for a finite tree with at most two branch
vertices. -/
theorem indepCount_logConcave_of_isTree_of_branchVerts_card_le_two (hG : G.IsTree)
    (hb : (branchVerts G).card ≤ 2) (j : ℕ) :
    indepCount G j * indepCount G (j + 2) ≤ indepCount G (j + 1) * indepCount G (j + 1) :=
  indepCount_logConcave_of_isC2Model
    (isC2Model_of_isTree_of_branchVerts_card_le_two hG hb) j

/-- **Log-concavity of the independence polynomial** for a finite tree with at most two
branch vertices. -/
theorem indepPoly_logConcave_of_isTree_of_branchVerts_card_le_two (hG : G.IsTree)
    (hb : (branchVerts G).card ≤ 2) (j : ℕ) :
    (indepPoly G).coeff j * (indepPoly G).coeff (j + 2)
      ≤ (indepPoly G).coeff (j + 1) * (indepPoly G).coeff (j + 1) :=
  indepPoly_logConcave_of_isC2Model
    (isC2Model_of_isTree_of_branchVerts_card_le_two hG hb) j

end ClanAudit
