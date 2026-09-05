import RequestProject.Main

/- Independent read-back of the requested statement, using fresh count definitions. -/
namespace FivefoldAudit

variable {V : Type*} [Fintype V]

noncomputable def count (G : SimpleGraph V) (d : ℕ) (extendable : Bool) : ℕ := by
  classical
  exact if d ≤ G.indepNum then
    (Finset.univ.filter fun S : Finset V =>
      S.card = G.indepNum - d ∧ G.IsIndepSet (S : Set V) ∧
        (if extendable then
          ∃ M : Finset V, G.IsMaximumIndepSet M ∧ S ⊆ M
        else ¬ ∃ M : Finset V, G.IsMaximumIndepSet M ∧ S ⊆ M)).card
  else 0

lemma count_false (G : SimpleGraph V) (d : ℕ) :
    count G d false = FivefoldForest.blockedCount G d := by
  classical
  simp [count, FivefoldForest.blockedCount, FivefoldForest.Extendable]

lemma count_true (G : SimpleGraph V) (d : ℕ) :
    count G d true = FivefoldForest.extendableCount G d := by
  classical
  unfold count FivefoldForest.extendableCount FivefoldForest.Extendable
  split
  · congr 1
    ext S
    simp
  · rfl

theorem requested_statement (G : SimpleGraph V)
    (hforest : G.IsAcyclic) (halpha : 5 ≤ G.indepNum)
    (hb2 : count G 2 false = 0) :
    5 * count G 3 false ≤ count G 3 true ∧
    5 * count G 4 false ≤ count G 4 true := by
  simp only [count_false] at hb2
  simp only [count_false, count_true]
  exact FivefoldForest.fivefold_of_forest G hforest halpha hb2

end FivefoldAudit

#check @FivefoldForest.fivefold_of_forest
#print axioms FivefoldAudit.requested_statement
#print axioms FivefoldForest.fivefold_of_forest
#print axioms FivefoldForest.fivefold
#print axioms FivefoldForest.good_or_starProf
#print axioms FivefoldForest.bD_one_eq_zero_of_bD_two
#print axioms FivefoldForest.forestLike_of_isAcyclic
#print axioms FivefoldForest.eD_split
#print axioms FivefoldForest.bD_le_of_split
#print axioms FivefoldForest.split_good
#print axioms FivefoldForest.star_good_or_star
#print axioms FivefoldForest.Pendant.eD_two_le
#print axioms FivefoldForest.Pendant.bD_two_le
#print axioms FivefoldForest.Pendant.eD_one
#print axioms FivefoldForest.Pendant.bD_one
#print axioms FivefoldForest.Pendant.alpha_core_one
#print axioms FivefoldForest.cone_eD_zero
#print axioms FivefoldForest.cone_eD_succ
#print axioms FivefoldForest.cone_bD_one
#print axioms FivefoldForest.card_le_two_alpha
#print axioms FivefoldForest.pendant_one_good
#print axioms FivefoldForest.pendant_three_good
#print axioms FivefoldForest.pendant_four_le_good
