import RequestProject.DepthThreeTargets

example : DepthThree.BlockedShadowTarget := DepthThree.blocked_shadow_of_forest
example : DepthThree.LowDensityTarget := DepthThree.low_density_of_forest

#print DepthThree.BlockedShadowTarget
#print DepthThree.LowDensityTarget
#print DepthThree.DepthThreeStrict
#print DepthThree.Extendable
#print DepthThree.e
#print DepthThree.b
#print axioms DepthThree.blocked_shadow_of_forest
#print axioms DepthThree.low_density_of_forest
#print axioms MatchingBag.exists_treeMatching_of_isAcyclic
#print axioms MatchingBag.TreeMatching.indepNum_eq_card_Bag
#print axioms MatchingBag.TreeMatching.isMaximumIndepSet_iff_mem_maxIndepSets
#print axioms MatchingBag.TreeMatching.card_extSets_fiber
#print axioms MatchingBag.TreeMatching.e_eq_erasure
#print axioms MatchingBag.TreeMatching.exists_small_blocked_subset
#print axioms MatchingBag.TreeMatching.blocked_shadow
#print axioms MatchingBag.TreeMatching.extendable_incidence
