import RequestProject.B1Zero.Target
import RequestProject.DepthThreeTargets

example : DepthThree.B1ZeroWindowTarget := DepthThree.b1_zero_window
example : DepthThree.BlockedShadowTarget := DepthThree.blocked_shadow_of_forest
example : DepthThree.LowDensityTarget := DepthThree.low_density_of_forest

#print DepthThree.B1ZeroWindowTarget
#print DepthThree.DepthThreeStrict
#print DepthThree.Extendable
#print DepthThree.e
#print DepthThree.b
#print axioms DepthThree.b1_zero_window
#print axioms DepthThree.blocked_shadow_of_forest
#print axioms DepthThree.low_density_of_forest
#print axioms MatchingBag.TreeMatching.density_window
#print axioms MatchingBag.TreeMatching.density_small
#print axioms MatchingBag.TreeMatching.density_four_aux
#print axioms MatchingBag.TreeMatching.flex_certificate
#print axioms MatchingBag.B1Zero.exists_forest_repr
#print axioms MatchingBag.B1Zero.CoronaData.cntRev_move_roots
#print axioms DepthThree.RootedForestCertificate.rooted_forest_simple_bound
