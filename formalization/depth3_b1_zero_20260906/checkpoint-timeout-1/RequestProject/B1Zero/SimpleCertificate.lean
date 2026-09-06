import RequestProject.RootedForestCertificate

/-!
# A second finite certificate for rooted forests of size ten

The blocked bound proved in `RequestProject/B1Zero/FourStructure.lean` produces the
coefficient vector `(73, 24, 3 ; 15, 6, 1)` rather than the sharper `(58, 21, 3 ; 30, 9, 1)`
of the written note: only the union bound over the four forbidden vertices is used, and no
concentration step.  The resulting numerator is still below one quarter of the denominator
for every plane rooted forest with ten nodes, which is what
`rooted_forest_simple_bound` records.  As with `RootedForestCertificate.all_codes_checked`
the check is a kernel evaluation over the 16796 enumerated codes.
-/

namespace DepthThree.RootedForestCertificate

/-- The numerator produced by the union-bound route. -/
def simpleNumerator (f : Forest) : ℕ :=
  let (p, q) := profile f
  73 * p.c0 + 24 * p.c1 + 3 * p.c2 + 15 * q.c0 + 6 * q.c1 + q.c2

/-- Four times the union-bound numerator is at most the denominator. -/
def simpleWithinBound (f : Forest) : Bool :=
  decide (4 * simpleNumerator f ≤ denominator f)

set_option maxRecDepth 100000 in
set_option maxHeartbeats 100000000 in
theorem all_simple_codes_checked : (enumerate 10).all simpleWithinBound = true := by
  decide +kernel

/-- **The finite certificate.**  For every plane rooted forest with ten nodes,
`4 · (73 p₀ + 24 p₁ + 3 p₂ + 15 q₀ + 6 q₁ + q₂) ≤ 126 p₀ + 84 p₁ + 36 p₂ + 9 p₃ + p₄`. -/
theorem rooted_forest_simple_bound (f : Forest) (h : size f = 10) :
    4 * simpleNumerator f ≤ denominator f :=
  of_decide_eq_true (List.all_eq_true.mp all_simple_codes_checked f (enumeration_complete h))

#print axioms rooted_forest_simple_bound

end DepthThree.RootedForestCertificate
