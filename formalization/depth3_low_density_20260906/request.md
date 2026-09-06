# Frozen Lean request: blocked-shadow bound and low-density depth-three theorem

Formalize the actual finite-graph statements in RequestProject/DepthThreeSpec.lean.
This file has been compiled locally with Lean 4.28.0; its target definitions are
SPECIFICATIONS, not theorems. Do not change their meanings or guards.

Required COMPLETE deliverables, all without extra assumptions:

1. `blocked_shadow_of_forest : DepthThree.BlockedShadowTarget`.
   Every finite forest with alpha >= 4 satisfies (alpha-4)b2 <= 6b3.
2. `low_density_of_forest : DepthThree.LowDensityTarget`.
   Every finite forest with alpha in {17,18,19} and
   3(alpha-3)b4 <= (alpha-7)e4 satisfies s2*s4 < s3^2.

Disconnected forests and isolated vertices are included. There is no restriction
on vertex count, matching deficiency, b1, or b2. Extending means containment in a
MAXIMUM independent set, never merely a maximal one. Counts at negative sizes
are zero, as encoded by the d <= indepNum guard. Count partial sets once, not once
per maximum completion. Equivalent internal definitions need proved bridges to
the exact specifications. This is not a request to solve Erdős #993.

## Written argument and existing library

See WRITTEN_PROOF.md Section 9 for the new argument. The earlier matching-poset
Lean source is supplied in RequestProject. It is legitimate to reuse its proved
declarations. They were built with the supplied pinned Lean/Mathlib versions.

IMPORTANT INTERFACE OBLIGATION: TreeCodeBridge.erasure_depth_three_reserve is
about `D.erasure`, a projected-code cardinality. It is NOT already a theorem
about `DepthThree.e`. Prove a finite-graph maximum matching exists, identify
`card D.Bag` with G.indepNum, identify D.maxIndepSets with actual maximum
independent sets, and prove the cardinality-preserving bijection between
extendable vertex sets and projected partial words (including flips and guard).
The bridge may use a simpler equivalent route, but must not be assumed.

For the blocked-shadow bound: the poset representation implies every blocked
independent vertex set contains a blocked subset of cardinality <= 2. A forbidden
constant-coordinate value is unary; otherwise nonextension of a partial ideal
indicator is witnessed by two specified values violating an order comparison.
Use the existing PosetCode.codeProj_idealCode characterization, together with the
graph/partial-assignment bridge. For a blocked (alpha-2)-set, at least alpha-4
one-vertex deletions preserve a chosen certificate. Each independent
(alpha-3)-set has at most six one-vertex extensions: it leaves three matching
bags empty, each of size <= 2. Double count the incidence relation. Do NOT
assume every deletion of a blocked set remains blocked.

For every graph, (alpha-3)e3 >= 4e4 by double counting extendable one-vertex
extensions. The code reserve transported to graph counts is
32(alpha-2)e2e4 <= 27(alpha-3)e3^2 for alpha >= 4.
For a=alpha, c=(5a+17)/(27(a-3)), h=(a-7)/(3(a-3)), v=b4/e4,
the full margin is at least
(c-v)e2e4 + (h-v)e4b2 + b3^2.
When a=17,18,19 and v<=h, the first coefficient is strictly positive:
c-h=4(20-a)/(27(a-3)); e2 and e4 are positive. Integer-cleared arithmetic
avoids denominators. The thresholds are respectively 5/21, 11/45, 1/4.

## Grading and proof integrity

G0 COMPLETE requires both exact target propositions proved. Supporting targets:
G1 actual graph/count↔code bridge and graph Pascal reserve;
G2 unary-or-pair certificate and blocked-shadow inequality;
G3 extendable incidence inequality and positivity;
G4 numerical assembly and final low-density theorem.
Expose the sublemmas as traceable declarations. Partial progress is welcome and
must be reported honestly with the unresolved bridge named. An exact refutation
or identified missing inference is a useful success mode. Do not silently weaken
the statements to profiles satisfying assumed recurrences or the desired bounds.

No sorry/admit, added axioms, opaque unproved hypotheses, implemented_by,
native_decide, or unchecked external certificates in an accepted proof. Use
kernel-checkable proof terms. A proof that only the scalar implication holds is
PARTIAL, not the graph theorem. Keep original dependency sources unchanged when
practical and add modules. Include final theorem modules in default Lake build,
a declaration map, precise completion status, and #print axioms results. The
local verifier will rebuild, inspect the exact statements, and audit axioms.

Source note SHA-256: 4d77890b084b1a7aec056217b84d2606f17c78b4ec5772aff2401426d300a8bb.
Repository HEAD: c1abc1ac7b4ebfedfc28d482ec7440e2da5a3bec (new source uncommitted).
