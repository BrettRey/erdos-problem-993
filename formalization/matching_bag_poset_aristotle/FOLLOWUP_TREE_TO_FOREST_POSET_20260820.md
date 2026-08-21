# Aristotle continuation: tree maximum matchings produce forest-poset codes

Continue the existing matching-bag poset project after `PascalBridge.lean`. Close the
remaining equation (2) bridge from an arbitrary finite tree with a maximum matching to
the already-formalized order-ideal code. Treat
`matching_bag_poset_reduction_2026-08-20.md` and
`matching_bag_csp_attack_2026-08-20.md` as proposed mathematics to audit.

## Final target

For every finite tree `T` and every maximum matching `M`, construct the binary code of
maximum independent sets (equivalently complements of minimum vertex covers) on the
matching bags. Prove that, after consistently orienting every matched edge by the
bipartition and propagating/removing forced variables, this code is isomorphic to

```text
idealCode P × {one constant word on c coordinates},
```

where `P` is a finite poset whose cover graph is a forest and `c` includes all unmatched
singleton bags and forced matching bags. Deduce that its projection polynomial is

```text
(1+t)^c * indepPoly P
```

and hence consumes the existing `PascalBridge` theorems.

## Required checkpoints

1. **Matching bags.** Contract or index the edges of `M`, retain unmatched singleton
   vertices, and prove every independent set chooses at most one endpoint per matched
   bag. Define the port-labelled quotient tree and prove the feasible bag-assignment
   bijection, including preservation of cardinality/number of nonempty bags.
2. **Maximum layer.** Prove the tree/bipartite König equality needed here rather than
   postulating it: a minimum cover has `|M|` vertices, hence exactly one endpoint from
   every matched edge and no unmatched vertex. Reuse `cover_card_eq_matching_card`.
3. **Inequality system.** With matched edge `L_i R_i` and Boolean variable `x_i`, prove
   every nonmatching edge `L_i R_j` gives `x_j ≤ x_i`; edges incident to unmatched
   vertices give unary forced values.
4. **Consistency and forest structure.** Prove the unary values are consistent from
   maximality of `M` (an inconsistency would yield an augmenting path). Prove the
   undirected comparison graph is a forest. After propagating and deleting forced
   variables, construct the remaining poset `P` and prove its cover graph is still a
   forest; do not merely assume acyclicity or antisymmetry.
5. **Code equivalence.** Give explicit mutually inverse maps between maximum independent
   sets/minimum covers and order ideals of `P` together with the constant coordinates.
6. **Polynomial consequence.** Derive the constant-coordinate factor and instantiate
   the existing poset/Pascal bridge, including the depth-eight, `M≤33`, and depth-three
   reserve statements for the extendable profile.

## Soundness and grading

- No `sorry`, `admit`, `axiom`, `implemented_by`, `native_decide`, weakened theorem,
  bounded enumeration, or hypotheses that simply postulate König equality, consistency,
  forest structure, or the desired code equivalence.
- The optional `ExhaustiveSmallCodes.lean` target remains separate and irrelevant.
- Add source-labelled files and update `FORMALIZATION_NOTES.md` with exact traceability.
- Run `lake build` and `#print axioms` for every principal bridge theorem.
- If full Mathlib graph infrastructure makes one checkpoint too large, return the
  strongest genuine prefix and the first exact missing theorem. If any claim is false,
  return the smallest explicit tree/matching witness and a repaired statement.
