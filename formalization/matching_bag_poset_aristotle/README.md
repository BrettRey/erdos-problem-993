# Matching-bag poset-layer Lean formalization

This project formalizes the abstract binary-code, order-ideal, bipartite-graph, constant
coordinate, antichain-transform, and matching-cover identities from
`matching_bag_poset_reduction_2026-08-20.md`.

- Aristotle project: `2a2d691b-91bb-4204-bb00-ad47a51861ba`.
- Independent local replay succeeded on 2026-08-20 (8,034 jobs), and the principal
  structural declarations passed a separate `#print axioms` audit.
- All principal structural declarations depend only on `propext`, `Classical.choice`,
  and `Quot.sound`.
- `ExhaustiveSmallCodes.lean` is deliberately separate: its length-four exhaustive
  census uses `native_decide` and therefore `Lean.ofReduceBool`. No structural theorem
  imports or depends on that file.

The tree-to-forest-poset construction and the blocked correction remain outside this
project; see `FORMALIZATION_NOTES.md`.

This project was edited by [Aristotle](https://aristotle.harmonic.fun).

To cite Aristotle:
- Tag @Aristotle-Harmonic on GitHub PRs/issues
- Add as co-author to commits:
```
Co-authored-by: Aristotle (Harmonic) <aristotle-harmonic@harmonic.fun>
```
