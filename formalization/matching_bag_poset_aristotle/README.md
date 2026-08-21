# Matching-bag poset-layer Lean formalization

This project formalizes the abstract order-ideal/graph-transform layer and connects it
to the universal Pascal-smoothing theorem.

- Aristotle project: `2a2d691b-91bb-4204-bb00-ad47a51861ba`.
- Independent local replay of the Pascal continuation succeeded through the structural
  bridge; all principal bridge declarations depend only on `propext`,
  `Classical.choice`, and `Quot.sound`.
- The separate optional length-four code census retains its documented `native_decide`
  dependency; no structural or Pascal-bridge theorem imports it.

The tree-to-forest-poset construction is the active continuation; the blocked
correction remains future work. See `FORMALIZATION_NOTES.md` for exact traceability.

This project was edited by [Aristotle](https://aristotle.harmonic.fun).

To cite Aristotle:
- Tag @Aristotle-Harmonic on GitHub PRs/issues
- Add as co-author to commits:
```
Co-authored-by: Aristotle (Harmonic) <aristotle-harmonic@harmonic.fun>
```
