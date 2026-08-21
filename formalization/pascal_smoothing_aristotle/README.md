# Pascal-smoothing Lean formalization

This project formalizes the universal Pascal-smoothing lemma and its stated
consequences from `pascal_smoothing_shadow_lemma_2026-08-20.md`.

- Aristotle project: `6519b82f-e80b-4c41-b633-990b72c7e3b9`.
- Independent local replay: `lake build` succeeded on 2026-08-20 (8,028 jobs).
- Local `#print axioms` checks on the principal declarations found only `propext`,
  `Classical.choice`, and `Quot.sound`.
- No `sorry`, `admit`, `axiom`, `implemented_by`, `native_decide`, or `unsafe`
  declaration occurs in `RequestProject/`.

The strict depth-three conclusion correctly assumes
`erasureShadow Δ 2 * erasureShadow Δ 4 > 0`; without it the claimed strict inequality
can degenerate to `0 = 0`.

This project was edited by [Aristotle](https://aristotle.harmonic.fun).

To cite Aristotle:
- Tag @Aristotle-Harmonic on GitHub PRs/issues
- Add as co-author to commits:
```
Co-authored-by: Aristotle (Harmonic) <aristotle-harmonic@harmonic.fun>
```
