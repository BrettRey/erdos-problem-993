This project was edited by [Aristotle](https://aristotle.harmonic.fun).

To cite Aristotle:
- Tag @Aristotle-Harmonic on GitHub PRs/issues
- Add as co-author to commits:
```
Co-authored-by: Aristotle (Harmonic) <aristotle-harmonic@harmonic.fun>
```

# Matching-bag forest-poset bridge

Formalization of `matching_bag_poset_reduction_2026-08-20.md` and its two continuations.
The declaration-by-declaration audit is in `FORMALIZATION_NOTES.md`.

The completed 2026-08-21 continuation is Aristotle task
`fcbe1949-5ec5-4c84-8e73-a085b2352297`. The new `KonigHall`, `TreeMatching`,
`ForestLemmas`, `BagPoset`, and `TreeCodeBridge` modules replay locally. Direct axiom
checks on König, equations (2) and (6), the code equivalence, the erasure bridge, and the
depth-three reserve report only `propext`, `Classical.choice`, and `Quot.sound`. The old,
unrelated `ExhaustiveSmallCodes.lean` target uses the already documented `native_decide`
enumeration and was not awaited during this replay.
