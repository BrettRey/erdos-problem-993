# Aristotle continuation: connect the poset transform to Pascal smoothing

Continue the existing matching-bag poset project. The uploaded
`PascalSmoothing.lean` is the independently replayed formalization of the companion
Pascal-smoothing theorem. Integrate or adapt it cleanly into this project and close
equations (7) and (7a) from `matching_bag_poset_reduction_2026-08-20.md`.

## Target

For any finite poset `P` and `c : ℕ`, let

```text
M = Fintype.card P + c
E d = coeff ((1 + X)^c * indepPoly P) (M - d).
```

Using the already-proved antichain coefficient formula, prove that `E` is precisely the
erasure shadow of the downward-closed family of antichains of `P`. In particular prove,
without assuming log-concavity:

1. the antichains of `P` form a downward-closed family;
2. `E d = Σ_j antichainCount P j * choose (M-j) d`;
3. the denominator-cleared Pascal-smoothing inequality for adjacent `E` values;
4. log-concavity through defect depth eight for arbitrary `M`;
5. log-concavity at every interior defect when `M ≤ 33`;
6. the quantitative depth-three reserve

```text
32 * (M-2) * E 2 * E 4 ≤ 27 * (M-3) * (E 3)^2
```

for `M ≥ 4`, plus the strict form under the necessary positivity hypothesis.

The imported Pascal development uses subsets of `Fin M`, whereas `P` is an arbitrary
finite type. Either prove an equivalence-invariance/transport lemma or generalize the
Pascal definitions to arbitrary finite ground types. Do not assume that `P = Fin r` or
postulate equality of the two profiles.

## Soundness and deliverable

- No `sorry`, `admit`, `axiom`, `implemented_by`, `native_decide`, or weakened target
  in the structural bridge.
- The optional `ExhaustiveSmallCodes.lean` computation may remain separate and may
  retain its documented `native_decide`; no new theorem may depend on it.
- Add a dedicated bridge file, update `FORMALIZATION_NOTES.md` and traceability, run
  `lake build`, and report `#print axioms` for the principal bridge theorems.
- If a definition mismatch makes any displayed identity false, return the smallest
  explicit mismatch and a repaired theorem rather than coercing the definitions.
