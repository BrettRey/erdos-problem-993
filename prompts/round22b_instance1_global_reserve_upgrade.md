# Round 22b (Instance 1): Upgrade Global Reserve Theorem After n=19 Break

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context update

Conventions unchanged: support-root, boundary-correct indexing, prefix `k < mode(I_new)`.

Facts through `n<=18` (old):
- `sum_err <= D + lambda0*(C10+C01+C11)` held with
  `lambda0 = 0.05201381704686925`.

New facts at `n=19` (exhaustive in the same regime):
- total `X<0`: `428,434`
- `xneg_step_not2 = 12` (all at step 3)
- `sum_err <= D + lambda0*R_shift` fails in `33` cases
- worst required lambda:
  `lambda_needed = 0.08144365672607116`
  witness:
  - `n=19`, `g6=R?????????????????C??w?A^_?_~?`
  - `root=0`, `step=2`, `k=5`, pair `(a,b)=(3,14)`

Definitions:
- `R_shift := C10 + C01 + C11`
- `C10,C01,C11` and `sum_err,D` as in prior rounds.

## Task

Replace the broken fixed-`lambda0` theorem with the minimal robust global statement.

1. Propose one upgraded theorem candidate of the form
   `sum_err <= D + Reserve_upgraded`.
2. `Reserve_upgraded` should be the smallest natural extension of `lambda*R_shift`, e.g.
   - larger constant,
   - class-dependent constant,
   - step-dependent constant,
   - or one extra channel family.
3. Give one missing lemma only.
4. Explain how it covers both:
   - legacy `n<=18` regime,
   - new `n=19` failures.

## Output format

1. `Upgraded theorem (exact inequality)`
2. `Why fixed lambda0 fails`
3. `Minimal extension choice`
4. `Single missing lemma`
5. `Finite falsification checklist`
