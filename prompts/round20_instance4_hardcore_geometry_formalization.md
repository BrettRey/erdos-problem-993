# Round 20 (Instance 4): Formalize the Hard-Core Geometry Criterion

You previously proposed that only a finite step-2 hard core matters and suggested an overlap-width / boundary mechanism.

Use only this prompt + your prior output. Do not claim fresh scans.

## Locked context

Through n<=18 (same conventions):

- all negatives at step `t=2`
- pair classes among negatives: 9 total
- nonzero lambda only for `(2,14),(3,13),(2,13)`
- pairwise lambda maxima:
  - `0.05201381704686925`, `0.04386927442810327`, `0.023760967407659456`
- all other observed classes have lambda `=0` under
  `sum_err <= D + lambda*(C10+C01+C11)`.

## Task

Turn your qualitative criterion into a **precise mathematical classifier** on step-2 instances.

1. Define explicit quantities (e.g., overlap width `W(k)`, boundary-touch indicators, mode distance) with exact formulas.
2. State a classifier condition `Hard(instance)` in these quantities.
3. Show analytically why `Hard(instance)` implies need for nonzero shifted reserve.
4. Show analytically why failure of `Hard(instance)` implies `sum_err <= D` (or identify one missing sublemma if needed).
5. Reduce this to a finite inequality family for `a=2` and `a=3` (the only hard first-child sizes in data).

## Output format

1. `Formal classifier definition`
2. `Implication to odd residual magnitude`
3. `Finite inequality family (a=2 and a=3)`
4. `Which part is proved vs conjectural`

No vague terms like “short overlap” without a formula.
