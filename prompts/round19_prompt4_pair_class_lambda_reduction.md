# Round 19, Prompt 4: Pair-Class Lambda Reduction (Step-2 Hard Core)

**Target model:** GPT 5.2 Pro or Codex 5.3 (hybrid theory + computation)

---

## Hard facts (locked)

All scans use boundary-correct indexing and support-step prefix regime (`k < mode(I_new)`, smallest mode index).

Through `n <= 18`:

- `X_k < 0` cases: `135,976`
- all such cases occur at `step t = 2`
- global shifted-reserve scalar holds:
  - `sum_err <= D + lambda*(C10 + C01 + C11)`
  - minimal observed `lambda* = 0.05201381704686925`

Step-2 pair-class profile (`(a,b)` = sizes of first two non-leaf child rooted subtrees):

- only 9 classes appear among all negatives
- only 3 classes need nonzero shifted reserve:
  - `(2,14)` with `lambda*_pair = 0.05201381704686925`
  - `(3,13)` with `lambda*_pair = 0.04386927442810327`
  - `(2,13)` with `lambda*_pair = 0.023760967407659456`
- all other classes satisfy `sum_err <= D` (i.e. lambda `= 0`)

This is strong evidence that closure can be reduced to a finite hard core.

---

## Task 1: Prove / explain why only the three pair classes are hard

Using step-2 normal form

- `A = (1+x)^ell * I_1`
- `B = E_1`
- `(P,Q) = (I_2,E_2)`

derive a structural criterion (in terms of child-size geometry, mode distance, or diagonal support overlap) that separates:

- easy classes where `sum_err <= D` should be automatic,
- hard classes `(2,14),(3,13),(2,13)` where odd reserve is required.

If full proof is not possible, provide a falsifiable finite criterion that matches the exhaustive data through `n<=18`.

---

## Task 2: Pair-specific reserve theorem candidate

Build a theorem of the form

`sum_err <= D + lambda_{a,b} * (C10 + C01 + C11)`

for step-2 negatives, with explicit constants by pair class.

Targets to justify:

- `lambda_{2,14} <= 0.0521`
- `lambda_{3,13} <= 0.044`
- `lambda_{2,13} <= 0.024`
- `lambda_{a,b}=0` for the other six observed classes.

Need either a derivation or a reduced finite inequality family that can be machine-checked quickly.

---

## Task 3: Theoretical compression

Express the odd-residual channel in a way that isolates one scalar obstruction per pair class, then show how it is bounded by the shifted channels.

Prefer one of:

1. diagonal-Abel parity split with a pair-class monotonicity lemma,
2. two-child kernel inequality with explicit extremal coefficients,
3. Sym2/Toeplitz minor inequality that implies the shifted reserve bound.

---

## Output format

Provide:

1. theorem candidate(s) with precise hypotheses,
2. proof skeleton (or exact blocker with one missing lemma),
3. per-pair constants and whether they are proven or empirical,
4. top-5 witnesses for each hard class if any bound fails.

Include any generated JSON/tables if computation is used.
