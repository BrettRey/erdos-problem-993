# Round 20 (Instance 1): Make the Transport Idea Concrete

You previously argued for a diagonal-Abel / transport-type separator and suggested an odd-residual bound by adjacent shifted channels.

Use only the context in this prompt + your prior output. Do not claim to rerun scans.

## Locked empirical context (through n<=18)

Conventions: boundary-correct indexing, support-root step-prefix regime (`k < mode(I_new)`, smallest mode index).

- Total prefix instances checked: `2,656,341`
- `X_k < 0`: `135,976` (`5.119%`)
- all `X_k < 0` occur at `step t=2`
- `D_k/|X_k|` on negative cases: minimum `3.863636`
- global odd/even diagnostic:
  - even-diagonal bound had 0 failures
  - global `sum_err <= D` fails in 85 cases
- shifted reserve works globally:
  - `sum_err <= D + lambda*(C10+C01+C11)`
  - minimal observed `lambda* = 0.05201381704686925`
- only step-2 pair classes needing nonzero lambda:
  - `(2,14): 0.05201381704686925`
  - `(3,13): 0.04386927442810327`
  - `(2,13): 0.023760967407659456`
  - other 6 observed classes have lambda `=0`

Definitions:

- `D = sum_i Lambda_old(i) P(k-i)Q(k-i)`
- `C10 = sum_i Lambda_old(i) P(k-i-1)Q(k-i)`
- `C01 = sum_i Lambda_old(i) P(k-i)Q(k-i-1)`
- `C11 = sum_i Lambda_old(i) P(k-i-1)Q(k-i-1)`
- `sum_err` = sum over diagonals of `err_s`, where for each diagonal,
  - `m` = first index with `D_i^(s) >= 0`
  - `err_s = sum_{i>=m, D_i^(s)>0} (u_m^(s)-u_i^(s))_+ * D_i^(s)`

## Task

Turn your transport narrative into a **single explicit theorem template** with one missing lemma only.

1. State the theorem in exact symbols above (no placeholders).
2. Split it into:
   - proved algebraic identities,
   - one missing inequality lemma.
3. The missing lemma must be local/diagonal and directly testable from finite coefficient arrays.
4. Explain why this lemma naturally predicts only the three hard pair classes.

## Output format

1. `Theorem candidate`
2. `Proof skeleton (done parts)`
3. `Single missing lemma (exact inequality)`
4. `Why hard core = {(2,14),(3,13),(2,13)}`
5. `Finite test checklist`

Do not give generic prose; produce math-ready statements.
