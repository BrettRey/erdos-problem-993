# Round 24 (Instance 1): Repair the Defect-to-R_shift Bound After Falsification

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

Conventions unchanged: support-root, boundary-correct indexing, prefix `k < mode(I_new)` (smallest maximizing index).

Constants:
- `lambda0 = 0.05201381704686925`
- `lambda19 = 0.08144365672607116`
- `delta = lambda19 - lambda0 = 0.02942983967920191`

Defect framework (already accepted):

`sum_err <= D + lambda0*R_shift + Defect`,

with `R_shift = C10 + C01 + C11`.

New exhaustive falsification on full `n<=19` X<0 corpus (`428,434` cases):
- Candidate `Defect <= delta*R_shift` fails in **1** case.
- Worst ratio `Defect/(delta*R_shift) = 1.0225408414008559`.
- Witness: `n=19`, `g6=R?????????????????C??w?A^_?_~?`, `root=0`, `step=2`, `k=5`, `(a,b)=(3,14)`.

The proposed local lemma from Round 23 output also fails badly:
- `d_t <= e_{t+1} + delta*Lambda_t*K_t` fails in **5,417** cases.
- Worst local ratio `lhs/rhs = 4.470112793366624`.

## Task

Produce a **repaired** defect theorem that survives these failures.

You must give one of:

1. Pure-R_shift repair:
   `Defect <= delta'*R_shift` with explicit `delta' > delta` and a sharp lower bound from the witness.

or

2. Mixed repair:
   `Defect <= delta*R_shift + c*Extra(k)` with explicit nonnegative `Extra(k)` built from `C_ab` channels.

Requirements:
1. Use the witness above to derive a nontrivial lower bound on any proposed constants.
2. Explain exactly why the old local lemma failed and what structural change your repaired lemma uses.
3. Give one single "missing lemma" whose proof would close the repaired theorem.
4. Keep formulas explicit; no abstract existence statements.

## Output format

1. `Repaired defect bound (exact formula)`
2. `Constant lower bounds forced by witness`
3. `Derivation to global theorem`
4. `Single missing lemma`
5. `Targeted falsification checklist`
