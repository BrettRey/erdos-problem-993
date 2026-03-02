# Round 23 (Instance 1): Prove a Nontrivial Defect-to-Channel Bound

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

Same conventions as prior rounds:
- support-root
- boundary-correct indexing
- prefix `k < mode(I_new)` (smallest maximizing index)

Current best theorem shape:

`sum_err <= D + lambda * R_shift`, where `R_shift = C10 + C01 + C11`.

Locked `n<=19` facts:
- `lambda0 = 0.05201381704686925` fails (33 failures)
- upgraded constant works: `lambda19 = 0.08144365672607116`
- all tested `n=19` negatives satisfy global bound with `lambda19`.

You previously reframed the gap via defect:
- `d_t := (err_{2t+1} - lambda0*Lambda_t*(10+01+11))_+`
- `e_t := (Lambda_t*P(k-t)Q(k-t) - err_{2t})_+`
- `Defect := max_t (sum_{u<=t} d_u - sum_{u<=t+1} e_u)_+`
- then `sum_err <= D + lambda0*R_shift + Defect`.

## Task

Propose and justify one **nontrivial** inequality of the form:

`Defect <= delta * R_shift`  or  `Defect <= c * Extra(k)`

where `Extra(k)` is an explicit nonnegative channel built from `Lambda_old` and shifted `(P,Q)` products.

Requirements:
1. Not tautological (cannot define Extra := Defect).
2. Explicit formula(s) only.
3. Must imply upgraded global theorem with concrete constant(s).
4. If not fully provable, isolate a single sharp missing lemma.

## Output format

1. `Candidate defect bound (exact formula)`
2. `Derivation to global theorem`
3. `Single missing lemma`
4. `Why this is strictly stronger than bookkeeping`
5. `Finite falsification checklist`
