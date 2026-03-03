# Round 25 (Instance 1): Defect Bound with Tiny Residual Extra

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

Conventions unchanged: support-root, boundary-correct indexing, prefix `k < mode(I_new)`.

Constants:
- `lambda0 = 0.05201381704686925`
- `lambda19 = 0.08144365672607116`
- `delta = lambda19 - lambda0 = 0.02942983967920191`

Known from Round 24:
- `Defect <= delta*R_shift` fails in exactly 1 n<=19 case.
- Worst ratio: `1.0225408414008559`.
- Residual can be absorbed by tiny second-shift coefficients.
- For `Extra = C20+C02+C21+C12`, needed residual coefficient is
  `c_resid = 0.0015447425039931744` on full n<=19 X<0 corpus.

## Task

Propose and justify a repaired defect theorem of the form

`Defect <= delta*R_shift + c_resid*Extra`,

with explicit `Extra` and explicit constants.

Requirements:
1. Use the single-failure witness to derive lower bounds on `c_resid`.
2. Explain why this route is strictly better than inflating `delta` alone.
3. Give one single missing lemma that would prove this form.
4. Keep all formulas explicit in `C_ab` channels.

## Output format

1. `Repaired defect theorem (exact constants)`
2. `Witness-forced lower bounds`
3. `Derivation to global reserve inequality`
4. `Single missing lemma`
5. `Finite falsification checklist`
