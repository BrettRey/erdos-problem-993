# Round 25 Prompt Prep (2026-03-02)

Purpose: prepare next prompts using Round 24 outputs + additional diagnostics.

## Locked carry-forward facts

- `lambda0 = 0.05201381704686925`
- `lambda19 = 0.08144365672607116`
- `delta = lambda19 - lambda0 = 0.02942983967920191`

From `results/round24_falsification_checks_n19.json`:
- `Defect <= delta*R_shift` fails in 1 case.
- Worst ratio: `1.0225408414008559` (same `(3,14)` witness).
- Local one-step lemma fails massively (`5417` failures).

From `results/round24_defect_extra_fit_n19_full.json`:
- Residual after `delta*R_shift` is tiny.
- For `Extra = C20+C02+C21+C12`, needed residual coefficient is
  `c = 0.0015447425039931744` on full `n<=19` X<0 corpus.

From `results/round24_lambda_c_tradeoff_n19.json`:
- Mixed-reserve tradeoff is explicit.
- With `Extra = C20+C02+C21+C12` and `c=0.06853082706728619`,
  needed `lambda` returns to `lambda0` on this corpus.

## New obstruction diagnostic

From `results/round24_split_domination_check_n19.json`:
- For `c=0.06853082706728619`, `Extra=C20+C02+C21+C12`, and all prefix instances:
  `max (c*Extra)/((1-lambda0)*R_shift) = Infinity`.
- Witness has `R_shift=0` but `Extra>0`.

Implication:
- Any proof route that requires global side condition
  `(1-lambda0)R_shift >= c*Extra` is invalid as a theorem skeleton.
- Split-reserve route must instead prove an **extended bookkeeping inequality** that extracts positive `Extra` directly in the lower bound for `Lambda_k`, or avoid this condition entirely.

## Bookkeeping validity check (xneg corpus)

Artifact:
- `results/round25_bookkeeping_validity_n19_xneg.json`

Tested on full `n<=19` X<0 corpus (`428,434` cases):

Candidate inequality (as previously assumed in some prompts):
`Lambda_k >= D + R_shift - sum_err`.

Results:
- Using `sum_err = sum_s err_s` (all diagonals): **9,223 failures**.
- Using `sum_err = sum_{s odd} err_s` (odd-only): **428,434 failures**.

Worst all-diagonal slack witness:
- `n=19`, `(a,b)=(2,15)`, `step=2`, `k=7`,
- slack = `Lambda_k - (D+R_shift-sum_err)` = `-56,361,669`.

Conclusion:
- This bookkeeping inequality cannot be treated as an accepted base identity.
- Future theorem variants must avoid relying on it directly.
