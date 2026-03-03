# Round 25 (Instance 2): Extended Bookkeeping that Extracts Extra Positivity

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

You cannot rely on global domination condition `(1-lambda0)R_shift >= c*Extra`:
- there are prefix instances with `R_shift=0` but `Extra>0`.

So split-reserve closure must avoid that side condition.

Current accepted bookkeeping shape for k>=2:

`Lambda_k >= D + R_shift - sum_err`.

## Task

Construct a stronger bookkeeping inequality that includes positive `Extra` directly, e.g.

`Lambda_k >= D + R_shift + eta*Extra - sum_err`

for explicit `eta>0` and explicit `Extra` in `C_ab` terms.

Then combine with a split-reserve upper bound on `sum_err` to get closure **without** any separate domination condition involving `R_shift`.

Requirements:
1. Extra must be explicit and nonnegative in the same channel family (`C_ab`).
2. Show how this avoids the `R_shift=0, Extra>0` obstruction.
3. Isolate one sharp unresolved lemma where proof difficulty concentrates.

## Output format

1. `Extended bookkeeping inequality (exact formula)`
2. `Split-reserve closure derivation`
3. `Why obstruction is removed`
4. `Single unresolved lemma`
5. `Targeted verification plan`
