# Round 21 (Instance 2): Design the Odd-Diagonal Compensation Mechanism

You gave the best odd decomposition (`midpoint + right-tail`) and showed synthetic falsification for naive local bounds.

Use only your prior output + this prompt.

## Locked facts

Through `n<=18` in the step-prefix regime:

- global bound holds:
  `sum_err <= D + lambda0*(C10+C01+C11)`, `lambda0=0.05201381704686925`
- local odd lemma with `Lambda_t` fails massively (`38,642` failures on odd checks)
- local odd lemma with `(Lambda_t+Lambda_{t+1})` still fails (`1,600` failures)

So any successful theory must include **cross-diagonal compensation**.

## Task

Starting from your exact odd split
`err_{2t+1} = err_mid + err_tail`,
derive a compensation identity/inequality of the form:

`sum_t err_{2t+1} <= sum_t F_t - sum_t G_t`

where:

- `F_t` is a local shifted reserve term (built from `Lambda_old`, `10`, `01`, `11`),
- `G_t` is a telescoping or monotone correction term that captures cancellation across neighboring odd diagonals.

Goal is to recover a global bound by summing, even though local inequalities fail.

## Output requirements

1. Exact compensated odd decomposition (not just heuristic)
2. Explicit formulas for `F_t` and `G_t`
3. One final global corollary in `C10,C01,C11`
4. Precise missing sublemma if any (single statement)
5. One tiny synthetic example showing compensation succeeds where local fails

No per-diagonal bound claims unless explicitly labeled false.
