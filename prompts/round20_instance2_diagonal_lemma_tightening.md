# Round 20 (Instance 2): Tighten H_odd to a Checkable Local Inequality

You previously provided the strongest formal diagonal-Abel decomposition and proposed a local odd-diagonal hypothesis.

Use only this prompt + your prior output. Do not claim fresh scans.

## Locked context

Same conventions as before (boundary-correct, support-root step-prefix).

Empirical facts through n<=18:

- `X_k<0` only at `t=2`
- even-diagonal channel bound: 0 failures
- global `sum_err <= D` has 85 failures
- global shifted reserve closes:
  - `sum_err <= D + lambda*(C10+C01+C11)`
  - minimal observed `lambda* = 0.05201381704686925`

Definitions are fixed:

- `err_s = sum_{i>=m, D_i^(s)>0} (u_m^(s)-u_i^(s))_+ D_i^(s)`, `m = first i with D_i^(s)>=0`
- `sum_err = sum_s err_s`
- `D, C10, C01, C11` as in prior rounds.

## Task

Starting from your prior decomposition, produce a **sharpened odd-diagonal inequality** that is strictly stronger than a vague `H_odd` and directly links to shifted channels.

Required structure:

1. Derive explicit decomposition of odd `err_{2t+1}` into
   - midpoint term,
   - right-tail term.
2. Give one local inequality that upper-bounds each part by a linear combination of local shifted terms (`10`,`01`,`11`) with coefficients involving `Lambda_old(t)` and optionally `Lambda_old(t+1)`.
3. Summation over odd diagonals must yield a global bound of the form
   - `sum_odd_err <= alpha*C10 + beta*C01 + gamma*C11 + (optional neighbor-Lambda terms)`.
4. Minimize number of free coefficients (target <=3).

## Output format

1. `Exact odd decomposition`
2. `Local bound lemma (final form)`
3. `Global corollary after summing`
4. `What remains unproved`
5. `How to falsify lemma on one explicit coefficient instance`

No new notation unless immediately mapped to `{u_i^(s), D_i^(s), err_s, Lambda_old, C10,C01,C11}`.
