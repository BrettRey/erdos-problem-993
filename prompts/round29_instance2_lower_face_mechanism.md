# Round 29 (Instance 2): Mechanism for the Lower Face (G_k Bound)

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

From Round 28 outputs:
- exact object: `G_k := Lambda_k - D + sum_all`
- exact identity: `G_k = X_k + sum_all`
- and channel split form: `G_k = R_shift + H_k` with sign-indefinite remainder.

From n<=21 alpha calibration:
- feasible lower-face constant is at least as low as
  `alpha_* = 0.21034113597068071`.

Target lower face:
`G_k >= alpha_* * R_shift`.

## Task

Build one concrete proof mechanism for the lower face:

1. Start from exact identity and remainder taxonomy (no inherited bookkeeping).
2. Propose one local/paired inequality family (e.g., odd-even paired diagonals, gap bands, or packet inequalities) whose summation implies `G_k >= alpha_* R_shift`.
3. Make all index shifts explicit in `C_ab` notation.
4. Identify exactly one bottleneck sublemma and why it is plausible in tree-realizable regime.

## Output format

1. `Exact lower-face target restatement`
2. `Proposed mechanism (local-to-global)`
3. `Explicit summed inequality implication`
4. `Single bottleneck sublemma`
5. `Minimal verification harness`
