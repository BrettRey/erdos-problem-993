# Round 28 (Instance 2): Exact Decomposition of G_k := Lambda_k - D + sum_all

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

- Inherited bookkeeping assumptions were wrong.
- We now work with `sum_all = sum_s err_s` (all diagonals), not odd-only.
- The calibrated target is of the form
  `G_k := Lambda_k - D + sum_all >= alpha*R_shift`.

Definitions:
- `R_shift = C10 + C01 + C11`.
- `C_ab(k) = sum_i Lambda_old(i) P(k-i-a)Q(k-i-b)`.

## Task

Identity-first rebuild:

1. Derive one exact decomposition of `G_k` into explicit channel/remainder pieces.
2. Separate terms by sign-certifiable vs sign-indefinite.
3. Identify exactly where previous wrong-sign collapse happened.
4. Propose one viable global inequality target that could imply
   `G_k >= alpha*R_shift` with `alpha≈0.2437`.
5. Provide one explicit non-viable target (and why it must fail).

Keep notation compatible with `C_ab`, `R_shift`, and optional radius-2 extras.

## Output format

1. `Exact identity for G_k`
2. `Sign-sensitive term taxonomy`
3. `Viable inequality target`
4. `Non-viable target and obstruction`
5. `Minimal check harness`
