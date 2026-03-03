# Round 26 (Instance 1): Re-derive Exact CB/Abel Identity After Bookkeeping Failure

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

Conventions unchanged: support-root, boundary-correct indexing, prefix `k < mode(I_new)`.

Critical correction from local scans (`n<=19`, X<0 corpus):
- The inequality `Lambda_k >= D + R_shift - sum_err` is **not valid** in general.
- Fail counts on X<0 cases:
  - using `sum_err = sum_s err_s`: `9,223` failures,
  - using odd-only `sum_err`: `428,434` failures.

So this route is dead and must not be reused.

## Task

Start over from the exact Cauchy–Binet expansion and derive a **correct exact identity** for

`Lambda_k - D`

in terms of diagonal contributions and explicit remainder terms.

Requirements:
1. No heuristic lower bounds unless you show where each inequality enters.
2. Isolate the exact term that was previously mis-lower-bounded by `R_shift - sum_err`.
3. Provide one corrected decomposition suitable for theorem use.
4. Clearly label identity vs inequality steps.

## Output format

1. `Exact identity (final algebraic form)`
2. `Where old bookkeeping failed`
3. `Corrected decomposition pieces`
4. `Single viable inequality target`
5. `Minimal falsification checklist`
