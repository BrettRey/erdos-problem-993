# Route-1 5.2 Status: Fixed-lambda vs Canonical (2026-02-25)

## Current status

- Proven: fixed-`lambda` finite-`r` obstruction in the unconstrained product model.
  - For any fixed finite `r`, `(d, lambda, mu_1, ..., mu_r)` does not determine `N=[x]P`.
  - See `notes/fixed_lambda_finite_r_obstruction_2026-02-25.md`.
- Not yet proven: canonical no-go or canonical injectivity for `K2 -> N`.
  - Canonical coupling through `I(T;x)=(1+2x)P(x)+(1+x)Q(x)` remains the unresolved bridge.
- Canonical witness now recorded that `P` does not determine canonical-derived `lambda`:
  - see `notes/canonical_sameP_different_lambda_pair_2026-02-25.md`.

## What is accepted

- `i_k = p_k + 2p_{k-1} + q_k + q_{k-1}` (with `p_{-1}=q_{-1}=0`).
- `p_0=1`, `q_0=0`, `q_1=1`.
- Therefore `i_1 = p_1 + 3` and `N = p_1 = i_1 - 3`.

## Exact bottleneck

- To conclude `K2=(d,m,lambda,mu1,mu2) -> N`, one must prove `K2 -> i_1`.
- That step is currently missing.

## Next 5.2 gate (strict)

Use this as the next instruction:

```text
Accepted: fixed-lambda unconstrained no-go is proved.

Now do exactly one:

(A) Canonical lift theorem:
Prove both bridges rigorously:
1) canonical-admissible infinite motif family replacement for stars,
2) lambda-coupled preservation under swap:
   keep canonical derived (m,lambda) from I(T)=(1+2x)P+(1+x)Q unchanged,
   with all canonical gates valid.

or

(B) Canonical blocked witness:
Give an explicit canonical pair (T1,T2) with gate outputs,
same target key (K2 or K_r), different i1/N,
including exact g6, P,Q,I(T), m, lambda, and mu-data.

No genericity language. Exact equations + gate-verified data only.
If neither is completed, output INSUFFICIENT.
```
