# Family Miner Demo (B Output) Validation (2026-02-28)

## Artifact imported
- `results/family_miner_kernel_scaffold_demo_n10_m4_box42_f10_adjkernel.json`

Source command reported by B:
- `python3 scripts/family_miner_kernel_scaffold.py --max-comp-n 10 --m-min 4 --max-kernel-total 4 --max-scaffold 2 --max-families 10 --require-adjacent-kernel --out results/family_miner_kernel_scaffold_demo_n10_m4_box42_f10_adjkernel.json`

## Schema notes
Artifact top-level keys:
- `all`, `library`, `params`, `results`, `runtime_s`, `top`

Validated headline fields:
- `library.rooted_types = 324`
- `results.candidates = 10`
- `results.any_mlr_split_in_box = false`
- `runtime_s = 0.31648475499991946` (same order as reported)

Top candidate structure:
- shared kernel pair across all top-10:
  - `K1`: `n=8`, `F=[1,8,21,21,6]`, `G=[1,7,16,13,2]`
  - `K2`: `n=9`, `F=[1,9,28,36,18,3]`, `G=[1,8,22,24,9,1]`
  - `delta_n = 1`
- each `top[i].hits` confirms no collision witness:
  - `ml_collision = null`
  - `mlr_collision = null`
  - (and kernel-line variants null as well)

## Interpretation
- In this demo box (`a+b<=4`, scaffold counts `<=2`, top-10 candidates), there is no same-`(m,lambda,rho)` / different-`N` witness.
- Diagnostics on top candidates report `sign_mixed = false`, aligning with the lock-collapse pattern seen in other family-specific no-go results.

## Caveat
- This is a bounded demo over selected candidates (`max-families=10`), not an exhaustive miner over all rooted-type combinations.
