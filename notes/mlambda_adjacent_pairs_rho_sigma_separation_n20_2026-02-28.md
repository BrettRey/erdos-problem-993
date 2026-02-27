# Adjacent `(m,lambda)` Pairs: `rho` / `sigma` Separation (`n<=20`)

## Inputs
- Adjacent same-`(m,lambda)` pairs from:
  - `results/mlambda_overlap_profile_minu_mge4_n20.json`
- Rebuilt canonical records via:
  - `scripts/canonical_kstar_split_scan_minu.py::rebuild_record`

## Command used
- See artifact generator output:
  - `results/mlambda_adjacent_pairs_rho_sigma_separation_n20.json`

## Summary
From `results/mlambda_adjacent_pairs_rho_sigma_separation_n20.json`:
- `pairs_total = 12`
- `rho_equal_count = 0`
- `sigma_equal_count = 0`

So for every currently observed adjacent odd/even overlap in `(m,lambda)` through `n<=20`, both `rho` and `sigma` already separate the pair.

## Interpretation
- This is consistent with the broader scans showing no `(m,lambda,rho)` or `(m,lambda,sigma)` splits in tested ranges.
- It does not prove global injectivity, but it narrows practical risk: known adjacent `(m,lambda)` overlaps in the small range are not problematic once `rho` (or `sigma`) is included.
