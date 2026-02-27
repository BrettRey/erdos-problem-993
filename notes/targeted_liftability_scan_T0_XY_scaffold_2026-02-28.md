# Targeted Liftability Scan: `T0 + X/Y kernel + E/P3 scaffold` (2026-02-28)

## What was validated
- Downloaded artifact:
  - `/Users/brettreynolds/Downloads/targeted_liftability_scan_T0_XY_scaffold_t20_e20_p3_20_mge4.json`
- Local reproducible script added:
  - `scripts/targeted_liftability_scan_t0_xy_scaffold.py`
- Local rerun artifact:
  - `results/targeted_liftability_scan_t0_xy_scaffold_t20_e20_p3_20_mge4_local.json`

## Command (local reproducible rerun)
- `nice -n 15 python3 scripts/targeted_liftability_scan_t0_xy_scaffold.py --max-t 20 --max-e 20 --max-p3 20 --m-min 4 --out results/targeted_liftability_scan_t0_xy_scaffold_t20_e20_p3_20_mge4_local.json`

## Exact totals (downloaded vs local rerun)
- `combos_scanned = 101871` (both)
- `combos_passing_m_gate = 101871` (both)
- `collision_count = 0` (both)
- `split_found = false` (both)
- Runtime differs by machine/load only (`~30.0s` downloaded vs `~21.3s` local).

## Interpretation
- Within this scaffolded family and bounds (`t<=20`, `cE<=20`, `dP3<=20`, `m>=4`):
  - no same-`(m,lambda)` collisions across distinct `b` at fixed `(t,cE,dP3)`,
  - so no lifted same-`(m,lambda,rho)` / different-`N` split was found.
- This is a bounded negative result (family-specific), not a global proof for the full canonical gated class.
