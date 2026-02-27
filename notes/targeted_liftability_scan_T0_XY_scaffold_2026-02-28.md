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

## Extended local bound
- Command:
  - `nice -n 18 python3 scripts/targeted_liftability_scan_t0_xy_scaffold.py --max-t 30 --max-e 30 --max-p3 30 --m-min 4 --out results/targeted_liftability_scan_t0_xy_scaffold_t30_e30_p3_30_mge4_local.json`
- Totals:
  - `combos_scanned = 476656`
  - `combos_passing_m_gate = 476656`
  - `collision_count = 0`
  - `split_found = false`
  - runtime `~221.34s`

## Interpretation
- Within this scaffolded family and bounds (`t<=20`, `cE<=20`, `dP3<=20`, `m>=4`):
  - no same-`(m,lambda)` collisions across distinct `b` at fixed `(t,cE,dP3)`,
  - so no lifted same-`(m,lambda,rho)` / different-`N` split was found.
- The extended run (`t,e,p3 <= 30`) preserves the same negative outcome.
- Further extended run (`t,e,p3 <= 40`) also preserves the same negative outcome:
  - artifact: `results/targeted_liftability_scan_t0_xy_scaffold_t40_e40_p3_40_mge4_local.json`
  - totals:
    - `combos_scanned = 1447341`
    - `combos_passing_m_gate = 1447341`
    - `collision_count = 0`
    - `split_found = false`
    - runtime `~1186.82s`
- This is a bounded negative result (family-specific), not a global proof for the full canonical gated class.
