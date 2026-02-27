# Adjacent n<=27 Consistency Check (2026-02-28)

## Scope
Post-reconciliation validation of:
- `results/adjacent_rho_split_scan_minu_mge4_n27_exact.json`
- `results/adjacent_rho_split_scan_minu_mge4_n27_odd_even_exact.json`

## Commands
- `python3 test_all.py`
- consistency script comparing full vs odd-even derived artifact fields.

## Results
- Test suite: `37` tests passed.
- Full adjacent projection:
  - `split_found=false`
  - `adjacent_split_found=false`
- Odd-even-only projection:
  - `split_found=false`
  - `adjacent_split_found=false`
- Field equality checks (full vs odd-even artifact):
  - `checked_total`: match
  - `full_unique_keys`: match
  - `full_collisions`: match

## Conclusion
Artifacts are internally consistent and support:
- no adjacent same-`(m,lambda,rho)` split through `n<=27`;
- therefore no odd->even adjacent split through `n<=27`.
