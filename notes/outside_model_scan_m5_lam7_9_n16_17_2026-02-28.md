# Outside-Model Scan for `(m,lambda)=(5,7/9)` on `n=16..17`

## Imported files
- Script: `scripts/outside_model_scan_m5_lam7_9.py`
- Artifact: `results/outside_model_scan_m5_lam7_9_n16_17.json`

## Reported command (from C)
- `python3 /mnt/data/outside_model_scan_m5_lam7_9.py --min-n 16 --max-n 17 --m-target 5 --lambda-target 7/9 --out /mnt/data/outside_model_scan_m5_lam7_9_n16_17.json`

## Decisive verdict
- `outside-model witness not found in tested bounds`

## Exact artifact contents
- Params:
  - `m_target=5`
  - `lambda_target=7/9`
  - `min_n=16`, `max_n=17`
- `n=16`:
  - `total_trees=19320`
  - `gated_ok=1692`
  - `lock_ok=1`
  - match: `N=13`, `isomorphic_to_A=true`, `isomorphic_to_B=false`
- `n=17`:
  - `total_trees=48629`
  - `gated_ok=3726`
  - `lock_ok=1`
  - match: `N=14`, `isomorphic_to_A=false`, `isomorphic_to_B=true`
- cross-N rho matches:
  - `rho_cross_matches=[]`
