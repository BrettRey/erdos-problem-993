# `(m,lambda,rho)` / `(m,lambda,sigma)` Exact Scan (n<=26, m>=4)

Command:

```bash
python3 scripts/canonical_projection_battery_minu.py \
  --min-n 3 --max-n 26 --m-min 4 \
  --projections 'm,lambda,rho;m,lambda,sigma' \
  --progress-every 0 \
  --out results/canonical_projection_battery_minu_mlambda_rho_sigma_mge4_n26_exact.json
```

Artifact:
- `results/canonical_projection_battery_minu_mlambda_rho_sigma_mge4_n26_exact.json`

Totals:
- `checked_total = 11,850,490`
- `full_unique_keys = 11,026,431`
- `full_collisions = 824,059`
- `skipped_m = 82`
- `elapsed_sec_scan = 16091.41`
- `elapsed_sec_total = 16139.65`

Projection outcomes:
- key `(m,lambda,rho)`:
  - `unique_keys = 11,026,431`
  - `collisions = 824,059`
  - `split_found = false`
- key `(m,lambda,sigma)`:
  - `unique_keys = 11,026,431`
  - `collisions = 824,059`
  - `split_found = false`

Conclusion:
- Through `n<=26` (with `m>=4`, min-`u` canonical tie-break), neither `(m,lambda,rho)` nor `(m,lambda,sigma)` exhibits a same-key/different-`N` split.
