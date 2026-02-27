# C1/C2 Coverage Check for `(m,lambda,rho)` (n<=24, m>=4)

Command:

```bash
python3 scripts/canonical_projection_battery_minu.py \
  --min-n 3 --max-n 24 --m-min 4 \
  --projections 'm,lambda,rho' \
  --analyze-c1c2 \
  --progress-every 0 \
  --out results/c1c2_coverage_mlambda_rho_mge4_n24_exact.json
```

Artifact:
- `results/c1c2_coverage_mlambda_rho_mge4_n24_exact.json`

Scan totals:
- `checked_total = 2,164,055`
- `full_unique_keys = 2,027,329`
- `full_collisions = 136,726`

Projection `(m,lambda,rho)`:
- `unique_keys = 2,027,329`
- `collisions = 136,726`
- `split_found = false`

C1/C2 criterion coverage:
- `keys_total = 2,027,329`
- `c1_true = 0`
- `c2_true = 0`
- `c1c2_true = 0`
- `c1c2_and_NeqL = 0`
- `c1c2_and_NneL = 0`

Interpretation:
- The current sufficient singleton criterion `(C1,C2)` is too strict to certify any
  observed key at this bound.
- This does **not** contradict no-split scans; it only shows the criterion has zero
  practical coverage on this dataset.
