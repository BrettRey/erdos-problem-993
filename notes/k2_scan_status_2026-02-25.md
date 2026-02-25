# K2 Scan Status (2026-02-25)

## Scripts added

- `scripts/canonical_k2_split_scan.py`
  - Exact canonical scan for K2 collisions/splits:
    - key: `(d,m,lambda,mu1,mu2)`
    - gate: `is_dleaf_le_1` + `bridge_decomposition(require_dleaf=True)`
    - split criterion: same key, different `i1` / `N=[x]P`.

- `scripts/canonical_same_p_diff_lambda_scan.py`
  - Finds canonical-valid pairs with same `P=dp0[u]` but different derived `lambda`.
  - Default enforces same `m`.

## Reproducible runs

1) K2 split scan via script

```bash
python3 scripts/canonical_k2_split_scan.py \
  --min-n 4 --max-n 22 \
  --out results/k2_split_scan_n22_via_script.json
```

Outcome:
- `checked_total=403400`
- `unique_keys=370408`
- `collisions=32992`
- `split_found=false`

This matches the previously committed exact artifact:
- `results/k2_split_scan_n22_exact.json`.

2) Same-P / different-lambda witness scan

```bash
python3 scripts/canonical_same_p_diff_lambda_scan.py \
  --min-n 9 --max-n 9 \
  --out results/canonical_same_p_diff_lambda_witness_n9.json
```

Outcome:
- `witness_found=true` (same `m`, same `P`, different canonical-derived `lambda`).
- Pair:
  - `H?AA@bg`: `P=[1,6,10,5]`, `Q=[0,1,5,8,4]`, `I=[1,9,28,38,22,4]`, `m=3`, `lambda=14/19`
  - `H?ABAag`: `P=[1,6,10,5]`, `Q=[0,1,5,6,2]`, `I=[1,9,28,36,18,2]`, `m=3`, `lambda=7/9`

3) K2+X scan (X = Q'(lambda_hat))

```bash
python3 - <<'PY'
# ad-hoc exact scan (canonical gate) for key:
# (d,m,lambda,mu1,mu2,X), X=Q'(lambda_hat)
# output: results/k2_plus_X_split_scan_n22_exact.json
PY
```

Outcome (`results/k2_plus_X_split_scan_n22_exact.json`):
- `checked_total=403400`
- `unique_keys=380446`
- `collisions=22954`
- `split_found=false`

Comparison vs K2-only through `n<=22`:
- K2 collisions: `32992`
- K2+X collisions: `22954`
- no observed split in either scan through this bound.

## Interpretation

- Empirical status through `n<=22` remains:
  - many K2 collisions, no observed K2 split in `i1/N`.
- Canonical bridge-2 obstruction is explicit:
  - preserving `P` does not preserve canonical-derived `lambda`
  - therefore any canonical lift from fixed-`lambda` no-go must include `Q`-sensitive control.
- Empirically, `X=Q'(lambda_hat)` is a plausible next key component:
  - it reduces collision count substantially at fixed scan bound,
  - but does not yet constitute a proof of injectivity (`K2+X -> i1/N` remains open).
