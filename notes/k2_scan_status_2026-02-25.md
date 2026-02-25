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

4) K2 + (Q'(lambda_hat), Q''(lambda_hat)) scan

```bash
python3 - <<'PY'
# ad-hoc exact scan for key:
# (d,m,lambda,mu1,mu2,Q'(lambda_hat),Q''(lambda_hat))
# output: results/k2_plus_Qjet12_split_scan_n22_exact.json
PY
```

Outcome (`results/k2_plus_Qjet12_split_scan_n22_exact.json`):
- `checked_total=403400`
- `unique_keys=380446`
- `collisions=22954`
- `split_found=false`

Notably, through `n<=22` this is identical to the `K2+Q'` scan
(same unique key count and same collision count), i.e. adding `Q''` gives
no extra observed discrimination at this bound.

5) Extended bound: K2 + Q' through `n<=23`

```bash
python3 scripts/canonical_k2_split_scan.py \
  --min-n 4 --max-n 23 \
  --q-jet-max-order 1 \
  --out results/k2_plus_qjet1_split_scan_n23_exact.json
```

Outcome (`results/k2_plus_qjet1_split_scan_n23_exact.json`):
- `checked_total=931596`
- `unique_keys=874807`
- `collisions=56789`
- `split_found=false`

Comparison at `n<=23`:
- K2-only collisions: `73516` (from earlier exact artifact)
- K2+Q' collisions: `56789`
- no observed split in either scan at this bound.

6) Extended bound: K2 + (Q',Q'') through `n<=23`

```bash
python3 scripts/canonical_k2_split_scan.py \
  --min-n 4 --max-n 23 \
  --q-jet-max-order 2 \
  --out results/k2_plus_qjet12_split_scan_n23_exact.json
```

Outcome (`results/k2_plus_qjet12_split_scan_n23_exact.json`):
- `checked_total=931596`
- `unique_keys=874807`
- `collisions=56789`
- `split_found=false`

This is exactly identical to the `n<=23` K2+Q' result; adding `Q''`
again gave no extra observed discrimination.

7) Extended bound: K2 + (Q',Q'') through `n<=24`

```bash
python3 scripts/canonical_k2_split_scan.py \
  --min-n 1 --max-n 24 \
  --q-jet-max-order 2 \
  --out results/k2_plus_qjet12_split_scan_n24_exact.json
```

Outcome (`results/k2_plus_qjet12_split_scan_n24_exact.json`):
- `checked_total=2164137`
- `unique_keys=2025984`
- `collisions=138153`
- `split_found=false`

No same-key/different-`i1` split was observed at this larger bound.

8) Scale-anchored key test: K2 + (rho,sigma) through `n<=23`

Key extension:
- `rho = Q(lambda_hat)/P(lambda_hat)`
- `sigma = lambda_hat*Q'(lambda_hat)/P(lambda_hat)`

```bash
python3 scripts/canonical_k2_split_scan.py \
  --min-n 1 --max-n 23 \
  --include-rho-sigma \
  --out results/k2_plus_rhosigma_split_scan_n23_exact.json
```

Outcome (`results/k2_plus_rhosigma_split_scan_n23_exact.json`):
- `checked_total=931596`
- `unique_keys=874807`
- `collisions=56789`
- `split_found=false`

This is exactly identical to the `n<=23` runs for:
- `K2 + Q'`
- `K2 + (Q',Q'')`

So, at this bound, adding `(rho,sigma)` did not provide extra observed separation
beyond the existing Q-sensitive extensions.

9) Scale-anchored key test: K2 + (rho,sigma) through `n<=24`

```bash
python3 scripts/canonical_k2_split_scan.py \
  --min-n 1 --max-n 24 \
  --include-rho-sigma \
  --out results/k2_plus_rhosigma_split_scan_n24_exact.json
```

Outcome (`results/k2_plus_rhosigma_split_scan_n24_exact.json`):
- `checked_total=2164137`
- `unique_keys=2025984`
- `collisions=138153`
- `split_found=false`

This is exactly identical to the `n<=24` run for `K2+(Q',Q'')`
(`results/k2_plus_qjet12_split_scan_n24_exact.json`), i.e. no additional
observed discrimination at this bound.

10) Component-level injectivity probe at fixed lambda

Added script:
- `scripts/component_aggregate_collision_scan.py`
  - Enumerates rooted component library (from rooted trees up to `--max-comp-n`)
  - Scans multiset aggregates at a fixed lambda for collisions on:
    - `d = sum deg(F_j)`
    - `mu1 = sum a_j`
    - `mu2 = sum b_j + 2*sum_{i<j} a_i a_j`
    - `rho = lambda * product r_j`
    - `sigma = rho * (1 + sum g_j)`
  - Reports split if same aggregate key but different `N = sum |V(C_j)|`.

Runs:

```bash
python3 scripts/component_aggregate_collision_scan.py \
  --lambda 14/19 --max-comp-n 8 --multiset-size 2 \
  --out results/component_collision_lambda14_19_n8_k2.json
```
- tested `19900`, unique `19577`, collisions `323`, `split_found=false`

```bash
python3 scripts/component_aggregate_collision_scan.py \
  --lambda 14/19 --max-comp-n 8 --multiset-size 3 \
  --out results/component_collision_lambda14_19_n8_k3.json
```
- tested `1333300`, unique `1264688`, collisions `68612`, `split_found=false`

```bash
python3 scripts/component_aggregate_collision_scan.py \
  --lambda 14/19 --max-comp-n 9 --multiset-size 3 \
  --out results/component_collision_lambda14_19_n9_k3.json
```
- tested `18779684`, unique `18087293`, collisions `692391`, `split_found=false`

```bash
python3 scripts/component_aggregate_collision_scan.py \
  --lambda 7/9 --max-comp-n 8 --multiset-size 3 \
  --out results/component_collision_lambda7_9_n8_k3.json
```
- tested `1333300`, unique `1264688`, collisions `68612`, `split_found=false`

```bash
python3 scripts/component_aggregate_collision_scan.py \
  --lambda 7/9 --max-comp-n 9 --multiset-size 3 \
  --out results/component_collision_lambda7_9_n9_k3.json
```
- tested `18779684`, unique `18087293`, collisions `692391`, `split_found=false`

```bash
python3 scripts/component_aggregate_collision_scan.py \
  --lambda 11/12 --max-comp-n 8 --multiset-size 3 \
  --out results/component_collision_lambda11_12_n8_k3.json
```
- tested `1333300`, unique `1264688`, collisions `68612`, `split_found=false`

```bash
python3 scripts/component_aggregate_collision_scan.py \
  --lambda 11/12 --max-comp-n 9 --multiset-size 3 \
  --out results/component_collision_lambda11_12_n9_k3.json
```
- tested `18779684`, unique `18087293`, collisions `692391`, `split_found=false`

Interpretation: on this rooted-component superset probe (fixed lambda, multiset size <=3,
component size <=9), no aggregate-key split in `N` was found.

11) Independent verification of 5.2 component-level BLOCKED artifacts (A and C)

Added verifier:
- `scripts/verify_component_blocked_artifacts.py`
  - Checks reported case A (`lambda=1`) and case C (`lambda=1/2`) exactly.
  - Verifies equality of `(d,mu1,mu2,rho,sigma)` and inequality of `N`.

Run:

```bash
python3 scripts/verify_component_blocked_artifacts.py \
  --case both \
  --out results/component_blocked_artifact_verify_A_C.json
```

Outcome:
- Case A: equal `(d,mu1,mu2,rho,sigma)` and `N(A)-N(B)=-16`
- Case C: equal `(d,mu1,mu2,rho,sigma)` and `N(A)-N(B)=-173070`

This confirms the A/C outputs are valid at the **component aggregate** level.
As noted by 5.2, these do not automatically lift to full canonical same-`K*`
tree pairs because mode/fugacity locking (`m, lambda_hat`) is still unresolved.

## Interpretation

- Empirical status through `n<=22` remains:
  - many K2 collisions, no observed K2 split in `i1/N`.
- Canonical bridge-2 obstruction is explicit:
  - preserving `P` does not preserve canonical-derived `lambda`
  - therefore any canonical lift from fixed-`lambda` no-go must include `Q`-sensitive control.
- Empirically, `X=Q'(lambda_hat)` is a plausible next key component:
  - it reduces collision count substantially at fixed scan bound,
  - but does not yet constitute a proof of injectivity (`K2+X -> i1/N` remains open).
- Empirically, adding `Q''(lambda_hat)` on top of `Q'(lambda_hat)` did not improve
  separation through `n<=22`.
- Extending to `n<=23`, `K2+Q'` remains split-free and continues to reduce
  collision count relative to K2-only.
- Through `n<=23`, `Q''` remains empirically redundant once `Q'` is included.
- Through `n<=24`, `K2+(Q',Q'')` remains split-free on the canonical gated scan.
- Through `n<=23`, `K2+(rho,sigma)` is also split-free and empirically
  indistinguishable (in key count/collision count) from `K2+Q'` and
  `K2+(Q',Q'')`.
- Through `n<=24`, `K2+(rho,sigma)` remains split-free and still matches
  `K2+(Q',Q'')` exactly in aggregate key/collision counts.
- Fixed-lambda component-level probes (size <=9, multiset size <=3) also found
  no `N`-split for the anchored aggregate tuple `(d,mu1,mu2,rho,sigma)` at
  `lambda in {14/19, 7/9, 11/12}`.
- Separately, two explicit large component-level collisions from 5.2 (A/C) are
  now independently verified; the open gap is liftability to full canonical
  trees with the same derived `(m,lambda_hat)`.
