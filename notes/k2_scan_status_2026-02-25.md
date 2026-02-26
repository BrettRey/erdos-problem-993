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

12) Targeted liftability search in a structured family

Added script:
- `scripts/star_component_kstar_scan.py`
  - Builds trees of the form `leaf-support-u` plus a multiset of rooted components
    attached to `u`.
  - Applies canonical gate via `bridge_decomposition(..., require_dleaf=True)`.
  - Scans for same
    `K* = (d,m,lambda_hat,mu1,mu2,rho,sigma)`
    with different `i1/N`.

Runs:

```bash
python3 scripts/star_component_kstar_scan.py \
  --min-comp-n 2 --max-comp-n 5 --multiset-size 4 \
  --out results/star_component_kstar_scan_n5_k4.json
```
- tested `3876`, passed_gate `330`, unique `330`, collisions `0`, `split_found=false`

```bash
python3 scripts/star_component_kstar_scan.py \
  --min-comp-n 2 --max-comp-n 6 --multiset-size 5 \
  --out results/star_component_kstar_scan_n6_k5.json
```
- tested `658008`, passed_gate `15504`, unique `13559`, collisions `1945`, `split_found=false`

```bash
python3 scripts/star_component_kstar_scan.py \
  --min-comp-n 2 --max-comp-n 7 --multiset-size 5 \
  --max-tested 1200000 \
  --out results/star_component_kstar_scan_n7_k5_cap1p2m.json
```
- tested `1200001`, passed_gate `35505`, unique `32305`, collisions `3200`, `split_found=false`

Interpretation: in this explicit constructive family, no lifted same-`K*`/different-`N`
pair has been found at current bounds.

Randomized coverage (same family):

```bash
python3 scripts/star_component_kstar_scan.py \
  --min-comp-n 2 --max-comp-n 7 --multiset-size 6 \
  --random-samples 500000 --seed 1 \
  --out results/star_component_kstar_scan_n7_k6_rand500k_s1.json
```
- tested `500000`, passed_gate `1826`, unique `1825`, collisions `1`, `split_found=false`

```bash
python3 scripts/star_component_kstar_scan.py \
  --min-comp-n 2 --max-comp-n 7 --multiset-size 6 \
  --random-samples 300000 --seed 2 \
  --out results/star_component_kstar_scan_n7_k6_rand300k_s2.json
```
- tested `300000`, passed_gate `1098`, unique `1098`, collisions `0`, `split_found=false`

```bash
python3 scripts/star_component_kstar_scan.py \
  --min-comp-n 2 --max-comp-n 7 --multiset-size 6 \
  --random-samples 300000 --seed 3 \
  --out results/star_component_kstar_scan_n7_k6_rand300k_s3.json
```
- tested `300000`, passed_gate `1115`, unique `1115`, collisions `0`, `split_found=false`

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
- Structured constructive liftability scans (star-of-components family) remain
  split-free at current tested bounds (both lexicographic and randomized samples).

Additional randomized liftability coverage:

```bash
python3 scripts/star_component_kstar_scan.py \
  --min-comp-n 2 --max-comp-n 8 --multiset-size 7 \
  --random-samples 150000 --seed 11 \
  --out results/star_component_kstar_scan_n8_k7_rand150k_s11.json
```
- tested `150000`, passed_gate `67`, unique `67`, collisions `0`, `split_found=false`
- elapsed `126.68s`

Interpretation: still no same-`K*`/different-`N` split in this broader randomized family probe.

Additional randomized liftability coverage (new batch):

```bash
python3 scripts/star_component_kstar_scan.py \
  --min-comp-n 2 --max-comp-n 8 --multiset-size 7 \
  --random-samples 150000 --seed 12 \
  --out results/star_component_kstar_scan_n8_k7_rand150k_s12.json
```
- tested `150000`, passed_gate `68`, unique `68`, collisions `0`, `split_found=false`
- elapsed `108.16s`

```bash
python3 scripts/star_component_kstar_scan.py \
  --min-comp-n 2 --max-comp-n 8 --multiset-size 7 \
  --random-samples 150000 --seed 13 \
  --out results/star_component_kstar_scan_n8_k7_rand150k_s13.json
```
- tested `150000`, passed_gate `82`, unique `82`, collisions `0`, `split_found=false`
- elapsed `99.98s`

```bash
python3 scripts/star_component_kstar_scan.py \
  --min-comp-n 2 --max-comp-n 8 --multiset-size 7 \
  --random-samples 150000 --seed 14 \
  --out results/star_component_kstar_scan_n8_k7_rand150k_s14.json
```
- tested `150000`, passed_gate `88`, unique `88`, collisions `0`, `split_found=false`
- elapsed `100.98s`

```bash
python3 scripts/star_component_kstar_scan.py \
  --min-comp-n 2 --max-comp-n 9 --multiset-size 7 \
  --random-samples 120000 --seed 21 \
  --out results/star_component_kstar_scan_n9_k7_rand120k_s21.json
```
- tested `120000`, passed_gate `27`, unique `27`, collisions `0`, `split_found=false`
- elapsed `100.15s`

Interpretation: no lifted same-`K*`/different-`N` split found in this additional randomized coverage.

Additional randomized liftability coverage (`max-comp-n=9`, `k=7`):

```bash
python3 scripts/star_component_kstar_scan.py \
  --min-comp-n 2 --max-comp-n 9 --multiset-size 7 \
  --random-samples 120000 --seed 22 \
  --out results/star_component_kstar_scan_n9_k7_rand120k_s22.json
```
- tested `120000`, passed_gate `36`, unique `36`, collisions `0`, `split_found=false`
- elapsed `102.65s`

```bash
python3 scripts/star_component_kstar_scan.py \
  --min-comp-n 2 --max-comp-n 9 --multiset-size 7 \
  --random-samples 120000 --seed 23 \
  --out results/star_component_kstar_scan_n9_k7_rand120k_s23.json
```
- tested `120000`, passed_gate `20`, unique `20`, collisions `0`, `split_found=false`
- elapsed `108.09s`

```bash
python3 scripts/star_component_kstar_scan.py \
  --min-comp-n 2 --max-comp-n 9 --multiset-size 7 \
  --random-samples 120000 --seed 24 \
  --out results/star_component_kstar_scan_n9_k7_rand120k_s24.json
```
- tested `120000`, passed_gate `25`, unique `25`, collisions `0`, `split_found=false`
- elapsed `106.63s`

Aggregate randomized `star_component_kstar_scan` coverage so far:
- result files: `11`
- tested multisets: `2,180,000`
- passed canonical gate: `4,452`
- unique `K*` keys: `4,451`
- collisions: `1`
- observed same-`K*`/different-`N` splits: `0`

## Gate-aware component filter added

Script update (`scripts/star_component_kstar_scan.py`):
- Added `--require-component-dleaf` option.
- When enabled, rooted component library keeps only components that satisfy
  `d_leaf<=1` after attaching a parent to the chosen root.

This improves canonical-gate hit rate in constructive scans.

### A/B benchmark (same seed)

```bash
python3 scripts/star_component_kstar_scan.py \
  --min-comp-n 2 --max-comp-n 9 --multiset-size 7 \
  --random-samples 50000 --seed 41 \
  --out results/star_component_kstar_scan_n9_k7_rand50k_s41_baseline.json
```
- tested `50000`, passed_gate `6`, unique `6`, collisions `0`, split `false`

```bash
python3 scripts/star_component_kstar_scan.py \
  --min-comp-n 2 --max-comp-n 9 --multiset-size 7 \
  --require-component-dleaf \
  --random-samples 50000 --seed 41 \
  --out results/star_component_kstar_scan_n9_k7_rand50k_s41_compdleaf.json
```
- tested `50000`, passed_gate `50000`, unique `50000`, collisions `0`, split `false`

### Additional scans with new option

```bash
python3 scripts/star_component_kstar_scan.py \
  --min-comp-n 2 --max-comp-n 6 --multiset-size 7 \
  --require-component-dleaf \
  --out results/star_component_kstar_scan_n6_k7_compdleaf_full.json
```
- tested `11440`, passed_gate `11440`, unique `11440`, collisions `0`, split `false`

```bash
python3 scripts/star_component_kstar_scan.py \
  --min-comp-n 2 --max-comp-n 7 --multiset-size 7 \
  --require-component-dleaf \
  --max-tested 350000 --progress-every 50000 \
  --out results/star_component_kstar_scan_n7_k7_compdleaf_cap350k.json
```
- tested `350001`, passed_gate `350000`, unique `333341`, collisions `16659`, split `false`

### Control check (root-min-deg=2 only, no component filter)

```bash
python3 scripts/star_component_kstar_scan.py \
  --min-comp-n 2 --max-comp-n 9 --root-min-deg 2 --multiset-size 7 \
  --random-samples 120000 --seed 31 \
  --out results/star_component_kstar_scan_n9_k7_deg2_rand120k_s31.json
```
- tested `120000`, passed_gate `10`, unique `10`, collisions `0`, split `false`

Interpretation:
- `root-min-deg=2` alone does not improve hit rate.
- `--require-component-dleaf` is effective and enables high-collision liftability probes.
- Even with `16659` same-`K*` collisions in the capped exhaustive regime, no same-`K*`/different-`N` split observed.

## Structural fact verification (new exact artifact)

Added script:
- `scripts/verify_canonical_structural_facts.py`

Command:
```bash
python3 scripts/verify_canonical_structural_facts.py \
  --min-n 3 --max-n 22 \
  --out results/verify_canonical_structural_facts_n22.json
```

Outcome (`results/verify_canonical_structural_facts_n22.json`):
- checked trees: `9,114,283`
- `d_leaf<=1` trees: `403,400`
- canonical-gated decompositions: `403,400`

Verified facts (exhaustive through `n<=22`):
1. `lambda_hat = i_{m-1}/i_m` is always strictly `< 1`
   - `max_lambda = 13513/13514 ≈ 0.9999260027`
2. Number of admissible degree-2 bridge triplets is never exactly `1`
   - `one_triplet_count = 0`
3. In all `m=2` canonical cases, the identity holds:
   - `i2 = (i1-1)(i1-2)/2`
   - observed `m2_cases = 10`, `m2_identity_failures = 0`

Interpretation:
- `lambda_hat=1` component-level collisions are non-liftable by definition.
- Any lifted same-key split witness must lie in `m>=3` regime.

## Canonical K* split scan with explicit min-u tie-break (new)

Added script:
- `scripts/canonical_kstar_split_scan_minu.py`

Triplet rule implemented explicitly:
- among admissible `(leaf,support,u)` with `leaf` leaf and `deg(support)=2`, choose minimum by `(u,support,leaf)`.

Sanity against accepted blocked witnesses:
- reproduces the same canonical triplets/lambdas for
  - `KpCHA?@?G?_@` vs `KpE?GD??G?_@`
  - `N??????_A?E?b?OoDS?` vs `N??????_A?E?B_aOKc?`

Exhaustive run:
```bash
python3 scripts/canonical_kstar_split_scan_minu.py \
  --min-n 3 --max-n 22 \
  --out results/canonical_kstar_split_scan_minu_n22_exact.json
```

Outcome (`results/canonical_kstar_split_scan_minu_n22_exact.json`):
- checked gated trees: `403,400`
- unique `K*` keys: `380,682`
- key collisions: `22,718`
- same-`K*`/different-`N` split found: `false` through `n<=22`
- elapsed: `262.94s`

## Canonical K* split scan with explicit min-u tie-break extended to n<=23

Command:
```bash
python3 scripts/canonical_kstar_split_scan_minu.py \
  --min-n 3 --max-n 23 --progress-every 2000000 \
  --out results/canonical_kstar_split_scan_minu_n23_exact.json
```

Outcome (`results/canonical_kstar_split_scan_minu_n23_exact.json`):
- checked gated trees: `931,596`
- unique `K*` keys: `875,604`
- key collisions: `55,992`
- same-`K*`/different-`N` split found: `false` through `n<=23`
- elapsed: `708.66s`

This tightens the min-u tie-break frontier from `n<=22` to `n<=23` with no observed split.

## Triplet-invariance scan (tie-break independent subclass)

Added script:
- `scripts/canonical_kstar_triplet_invariance_scan.py`

Definition used by script:
- A tree is `triplet-invariant` if all admissible triplets `(leaf,support,u)` yield the same
  `K*=(d,m,lambda,mu1,mu2,rho,sigma)`.

### m>=3 targeted run through n<=22

```bash
python3 scripts/canonical_kstar_triplet_invariance_scan.py \
  --min-n 3 --max-n 22 --m-min 3 --progress-every 2000000 \
  --out results/canonical_kstar_triplet_invariance_kstar_mge3_n22_exact.json
```

Outcome:
- dleaf trees: `403,400`
- gated trees: `403,400`
- triplet-invariant trees: `659`
- checked (triplet-invariant and `m>=3`): `651`
- unique keys: `651`
- collisions: `0`
- split_found: `false`

### m>=3 targeted run through n<=23

```bash
python3 scripts/canonical_kstar_triplet_invariance_scan.py \
  --min-n 3 --max-n 23 --m-min 3 --progress-every 2000000 \
  --out results/canonical_kstar_triplet_invariance_kstar_mge3_n23_exact.json
```

Outcome:
- dleaf trees: `931,596`
- gated trees: `931,596`
- triplet-invariant trees: `783`
- checked (triplet-invariant and `m>=3`): `775`
- unique keys: `775`
- collisions: `0`
- split_found: `false`

Interpretation:
- On the tie-break-independent `triplet-invariant` subclass, there are currently no
  observed `K*` collisions at all in the `m>=3` regime through `n<=23`.

## m>=4 targeted min-u canonical K* scan (new)

Updated script:
- `scripts/canonical_kstar_split_scan_minu.py` now supports `--m-min`.

Command:
```bash
python3 scripts/canonical_kstar_split_scan_minu.py \
  --min-n 3 --max-n 23 --m-min 4 --progress-every 2000000 \
  --out results/canonical_kstar_split_scan_minu_mge4_n23_exact.json
```

Outcome (`results/canonical_kstar_split_scan_minu_mge4_n23_exact.json`):
- checked (`m>=4`): `931,514`
- skipped by mode gate (`m<4`): `82`
- unique `K*` keys: `875,523`
- key collisions: `55,991`
- same-`K*`/different-`N` split found: `false` through `n<=23`

Interpretation:
- Restricting to the only unresolved mode regime (`m>=4`) leaves the no-split result unchanged.

## Extended min-u canonical K* scan to n<=24 (m>=4)

Command:
```bash
python3 scripts/canonical_kstar_split_scan_minu.py \
  --min-n 3 --max-n 24 --m-min 4 --progress-every 5000000 \
  --out results/canonical_kstar_split_scan_minu_mge4_n24_exact.json
```

Outcome (`results/canonical_kstar_split_scan_minu_mge4_n24_exact.json`):
- checked (`m>=4`): `2,164,055`
- unique `K*` keys: `2,027,329`
- key collisions: `136,726`
- same-`K*`/different-`N` split found: `false` through `n<=24`
- elapsed: `2286.63s`

Per-layer highlight at `n=24`:
- checked: `1,232,541`
- unique keys cumulative: `2,027,329`
- collisions cumulative: `136,726`

## Triplet-invariant K* scan with m-histograms (m>=4 gate)

Updated script:
- `scripts/canonical_kstar_triplet_invariance_scan.py` now records
  - `m_hist_invariant`
  - `m_hist_checked`

Command:
```bash
python3 scripts/canonical_kstar_triplet_invariance_scan.py \
  --min-n 3 --max-n 23 --m-min 4 --progress-every 2000000 \
  --out results/canonical_kstar_triplet_invariance_kstar_mge4_n23_exact.json
```

Outcome (`results/canonical_kstar_triplet_invariance_kstar_mge4_n23_exact.json`):
- dleaf trees: `931,596`
- gated trees: `931,596`
- triplet-invariant trees: `783`
- checked (`m>=4`): `756`
- unique keys: `756`
- collisions: `0`
- split_found: `false`

Histogram on triplet-invariant trees (all m):
- `m_hist_invariant = {1:1, 2:7, 3:19, 4:50, 5:122, 6:271, 7:313}`

Histogram in checked slice (`m>=4`):
- `m_hist_checked = {4:50, 5:122, 6:271, 7:313}`

Interpretation:
- The collision-free result in the triplet-invariant subclass is not just a low-mode artifact;
  it remains collision-free on the substantive `m>=4` slice through `n<=23`.

## Projected key scan: (m, lambda) only (min-u, m>=4)

Added script:
- `scripts/canonical_mlambda_split_scan_minu.py`

Command:
```bash
python3 scripts/canonical_mlambda_split_scan_minu.py \
  --min-n 3 --max-n 24 --m-min 4 --progress-every 5000000 \
  --out results/canonical_mlambda_split_scan_minu_mge4_n24_exact.json
```

Outcome (`results/canonical_mlambda_split_scan_minu_mge4_n24_exact.json`):
- run terminated early at first split (`max_n reached in run: 17`)
- checked (`m>=4`): `6,791`
- skipped by mode gate (`m<4`): `82`
- unique `(m,lambda)` keys: `1,417`
- key collisions: `5,375`
- split_found: `true`

First split witness (same `(m,lambda)`, different `N`):
- A: `n=16`, `g6=O??????_A?C?E?@_WG@j?`, `N=13`, `(m,lambda)=(5,7/9)`
- B: `n=17`, `g6=P????A?OD?E?E?B??o?E?OO?`, `N=14`, `(m,lambda)=(5,7/9)`

Interpretation:
- `(m,lambda) -> N` is false on the gated class under min-u canonicalization.
- So any proof route must use additional anchored coordinates (e.g. `d, mu1, mu2, rho, sigma`),
  not mode-ratio data alone.

## One-pass projection battery (single enumeration) through n<=24, m>=4

Added script:
- `scripts/canonical_projection_battery_minu.py`

Command:
```bash
python3 scripts/canonical_projection_battery_minu.py \
  --min-n 3 --max-n 24 --m-min 4 \
  --out results/canonical_projection_battery_minu_mge4_n24_exact.json
```

Full-scan stage:
- checked: `2,164,055`
- unique full `K*` keys: `2,027,329`
- full-key collisions: `136,726`

Projection results:
- `(m,lambda,d)`: `split_found=true`
  - confirms expected split on the known `(m,lambda)=(5,7/9)` witness pair with shared `d=7`.
- `(m,lambda,mu1)`: `split_found=false`
  - `unique=1,993,533`, `collisions=170,522`
- `(m,lambda,rho)`: `split_found=false`
  - `unique=2,027,329`, `collisions=136,726`
- `(m,lambda,sigma)`: `split_found=false`
  - `unique=2,027,329`, `collisions=136,726`
- `(m,lambda,rho,sigma)`: `split_found=false`
  - `unique=2,027,329`, `collisions=136,726`
- `(d,m,lambda,mu1,mu2)`: `split_found=false`
  - `unique=1,993,533`, `collisions=170,522`
- `(d,m,lambda,rho,sigma)`: `split_found=false`
  - `unique=2,027,329`, `collisions=136,726`

Notes:
- Through `n<=24`, any projection including `(m,lambda)` plus either `rho` or `sigma`
  already matches full-`K*` key cardinality in this regime.
- The key `(m,lambda,d)` fails quickly, while all tested Q-sensitive augmentations above remain split-free.
