# Mixed Spider Combined Bound (2026-02-18)

Target family: `T = S(2^k,1^j)` with

`I_T(x) = (1+2x)^k(1+x)^j + x(1+x)^k`.

At `m = mode(I_T)` (leftmost mode at `lambda=1`) and
`lambda_m = i_{m-1}/i_m`, define

`margin = mu(T,lambda_m) - (m-1)`.

For a leaf decomposition `I_T = I_A + x I_B`:

- `c1 = mu(A,lambda_m) - (m-1)`
- `c2 = mu(B,lambda_m) - (m-2)`
- `w1 = I_A(lambda_m)/I_T(lambda_m)` (weight on `c1`)
- `w2 = lambda_m I_B(lambda_m)/I_T(lambda_m)` (weight on `c2`)
- `combined = w1*c1 + w2*c2`

Identity: `combined = margin` (checked numerically to ~`1e-13`).

## Closed-form pieces used

For any mixed spider `S(2^a,1^b)`:

- `I = (1+2x)^a(1+x)^b + x(1+x)^a`
- with `r = x(1+x)^(a-b)/(1+2x)^a`,
  `mu = [a*2x/(1+2x) + b*x/(1+x) + r*(1 + a*x/(1+x))]/(1+r)`.

Leaf types:

1. Tip leaf (from a length-2 arm):
- `A = S(2^(k-1),1^(j+1))`
- `B = S(2^(k-1),1^j)`

2. Unit leaf (adjacent to hub):
- `A = S(2^k,1^(j-1))`
- `B` has polynomial `(1+2x)^k(1+x)^(j-1)` (product graph after deleting hub + unit leaf).

## Script

- `conjecture_a_mixed_spider_combined_scan.py`

## Scan results

### Grid A: `k=1..400`, `j=0..400`

Artifact:
- `results/whnc_mixed_spider_combined_k400_j400.json`

Summary:
- `combined_fail = 0` (no `combined < 0` cases),
- tip-leaf:
  - `c1_neg_count = 0`,
  - `min combined = 0.33416458852869274` at `(k,j)=(400,0)`,
- unit-leaf:
  - `c1_neg_count = 157,802 / 160,000`,
  - `min c1 = -0.1653177570093476`,
  - `min combined = 0.33421626629649626` at `(k,j)=(399,1)`.

### Grid B: `k=1..3000`, `j=0..40`

Artifact:
- `results/whnc_mixed_spider_combined_k3000_j40.json`

Summary:
- `combined_fail = 0`,
- global minimum combined:
  - `0.33344448149394906` at `(k,j)=(2998,0)` (tip case),
- tip-leaf:
  - `c1_neg_count = 0`,
  - `min c1 = 0.1667639317299745`,
- unit-leaf:
  - `c1_neg_count = 119,782 / 120,000`,
  - `min c1 = -0.1664862649247425` at `(k,j)=(3000,1)`,
  - `min combined = 0.3334513053243715`.

## Structural takeaway

Empirically on this family:

1. Tip-leaf branch is easy (`c1` stayed positive in all scanned ranges).
2. Unit-leaf branch is the difficult regime (`c1 < 0` is common), but
   combined stays strictly positive with a large buffer (`~1/3` asymptotically).
3. The worst combined cases occur near `j=0`/`j=1` and large `k`, matching the
   mixed-spider extremal pattern seen elsewhere.

## Per-k minimizer structure (new, stronger)

We extended with per-`k` minimizer tracking:

- `results/whnc_mixed_spider_combined_k3000_j40_v2.json`
- `results/whnc_mixed_spider_combined_k8000_j6.json`

### `k<=3000`, `j<=40`

- `combined_fail=0`.
- For each residue class `k mod 3`, per-`k` minima are monotone decreasing from
  `k>=20`.
- Dominant minimizer signatures by residue:
  - `k ≡ 1 (mod 3)`: almost always `j=0`, tip leaf (`999/1000` in this run).
  - `k ≡ 0 (mod 3)`: almost always `j=1` (tip or unit).
  - `k ≡ 2 (mod 3)`: almost always `j=1` (tip or unit).

### `k<=8000`, `j<=6`

- `combined_fail=0`.
- Same residue pattern persists:
  - `k ≡ 1`: `j=0`, tip in `2666/2667` cases.
  - `k ≡ 0`: `j=1` (tip/unit split `1337/1329`).
  - `k ≡ 2`: `j=1` (tip/unit split `1348/1318`).
- No monotonicity violations for per-`k` minima from `k>=20` in any residue
  class.

This gives a compressed target:

1. prove the minimizer reduction to `j in {0,1}` by residue class, then
2. prove explicit positivity for the reduced branches:
   - `j=0` tip (balanced spider lane),
   - `j=1` tip/unit (mixed lane).
