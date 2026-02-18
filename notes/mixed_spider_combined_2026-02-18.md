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

Artifacts:
- `results/whnc_mixed_spider_combined_k3000_j40.json`
- `results/whnc_mixed_spider_combined_k3000_j40_v2.json`

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

## New: analytic positivity templates for reduced branches

Script:

- `prove_mixed_spider_j0_j1_branches.py`

This script checks exact finite bases and validates the closed-form inequalities
used below.

### Branch 1: `j=0` (`T = S(2^k)`)

Let `m=(2k+1)//3`, `lambda=lambda_m(T)`, and write

- `margin = mu(T,lambda)-(m-1)`,
- `margin = (A + rB)/(1+r)`,
- `A = 2k*lambda/(1+2lambda) - (m-1)`,
- `B = 1 + k*lambda/(1+lambda) - (m-1)`,
- `r = lambda * ((1+lambda)/(1+2lambda))^k`.

Using coefficient split `i_t = U_t + V_t` with
`U_t = C(k,t)2^t`, `V_t = C(k,t-1)`, mediant gives
`lambda >= min(U_{m-1}/U_m, V_{m-1}/V_m)`.
For `k>=5`, `U_{m-1}/U_m <= V_{m-1}/V_m`, so

`lambda >= m/(2(k-m+1))`.

Hence

`A >= 1 - m/(k+1) >= (k+2)/(3(k+1))`.

Also:

- `B >= 2-m`,
- `r <= (2/3)^k`.

So, for `k>=5`:

`(1+r)margin = A + rB >= (k+2)/(3(k+1)) - (m-2)(2/3)^k`
`>= (1/3)*((k+2)/(k+1) - (2k-5)(2/3)^k) > 0`.

(`(2k-5)(2/3)^k` decreases for `k>=5` and is `<1` at `k=5`.)

Finite exact checks for `k=2,3,4` are positive.

### Branch 2: `j=1` (`T = S(2^k,1)`)

Let `m=(2k+3)//3`, same decomposition
`margin=(A+rB)/(1+r)` with

- `A = 2k*lambda/(1+2lambda) + lambda/(1+lambda) - (m-1)`,
- `B = 1 + k*lambda/(1+lambda) - (m-1)`,
- `r = lambda * ((1+lambda)/(1+2lambda))^(k-1)`.

Split coefficients as `i_t = F_t + Q_t`, where
`F_t` comes from `(1+2x)^k(1+x)` and `Q_t=C(k,t-1)`.
Mediant gives `lambda >= min(F_{m-1}/F_m, Q_{m-1}/Q_m)`.
For `k>=3`, `F_{m-1}/F_m <= Q_{m-1}/Q_m`, so

`lambda >= r_F := F_{m-1}/F_m`.

Since `A` increases in `lambda`,
`A >= A(r_F)`.
And `A(r_F)-1/3` has explicit positive residue formulas:

- `k=3t`:
  `(68 t^3 + 102 t^2 + 44 t + 6) / (3*(192 t^4 + 424 t^3 + 330 t^2 + 106 t + 12))`
- `k=3t+1`:
  `(196 t^3 + 582 t^2 + 528 t + 152) / (3*(192 t^4 + 776 t^3 + 1134 t^2 + 708 t + 160))`
- `k=3t+2`:
  `(132 t^3 + 492 t^2 + 588 t + 228) / (3*(192 t^4 + 984 t^3 + 1860 t^2 + 1536 t + 468))`

So `A(r_F) > 1/3`.

Also:

- `B >= 2-m`,
- `r <= (2/3)^(k-1)`.

Thus for `k>=7`:

`(1+r)margin = A + rB >= A(r_F) - (m-2)(2/3)^(k-1)`
`> 1/3 - ((2k-3)/3)(2/3)^(k-1) > 0`.

(`((2k-3)/3)(2/3)^(k-1)` decreases for `k>=3` and at `k=7` is `704/2187 < 1/3`.)

Finite exact checks for `k=1..6` are positive.

### Status of this proof attempt

- Positivity is now established for the reduced branches `j=0` and `j=1`
  under the above templates.
- Remaining gap is still the minimizer-reduction theorem
  (`j in {0,1}` extremality by residue), which is strongly supported by
  scans but not yet proved.

---

## Unit-leaf `c2 >= 0`: algebraic proof (2026-02-18)

See full write-up: `notes/unit_leaf_c2_algebra_2026-02-18.md`

For the unit-leaf decomposition `I_T = (1+x)I_B + x(1+x)^k` with
`B = (1+2x)^k(1+x)^{j-1}`:

**Step 1 (lower bound on `lambda_m^T`)**: Using log-concavity of `B` (product of
linear factors) and the elementary-symmetric bound `e_r/e_{r-1} >= (n-r+1)/r`:

`lambda_m^T >= tau = b_{m-2}/b_{m-1}`.

Verified: `k<=3000`, `j<=40`, `1,380,000` pairs, 0 violations.
Minimum gap: `lambda_m^T - tau = 0.000738` at `(k=2999, j=40, m=2019)`.

**Step 2 (mean at tau)**: Since `mu_B(lambda)` is increasing, `mu_B(lambda_m^T) >= mu_B(tau)`.
At `lambda=tau`, the distribution `w_t = b_t tau^t` is log-concave with mode at `m-1`.
By **Darroch's mode-mean inequality** (`|mode - mean| < 1`):
`mu_B(tau) > (m-1)-1 = m-2`.

Therefore `c2 = mu_B(lambda_m^T) - (m-2) > 0`. QED.

**Corollary (alternative formula)**: The algebraic identity
`c2 = 1/(1+lambda) + (1+r)*margin - r*(k*lambda/(1+lambda) - m + 2)`
(where `r = lambda*(1+lambda)^{k-j}/(1+2lambda)^k`) is verified exactly.
Given `margin >= 1/3`, this gives `c2 >= 7/18 > 0` (using `k*(2/3)^k <= 8/9` for all `k>=1`).

---

## Minimizer-reduction: structural analysis (2026-02-18)

New verifier script:
- `verify_mixed_spider_minimizer_reduction.py`

Artifacts:
- `results/whnc_mixed_spider_minreduction_k1_5_j80.json`
- `results/whnc_mixed_spider_minreduction_k6_3000_j80.json`
- `results/whnc_mixed_spider_minreduction_k3001_5000_j80.json`
- `results/whnc_mixed_spider_minreduction_k5001_8000_j80.json`
- `results/whnc_mixed_spider_minreduction_k1_8000_j80_aggregate.json`

Combined scan size: `648,000` pairs (`k=1..8000`, `j=0..80`).

### Verified structural claims (current strongest form)

1. **Parity-tail monotonicity (A-even/A-odd)**  
   For `k >= 6`:
   - `margin(k, j+2) >= margin(k, j)` for all even `j >= 4`,
   - `margin(k, j+2) >= margin(k, j)` for all odd `j >= 5`.
   No failures in `k=6..8000`, `j<=80`.

2. **Min over `j>=2` is at `j in {2,3}` (A')**  
   No failures for `k=4..8000`, `j<=80`.  
   Only exceptions are tiny: `k=1,2,3`.

3. **Residue comparison `j=0` vs `j=1`**
   - `k ≡ 1 (mod 3)`: `margin(k,0) <= margin(k,1)` for all checked `k=1..8000` (strict for `k>=4`).  
     Sub-claim B is proved algebraically; see `notes/subclaim_B_mod1_algebra_2026-02-18.md`.
   - `k ≡ 0,2 (mod 3)`: `margin(k,1) <= margin(k,0)` for all checked `k=3..8000`; only exception is `k=2`.  
     Sub-claim C is now proved; see `notes/subclaim_c_complete_proof_2026-02-18.md`.

4. **Global minimizer over scanned j-range**
   - For `k>=2`, `argmin_{0<=j<=80} margin(k,j)` is always in `{0,1}`.
   - Only exception in `k=1..8000` is `k=1` (minimum at `j=6`).

### Consequence for reduction strategy

Given the proved reduced-branch positivity (`j=0`, `j=1`) and the above structure,
the remaining analytic gap is now sharply concentrated in one place:

1. prove the parity-tail monotonicity claim for `k>=6`, and

once this is closed, minimizer reduction to `j in {0,1}` follows (with the tiny
`k=1,2` cases handled explicitly).

## Sub-claim C Closed (2026-02-18)

Sub-claim C (`k ≡ 0,2`: `margin(k,1) < margin(k,0)`) is now proved.

Proof package:
- `notes/subclaim_c_E1_le0_proof_2026-02-18.md` (proved `E_1 <= 0`),
- `notes/subclaim_c_complete_proof_2026-02-18.md` (new complete closure),
- `prove_subclaim_c_complete.py`,
- `results/whnc_subclaim_c_complete_proof_check.txt`.
