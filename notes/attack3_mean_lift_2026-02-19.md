# Attack 3: Mean-Lift from `mu_P` to `mu_T` and Mode Consequences (2026-02-19)

## Goal

Analyze the bridge decomposition at mode tie fugacity and test whether controlling
`mu_T - mu_P` can force `mode(P) >= m-1`.

Setup (same notation as current Route C notes):

- `I(T) = (1+2x)P + (1+x)Q`
- `I(B) = P + Q`
- `m = mode(I(T))` (leftmost)
- `lambda = lambda_m(T) = i_{m-1}(T)/i_m(T)`
- `mu_F(lambda) = lambda F'(lambda)/F(lambda)`
- `Z_F = F(lambda)`
- `p_u = Z_Q/Z_B`, `D = mu_B - mu_P`

## 1. Exact lift formulas

From

`I_T'(x) = 2P + (1+2x)P' + Q + (1+x)Q'`,

at `x=lambda`:

`lambda I_T' = lambda(2Z_P+Z_Q) + (1+2lambda)mu_P Z_P + (1+lambda)mu_Q Z_Q`.

Hence

`mu_T - mu_P`
`= [lambda(2Z_P+Z_Q) + (1+lambda)Z_Q(mu_Q-mu_P)] / Z_T`,

where

`Z_T = (1+2lambda)Z_P + (1+lambda)Z_Q`.

Using `p_u = Z_Q/(Z_P+Z_Q)` and `D = mu_B-mu_P`:

`mu_T - mu_P = [lambda(2-p_u) + (1+lambda)D] / [1+2lambda-lambda p_u]`.

This gives the useful normalized form:

`(mu_T - mu_P) - 1 = ((1+lambda)D - 1) / (1+2lambda-lambda p_u)`.

So:

- `mu_T - mu_P < 1` iff `D < 1/(1+lambda)`.
- `mu_T - mu_P > 1` iff `D > 1/(1+lambda)`.

Also from `I(T)=I(A)+xI(B)` with `I(A)=(1+x)P+Q`:

`mu_T = (Z_A/Z_T)mu_A + (lambda Z_B/Z_T)(1+mu_B)`.

All three identities were checked numerically.

## 2. New verification scripts

- `attack3_mean_lift_scan_2026_02_19.py`
- `attack3_product_mode_mean_probe_2026_02_19.py`

## 3. Exhaustive scan results (`d_leaf<=1`, `n<=23`)

### 3.1 Canonical leaf (one decomposition per tree)

Command:

```bash
python3 attack3_mean_lift_scan_2026_02_19.py \
  --min-n 4 --max-n 23 \
  --leaf-mode canonical \
  --out results/attack3_mean_lift_scan_n23_canonical_2026_02_19.json \
  --no-resume
```

Summary:

- trees checked: `931,596`
- decompositions checked: `931,596`
- `max(mu_T-mu_P) = 1.0040727976730937`
- `min(1-(mu_T-mu_P)) = -0.004072797673093653`
- `mu_T-mu_P >= 1` failures of strict-`<1`: `4`
- `mode(P) < m-1` failures: `0`
- `mode(P) < floor(mu_P(lambda_m(T)))` failures: `0`
- identity errors:
  - max Z-form error: `6.217248937900877e-15`
  - max `(p_u,D)`-form error: `6.217248937900877e-15`
  - max A/B decomposition error: `5.329070518200751e-15`

Global max-lift witness:

- `n=23`, `g6=V?????????????_?G?@??C??G??G??E???o?oB?@|S??`
- `m=8`, `lambda=0.9972305185989683`
- `mu_T=7.434366273337189`
- `mu_P=6.4302934756640955`
- `mu_T-mu_P=1.0040727976730937`
- `D=0.5058006710219001`, threshold `1/(1+lambda)=0.5006933304331276`
- `D - 1/(1+lambda) = 0.005107340588772491` (so lift exceeds `1`)

### 3.2 All degree-2-support leaves (all valid decompositions)

Command:

```bash
python3 attack3_mean_lift_scan_2026_02_19.py \
  --min-n 4 --max-n 23 \
  --leaf-mode all_deg2 \
  --out results/attack3_mean_lift_scan_n23_all_deg2_2026_02_19.json \
  --no-resume
```

Summary:

- trees checked: `931,596`
- decompositions checked: `4,543,370`
- `max(mu_T-mu_P) = 1.0040727976730937` (same global witness)
- strict `mu_T-mu_P < 1` failures: `14`
- `mode(P) < m-1` failures: `0`
- `mode(P) < floor(mu_P(lambda_m(T)))` failures: `0`
- identity errors:
  - max Z-form error: `6.439293542825908e-15`
  - max `(p_u,D)`-form error: `6.5503158452884236e-15`
  - max A/B decomposition error: `5.329070518200751e-15`

So the strict `<1` hypothesis is false even in canonical mode, and remains false in all-leaf mode.

## 4. Mean-mode facts seen in the scan

- `mu_T(lambda_m)` always remained in `[m-1, m]` on this frontier:
  - `mu_T < m-1` failures: `0`
  - `mu_T > m` failures: `0`
- `mu_P >= m-1` is **not** true (fails almost always):
  - canonical failures: `931,591 / 931,596`
  - all-leaf failures: `4,543,328 / 4,543,370`
- Therefore, proving `mode(P) >= m-1` cannot proceed by trying to force `mu_P >= m-1`
  from this lift bound alone.

## 5. Step-3 product probe (`mode(P) >= floor(mu_P(lambda))`)

Command:

```bash
python3 attack3_product_mode_mean_probe_2026_02_19.py \
  --max-tree-n 10 --max-factors 4 --max-total-vertices 14 \
  --check-ties --tie-lambda-max 1.0 \
  --out results/attack3_product_mode_mean_probe_n10f4v14_2026_02_19.json
```

Sampled product class:

- unique factor polynomials: `186`
- product combos tested: `4,216`

Result:

- For `lambda in {0.25, 0.5, 0.75, 1.0}`: `0` failures.
- For tie fugacities of those products with `lambda_tie <= 1`: `18,394` checks, `0` failures.
- For `lambda > 1`, failures appear (`54` at `1.25`, `941` at `1.5`).

Interpretation:

- In the fugacity regime relevant to tie points here (`lambda <= 1`), this relation
  looks robust in tested product families.
- It is not a universal all-`lambda` statement.

## 6. Status vs converging target

What we now have:

1. Exact formulas for `mu_T-mu_P` with machine-precision verification on full frontier.
2. Exhaustive evidence that `mode(P) >= m-1` still holds (`0` failures) in both canonical and all-leaf scans.
3. Exhaustive disproof of the stronger shortcut `mu_T-mu_P < 1` (counterexamples exist).
4. Strong evidence (not proof) for `mode(P) >= floor(mu_P(lambda))` in relevant `lambda<=1` product regime.

So the promising route remains:

- prove `mode(P) >= floor(mu_P(lambda_m(T)))` (or directly `mode(P) >= m-1`)
  using structural properties of this bridge/product class,
- not via a global strict `<1` lift bound.
