# Path 1 Status: Algebraic `Phi_{m-1}(P) >= 0` (2026-02-19)

This note tracks the Route-1 target from the degree-2 bridge decomposition:

- Choose leaf `l` with minimum support degree, support `s`, and `deg(s)=2`.
- Let `u` be the other neighbor of `s`.
- Let `B = T - {l,s}` and `P = dp_B[u][0]`.
- Let `m = mode(I(T))` and `lambda = lambda_m(T) = i_{m-1}(T)/i_m(T)`.

Target:

`Phi_{m-1}(P;lambda) = lambda P'(lambda) - (m-2)P(lambda) >= 0`

equivalently:

`mu_P(lambda) >= m-2`.

---

## 1) What was run

Scanner:

- `conjecture_a_phiP_scan.py`

Primary outputs:

- `results/whnc_phiP_scan_n23.json` (full run with partial `n=23`, retained for provenance)
- `results/whnc_phiP_scan_n23_mod8_r0.json` ... `results/whnc_phiP_scan_n23_mod8_r7.json` (full `n=23` in 8 disjoint chunks)
- merged:
  `results/whnc_phiP_scan_n23_merged.json`

The merged file uses `n<=22` from the full run and replaces `n=23` by exact 8-way chunk merge.

---

## 2) Full-frontier facts (d_leaf<=1, n<=23, 931,596 trees)

From `results/whnc_phiP_scan_n23_merged.json`:

- checked trees: `931,596`
- `Phi_{m-1}(P;lambda) < 0` failures: `0`
- `mode(P) < m-1` failures: `0`
- `P` not ULC failures: `0`

Distributions:

- `mode(P) - (m-1)`:
  - `0`: `340,021`
  - `1`: `63,682`

- `deg(P) - (m-1)`:
  - `0`: `1`
  - `1`: `17`
  - `2`: `556`
  - `3`: `14,689`
  - `4`: `241,696`
  - `5`: `140,237`
  - `6`: `6,421`
  - `7`: `86`

Extremal margins:

- min `mu_P(lambda) - (m-2)`:
  `0.38345129375473075`
- min `mu_P(tau_P) - (m-2)`, `tau_P = p_{m-2}/p_{m-1}`:
  `0.36287712223494495`

So Route-1 has substantial positive slack on the full frontier.

---

## 3) Tie-point comparison with `tau_P = p_{m-2}/p_{m-1}`

Observed:

- `lambda - tau_P` failures: `1`
- minimum `lambda - tau_P`:
  `-0.00013108452027321693`

Exact witness (fraction check):

- `lambda = 37595/37851`
- `tau_P = 1348/1357`
- `lambda - tau_P = -6733/51363807`

For that witness, Route-1 is still far from failing:

- `mu_P(lambda) - (m-2) = 0.45768794816449`
- `mu_P(tau_P) - (m-2) = 0.4580223988607841`

So the strict inequality `lambda >= tau_P` is almost always true but not universal.

---

## 4) Current algebraic reduction for Path 1

Path 1 now reduces to proving a direct lower bound for `mu_P(lambda)` rather than relying
on `lambda >= tau_P`.

A robust split is:

1. prove a structural inequality giving `mu_P(tau_P) >= m-2 + c` for some absolute `c>0`
   in this bridge regime;
2. control drift from `tau_P` to `lambda_m(T)` using
   `d mu_P / d lambda = Var_lambda(X_P)/lambda`.

Given observed slack (`>= 0.3628`) and observed worst tie offset (`1.31e-4`), this route
is numerically very stable; the remaining work is symbolic, not computational.

---

## 5) Transfer-constant calibration (Route-2 -> Route-1)

A separate calibration was run to test whether the exact Route-2 threshold

`mu_B >= m - 1 - lambda/(1+lambda)`

implies Route-1 directly by bounding

`D := mu_B - mu_P`.

Define

`exact_excess := D - (1 - lambda/(1+lambda))`.

If `exact_excess <= 0` always, Route-2 exact would imply Route-1 immediately.

Full `d_leaf<=1` frontier through `n<=23` (canonical leaf) gives:

- max `exact_excess` = `0.005107340588772491`,
- positive `exact_excess` cases: `4 / 931,596`.

So the direct implication is not exact, but only misses by a very small additive constant.
Details/artifacts are in:

- `notes/path1_transfer_constant_2026-02-19.md`
- `results/whnc_route1_transfer_scan_n23_merged.json`

---

## 6) Commands used (exact)

Full run:

```bash
python3 conjecture_a_phiP_scan.py --min-n 4 --max-n 23 --out results/whnc_phiP_scan_n23.json
```

Chunked `n=23`:

```bash
for r in 0 1 2 3 4 5 6 7; do
  python3 conjecture_a_phiP_scan.py \
    --min-n 23 --max-n 23 \
    --res "$r" --mod 8 \
    --no-resume \
    --out "results/whnc_phiP_scan_n23_mod8_r${r}.json"
done
```

Merged summary creation:

```bash
python3 - <<'PY'
# merge script used in session (writes results/whnc_phiP_scan_n23_merged.json)
PY
```
