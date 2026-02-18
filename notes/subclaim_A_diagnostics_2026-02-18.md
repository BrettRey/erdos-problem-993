# Sub-claim A diagnostics package (2026-02-18)

Target:

- `margin(k,j+2) >= margin(k,j)` for all `k>=6`, `j>=0`,
  where `margin(k,j)` is the tie-fugacity margin on `S(2^k,1^j)`.

This note records a diagnostics package to support the remaining analytic proof.

---

## 1) New script

- `analyze_subclaim_A_components.py`

It scans `(k,j)` and reports:

1. direct Sub-claim A failures (`full_delta < 0`),
2. mode-step failures (`m(k,j+2)-m(k,j) != 1`),
3. tie-fugacity monotonicity failures (`lambda_{j+2} < lambda_j`),
4. decomposition terms:
   - `full_delta = margin(k,j+2)-margin(k,j)`,
   - `core_delta = mu(k,j+2;lambda_j)-mu(k,j;lambda_j)-(m_{j+2}-m_j)`,
   - `lift_delta = full_delta-core_delta`.

Optional exact checks (`conjecture_a_mixed_spider_exact_margin.py`) are run on the
tightest floating witnesses.

---

## 2) Run used

`python3 analyze_subclaim_A_components.py --k-min 6 --k-max 1200 --j-min 0 --j-max 120 --exact-top 6 --out results/whnc_subclaim_A_components_k6_1200_j120.json`

Observed:

- `Sub-claim A failures = 0`,
- `mode-step failures = 0`,
- `lambda-monotonic failures = 0`,
- `min full delta = 4.6647404176e-05` at `(k,j)=(8,118)`,
- `min core delta = -0.2290221646` at `(k,j)=(8,0)`,
- `min lift delta = 2.3201967456e-04` at `(k,j)=(6,117)`.

Exact checks on tight witnesses stayed positive:

- `(k,j)=(8,118)`: exact `full_delta = 4.6647404168e-05`,
- `(k,j)=(8,116)`: exact `full_delta = 4.8201911896e-05`,
- `(k,j)=(7,117)`: exact `full_delta = 4.9387128876e-05`,
- `(k,j)=(8,114)`: exact `full_delta = 4.9835351140e-05`,
- `(k,j)=(6,118)`: exact `full_delta = 5.0758621691e-05`,
- `(k,j)=(7,115)`: exact `full_delta = 5.1039868636e-05`.

---

## 3) Mechanistic implication

The naive fixed-`lambda_j` term can be substantially negative (`core_delta < 0`), so
Sub-claim A cannot be proved by the fixed-lambda increment alone.

The data indicates the compensation mechanism is:

1. `lambda_{j+2} >= lambda_j` (no failures in this range), and
2. the lambda-lift term `mu(k,j+2;lambda_{j+2})-mu(k,j+2;lambda_j)` dominates the
   negative part of `core_delta`.

So a plausible proof lane is to bound these two pieces separately.
