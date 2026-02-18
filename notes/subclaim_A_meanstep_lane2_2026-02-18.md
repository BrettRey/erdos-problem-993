# Sub-claim A via fixed-`lambda` identity + lane 2 (2026-02-18)

Target (`k>=6`, `j>=0`):

- `margin(k,j+2) >= margin(k,j)`.

When `mode(k,j+2)=mode(k,j)+1`, this is equivalent to:

- `mu_{k,j+2}(lambda_{j+2}) - mu_{k,j}(lambda_j) >= 1`.

---

## 1) Fixed-`lambda` identity

For fixed `k,j` and fixed `lambda>0`, define:

- `u = lambda/(1+lambda)`,
- `q = 2lambda/(1+2lambda)`,
- `a = k q`,
- `b = 1 + k u`,
- `r_j = lambda (1+lambda)^(k-j)/(1+2lambda)^k`,
- `s = (1+lambda)^2`,
- `Y_j = (a + j u) - b`.

Then

- `mu_{k,j}(lambda) = (a + j u + r_j b)/(1+r_j)`,
- `mu_{k,j+2}(lambda) = (a + (j+2)u + (r_j/s)b)/(1+r_j/s)`.

Subtracting gives the exact identity:

`mu_{k,j+2}(lambda)-mu_{k,j}(lambda)`
`= [Y_j r_j (s-1) + 2 s u (1+r_j)] / ((1+r_j)(r_j+s))`.

This is the fixed-`lambda` mean-step used below.

---

## 2) Total-step decomposition at tie fugacities

Let `lambda_j` be the tie fugacity for `(k,j)` and `lambda_{j+2}` for `(k,j+2)`.
Then:

`Delta_total := mu_{k,j+2}(lambda_{j+2}) - mu_{k,j}(lambda_j)`
`= Delta_fixed + Delta_gain`,

with

- `Delta_fixed = mu_{k,j+2}(lambda_j) - mu_{k,j}(lambda_j)`,
- `Delta_gain = mu_{k,j+2}(lambda_{j+2}) - mu_{k,j+2}(lambda_j)`.

Lane 2 gives `lambda_{j+2} >= lambda_j`, hence `Delta_gain >= 0` by monotonicity
of `mu` in `lambda`.

So this is exactly the requested lane:

- fixed-`lambda` step + lane-2 gain.

---

## 3) New verifier for this lane

Script:

- `prove_subclaim_A_meanstep_lane2.py`

What it checks:

1. `mode(k,j+2)-mode(k,j)=1`,
2. lane-2 monotonicity `lambda_{j+2}>=lambda_j`,
3. identity `Delta_total = Delta_fixed + Delta_gain`,
4. target `Delta_total >= 1`,
5. helper lower bound using lane-2 proxy
   `u2 = a'_{m_j}/a'_{m_j+1}` from `(1+2x)^k(1+x)^(j+2)`:
   `lb_u2 := mu_{k,j+2}(u2)-mu_{k,j}(lambda_j)`.

---

## 4) Computational outcomes

### Run A (with exact spot checks)

`python3 prove_subclaim_A_meanstep_lane2.py --k-min 6 --k-max 1200 --j-max 120 --exact-top 8 --out results/whnc_subclaim_A_meanstep_lane2_k6_1200_j120_v2.json`

Results:

- `mode-step` failures: `0`,
- lane-2 failures: `0`,
- `Delta_total >= 1` failures: `0`,
- minimum `Delta_total`: `1.000045166827945` at `(k,j)=(8,120)`,
- identity residual max: `0` (to floating precision shown).

Exact-rational checks on tight witnesses are all `>1` (see output JSON).

### Run B (larger range)

`python3 prove_subclaim_A_meanstep_lane2.py --k-min 6 --k-max 4000 --j-max 80 --exact-top 0 --out results/whnc_subclaim_A_meanstep_lane2_k6_4000_j80.json`

Results:

- `mode-step` failures: `0`,
- lane-2 failures: `0`,
- `Delta_total >= 1` failures: `0`,
- minimum `Delta_total`: `1.000089495361863` at `(k,j)=(3998,80)`.

---

## 5) About the `u2` helper bound

The surrogate lower bound `lb_u2 >= 1` is false only on five small pairs:

- `(k,j) = (6,0), (6,2), (7,1), (8,0), (8,2)`.

These are checked exactly in the script output; despite `lb_u2<1`, the actual
`Delta_total` is still strictly `>1` on each of them.

So the fixed-`lambda` + lane-2 decomposition is consistent with Sub-claim A on
all scanned ranges, and the only observed obstruction is to the stronger helper
surrogate, not to `Delta_total >= 1`.
