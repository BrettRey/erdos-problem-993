# Sub-claim C: algebraic attempt (2026-02-18)

Target:

- For `k >= 3`, `k ≡ 0,2 (mod 3)`, show
  `margin(k,1) < margin(k,0)`
  for mixed spiders `S(2^k,1^j)`.

## 1) Decomposition used

For branch `j in {0,1}` at its own tie fugacity `lambda_j` and mode `m_j`:

- `margin(k,j) = A_j + E_j`,
- `E_j = r_j (B_j - A_j)/(1+r_j)`,
- `r_j = lambda_j (1+lambda_j)^(k-j) / (1+2lambda_j)^k`.

With

- `A_0 = 2k*lambda_0/(1+2lambda_0) - (m_0-1)`,
- `B_0 = 1 + k*lambda_0/(1+lambda_0) - (m_0-1)`,
- `A_1 = 2k*lambda_1/(1+2lambda_1) + lambda_1/(1+lambda_1) - (m_1-1)`,
- `B_1 = 1 + k*lambda_1/(1+lambda_1) - (m_1-1)`.

Let `A_gap = A_0 - A_1`. Then

`margin(k,1)-margin(k,0) = -A_gap + (E_1 - E_0)`.

Hence a sufficient condition is:

1. `E_1 <= 0`, and
2. `A_gap > -E_0`.

Indeed then
`margin(k,1)-margin(k,0) <= -A_gap - E_0 < 0`.

## 2) Closed forms for tie fugacities

`lambda_0 = i_{m_0-1}/i_{m_0}`, `lambda_1 = i_{m_1-1}/i_{m_1}`.

For `k=3t`:

- `lambda_0 = [4^t/2 + (2t-1)/(t+2)] / [1 + 4^t (t+1)/(2t)]`,
- `lambda_1 = (2t+1)[4^t(2t+1)+2t] / [(t+1)(4^t(4t+1)+2t+1)]`.

For `k=3t+2`:

- `lambda_0 = [4^t + 2t/(t+3)] / [1 + 2*4^t (t+2)/(2t+1)]`,
- `lambda_1 = [4^t(4t+5)+(2t+1)] / [(t+2)(4*4^t+1)]`.

These formulas are validated exactly against coefficient ratios in the script.

## 3) Current status (computer-assisted algebraic checks)

Script:

- `prove_subclaim_c_algebra.py`

Run:

- `python3 prove_subclaim_c_algebra.py --k-max 500 --validate-formulas 180`

Exact-rational outcomes:

- `margin(k,1) < margin(k,0)` for all checked `k ≡ 0,2 (mod 3)`, `3 <= k <= 500`.
- `B_1 - A_1 <= 0` for all checked `k >= 5`.
- `A_gap >= 1/(4k)` for all checked `k >= 6`.
- `A_gap > -E_0` for all checked `k >= 6`.

Finite base witnesses (exact):

- `k=3`: `diff = -709603816133 / 6108986807044`.
- `k=5`: `diff = -871270001341 / 10827617187390`.
- `k=6`: `diff = -21495253443141677102383384978632 / 351455149421499530933336519680607`.
- `k=8`: `diff = -24550857995345285589064229633421759052 / 501196238048458555120781129690463817441`.

## 4) Remaining gap to a full symbolic proof

The checks above support the sufficient lane, but a complete paper-proof still needs
closed-form inequalities for all `t` proving:

- `E_1 <= 0` (equivalently `B_1 <= A_1`) for all relevant `k`, and
- `A_gap > -E_0` globally.

So this is a strong algebraic reduction + exact verification, not yet a fully closed
all-`k` symbolic proof.
