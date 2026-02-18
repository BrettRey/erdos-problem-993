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

## 4) Status update (2026-02-18)

**`E_1 <= 0` is now PROVED** (see `notes/subclaim_c_E1_le0_proof_2026-02-18.md`):

Key algebraic identity: `B_1 - A_1 = [1 - (k-2)lambda_1] / [(1+lambda_1)(1+2lambda_1)]`,
so `E_1 <= 0 iff lambda_1 >= 1/(k-2)`.

The b-bound from `unit_leaf_c2_algebra_2026-02-18.md` gives `lambda_1 >= tau =
(m_1-1)/[2(k-m_1+2)]`. Then `(k-2)*tau >= 1` reduces to:
- `k = 3t+2`, `t >= 1`: `6t^2+t-4 >= 0` (true for all `t >= 1`, min value 3 at `t=1`).
- `k = 3t`,   `t >= 2`: `3t^2-3t-1 >= 0` (true for all `t >= 2`, min value 5 at `t=2`).

Verified exactly: `verify_subclaim_c_E1_le0.py`, all 5 checks PASS (k <= 200).

**Remaining gap for a full paper-proof of Sub-claim C:**

`A_gap > -E_0` globally. This is verified for `k <= 500` in `prove_subclaim_c_algebra.py`
(0 failures). For large `k`, `A_gap >= 1/(4k)` while `|E_0|` decays exponentially. A
closed-form proof of `A_gap >= 1/(4k)` for all `k` would close Sub-claim C completely.
