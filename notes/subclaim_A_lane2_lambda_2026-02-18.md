# Sub-claim A, Lane #2 (lambda monotonicity): progress note

Goal for lane #2:

- Prove `lambda(k,j+2) >= lambda(k,j)` for mixed spiders `S(2^k,1^j)`.

Here `lambda(k,j) = i_{m-1}/i_m` at `m = mode(I_{k,j})`.

## New script and artifact

- Script: `prove_subclaim_A_lambda_monotone.py`
- Run:
  `python3 prove_subclaim_A_lambda_monotone.py --k-max 3000 --j-max 120 --out results/subclaimA_lane2_lambda_monotone_k3000_j120.json`
- Output artifact:
  `results/subclaimA_lane2_lambda_monotone_k3000_j120.json`

Checked `362,395` `(k,j)` pairs.

Observed failures: **0** for all of:

1. `mode(k,j+2) = mode(k,j) + 1` (in this scan range).
2. `lambda(k,j+2) >= lambda(k,j)`.
3. Stronger sufficient condition `u2 >= lambda(k,j)` (defined below).

## Algebraic reduction used

Fix `(k,j)`, let `m = mode(k,j)` and write:

- `a_t = [x^t](1+2x)^k(1+x)^j`.
- `a'_t = [x^t](1+2x)^k(1+x)^(j+2)`.
- `g_t = [x^t]x(1+x)^k = C(k,t-1)`.

When `mode(k,j+2)=m+1`,

- `lambda(k,j)   = (a_{m-1}+g_{m-1})/(a_m+g_m)`.
- `lambda(k,j+2) = (a'_m+g_m)/(a'_{m+1}+g_{m+1})`.

Define

- `u2 = a'_m/a'_{m+1}`.

Then `lambda(k,j+2)` is a weighted average of `u2` and `g_m/g_{m+1}`, so
proving `u2 >= lambda(k,j)` is sufficient for lane #2.

Equivalent cross-product form:

`u2 >= lambda(k,j)` iff

`F := (A+2B+C)(A+G) - (D+2A+B)(B+H) >= 0`,

with

- `A=a_m, B=a_{m-1}, C=a_{m-2}, D=a_{m+1}, G=g_m, H=g_{m-1}`.

Decomposition:

- `F = T1 + T2`, where
- `T1 = (A^2 - BD) + (AC - B^2)`,
- `T2 = G(A+2B+C) - H(2A+B+D)`.

Scan diagnostics (`k<=3000, j<=120`):

- minimum `F` witness at `(k,j)=(6,1)`: `F=10644`.
- minimum `T1` witness at `(k,j)=(6,0)`: `T1=15680`.
- `T2` was always negative in the tested range (`0` positive, `135,187` negative in a focused decomposition scan).

## Current proof status

- Computationally very strong support for lane #2, including the stronger `u2 >= lambda_j` condition.
- Remaining symbolic step: close `F>=0` algebraically for all `k>=6,j>=0`.

This reduction should help the full Sub-claim A proof effort: once `F>=0` is proved
symbolically, lane #2 is done.

---

## Update: envelope lane (same date)

A stronger algebraic envelope step was added in:

- `notes/subclaim_A_lane2_envelope_2026-02-18.md`
- script `prove_subclaim_A_lane2_envelope.py`

Key idea: in the hard `T2<0` regime (`F=T1+T2`), write
`F/A^2 = T1/A^2 + s*(T2/(AG))` with `s=G/A`, then bound `s` above by a closed
coefficient lower bound on `A`. This yields an explicit lower envelope for `F`
that was verified with zero failures on large scans (up to `k<=3000, j<=120`).
