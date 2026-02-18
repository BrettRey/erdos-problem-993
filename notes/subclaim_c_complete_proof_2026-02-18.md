# Sub-Claim C Closed (2026-02-18)

Target:

- For all `k >= 3` with `k ≡ 0,2 (mod 3)`,
  `margin(k,1) < margin(k,0)` for mixed spiders `S(2^k,1^j)`.

This note closes Sub-claim C by combining:

1. the previously proved lemma `E1 <= 0` (see `notes/subclaim_c_E1_le0_proof_2026-02-18.md`);
2. a new analytic lower bound `A_gap > 1/(4k)` for `k >= 14`;
3. a new analytic upper bound `-E0 < 1/(4k)` for `k >= 14`;
4. exact finite checks for the remaining small `k`.

---

## 1) Decomposition

For `j in {0,1}`:

- `margin_j = A_j + E_j`,
- `A_gap = A_0 - A_1`,
- `diff := margin(k,1) - margin(k,0) = -A_gap + (E1 - E0)`.

With `E1 <= 0`, we have:

`diff <= -A_gap - E0`.

So it suffices to prove `A_gap > -E0`.

---

## 2) New bound: `A_gap > 1/(4k)` for `k >= 14`

Write `k=3t` or `k=3t+2`, `p=4^t`.

Using the closed-form tie fugacities (`lambda_j0_closed`, `lambda_j1_closed`), the exact
quantity

`D(k) := A_gap - 1/(4k)`

has exact rational forms:

- For `k=3t`:
  `D = N0 / Den0`, with
  `N0 = p^3*c3 + p^2*c2 + p*c1 + c0`
  and explicit coefficient polynomials `c3,c2,c1,c0`.
- For `k=3t+2`:
  `D = N2 / Den2`, with
  `N2 = p^3*d3 + p^2*d2 + p*d1 + d0`
  and explicit coefficient polynomials `d3,d2,d1,d0`.

The exact coefficient formulas are implemented in
`prove_subclaim_c_complete.py` and validated against direct exact-fraction computation.

For `t >= 4` (`k >= 12` in residue-0, `k >= 14` in residue-2):

- residue 0: `c2>0`, `c1>0`,
- residue 2: `d2>0`, `d1>0`,

so

- `N0 > p^3*c3 + c0`,
- `N2 > p^3*d3 + d0`.

Then coarse bounds give:

- residue 0: `N0 > 72 t^6 (64^t - 138 t)`,
- residue 2: `N2 > 144 t^4 (64^t - 69 t)`.

Since `64^t/t` is increasing and huge at `t=4`, both are strictly positive.
Hence `D>0`, i.e. `A_gap > 1/(4k)` for `k >= 14`.

Small `t` bridge checks (exact):

- residue 0: `t=2,3` (`k=6,9`) pass;
- residue 2: `t=2,3` (`k=8,11`) pass.

---

## 3) New bound: `-E0 < 1/(4k)` for `k >= 14`

For `j=0` branch:

`-E0 = r0*(A0-B0)/(1+r0) < r0*(k/6)`

because `lambda<=1` at tie fugacity and
`A0-B0 = k*lambda/((1+lambda)(1+2lambda)) - 1 < k/6`.

Also:

- `lambda0 >= a_k`,
  where
  - `a_k = k/(k+3)` for `k ≡ 0`,
  - `a_k = (k-1)/(k+2)` for `k ≡ 2`.

Thus

`r0 = lambda0 * ((1+lambda0)/(1+2lambda0))^k <= s_k^k`,

with

- `s_k = (2k+3)/(3k+3)` for `k ≡ 0`,
- `s_k = (2k+1)/(3k)` for `k ≡ 2`.

For all `k>=14` in these residues, `s_k <= 29/42`, so

`-E0 < (k/6)*(29/42)^k`.

Set `f(k)=k^2*(29/42)^k`. Then

`f(k+1)/f(k) = (29/42)*((k+1)/k)^2 <= (29/42)*(15/14)^2 < 1` for `k>=14`,

so `f` decreases on `k>=14`. With exact base:

`f(14) < 3/2`,

we get `f(k) < 3/2` for all `k>=14`, hence

`-E0 < (1/(6k))*(3/2) = 1/(4k)`.

---

## 4) Finish Sub-claim C

For `k>=14`, combine:

- `A_gap > 1/(4k)`,
- `-E0 < 1/(4k)`,
- `E1 <= 0` (proved separately),

to get

`diff = -A_gap + E1 - E0 < -A_gap - E0 < 0`.

Remaining small cases are exact-checked:

`k in {3,5,6,8,9,11,12}` all satisfy `diff < 0`.

Therefore Sub-claim C holds for all `k>=3`, `k ≡ 0,2 (mod 3)`.

---

## Reproducibility

Script:

- `prove_subclaim_c_complete.py`

Output artifact:

- `results/whnc_subclaim_c_complete_proof_check.txt`

The script performs:

1. exact identity validation for the transcribed `A_gap - 1/(4k)` rational forms,
2. analytic inequality checks for the large-`k` lane,
3. exact finite-base checks,
4. a direct sanity scan up to `k=500`.
