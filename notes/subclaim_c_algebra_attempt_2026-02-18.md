# Sub-claim C closed (finite base + analytic tail)

Target (Sub-claim C): for `k >= 3`, `k ≡ 0,2 (mod 3)`,

`margin(k,1) < margin(k,0)`.

Here `margin(k,j)` is the tie-fugacity margin of `S(2^k,1^j)`.

---

## 1) Core decomposition

For `j in {0,1}` at its own tie-fugacity `lambda_j` and mode `m_j`:

- `margin_j = A_j + E_j`,
- `E_j = r_j (B_j - A_j)/(1+r_j)`,
- `r_j = lambda_j (1+lambda_j)^(k-j)/(1+2lambda_j)^k`.

So

`margin_1 - margin_0 = -(A_gap) + (E_1 - E_0)`, where `A_gap := A_0 - A_1`.

A sufficient condition is:

1. `E_1 <= 0`,
2. `A_gap > -E_0`.

Then `margin_1 - margin_0 <= -(A_gap) - E_0 < 0`.

---

## 2) `E_1 <= 0` (all `k>=5`)

For `j=1` one gets an exact simplification:

`B_1 - A_1 = [1 - (k-2)lambda_1] / [(1+lambda_1)(1+2lambda_1)]`.

So it is enough to show `lambda_1 > 1/(k-2)`.

From closed forms (with `p=4^t`):

- if `k=3t`:
  `lambda_1 = (2t+1)(p(2t+1)+2t) / ((t+1)(p(4t+1)+2t+1))`.
  Then
  `2*Num - Den = p(4t^2+3t+1) + (2t+1)(3t-1) > 0`, so `lambda_1 > 1/2`.

- if `k=3t+2`:
  `lambda_1 = (p(4t+5)+(2t+1))/((t+2)(4p+1))`.
  Then
  `2*Num - Den = p(4t+2) + 3t > 0`, so `lambda_1 > 1/2`.

Hence for `k>=5`, `1/(k-2) <= 1/3 < 1/2 < lambda_1`, so `B_1-A_1<0`, hence `E_1<0`.

---

## 3) Lower bound on `A_gap`

Use explicit lower benchmark `G(k)` (the `p->infinity` comparison value):

- `k=3t`:
  `G = (t+1)^2(10t+3)/((3t+1)(4t+3)(8t^2+9t+2))`.
- `k=3t+2`:
  `G = 5(t+2)^2/(3(t+1)(2t+3)(8t+13))`.

Exact algebra gives

- `G - 1/(4k) > 0` in both residues:
  - `k=3t`: numerator `24t^4+64t^3+27t^2-17t-6 > 0` for `t>=1`.
  - `k=3t+2`: numerator `12t^3+82t^2+133t+43 > 0` for `t>=0`.

So `G > 1/(4k)`.

Next, `A_gap - G` is a rational function in `t,p` with denominator positive and numerator quadratic in `p`.

For `k=3t` (`t>=4`), numerator has form `a2 p^2 + a1 p + a0` with:

- `a2 >= 27648 t^11`,
- `a1 >= t^10(20736 t - 67263) > 0`,
- `|a0| <= 54000 t^10`.

Since `p=4^t`,

`a2 p^2 + a0 >= 27648*16^t*t^11 - 54000*t^10 > 0`.

Hence `A_gap >= G` for `t>=4` (`k>=12`, residue 0).

For `k=3t+2` (`t>=4`), numerator has form `b2 p^2 + b1 p + b0` with:

- `b2 >= 10500 t^8`,
- `b1 > 0`,
- `|b0| <= 14000 t^7`.

Thus

`b2 p^2 + b0 >= 10500*16^t*t^8 - 14000*t^7 > 0`,

so `A_gap >= G` for `t>=4` (`k>=14`, residue 2).

Therefore on those tails:

`A_gap >= G > 1/(4k)`.

---

## 4) Upper bound on `-E_0`

`E_0 = r_0(B_0-A_0)/(1+r_0)`, so

`-E_0 < r_0(A_0-B_0)`.

Let

- `f(lambda) := A_0-B_0 = [-1 + (k-3)lambda - 2lambda^2]/[(1+lambda)(1+2lambda)]`,
- `q(lambda) := (1+lambda)/(1+2lambda)`.

Then `r_0 = lambda q(lambda)^k`, hence

`-E_0 < q(lambda_0)^k f(lambda_0)` (using `lambda_0<1`).

Also `lambda_0 >= L` from mediant lower bound (`L = U_{m-1}/U_m`):

- residue 0 (`k=3t`): `L = t/(t+1)`;
- residue 2 (`k=3t+2`): `L = (2t+1)/(2(t+2))`.

And `f'(lambda) = -k(2lambda^2-1)/((1+lambda)^2(2lambda+1)^2)`,
so for `lambda >= 1/sqrt(2)`, `f` is decreasing.

### Residue 0 tail (`k=3t`, `t>=5`)

Here `L >= 5/6 > 1/sqrt(2)`, so

`f(lambda_0) <= f(L)`, `q(lambda_0) <= q(L)`.

Compute:

- `f(L) = (3t^3-3t^2-5t-1)/((2t+1)(3t+1)) < t/2`,
- `q(L) = (2t+1)/(3t+1) <= 11/16` for `t>=5`.

Hence

`-E_0 < (t/2)*(11/16)^(3t)`.

So it suffices to show `6t^2*(11/16)^(3t) < 1`.
At `t=5`, this is `0.536... < 1`, and

`S_{t+1}/S_t = ((t+1)/t)^2*(11/16)^3 <= (6/5)^2*(11/16)^3 < 0.47 < 1`,

so true for all `t>=5`.
Thus `-E_0 < 1/(12t) = 1/(4k)` for residue 0 tail.

### Residue 2 tail (`k=3t+2`, `t>=4`)

Here `L >= 3/4 > 1/sqrt(2)`, so again

`f(lambda_0) <= f(L)`, `q(lambda_0) <= q(L)`.

Compute:

- `f(L) = (6t^3+7t^2-11t-11)/(3(t+1)(4t+5)) < t/2`,
- `q(L) = (4t+5)/(6(t+1)) <= 7/10` for `t>=4`.

Hence

`-E_0 < (t/2)*(7/10)^(3t+2)`.

Need `2t(3t+2)*(7/10)^(3t+2) < 1`.
At `t=4`, value is `0.759... < 1`, and

`R_{t+1}/R_t = ((t+1)(3t+5)/(t(3t+2)))*(7/10)^3 <= (5*17)/(4*14)*(7/10)^3 < 0.52 < 1`,

so true for all `t>=4`.
Thus `-E_0 < 1/(4(3t+2)) = 1/(4k)` for residue 2 tail.

---

## 5) Finish + finite base

From Sections 2–4, on analytic tails:

- residue 0, `t>=5` (`k>=15`): `E_1<=0`, `A_gap>1/(4k)>-E_0`;
- residue 2, `t>=4` (`k>=14`): `E_1<=0`, `A_gap>1/(4k)>-E_0`.

So `margin(k,1)-margin(k,0) < 0` on these tails.

Remaining small `k` in Sub-claim C:

`k in {3,5,6,8,9,11,12}`.

These are checked exactly by `prove_subclaim_c_algebra.py` and all satisfy
`margin(k,1) < margin(k,0)`.

Therefore Sub-claim C holds for all `k>=3`, `k ≡ 0,2 (mod 3)`.

---

## Repro command

`python3 prove_subclaim_c_algebra.py --k-max 500 --validate-formulas 180`

Observed in this run:

- `margin(k,1) < margin(k,0)` failures: `0`
- `B1-A1 <= 0` failures (`k>=5`): `0`
- `A_gap >= 1/(4k)` failures (`k>=6`): `0`
- `-E0 < 1/(4k)` failures (`k>=6`): `0`
- `A_gap > -E0` failures (`k>=6`): `0`
