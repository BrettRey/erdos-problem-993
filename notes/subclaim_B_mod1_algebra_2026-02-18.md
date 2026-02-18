# Sub-Claim B (`k ≡ 1 mod 3`) Algebraic Lane (2026-02-18)

Target (mixed-spider minimizer reduction):

For `k >= 4`, `k ≡ 1 (mod 3)`,

`margin(k,0) < margin(k,1)`,

where `margin(k,j) = mu(S(2^k,1^j), lambda_m) - (m-1)` at the tie fugacity.

Set `k = 3t+1` (`t>=1`). Then both branches have mode `m = 2t+1`.

---

## 1) Decompose both margins

For branch `j` (`j=0` or `1`):

`margin_j = A_j(lambda_j) - p_j * g_j`,

with:
- `A_0(lambda) = k*2lambda/(1+2lambda) - (m-1)`,
- `A_1(lambda) = k*2lambda/(1+2lambda) + lambda/(1+lambda) - (m-1)`,
- `p_j = r_j/(1+r_j)`, `r_j = lambda_j * ((1+lambda_j)/(1+2lambda_j))^(k-j)`,
- `g_j = k*lambda_j/((1+lambda_j)(1+2lambda_j)) - 1`.

So

`margin_1 - margin_0 = (A_1-A_0) + p_0 g_0 - p_1 g_1`.

For `t>=2` (`k>=7`), `lambda_0 <= 1` and hence `g_0 >= k/6 - 1 > 0`, so

`margin_1 - margin_0 >= (A_1-A_0) - p_1 g_1`.

---

## 2) Main positive gap term

Define:
- `a = (2t+1)/(2t+2)` (dominant-step ratio for `j=0`),
- `r_F = 2(2t+1)(t+1)/((t+2)(4t+3))` (dominant-step ratio for `j=1`).

Using mediant:
- `lambda_1 >= r_F`,
- `A_1(lambda_1) >= A_1(r_F)`.

Also:

`A_0(lambda_0) = A_0(a) + delta_0`, `delta_0 >= 0`.

Hence:

`A_1(lambda_1) - A_0(lambda_0) >= G(t) - delta_0`,

where

`G(t) := A_1(r_F) - A_0(a)`

and exact simplification gives:

`G(t) = (66 t^4 + 227 t^3 + 269 t^2 + 134 t + 24) /`
`       (288 t^5 + 1356 t^4 + 2477 t^3 + 2196 t^2 + 948 t + 160)`.

All coefficients are positive.

Useful lower bound (for `t>=2`):

`G(t) > 1/(12t)`.

Equivalent polynomial after clearing denominators:

`12t*G(t)-1 = P(t)/[t*(...positive...)]`,

`P(t)=504 t^5 + 1368 t^4 + 751 t^3 - 588 t^2 - 660 t - 160 > 0` for `t>=2`.

---

## 3) Upper bounds on error terms

### 3a) `delta_0` bound

Write `j=0` coefficients as `U_r + V_r` with:
- `U_r = 2^r C(k,r)`,
- `V_r = C(k,r-1)`.

Then

`lambda_0 = (a + b c)/(1+c)`,

with
- `b = 2t/(t+2)`,
- `c = (2t+1)/(2^(2t+1)(t+1))`.

Using monotonicity of `A_0` and mean-value:

`delta_0 <= A_0'(a) * (lambda_0-a)`,

`A_0'(a)=2(3t+1)(t+1)^2/(3t+2)^2`,

`lambda_0-a <= c(b-a)`,

which yields

`delta_0 <= ((3t+1)(2t^2-t-2)(2t+1)) / ((3t+2)^2 (t+2) 2^(2t+1))`.

Crude simplification (for `t>=4`):

`delta_0 < (4/3) * t / 4^t <= 4/(3 t^3)`.

### 3b) `p_1 g_1` bound

`p_1 <= r_1 <= (2/3)^(k-1) = (2/3)^(3t)`,

and

`g_1 < k/3 - 1 = (3t-2)/3`.

So

`p_1 g_1 <= ((3t-2)/3) * (2/3)^(3t)`.

For `t>=6`, use `t*(8/27)^t <= 1/t^3`:

`p_1 g_1 <= 1/t^3`.

---

## 4) Final inequality for large `t`

For `t>=6`:

`margin_1 - margin_0`
`>= G(t) - delta_0 - p_1 g_1`
`> 1/(12t) - 4/(3t^3) - 1/t^3`
`= (t^2 - 28)/(12 t^3) > 0`.

So Sub-claim B is proved for `k=3t+1 >= 19`.

---

## 5) Finite base cases

Remaining `t=1..5` (`k=4,7,10,13,16`) were checked exactly:

- `k=4`: `margin(1)-margin(0) = 0.08144364258... > 0`
- `k=7`: `margin(1)-margin(0) = 0.06679986417... > 0`
- `k=10`: `margin(1)-margin(0) = 0.05245910462... > 0`
- `k=13`: `margin(1)-margin(0) = 0.04303670755... > 0`
- `k=16`: `margin(1)-margin(0) = 0.03639627231... > 0`

Therefore:

`margin(k,0) < margin(k,1)` for all `k>=4`, `k ≡ 1 (mod 3)`.

---

## Reproducibility

Script:

- `prove_subclaim_B_mod1.py`

It prints:
- exact base-case gaps (`t=1..5`),
- positivity of the analytic lower bound for `t=6..500`,
- helper inequalities used in the large-`t` bound.
