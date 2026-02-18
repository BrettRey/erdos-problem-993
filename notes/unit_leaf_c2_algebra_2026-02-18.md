# Unit-Leaf `c2 >= 0` Algebraic Proof Draft (2026-02-18)

Goal (unit-leaf branch for `T = S(2^k,1^j)`, `j >= 1`):

`c2 = mu_B(lambda_m^T) - (m-2) >= 0`,

where:
- `m = mode(I_T)` (leftmost mode at `lambda=1`),
- `lambda_m^T = i_{m-1}(T)/i_m(T)`,
- `B` is the unit-leaf deletion graph with
  `I_B(x) = (1+2x)^k (1+x)^(j-1) = sum_t b_t x^t`,
- `mu_B(lambda) = k*2lambda/(1+2lambda) + (j-1)*lambda/(1+lambda)`.

---

## 1) Reduce to a tie-fugacity comparison

Write

`I_T(x) = (1+x) I_B(x) + x(1+x)^k`.

Let `g_t = C(k,t-1)` be coefficients of `x(1+x)^k`, so

`i_t(T) = b_t + b_{t-1} + g_t`.

For `m >= 2`, set `r = m-1` and `tau = b_{r-1}/b_r` (the `B` tie-point between `r-1` and `r`).

Then

`lambda_m^T = (b_r + b_{r-1} + g_r) / (b_{r+1} + b_r + g_{r+1})`.

To prove `lambda_m^T >= tau`, clear denominators:

`(b_r + b_{r-1} + g_r)b_r - (b_{r+1} + b_r + g_{r+1})b_{r-1}`

`= (b_r^2 - b_{r+1}b_{r-1}) + (g_r b_r - g_{r+1}b_{r-1})`.

So it is enough to show both terms are nonnegative.

- First term: `b_r^2 - b_{r+1}b_{r-1} >= 0` because `b_t` is log-concave (`I_B` is a product of linear factors with negative real roots).
- Second term: `g_r b_r - g_{r+1}b_{r-1} >= 0`.

For the second term:
- If `r > k+1`, then `g_r = g_{r+1} = 0`.
- If `r = k+1`, then `g_r = 1`, `g_{r+1}=0`, so term is `b_r > 0`.
- If `1 <= r <= k`, then this is equivalent to
  `b_r/b_{r-1} >= (k-r+1)/r`.

Now `I_B(x) = prod_{i=1}^n (1 + q_i x)` with `n = k+j-1`, and `q_i in {2,1}`.
Thus `b_t = e_t(q_1,...,q_n)` (elementary symmetric polynomials). Standard identity:

`e_r/e_{r-1} = ((n-r+1)/r) * (weighted average of q_i) >= (n-r+1)/r`.

Since `n-r+1 = k+j-r >= k-r+1` (`j>=1`),

`b_r/b_{r-1} >= (k-r+1)/r`,

hence `g_r b_r - g_{r+1}b_{r-1} >= 0`.

Therefore:

`lambda_m^T >= tau = b_{m-2}/b_{m-1}`.

---

## 2) Convert to the mean inequality

Because

`mu_B(lambda) = k*2lambda/(1+2lambda) + (j-1)*lambda/(1+lambda)`

is strictly increasing for `lambda > 0`, we get

`mu_B(lambda_m^T) >= mu_B(tau)`.

So it remains to show `mu_B(tau) >= m-2`.

At `lambda=tau`, weighted coefficients are `w_t = b_t tau^t` and satisfy `w_{r-1}=w_r`.
Since `b_t` is log-concave, `w_t` is log-concave, so ratio `w_t/w_{t-1}` is nonincreasing.
With `w_r/w_{r-1}=1`, index `r` is a mode of `w`.

Now interpret `w` as a Poisson-binomial distribution (independent Bernoulli sum with odds `tau q_i`).
By Darroch's mode-mean bound (`|mode - mean| < 1`), if `r` is a mode then

`mu_B(tau) > r-1 = m-2`.

Hence

`mu_B(lambda_m^T) >= mu_B(tau) >= m-2`,

which is exactly `c2 >= 0`.

---

## 3) Sanity checks run in this workspace

Numerical scans (fast coefficient engine) on mixed spiders:
- Range `k=1..3000`, `j=1..40`.
- No violations of `lambda_m^T >= b_{m-2}/b_{m-1}`.
- Minimum observed `lambda` gap:
  `lambda_m^T - b_{m-2}/b_{m-1} = 0.000738015...` at `(k,j,m)=(2999,40,2019)`.
- Minimum observed `mu_B(b_{m-2}/b_{m-1}) - (m-2) = 0.333555...`.
- Minimum observed `c2 = mu_B(lambda_m^T) - (m-2) = 0.833513...`.

So the proof skeleton matches large-range numerics.
