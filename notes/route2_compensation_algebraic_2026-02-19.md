# Route-2 Compensation: Algebraic Reduction and Slope-Closure (2026-02-19)

## 0) Setup

Fix a `d_leaf <= 1` tree `T`, choose a bridge leaf `l` whose support `s` has
`deg(s)=2`, let `u` be the other neighbor of `s`, and define

- `B = T - {l,s}`,
- `m = mode(I(T))` (leftmost),
- `lambda = lambda_m(T) = i_{m-1}(T)/i_m(T)`,
- `tau = b_{m-2}/b_{m-1}` where `I(B)=sum_k b_k x^k`.

Route-2 target:

`mu_B(lambda) >= m - 3/2`.

Stronger pendant-bonus threshold:

`mu_B(lambda) >= m - 1 - lambda/(1+lambda)`.

---

## 1) Exact deficit/lift identity

Define

- `deficit_tau := (m-3/2) - mu_B(tau)`,
- `gain := mu_B(lambda) - mu_B(tau)`.

Then exactly:

`mu_B(lambda) - (m-3/2) = gain - deficit_tau`.

So Route-2 is equivalent to

`gain >= deficit_tau`.

Also

`d mu_B / d t = Var_B(X_t)/t` for `t>0`, where `X_t` is IS size under fugacity `t`.
Hence

`gain = int_{tau}^{lambda} Var_B(X_t)/t dt`.

Therefore Route-2 is exactly:

`int_{tau}^{lambda} Var_B(X_t)/t dt >= deficit_tau`.

This is the clean compensation form.

---

## 2) Exact tie-gap formula (`lambda - tau`)

From bridge coefficients

- `i_{m-1} = b_{m-1} + b_{m-2} + p_{m-2}`,
- `i_m     = b_m     + b_{m-1} + p_{m-1}`

(where `P=dp_B[u][0]=sum p_k x^k`), we get

`lambda = (b_{m-1}+b_{m-2}+p_{m-2})/(b_m+b_{m-1}+p_{m-1})`.

Subtract `tau=b_{m-2}/b_{m-1}` and clear denominators:

`lambda - tau = [ (b_{m-1}^2 - b_m b_{m-2}) + (p_{m-2} b_{m-1} - p_{m-1} b_{m-2}) ] / [ b_{m-1}(b_m+b_{m-1}+p_{m-1}) ]`.

So the tie-gap numerator is exactly the STRONG-C2 determinant split
`lc_surplus + mismatch`.

---

## 3) Concavity + endpoint slope gives a one-line closure

For `mu_B` concave on `[tau,lambda]`, the derivative is decreasing, so

`gain = int_{tau}^{lambda} mu_B'(t) dt >= (lambda-tau) * mu_B'(lambda)`.

Using `mu_B'(lambda)=Var_B(X_lambda)/lambda`:

`gain >= (lambda-tau) * Var_B(X_lambda)/lambda`.

Hence a sufficient condition for Route-2 is

`(Var_B(X_lambda)/lambda) * (lambda-tau) >= deficit_tau`.   (C1)

This turns Route-2 into a purely endpoint inequality once concavity is known.

---

## 4) New verification scripts and outputs

New script 1:

- `verify_route2_slope_compensation.py`

What it checks:

- Route-2 and stronger inequality slacks,
- exact compensation `gain - deficit_tau`,
- endpoint margins (deficit cases):
  - `M_tau = (Var_B(tau)/tau)*(lambda-tau) - deficit_tau`,
  - `M_lam = (Var_B(lambda)/lambda)*(lambda-tau) - deficit_tau`,
- optional sampled concavity of `mu_B` on `[tau,lambda]`.

Outputs produced:

- `results/route2_slope_compensation_n20.json`
- `results/route2_slope_compensation_n21.json`
- `results/route2_slope_compensation_n22.json`
- `results/route2_slope_compensation_n23.json`
- merged:
  `results/route2_slope_compensation_n23_merged.json`

Merged summary (`n<=23`, canonical leaf, `d_leaf<=1`):

- checked leaves: `931,595`
- route-2 failures: `0`
- stronger failures: `0`
- deficit-at-`tau` cases: `535,095`
- minimum route-2 slack: `0.2053400949297961`
- minimum stronger slack: `0.1913484628930444`
- minimum exact compensation margin `gain-deficit_tau`: `0.2053400949297961`
- minimum endpoint margins:
  - `min M_tau = 0.23382128170267158`
  - `min M_lam = 0.17929045776096725`
- worst `M_lam` witness is the same `n=20` star-like extremal route-2 witness
  (`g6 = S???????C?G?G?C?@??G??_?@??@?F~_?`).

Concavity sampling (deficit cases, grid=24) was enabled for `n<=21`:

- concavity checks: `39,903`
- concavity failures: `0`
- maximum sampled second difference: `0.0`.

---

New script 2:

- `verify_mu_lambda_concavity.py`

Purpose: direct sampled concavity/convexity profile of `mu(lambda)` for tree
IS polynomials (whole trees, not just bridge `B`).

Output:

- `results/mu_lambda_concavity_n16.json`

Result (`n<=16`, all trees, lambda-grid `0.1..1.9` with 19 points):

- checked trees: `32,508`
- sampled concavity failures: `0`
- sampled convexity failures: `32,508`
- all sampled second differences are strictly negative (`max_d2 < 0`).

So on this grid, `mu` is uniformly concave and never convex.

---

## 5) Route-2 status after this reduction

What is now rigorous algebra:

1. Route-2 is exactly `gain >= deficit_tau` with
   `gain = int_{tau}^{lambda} Var_B(X_t)/t dt`.
2. `lambda-tau` has an exact determinant split formula.
3. If concavity of `mu_B` on `[tau,lambda]` is available, endpoint condition
   `(C1)` is a sufficient one-line closure.

What remains to finish a full analytic proof:

- prove concavity of `mu_B` on the bridge interval `[tau,lambda]`
  (or an equivalent monotonicity of `Var_B(t)/t` strong enough for endpoint
  lower bounds), and
- prove endpoint slope condition `(C1)` structurally (without scan).

Empirically through the tested frontier, both ingredients are strongly supported:
- sampled concavity has no failures where tested (`n<=21`),
- endpoint slope margin `M_lam` stays uniformly positive with substantial slack
  on full canonical `n<=23` (`min M_lam = 0.17929045776096725`).

---

## 6) Relation to existing full-frontier compensation artifact (`n<=23`)

Existing canonical full frontier result:

- `results/whnc_route2_compensation_scan_n23_canonical.json`

already gives (through full `n<=23`):

- route-2 failures: `0`,
- stronger failures: `0`,
- max deficit at `tau`: `0.1000959551`,
- min gain in deficit cases: `0.2901109960`.

So the new slope/concavity reduction is consistent with the known full-frontier
compensation behavior and identifies a concrete analytic closure target.
