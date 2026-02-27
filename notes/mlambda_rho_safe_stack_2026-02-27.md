# Safe Stack for `(m,lambda,rho)` (Audited) — 2026-02-27

## Scope
Min-`u` canonical gated trees:
- `T` tree,
- `d_leaf(T) <= 1`,
- canonical admissible triplet `(ell,s,u)` chosen by min-`u` rule,
- `B := T - {ell,s}` rooted at `u`,
- `P:=dp0[u]`, `Q:=dp1[u]`,
- `I(T;x)=(1+2x)P(x)+(1+x)Q(x) = sum_k i_k x^k`,
- `m` leftmost mode index, `lambda:=i_{m-1}/i_m` (so `m>=1`),
- `rho:=Q(lambda)/P(lambda)`,
- `N:=[x]P=i1-3`.

## A. Well-definedness
- `P(x)` has nonnegative coefficients and constant term `1`.
- Therefore for `lambda>=0`, `P(lambda) >= 1`, so `rho` is well-defined.
- In canonical uses, `m>=1` and `i_m>0`, hence `lambda` is well-defined and `lambda>0`.

## B. Structural lemmas
### B1) `lambda < 1` from leftmost mode
If `m` is leftmost mode, then `i_{m-1} < i_m`, so `0 < lambda < 1`.

### B2) Admissible support count under `d_leaf<=1`
Define admissible supports as vertices `s` with `deg(s)=2` adjacent to a leaf.
Then the count is never exactly `1` (it is `0` or `>=2`).

### B3) Tree identity for `i2`
For any tree: `i2 = (i1-1)(i1-2)/2`.

## C. Rigidity at low modes
### C1) `m=2`
`lambda = i1/i2 = 2*i1/((i1-1)(i1-2))` and RHS is strictly decreasing in integer `i1>=3`.
So `(m,lambda)` uniquely determines `i1`, hence `N=i1-3`.

### C2) `m=3` with `d_leaf<=1`
Using `i3 = C(n,3) - (n-1)(n-2) + W(T)` and sharp `W` bounds under `d_leaf<=1`,
intervals for feasible `lambda=i2/i3` are disjoint across `n`.
So `(m,lambda)` uniquely determines `i1=n`, hence `N`.

Corollary: `(m,lambda)->N` for all `m<=3` on the canonical gated class.

## D. Exact message identities for `(m>=4)` lane
Define messages on rooted `B` by `t_v(x):=dp1[v](x)/dp0[v](x)`.
For a vertex `v` with children `c_j`:
` t_v(x) = x / prod_j (1 + t_{c_j}(x)) `.
Hence at `x=lambda`:
` 0 <= t_v(lambda) <= lambda `.
At root `u`:
` rho = t_u(lambda) = lambda / prod_j (1+t_{v_j}(lambda)) `.

If root has `r` children in `B` then:
- `rho >= lambda/(1+lambda)^r`,
- and because deleting `{ell,s}` preserves leaf-status of neighbors `v!=s` of `u`,
  `d_leaf<=1` implies at most one leaf child of `u` in `B`, thus `N >= 2r-1`.

Therefore:
` N >= 2*r_min - 1 `,
where `r_min` is least integer with `(1+lambda)^r >= lambda/rho`.

## E. General coefficient-ratio bound
For any graph with `n` vertices and independence coefficients `i_k`:
` i_{k-1}/i_k >= k/(n-k+1) ` (for `i_k>0`).
Applying at `k=m`, `n=i1=N+3`:
` N >= ceil(m - 4 + m/lambda) `.

## F. Proven finite upper bound from `(m,lambda,rho)`
Let `C_{lambda,rho}:=(1+2lambda)+(1+lambda)rho`.
Then:
` (1-lambda) * C_{lambda,rho} * (1+lambda)^{ceil(N/2)} <= binom(N+3,m) `.
So `N <= F(m,lambda,rho)` for explicit finite `F`.

## G. Current blocker (precise)
We have lower bounds and finite upper bounds on `N` from `(m,lambda,rho)`, but no
separation theorem showing disjoint feasible `rho`-sets across different `N` at fixed `(m,lambda)`.

Missing bridge:
For fixed `(m,lambda)` with `m>=4`, prove or refute
`R_{m,lambda}(N) ∩ R_{m,lambda}(N') = empty` for `N!=N'`.
This is equivalent to injectivity `(m,lambda,rho)->N`.
