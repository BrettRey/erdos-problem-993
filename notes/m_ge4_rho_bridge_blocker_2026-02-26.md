# `m >= 4` Bridge Blocker for `(m,lambda,rho) -> N`

## Context
For min-`u` canonical gated trees (`d_leaf <= 1`), we already have:
- `(m,lambda) -> N` for `m <= 3` (proved).
- no split through `n <= 25` for keys `(m,lambda,rho)` and `(m,lambda,sigma)` (empirical).

The unresolved step is a theorem for `m >= 4`.

## Exact missing bridge
Let `T` be in the min-`u` canonical gated class with:
- `I(T;x)=(1+2x)P(x)+(1+x)Q(x)`
- `m` leftmost mode, `lambda=i_{m-1}/i_m`, `rho=Q(lambda)/P(lambda)`
- `N=[x]P=i1-3`.

The single hard bridge can be stated as:

> **Root-message size separation lemma (missing).**  
> For fixed `(m,lambda)` with `m>=4`, the feasible `rho` values for different `N`
> are disjoint:
> `R_{m,lambda}(N) ∩ R_{m,lambda}(N') = empty` for `N != N'`,
> where `R_{m,lambda}(N)` is the set of `rho(T)` over canonical trees with those
> fixed `(m,lambda)` and size `N`.

This is equivalent to injectivity of `rho` in `N` at fixed `(m,lambda)`, i.e.
`(m,lambda,rho) -> N`.

## What is already rigorous
From DP messages `t_v = dp1[v]/dp0[v]`, we have
`t_u(lambda)=rho=lambda/prod_j (1+t_{v_j}(lambda))`, and `0 <= t_{v_j}(lambda) <= lambda`.
Hence:
- `rho <= lambda`
- with root child count `r`: `rho >= lambda/(1+lambda)^r`
- under `d_leaf<=1`: at most one size-1 child component at root, so `N >= 2r-1`.

Therefore:
`N >= 2*ceil(log(lambda/rho)/log(1+lambda)) - 1`.

Also from the general coefficient-ratio bound:
`i_{k-1}/i_k >= k/(n-k+1)` and `n=i1=N+3`, so at `k=m`:
`N >= m - 4 + m/lambda`.

These are clean lower bounds from `(m,lambda,rho)`, but they do not yet force
uniqueness of `N`; an upper/separation mechanism is exactly what is missing.

## Audit note
The only structural lemma requiring careful wording is the admissible triplet count:
it must be explicitly defined as counting degree-2 supports adjacent to leaves.
