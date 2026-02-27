# Even-N Adjacent Separation at Fixed `(m,lambda,rho)` (2026-02-27)

## Scope
Min-`u` canonical gated class:
- `T` tree, `d_leaf(T)<=1`, canonical triplet `(ell,s,u)` by min-`u`,
- `B=T-{ell,s}` rooted at `u`,
- `P=dp0[u]`, `Q=dp1[u]`,
- `I(T;x)=(1+2x)P(x)+(1+x)Q(x)=sum i_k x^k`,
- leftmost mode index `m>=4`, `lambda=i_{m-1}/i_m in (0,1)`,
- `N=[x]P=i1-3`, `rho=Q(lambda)/P(lambda)`.

For fixed `(m,lambda)`, define `R_{m,lambda}(N)` as feasible `rho` values at size `N`.

## 1) Strict incidence bound used at `k=m>=4`

For connected graph `G` on `n` vertices and `k>=2` with `i_k>0`:

`i_{k-1}/i_k >= k/(n-k)`.

Reason:
- double-count incidence between independent `(k-1)`-sets and `k`-sets,
- each `k`-set contributes `k` incidences,
- each independent `(k-1)`-set has at most `n-k` extensions since it is nonempty and
  connectedness forces at least one outside vertex adjacent to it.

Apply to tree `T` (`n=N+3`) at `k=m`:

`lambda = i_{m-1}/i_m > m/(N+4-m)`.

(The strict `>` follows since `m/(N+3-m) > m/(N+4-m)`.)

## 2) Mode/geometric + canonical decomposition inequality

Let
`C_{lambda,rho} := (1+2lambda) + (1+lambda)rho`.

Then
`(1-lambda) * C_{lambda,rho} * (1+lambda)^{ceil(N/2)} <= binom(N+3,m)`.

This is the established chain:
- `(1-lambda) I(lambda) <= i_m <= binom(N+3,m)`,
- `I(lambda)=C_{lambda,rho} P(lambda)`,
- `P(lambda)>= (1+lambda)^{ceil(N/2)}` for forest on `N` vertices.

## 3) Even-N adjacent separation theorem

Assume `N` is even and `rho` is feasible at both `N` and `N+1`:

`rho in R_{m,lambda}(N) cap R_{m,lambda}(N+1)`.

From the inequality above at `N=2k` and `N+1=2k+1`:
- `ceil(N/2)=k`, `ceil((N+1)/2)=k+1`.
- Dividing gives
  `1+lambda <= binom(N+4,m)/binom(N+3,m) = (N+4)/(N+4-m)`.
- Hence
  `lambda <= m/(N+4-m)`.

But strict incidence at `k=m>=4` gives
`lambda > m/(N+4-m)`.

Contradiction.

Therefore for every even `N`:

`R_{m,lambda}(N) cap R_{m,lambda}(N+1) = empty`.

## 4) Consequence (sharper blocker)

Any adjacent overlap, if it exists, must be odd-to-even:

`R_{m,lambda}(2k-1) cap R_{m,lambda}(2k)`.

So the unresolved adjacent-separation problem reduces to this parity half only.

## 5) Caveat

The strict incidence bound requires `k>=2`. It is false at `k=1`.
That is irrelevant here because the mode gate is `m>=4`.
