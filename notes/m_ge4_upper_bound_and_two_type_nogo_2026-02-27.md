# `m >= 4` Progress: Upper Bound from `(m,lambda,rho)` and Two-Type No-Go

## 1) Proven partial theorem: explicit finite upper bound on `N`

In min-`u` canonical gated trees (`d_leaf<=1`) with

- `I(T;x)=(1+2x)P(x)+(1+x)Q(x)`,
- `m` leftmost mode, `lambda=i_{m-1}/i_m`,
- `rho=Q(lambda)/P(lambda)`,
- `N=[x]P=i1-3`,

define

`C_{lambda,rho} := (1+2*lambda) + (1+lambda)*rho`.

Then

`(1-lambda) * C_{lambda,rho} * (1+lambda)^{ceil(N/2)} <= binom(N+3, m)`.  (★)

Sketch:
- `I(lambda)=C_{lambda,rho}*P(lambda)`.
- Mode bound: `i_m >= (1-lambda) I(lambda)`.
- Forest lower bound at fixed `N`: `P(lambda) >= (1+lambda)^{ceil(N/2)}` (bipartition argument).
- Trivial combinatorial upper bound: `i_m <= binom(N+3,m)`.

Hence (★). This gives

`N <= F(m,lambda,rho)` for an explicit finite `F`, i.e. only finitely many `N`
can match fixed `(m,lambda,rho)`.

Status: this is not injectivity; it is an upper bound plus finiteness.

## 2) Two-type heterogeneous root template: exact kernel criterion

For a two-type root decomposition at fixed `lambda`:

- `P = F1^a F2^b`,
- `Q = x G1^a G2^b`,
- `r1=G1(lambda)/F1(lambda)`, `r2=G2(lambda)/F2(lambda)`,
- `rho = lambda * r1^a * r2^b`.

Equal `(lambda,rho)` between `(a,b)` and `(a',b')` is equivalent to:

`r1^(a-a') * r2^(b-b') = 1`.

Define kernel

`H(r1,r2) = {(x,y) in Z^2 : r1^x r2^y = 1}`.

Then:
- uniqueness of `(a,b)` from `(lambda,rho)` iff `H={(0,0)}`,
- if `H=Z*(p,q)` rank-1, then all collisions are `(a',b')=(a,b)+k(p,q)`,
- and `Delta N = k*(n1*p + n2*q)` for type sizes `n1,n2`.

So two-type splits are possible only if there is a nontrivial kernel direction
with nonzero `n1*p+n2*q`.

## 3) Smallest explicit two-type canonical family and no-go result

Use child types:
- Type E: rooted edge endpoint, `F_E=1+2x`, `G_E=1+x`, size `n1=2`.
- Type P: rooted path-3 endpoint, `F_P=1+3x+x^2`, `G_P=1+2x`, size `n2=3`.

Construct min-`u` canonical gated trees by attaching `t` copies of E and `b` copies of P
at root `u`, plus one canonical support-leaf gadget `u-s-l`.
Then

`P=(1+2x)^t(1+3x+x^2)^b`,
`Q=x(1+x)^t(1+2x)^b`,
`N=2t+3b`,
`rho/lambda = ((1+lambda)/(1+2lambda))^t * ((1+2lambda)/(1+3lambda+lambda^2))^b`.

For rational `lambda in (0,1)`, equality of `(lambda,rho)` for two pairs `(t,b)!=(t',b')`
would force a monic integer polynomial identity in `lambda` with rational root `lambda`,
hence integer `lambda` (rational-root theorem), contradiction.

Therefore this smallest heterogeneous template cannot produce same `(lambda,rho)` with different `N`.

## 4) Frontier implication

We now have:
- lower bounds from `(m,lambda,rho)`,
- an explicit finite upper bound from `(m,lambda,rho)`,
- and a no-go for the smallest heterogeneous two-type family.

Still missing for closure:
- a global separation/injectivity lemma that forbids overlap of feasible `rho`
  across different `N` at fixed `(m,lambda)` in the full canonical gated class.

## 5) Audit-checked assumptions for the upper-bound chain

For the inequality

`(1-lambda) * C_{lambda,rho} * (1+lambda)^{ceil(N/2)} <= binom(N+3,m)`,

the assumptions that must be explicit are:

- `lambda` is derived from a leftmost mode, so `0 < lambda < 1`.
- `P` is the canonical `dp0[u]` from a rooted forest remainder, not an arbitrary polynomial.
- `P(lambda) > 0` (true since `P` has nonnegative coefficients and constant term `1`).
- `rho = Q(lambda)/P(lambda)` at the same `lambda`.
- `i1 = N+3` (canonical coefficient identity), so `|V(T)| = N+3`.

With these stated, every step is direct:

- `I(lambda) = C_{lambda,rho} * P(lambda)`,
  where `C_{lambda,rho} = (1+2*lambda) + (1+lambda)*rho`.
- Forest lower bound: `P(lambda) >= (1+lambda)^{ceil(N/2)}` (bipartition subset count).
- Mode upper-series bound: `(1-lambda) * I(lambda) <= i_m`.
- Counting upper bound: `i_m <= binom(N+3,m)`.

Combining gives the claimed inequality.

## 6) Two strict tightenings (same hypotheses)

### 6.1 Degree-based strengthening

Let `d := deg(P) = alpha(F)` for the forest `F` behind `P`.
Then `P(lambda) >= (1+lambda)^d`, hence:

`(1-lambda) * C_{lambda,rho} * (1+lambda)^d <= binom(N+3,m)`.

This is at least as strong as the `ceil(N/2)` version since `d >= ceil(N/2)` for bipartite `F`.

### 6.2 Two-disjoint-admissible-edge strengthening for `m>=4`

In the gated class with canonical bridge existence, admissible support count is never `1`;
when nonzero it is at least `2`, giving two vertex-disjoint leaf-support edges.
Excluding all `m`-sets containing either edge:

`i_m <= binom(N+3,m) - 2*binom(N+1,m-2) + binom(N-1,m-4)`  (for `m>=4`).

Substitute this in place of `binom(N+3,m)` to get:

`(1-lambda) * C_{lambda,rho} * (1+lambda)^{ceil(N/2)}`
`<= binom(N+3,m) - 2*binom(N+1,m-2) + binom(N-1,m-4)`.

This yields a sharper finite upper bound on feasible `N` at fixed `(m,lambda,rho)`.

## 7) Conditional singleton-interval criterion (from B's derivation)

Define lower endpoint

- `L1 = ceil(m - 4 + m/lambda)`,
- `L2 = max(0, 2*r_min - 1)`, where `r_min` is least integer with
  `(1+lambda)^r >= lambda/rho`,
- `L = max(L1, L2)`.

Let `A = (1-lambda)*((1+2*lambda)+(1+lambda)*rho)` and
`E = smallest even integer > L`.
Define

`g(N) = binom(N+3,m)/(1+lambda)^{N/2}` for even `N`.

If both conditions hold:

1. `((E+5)(E+4))/((E+5-m)(E+4-m)) <= 1+lambda`,
2. `A*(1+lambda)^{E/2} > binom(E+3,m)`,

then no feasible `N >= L+1` exists, so necessarily `N = L`.

This gives a concrete sufficient condition for `(m,lambda,rho) -> N`
without a full global separation theorem.
