# Adjacent-Separation Blocker at Fixed `(m,lambda)` (2026-02-27)

## Setup
In the min-`u` canonical gated class, define
- `N=[x]P`,
- `rho=Q(lambda)/P(lambda)`,
- feasible set `R_{m,lambda}(N)` at fixed `m,lambda`.

## Key reduction
A nontrivial finite-window separation claim

`exists K>=1` such that
`R_{m,lambda}(N) ∩ R_{m,lambda}(N') = empty` whenever `0<|N-N'|<=K`

is equivalent to proving adjacent disjointness as the minimal requirement:

- necessity: any such `K>=1` implies
  `R_{m,lambda}(N) ∩ R_{m,lambda}(N+1) = empty` for all `N`;
- sufficiency: adjacent disjointness implies the claim with `K=1`.

So the exact blocker is:

> **Adjacent `(m,lambda,rho)` separation:** no two canonical gated trees with same
> `(m,lambda,rho)` and `|Delta N|=1`.

## Why current envelopes are insufficient
Current universal bounds provide:
- lower bounds on `N` from incidence and root-message arity,
- finite upper bounds on `N` from exponential-vs-binomial inequalities.

But the parity exponents used in these envelopes coincide for adjacent pair
`(N,N+1)=(2k-1,2k)`:
- `floor((N+1)/2)=k`,
- `ceil(N/2)=k`, and same for `N+1`.

Hence these inequalities are structurally blind to the critical adjacency
`2k-1 <-> 2k`; they do not by themselves prove any `K>=1` separation.

## Practical implication
For proof progress, the highest-value next theorem is exactly adjacent separation.
For computation, the highest-value next scan is exactly adjacent witness search:
- same `(m,lambda,rho)`,
- `|Delta N|=1`.
