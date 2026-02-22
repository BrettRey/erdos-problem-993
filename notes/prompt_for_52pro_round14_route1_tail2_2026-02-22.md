# Round 14: Route-1 via Tail Certificate v2 (decisive pass)

Context: Erdos #993 proof-closure project, canonical degree-2 bridge setup.

## Objective
Prove `E_route1` non-circularly:

- `E_route1: exact_excess_D <= exact_slack_B`.
- Trusted identity: `mu_P - (m-2) = exact_slack_B - exact_excess_D`.

You must work only in packet symbols and trusted identities; do **not** introduce untrusted definitions of `exact_excess_D` or `exact_slack_B`.

## Fixed setup (trusted)
- `m = leftmost mode(I(T))`, `lambda = i_{m-1}(T)/i_m(T)`.
- `P = dp_B[u][0]` in canonical bridge decomposition.
- `TailDef := sum_{k=0}^{m-3} (m-2-k) * p_k * lambda^k`.
- Coefficients are nonnegative (`p_k >= 0`, `lambda > 0`, `P(lambda) > 0`).

Exact algebraic identity you may use:

```
lambda*P'(lambda) - (m-2)*P(lambda)
= -TailDef + p_{m-1}*lambda^(m-1) + 2*p_m*lambda^m
  + sum_{k>=m+1} (k-(m-2))*p_k*lambda^k.
```

(So no `p_k=0 for k>m` assumption is allowed.)

## New candidate to attack

`R1_tail2`:

```
TailDef <= p_{m-1}*lambda^(m-1) + 2*p_m*lambda^m.
```

Empirical status on canonical `d_leaf<=1`, `n<=20`: `0/77,141` failures.

## Requirements
1. Either produce a **complete symbolic derivation**
   `R1_tail2 => mu_P(lambda)>=m-2 => E_route1`,
   with every step explicit and no hidden assumptions,
   **and** give primitive one-step DP obligations sufficient to prove `R1_tail2`.

2. Or produce a **small no-go theorem**:
   show why `R1_tail2` cannot be proved from allowed local data/one-step obligations under these constraints.

3. You must avoid:
- any definition-equivalent restatement of `E_route1` as a "new lemma",
- global geometric-series tail bounds (`1/(1-lambda)` patterns),
- assuming `p_k=0` above `m`.

## Output format (strict)
- `1) Candidate statement`
- `2) Full closure chain`
- `3) Primitive one-step obligations (<=4)`
- `4) Binary verdict: SUCCESS or BLOCKED`
- If BLOCKED, provide exactly one minimal obstruction lemma needed.
