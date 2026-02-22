# Round 15: Repair `R1_tail2` proof obligations (decisive)

Context: Erdos #993 proof-closure project, canonical degree-2 bridge setup.

## Target
Prove Route-1:

- `E_route1: exact_excess_D <= exact_slack_B`.
- Trusted bridge: `mu_P - (m-2) = exact_slack_B - exact_excess_D`.

## Fixed trusted setup
- `m = leftmost mode(I(T))`, `lambda = i_{m-1}(T)/i_m(T)`.
- `P = dp_B[u][0]` in canonical bridge decomposition.
- `TailDef := sum_{k=0}^{m-3} (m-2-k) * p_k * lambda^k`.
- Coefficient nonnegativity and `P(lambda)>0` hold.
- Exact identity:

```
lambda*P'(lambda) - (m-2)*P(lambda)
= -TailDef + p_{m-1}*lambda^(m-1) + 2*p_m*lambda^m
  + sum_{k>=m+1} (k-(m-2))*p_k*lambda^k.
```

## Candidate retained
`R1_tail2`:

```
TailDef <= p_{m-1}*lambda^(m-1) + 2*p_m*lambda^m.
```

Empirically on canonical `d_leaf<=1`, `n<=20`: `0/77,141` failures.

## What failed in round 14
Do **not** reuse this obligation; it is false almost everywhere:

```
For all 0<=k<=m-2: p_{k+2}*lambda^(k+2) >= p_{k+1}*lambda^(k+1) + p_k*lambda^k.
```

Counterexample: `n=5`, `g6=DQo`, `m=2`, `lambda=5/6`, `p0=1,p1=2,p2=0` gives
`0 < 2*(5/6)+1`.

## Required output (strict)
1) `Candidate statement` (either keep `R1_tail2` or replace with strictly stronger non-circular local candidate).
2) `Full closure chain` to `E_route1` (all steps explicit).
3) `Primitive one-step obligations (<=4)` that are actually plausible under DP recurrences.
4) `Binary verdict: SUCCESS or BLOCKED`.
5) If BLOCKED: provide exactly one minimal obstruction lemma needed.

## Mandatory validation gates (no exceptions)
You may output `SUCCESS` **only if all gates pass**:

Gate A (non-circularity):
- Show each primitive obligation is not a restatement of `R1_tail2` or `E_route1`.
- Give one-line dependency arrows: `obligations -> R1_tail2 -> mu_P>=m-2 -> E_route1`.

Gate B (explicit falsification test on known witness):
- Test every proposed primitive obligation on:
  - `n=5`, `g6=DQo`, `m=2`, `lambda=5/6`, `p0=1`, `p1=2`, `p2=0`.
- For each obligation, print LHS and RHS numerically on this witness.
- If any obligation fails on this witness, verdict must be `BLOCKED`.

Gate C (quantifier discipline):
- For each obligation, state exact quantifiers (`for all k in ...` or single-index).
- If your proof needs extra side conditions (e.g. `m>=3`, `p_{k+2}>0`), state them explicitly.
- If those side conditions are not guaranteed by trusted setup, verdict must be `BLOCKED`.

Gate D (proof locality):
- Your obligations must be one-step DP-local (single child update or single coefficient-slice inequality).
- If a step silently uses global tails or `1/(1-lambda)` style bounds, verdict must be `BLOCKED`.

## Hard constraints
- No untrusted redefinition of `exact_slack_B` / `exact_excess_D`.
- No global geometric-tail (`1/(1-lambda)`) argument.
- No assumption `p_k=0` for `k>m`.
- No obligations already falsified above.
