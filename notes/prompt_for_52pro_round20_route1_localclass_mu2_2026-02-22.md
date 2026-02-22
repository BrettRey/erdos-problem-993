# Round 20: Route-1 local class with second moment (decisive)

Context: Erdos #993 proof-closure project, canonical degree-2 bridge setup.

## Fixed facts
- Route-1 target: `E_route1: exact_excess_D <= exact_slack_B`.
- Trusted bridge: `mu_P - (m-2) = exact_slack_B - exact_excess_D`.
- Accepted identity:

```
lambda*P'(lambda) - (m-2)*P(lambda)
= -TailDef + p_{m-1}*lambda^(m-1) + 2*p_m*lambda^m
  + sum_{k>=m+1} (k-(m-2))*p_k*lambda^k.
```

- Candidate retained:

```
R1_tail2: TailDef <= p_{m-1}*lambda^(m-1) + 2*p_m*lambda^m
```
(empirically 0/77,141 failures on canonical `d_leaf<=1`, `n<=20`).

- Accepted one-step decrement identity:

```
Tail_m(A*F) = Tail_m(A)F(lambda)
 - sum_{r=0}^{m-3} (sum_{i=0}^{m-3-r} a_i lambda^i)(sum_{j>=r+1} f_j lambda^j).
```

- Established no-go: class `C_{deg,mu}` (using only one-step identity + nonnegativity + `deg(A)` + `mu_A(lambda)`) cannot certify `R1_tail2` on W3.

## New class to resolve
Analyze class `C_{deg,mu,mu2}` that may use only:
- one-step identity,
- nonnegativity,
- `deg(A)`,
- `mu_A(lambda)`,
- `mu2_A(lambda) := (sum_i i(i-1)a_i lambda^i)/A(lambda)`,
and no other A-shape data.

## Required task
Produce exactly one output:

A) `SUCCESS`: a non-tautological local lemma in `C_{deg,mu,mu2}` that is sufficient to imply `R1_tail2` (and hence Route-1), with full derivation.

or

B) `BLOCKED`: a formal no-go theorem for `C_{deg,mu,mu2}` on the witness pack, with explicit bound-vs-threshold gap.

## Mandatory witness pack (use these exact strings)
- W1 g6: DQo
- W2 g6: G?`@F_
- W3 g6: S???????????_?O?C??o?@_?@_??oFig?

For each witness, first print computed `(m, i_{m-1}, i_m, lambda)` from the g6 decode before any further numbers.

## Mandatory non-tautology gate
Any proposed lemma must not mention final-target symbols (`TailDef`, `R1_tail2`, `p_{m-1}`, `p_m`, `exact_*`) in its statement.

## Output format (strict)
1) `Class/lemma statement`
2) `Derivation`
3) `Witness-pack numeric table`
4) `Binary verdict: SUCCESS or BLOCKED`
5) `If BLOCKED: one minimal next lemma outside C_{deg,mu,mu2}`

## Hard constraints
- No geometric-series tail bound (`1/(1-lambda)` patterns).
- No `p_k=0 for k>m` assumption.
- No reuse of already rejected lemmas/classes.
