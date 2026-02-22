# Round 19: Close local-class no-go and force a strictly stronger class

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

## Already rejected
- k+2 >= k+1+k recurrence (false massively)
- Drift3 (false massively)
- tail2 partial-product lemma M (false on 11,920/77,141 cases)
- tautological rearrangements equivalent to `R1_tail2`
- C_window-only bound insufficient on W3
- Markov-prefix (first moment only) insufficient on W3:
  max aggregate bound ≈ 2.9195 < threshold ≈ 4.3472

## Required task
Produce exactly one output:

A) `SUCCESS`: a new one-step local lemma outside the rejected classes that provably implies `R1_tail2`, or

B) `BLOCKED`: a formal no-go theorem for the class
`C_{deg,mu}` = bounds that use only
- one-step identity,
- nonnegativity,
- `deg(A)`,
- `mu_A(lambda)`
(and no other A-shape data),
showing this class cannot certify `R1_tail2` universally.

If you output A, your lemma must explicitly use extra local shape info beyond `(deg,mu)` and must pass witness checks.

## Witness pack (mandatory)
- W1: `n=5`, `g6=DQo`
- W2: n=8, g6 string = G?`@F_
- W3: `n=20`, `g6=S???????????_?O?C??o?@_?@_??oFig?`

For each witness, print numerical values relevant to your lemma/no-go.

## Output format (strict)
1) `Class statement / lemma statement`
2) `Derivation`
3) `Witness-pack numeric table`
4) `Binary verdict: SUCCESS or BLOCKED`
5) `If BLOCKED: one minimal next lemma outside C_{deg,mu}`

## Hard constraints
- No geometric-series tail bound (`1/(1-lambda)` patterns)
- No `p_k=0 for k>m` assumption
- No reuse of rejected lemmas/classes above
