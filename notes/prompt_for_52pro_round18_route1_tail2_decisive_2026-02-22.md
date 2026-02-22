# Round 18: Route-1 tail2 decisive pass (non-tautological or no-go)

Context: Erdos #993 proof-closure project, canonical degree-2 bridge setup.

## Fixed facts
- Target: `E_route1: exact_excess_D <= exact_slack_B`.
- Trusted bridge: `mu_P - (m-2) = exact_slack_B - exact_excess_D`.
- Accepted identity:

```
lambda*P'(lambda) - (m-2)*P(lambda)
= -TailDef + p_{m-1}*lambda^(m-1) + 2*p_m*lambda^m
  + sum_{k>=m+1} (k-(m-2))*p_k*lambda^k.
```

- Candidate retained (empirically 0/77,141 failures on canonical `d_leaf<=1`, `n<=20`):

```
R1_tail2: TailDef <= p_{m-1}*lambda^(m-1) + 2*p_m*lambda^m.
```

- Accepted one-step decrement identity:

```
Tail_m(A*F) = Tail_m(A)F(lambda)
 - sum_{r=0}^{m-3} (sum_{i=0}^{m-3-r} a_i lambda^i)(sum_{j>=r+1} f_j lambda^j).
```

## Already rejected (must not reuse)
1) `p_{k+2}lambda^(k+2) >= p_{k+1}lambda^(k+1)+p_k lambda^k` (false massively).
2) `Drift3` obstruction (false massively).
3) Tail2 partial-product lemma `M`: `sum_{i>=t} a_i lambda^i <= a_t lambda^t + a_{t+1} lambda^{t+1}`
   (fails: 11,920/77,141 cases; 12,295/106,706 partial products).
4) Any rearranged-target obstruction equivalent to `R1_tail2` (tautological).

## Required task
Produce exactly one of:

A) `SUCCESS`: a genuinely non-tautological local lemma `L*` (one-step DP level) that implies a lower bound `S>=B` and then implies `R1_tail2`.

B) `BLOCKED`: a formal no-go statement for the current local approach class, with explicit class definition and counterexample mechanism.

### If attempting A (`SUCCESS`)
- `L*` must avoid final-target symbols: no `Tail_m(P)`, no `R1_tail2`, no `p_{m-1},p_m`, no `exact_*`.
- `L*` must be step-local in `(A,F,m,lambda)` and derive `S>=B` via telescoping.
- You must prove separately that `B >= (m-2) - (p_{m-1}lambda^(m-1)+2p_m lambda^m)/P(lambda)`.

### If attempting B (`BLOCKED`)
- Define the local class explicitly (e.g., bounds depending only on degree-window or bounded prefix/tail summaries).
- Show why this class cannot certify `R1_tail2` universally.
- Give one concrete witness from the pack below where the class bound is insufficient.

## Mandatory witness-pack evaluation
Evaluate your proposed local lemma/bound numerically on:
- W1: `n=5`, `g6=DQo`
- W2: n=8, g6 string = G?`@F_
- W3: `n=20`, `g6=S???????????_?O?C??o?@_?@_??oFig?`

If any witness violates your proposed lemma, verdict must be `BLOCKED`.

## Output format (strict)
1) `Local lemma or no-go class statement`
2) `Derivation`
3) `Witness-pack table`
4) `Binary verdict: SUCCESS or BLOCKED`
5) `If BLOCKED: one minimal next lemma outside the rejected class`

## Hard constraints
- No geometric-series tail bound (`1/(1-lambda)` patterns).
- No assumption `p_k=0` for `k>m`.
- No already-rejected lemmas above.
