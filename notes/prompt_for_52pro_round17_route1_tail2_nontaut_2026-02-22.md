# Round 17: Route-1 tail2 with NON-TAUTOLOGICAL decrement bound

Context: Erdos #993 proof-closure project, canonical degree-2 bridge setup.

## Status
- `R1_tail2` remains empirically true on canonical `d_leaf<=1`, `n<=20` (`0/77,141` failures):

```
TailDef <= p_{m-1}*lambda^(m-1) + 2*p_m*lambda^m.
```

- The one-step identity from round 16 is accepted:

```
Tail_m(A*F)
= Tail_m(A)F(lambda)
  - sum_{r=0}^{m-3} (sum_{i=0}^{m-3-r} a_i lambda^i)(sum_{j>=r+1} f_j lambda^j).
```

- The round-16 obstruction lemma was rejected as tautological (algebraically equivalent to target after rearrangement).

## Required task
Work from the telescoped decrement representation and produce a **genuine, non-tautological** lower bound on decrement.

Let `A_0=1`, `A_t=A_{t-1}F_t`, final `P=A_d`, and define

```
S := sum_{t=1}^d sum_{r=0}^{m-3}
     (sum_{i=0}^{m-3-r} [A_{t-1}]_i lambda^i / A_{t-1}(lambda))
     (sum_{j>=r+1} [F_t]_j lambda^j / F_t(lambda)).
```

Then the telescoped identity is

```
Tail_m(P) = (m-2)P(lambda) - P(lambda)*S.
```

Your job is to provide a lemma `L` that yields a provable lower bound `S >= B` where:
- `L` is one-step DP-local,
- `B` is explicit and computable from step-local quantities,
- `L` does NOT mention `Tail_m(P)`, `R1_tail2`, `p_{m-1}`, `p_m`, or `E_route1`.

Then show how `S>=B` implies `R1_tail2`.

## Non-tautology gate (mandatory)
You may output `SUCCESS` only if all hold:
1. `L` contains no final-target symbols (`Tail_m(P)`, `R1_tail2`, `p_{m-1}`, `p_m`, `exact_*`).
2. Derivation of `S>=B` from `L` uses only one-step DP identities and nonnegativity.
3. The implication `B >= (m-2) - (p_{m-1} lambda^(m-1) + 2 p_m lambda^m)/P(lambda)` is proved separately and explicitly.
4. Witness-pack table (below) has no violated row for `L`.

If any gate fails, verdict must be `BLOCKED`.

## Witness-pack (must evaluate proposed `L` numerically)
- W1: `n=5`, `g6=DQo`
- W2: n=8, g6 string = G?`@F_
- W3: `n=20`, `g6=S???????????_?O?C??o?@_?@_??oFig?`

For each witness, print numeric LHS/RHS of `L`.

## Output format (strict)
1) `Local lemma L`
2) `Derive S>=B`
3) `Show B implies R1_tail2`
4) `Witness-pack table`
5) `Binary verdict: SUCCESS or BLOCKED`
6) If `BLOCKED`: exactly one minimal additional local lemma needed.

## Hard constraints
- No geometric-series tail bound (`1/(1-lambda)` patterns).
- No assumption `p_k=0` for `k>m`.
- No restatement-equivalent or rearranged-target "obstruction lemma".
