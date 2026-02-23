# Round 23: Beyond mu3 no-go (mu4 or non-moment local shape)

Context: Erdos #993 proof-closure project, canonical degree-2 bridge setup.

## Established state
- Route-1 target: `E_route1: exact_excess_D <= exact_slack_B`.
- Trusted bridge: `mu_P - (m-2) = exact_slack_B - exact_excess_D`.
- Candidate retained:

```
R1_tail2: TailDef <= p_{m-1}*lambda^(m-1) + 2*p_m*lambda^m.
```

- Accepted one-step decrement identity:

```
Tail_m(A*F) = Tail_m(A)F(lambda)
 - sum_{r=0}^{m-3} (sum_{i=0}^{m-3-r} a_i lambda^i)(sum_{j>=r+1} f_j lambda^j).
```

- No-go classes already established:
  - `C_{deg,mu}` insufficient
  - `C_{deg,mu,mu2}` insufficient
  - `C_{deg,mu,mu2,mu3}` insufficient globally (explicit failing witness below)

## Hard failing witness for mu3 class
- W4-fail (exact string): `S???????C?G?G?C?@??G??_?@??@?F~_?`
- canonical data: `m=7`, `i_6=9534`, `i_7=10068`, `lambda=9534/10068≈0.9469606675`
- mu3-class bound gap: `B_max - Threshold ≈ -0.0239342788` (fails)

## Required task
Produce exactly one of:

A) `SUCCESS`: a strictly stronger local class than mu3 (either `mu4` moment class or explicit non-moment local shape descriptor), with a complete non-tautological chain to `R1_tail2`.

or

B) `BLOCKED`: a formal no-go for your proposed stronger class on the mandatory witness set.

## If choosing A
- Define class precisely.
- Give one-step local lemma `L*` in that class.
- Show telescoping lower bound `S>=B` and prove `B>=Threshold`.
- Must pass all mandatory witnesses numerically.

## Mandatory witness set (exact strings only)
- W1: `DQo`
- W2: G?`@F_
- W3-old: `S???????????_?O?C??o?@_?@_??oFig?`
- W4-fail: `S???????C?G?G?C?@??G??_?@??@?F~_?`

For each witness, print decoded `(m, i_{m-1}, i_m, lambda)` first.
Reject any variant where `_` is replaced by `*`.

## Output format (strict)
1) `Class statement`
2) `Local lemma`
3) `Derivation to R1_tail2`
4) `Witness table`
5) `Binary verdict: SUCCESS or BLOCKED`
6) `If BLOCKED: one minimal next class`

## Hard constraints
- No geometric-series tail bounds (`1/(1-lambda)` patterns).
- No `p_k=0 for k>m` assumption.
- No tautological rearrangements equivalent to `R1_tail2`.
- Use exact g6 strings as listed.
