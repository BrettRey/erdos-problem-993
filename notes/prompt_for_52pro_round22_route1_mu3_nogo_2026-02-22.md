# Round 22: mu3-class no-go with explicit failing witness (decisive)

Context: Erdos #993 proof-closure project, canonical degree-2 bridge setup.

## Fixed facts
- Route-1 target: `E_route1: exact_excess_D <= exact_slack_B`.
- Trusted bridge: `mu_P - (m-2) = exact_slack_B - exact_excess_D`.
- Retained candidate:

```
R1_tail2: TailDef <= p_{m-1}*lambda^(m-1) + 2*p_m*lambda^m.
```

- Accepted one-step decrement identity:

```
Tail_m(A*F) = Tail_m(A)F(lambda)
 - sum_{r=0}^{m-3} (sum_{i=0}^{m-3-r} a_i lambda^i)(sum_{j>=r+1} f_j lambda^j).
```

- Established no-go classes already:
  - `C_{deg,mu}` insufficient.
  - `C_{deg,mu,mu2}` insufficient.

## New hard evidence (from full canonical n<=20 scan)
Class `C_{deg,mu,mu2,mu3}` is **not universally sufficient**.

Counterexample witness:
- g6: `S???????C?G?G?C?@??G??_?@??@?F~_?`
- canonical setup gives: `m=7`, `lambda≈0.9469606675`, root-child count `d=9`
- threshold for `R1_tail2` implication:
  `Threshold ≈ 4.3286125513`
- class-optimal bound over all 9! child-factor orders:
  `B_max(C_{deg,mu,mu2,mu3}) ≈ 4.3046782725`
- gap: `B_max - Threshold ≈ -0.0239342788` (fails).

So any `SUCCESS` claim for C_{deg,mu,mu2,mu3} is invalid globally.

## Required task
Produce exactly one output:

A) `BLOCKED`: a formal no-go statement for `C_{deg,mu,mu2,mu3}` using the witness above and the class definition.

or

B) `SUCCESS`: only if you provide a strictly stronger local class (outside mu3), and a full derivation to R1_tail2 **plus** witness checks including the failing witness above.

## If you choose A (expected)
You must provide:
1) clear class definition for `C_{deg,mu,mu2,mu3}`,
2) no-go derivation path (`best class bound < threshold => cannot certify R1_tail2`),
3) witness numeric table with decoded `(m, i_{m-1}, i_m, lambda)` and `(B_max, threshold, gap)`.

## If you choose B
You must name a strictly stronger local class (e.g. adding mu4 or non-moment local shape), and prove why it escapes the mu3 no-go.

## Mandatory witness set
- W1: `DQo`
- W2: G?`@F_
- W3-old: `S???????????_?O?C??o?@_?@_??oFig?`
- W4-new-fail: `S???????C?G?G?C?@??G??_?@??@?F~_?`

Use exact strings as above.

## Hard constraints
- No geometric-series tail bounds (`1/(1-lambda)` patterns).
- No `p_k=0 for k>m` assumption.
- No tautological rearrangements equivalent to `R1_tail2`.
