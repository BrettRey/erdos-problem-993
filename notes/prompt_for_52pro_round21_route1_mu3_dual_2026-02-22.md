# Round 21: Route-1 via mu3-class LP-dual local bound (decisive)

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

- Candidate retained:

```
R1_tail2: TailDef <= p_{m-1}*lambda^(m-1) + 2*p_m*lambda^m
```

- Established no-go: class `C_{deg,mu}` is insufficient on W3.
- Established no-go: class `C_{deg,mu,mu2}` is insufficient on W3.
- Verified on corrected W3 tree (`g6=S???????????_?O?C??o?@_?@_??oFig?`):
  - threshold for `R1_tail2` implication: `≈ 4.3472138532`
  - class-optimal `B_max(C_{deg,mu,mu2,mu3})` on that witness/order family: `≈ 4.3923839874` (above threshold).

## Required task
Resolve class `C_{deg,mu,mu2,mu3}` with a concrete non-tautological local lemma.

At one step, with tilted distribution `w_i = a_i lambda^i / A(lambda)`,

```
d_m(A,F) = sum_i beta_i(F,m,lambda) * w_i
```

where `beta_i` is determined by local `F` and `(m,lambda)`.

Use LP-dual style: find coefficients `(c0,c1,c2,c3)` such that

```
c0 + c1*i + c2*i(i-1) + c3*i(i-1)(i-2) <= beta_i    for all i=0..deg(A).
```

Then derive the local bound

```
d_m(A,F) >= c0 + c1*mu_A + c2*mu2_A + c3*mu3_A.
```

You must show how this step bound telescopes and implies `R1_tail2` (or fails).

## Output requirements (strict)
1) `Local lemma L*` in class `C_{deg,mu,mu2,mu3}` (non-tautological).
2) `Dual certificate derivation` (explicit inequality on beta_i and moment reduction).
3) `Telescoping chain to R1_tail2`.
4) `Witness-pack table` with exact decoded `(m, i_{m-1}, i_m, lambda)` first, then bound-vs-threshold:
   - W1 g6: DQo
   - W2 g6: G?`@F_
   - W3 g6: S???????????_?O?C??o?@_?@_??oFig?
5) `Binary verdict: SUCCESS or BLOCKED`.
6) If BLOCKED: one minimal next lemma outside `C_{deg,mu,mu2,mu3}`.

## Hard constraints
- No geometric-series tail bounds (`1/(1-lambda)` patterns).
- No assumption `p_k=0 for k>m`.
- No tautological rearrangement equivalent to `R1_tail2`.
- No corrupted witness strings: must use exactly the three g6 strings above.
