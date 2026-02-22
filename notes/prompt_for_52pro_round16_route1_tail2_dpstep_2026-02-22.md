# Round 16: Route-1 tail2 via exact one-child DP step (decisive)

Context: Erdos #993 proof-closure project, canonical degree-2 bridge setup.

## Current state (trusted)
- Target: `E_route1: exact_excess_D <= exact_slack_B`.
- Trusted bridge: `mu_P - (m-2) = exact_slack_B - exact_excess_D`.
- Exact identity:

```
lambda*P'(lambda) - (m-2)*P(lambda)
= -TailDef + p_{m-1}*lambda^(m-1) + 2*p_m*lambda^m
  + sum_{k>=m+1} (k-(m-2))*p_k*lambda^k.
```

- Candidate retained:

```
R1_tail2: TailDef <= p_{m-1}*lambda^(m-1) + 2*p_m*lambda^m.
```

- Empirical status: `R1_tail2` has `0/77,141` failures on canonical `d_leaf<=1`, `n<=20`.

## What already failed (do NOT reuse)
1) False obligation from round 14:
`p_{k+2}*lambda^(k+2) >= p_{k+1}*lambda^(k+1) + p_k*lambda^k`.
Counterexample: `n=5`, `g6=DQo`, `m=2`, `lambda=5/6`, `p0=1,p1=2,p2=0`.

2) False "obstruction" from round 15 (`Drift3`):
`(m-2-k)p_k lambda^k <= 2 p_{k+3} lambda^{k+3} - p_{k+2} lambda^{k+2} - p_{k+1} lambda^{k+1}`.
Fails massively (77,130/77,141 cases).

## Required method (new)
You must work from the exact one-child DP product step, not guessed recurrences.

Let the partial product before adding one child be
`A(x)=sum_i a_i x^i`, and child factor be
`F(x)=sum_j f_j x^j` (with nonnegative coefficients), so after one step:

```
A_new(x) = A(x) * F(x),
[a_new]_k = sum_{i+j=k} a_i f_j.
```

Construct the tail functional at level `m`:

```
Tail_m(A) := sum_{k=0}^{m-3} (m-2-k) a_k lambda^k.
```

Your job is to derive a valid one-step inequality of the form:

```
Tail_m(A*F) <= Phi_m(A,F;lambda)
```

that telescopes over child attachments to imply `R1_tail2` at the root, with `Phi_m` reducing exactly to
`p_{m-1} lambda^{m-1} + 2 p_m lambda^m` at the final step (plus only nonpositive disposable terms).

## Output format (strict)
1) `Candidate one-step inequality` (exact formula in symbols `a_i,f_j,lambda,m`).
2) `Telescoping derivation` to `R1_tail2` (fully explicit index manipulations).
3) `Closure chain` `R1_tail2 => mu_P>=m-2 => E_route1` (brief, exact).
4) `Witness-pack falsification table` for the proposed one-step inequality:
   - W1: `n=5`, `g6=DQo`
   - W2: n=8, g6 string = G?`@F_
   - W3: `n=20`, `g6=S???????????_?O?C??o?@_?@_??oFig?`
   For each witness print concrete numeric LHS and RHS.
5) `Binary verdict: SUCCESS or BLOCKED`.
6) If `BLOCKED`: provide exactly one minimal obstruction lemma needed at the one-step DP level.

## Hard constraints
- No untrusted redefinition of `exact_slack_B` / `exact_excess_D`.
- No geometric-series tail bound (`1/(1-lambda)` patterns).
- No assumption `p_k=0` for `k>m`.
- No reuse of already falsified obligations above.
- If any witness-pack row violates your proposed one-step inequality, verdict must be `BLOCKED`.
