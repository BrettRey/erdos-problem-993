# Round 24: TP2 kernel route (instantiate-or-no-go)

Context: Erdos #993 closure project (trees, independence polynomial), canonical degree-2 bridge framework.

## Why this round exists
A generic TP2/Cauchy-Binet claim is valid for true kernel composition:

- local step: `x^+ = K x`, `K_{i,j}>=0`
- TP2 minors: `Delta_K(i1,i2;j1,j2) >= 0`
- composition: `M = L K`
- `Delta_M = sum_{m1<m2} Delta_L * Delta_K` (nonnegative sum)

But in our project we previously hit a mismatch:
- `prove_strong_c2_tp2_closure.py` shows the **shifted-pair cross determinant** route fails for 3-factor products.
- So this round must either:
  1) instantiate a **true kernel state** where TP2 applies directly to our target quantity, or
  2) prove a no-go for plausible kernel-state shapes and drop this route.

## Established state (trusted)
- Route-1 target: `E_route1: exact_excess_D <= exact_slack_B`.
- Bridge: `mu_P - (m-2) = exact_slack_B - exact_excess_D`.
- Candidate sufficient inequality retained:
  `R1_tail2: TailDef <= p_{m-1}*lambda^(m-1) + 2*p_m*lambda^m`.
- Telescoped identity:
  `Tail_m(P) = (m-2)P(lambda) - P(lambda)S`,
  where `S` is the sum of one-step decrement terms.
- Known no-go classes already done:
  `C_{deg,mu}`, `C_{deg,mu,mu2}`, `C_{deg,mu,mu2,mu3}`.

## Required deliverable (binary)
Return exactly one:

A) `SUCCESS`: explicit TP2-kernel instantiation that is non-circular and closes a target; or
B) `BLOCKED`: formal no-go for the plausible TP2-kernel shapes attempted.

## Strict success criteria (must satisfy all)
1. Give explicit state vector `x` and one-step matrix `K(F,lambda,m)` such that each child attachment is exactly `x^+ = K x`.
2. Define explicit 2x2 minors of `x` or `K` whose nonnegativity implies the target inequality (Route-1 chain or hard-side chain).
3. Use only true matrix-kernel TP2 closure (Cauchy-Binet minors), not shifted-pair substitutions.
4. Show the target follows by algebraic chain to either:
   - `R1_tail2`, then `mu_P>=m-2`, then `E_route1`; or
   - `E_hard` then hard-ratio.
5. Provide witness verification table with mandatory strings (below).

If any criterion fails, verdict must be `BLOCKED`.

## Lift constructions to test first (required)
Try these in order and report exact formulas.

### Lift A: phase/parity lift
- Use lifted index set `Ihat = I x {0,1}` with lex order.
- State `y in R^Ihat`.
- One-step map must be literal linear kernel:
  - `y^+ = Khat y`
  - `Khat = [[K00, K01],[K10, K11]]` with nonnegative entries.
- TP2 minors are ordinary 2x2 minors on `Khat` with the lex order.
- Composition must be literal matrix multiplication on `Khat` (no hidden projection between steps).

### Lift B: order-2 companion lift
- If the local update is effectively a 2-term recurrence, encode
  - `x_{t+1} = A x_t + B x_{t-1}`
  - `y_t = (x_t, x_{t-1})`
  - `y_{t+1} = Khat y_t`, `Khat = [[A, B],[I, 0]]`.
- Again, only ordinary TP2 minors on `Khat`, with literal composition.

### Go/no-go test for either lift
To count as valid TP2 closure route, all must hold:
1. `Khat` is explicit and independent of the current state values.
2. `Khat >= 0` entrywise.
3. Multi-step closure is exactly product of these `Khat` matrices.
4. Any truncation/projection is only final row/column deletion (order-preserving), not interposed each step.
5. The target shifted-pair inequality is recovered as ordinary minors or a monotone consequence of ordinary minors on the lifted state.

If not, return `BLOCKED` and identify whether failure is:
- phase mismatch,
- lag mismatch,
- or repeated projection/trimming mismatch.

## Mandatory witness set (exact strings)
- W1: 'DQo'
- W2: 'G?`@F_'
- W3-old: 'S???????????_?O?C??o?@_?@_??oFig?'
- W4-fail: 'S???????C?G?G?C?@??G??_?@??@?F~_?'

For each witness, print decoded `(m, i_{m-1}, i_m, lambda)` first.
Reject any variant where `_` is replaced by `*`.

## Hard constraints
- No geometric-series tail bounds (`1/(1-lambda)` patterns).
- No `p_k=0 for k>m` truncation.
- No tautological restatement of `R1_tail2`/`E_route1` as a "lemma".
- No TP2 closure claim unless it is genuine matrix-kernel composition on your explicit state.

## Output format (strict)
1) `Kernel/state definition`
2) `One-step TP2 minor family`
3) `Composition rule and closure chain`
4) `Witness table`
5) `Binary verdict: SUCCESS or BLOCKED`
6) `If BLOCKED: smallest next non-TP2 class`
