# Round 32 (Web 5.2): Post-INSUFFICIENT Decisive Artifact

Context: You already produced an `INSUFFICIENT` diagnosis with the key identities:
- one-step decrement as expectation of a clipped variable,
- telescoped `S = E[min(K,m-2)]`,
- and the gap source: mu4 data does not obviously pin down near-mode endpoint mass terms in `Threshold`.

This round must produce a **decisive next artifact**, not another generic insufficiency.

## Authoritative packet (same as round31)
Use only these definitions/facts:
- Canonical tree class: `d_leaf<=1`, canonical `(leaf,support,u)` with `deg(support)=2`, `B=T-{leaf,support}` rooted at `u`.
- Rooted DP: `P=dp0[u]=prod_c(dp0[c]+dp1[c])`, `Q=dp1[u]=x prod_c dp0[c]`.
- `m` leftmost mode of `I(T)`, `lambda=i_{m-1}/i_m`.
- `TailDef = sum_{k=0}^{m-3}(m-2-k)p_k lambda^k`.
- `R1_tail2: TailDef <= p_{m-1} lambda^(m-1) + 2 p_m lambda^m`.
- `Threshold = (m-2) - (p_{m-1} lambda^(m-1)+2p_m lambda^m)/P(lambda)` and `R1_tail2 <=> S>=Threshold`.
- Mu4 class `C_{deg,mu,mu2,mu3,mu4}` with one-step LP-dual bound using `beta_i(F,m,lambda)`.
- Frontier facts: canonical `n<=23`: `931,596` checked, `0` failures, min gap `0.0002818864`.
- Collision facts: same `(deg(P),m,lambda,mu1..mu4)` has observed threshold-split count `0` (`n<=20`), and observed same-key/different-mu5 count `0` (`n<=21`).

## Required output (strict)
Return exactly one of:

A) `SUCCESS-ARTIFACT`
- Give one explicit **new canonical invariant** `J` (not just mu-order label), computable from packet objects, and a one-step lemma template:
  `d_m(A,F) >= L(mu1..mu4, J(A), F, m, lambda)`
  with clear dependency arrows showing how telescoping yields `S>=Threshold`.
- Must include a concrete formula for `J` (or a finite tuple), not prose only.

B) `BLOCKED-ARTIFACT`
- Give a canonical no-go statement that is strictly stronger than round31:
  prove from packet identities that any class using only `(deg,mu1..mu4)` cannot uniformly lower-bound the step-function expectation needed for `Threshold` control.
- Then specify the **minimal additional canonical invariant** as an explicit function/formula (not just “mu5”), and a one-step lemma template using it.

C) `INSUFFICIENT` is allowed only if you provide exactly one missing invariant as a formula AND explain why the collision facts above do not already imply it is redundant.

## Hard constraints
- No geometric tails, no truncation `p_k=0 for k>m`, no fabricated runtime outputs.
- No unrestricted perturbation arguments that are not mapped back into canonical DP symbols.
- If you propose moment-order extension, you must justify minimality against the packet fact: observed same-mu4-key/different-mu5 count is 0 through n<=21.

## Output format
1) `Canonical inputs used`
2) `Decisive artifact (formula-level)`
3) `Derivation`
4) `Binary verdict: SUCCESS-ARTIFACT or BLOCKED-ARTIFACT or INSUFFICIENT`
5) `Next step in one sentence`
