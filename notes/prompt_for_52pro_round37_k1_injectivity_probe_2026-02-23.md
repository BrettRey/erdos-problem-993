# Round 37 (Web 5.2): K1 Injectivity Probe

Context:
- Round36/35 established: child-factor multiset uniqueness is false and not the right target.
- Current valid target remains polynomial-level injectivity (`K* -> P`), then `K* -> Threshold`.

## Authoritative packet (same canonical definitions)
- canonical class: `d_leaf<=1`, canonical `(leaf,support,u)`, `deg(support)=2`, `B=T-{leaf,support}` rooted at `u`.
- rooted DP: `P=dp0[u]=prod_c f_c`, `f_c=dp0[c]+dp1[c]`.
- `m` leftmost mode of `I(T)`, `lambda=i_{m-1}/i_m`.
- `R1_tail2`, `Threshold`, `B_max`, `S` as before.

Moment identity:
- `mu_r = lambda^r P^{(r)}(lambda)/P(lambda)`.

## New audited facts
1) `K4 -> P` has no observed split on full canonical frontier `n<=23`.
2) Multiset uniqueness is false (explicit same-K4 same-P different-factor-multiset examples).
3) Additional empirical fact:
   On canonical `n<=21` (`175,722` checked), even
   `K1 := (deg(P), m, lambda, mu1)`
   has no observed split to different `P`.
4) `K0 := (deg(P), m, lambda)` splits heavily (many different `P`).

## Required output (binary)
Return exactly one:

A) `SUCCESS-ARTIFACT`
- Provide a constructive symbolic route that explains why adding `mu1` to `K0` could force injectivity on canonical DP-image (`K1 -> P`, or at least `K1 -> Threshold`).
- Must explicitly avoid multiset uniqueness assumptions.
- Include lemma graph and dependency chain to `E_route1`.

B) `BLOCKED-ARTIFACT`
- Provide an explicit canonical pair with same `K1` but different `P` (or different `Threshold`).
- Must include concrete g6 strings and exact key/target values.

If neither is achieved, output `INSUFFICIENT` and name the single missing theorem.

## Hard constraints
- No unrestricted perturbation arguments detached from canonical DP.
- No geometric tails / no truncation assumptions.
- Must use the counterexample insight correctly: factor multiset nonuniqueness does not imply P nonuniqueness.

## Output format
1) `Canonical inputs used`
2) `Artifact`
3) `Derivation`
4) `Binary verdict: SUCCESS-ARTIFACT or BLOCKED-ARTIFACT or INSUFFICIENT`
5) `Single missing theorem`
