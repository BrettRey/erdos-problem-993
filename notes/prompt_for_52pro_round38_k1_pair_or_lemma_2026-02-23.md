# Round 38 (Web 5.2): K1 Pair-Or-Lemma Gate

Context:
- Round37 output ended at `INSUFFICIENT` with missing theorem:
  injectivity of `P -> (deg(P),mu1(P))` at fixed `(m,lambda)` on canonical DP image.
- We now force a concrete next artifact.

## Authoritative packet (unchanged)
Canonical definitions:
- canonical class `d_leaf<=1`, canonical `(leaf,support,u)`, `deg(support)=2`, `B=T-{leaf,support}` rooted at `u`
- `P=dp0[u]=prod_c f_c`, `f_c=dp0[c]+dp1[c]`
- `m` leftmost mode of `I(T)`, `lambda=i_{m-1}/i_m`
- `K1=(deg(P),m,lambda,mu1)` where `mu1=lambda P'(lambda)/P(lambda)`
- `Threshold`, `B_max`, `R1_tail2`, `S` as before

Known facts:
- `K0=(deg,m,lambda)` splits heavily (different `P` common)
- no observed `K4->P` split through canonical `n<=23`
- multiset uniqueness of child factors is false (same K4/same P/different factor multisets exists)

## Required output (binary)
Return exactly one:

A) `SUCCESS-ARTIFACT`
- Provide one explicit lemma `L1` toward `K1->P` injectivity that is formula-level and canonical.
- `L1` must be nontrivial (not restating injectivity) and must reduce the space of possible `P` from `K1` data.
- Then provide a 2-3 lemma dependency sketch showing how `L1` could lead to `K1->P` (or `K1->Threshold`).

B) `BLOCKED-ARTIFACT`
- Provide an explicit canonical pair `(T1,T2)` with same `K1` and different `P` (or different `Threshold`).
- Include concrete g6 strings and exact equality/inequality checks.

If neither is possible from packet data, output `INSUFFICIENT` and give a single sharper missing theorem than last round (must include a checkable formula).

## Hard constraints
- No unrestricted perturbation arguments detached from canonical DP.
- No multiset-uniqueness assumptions.
- No geometric tails or truncation assumptions.
- No fabricated runtime outputs.

## Output format
1) `Canonical inputs used`
2) `Artifact`
3) `Derivation`
4) `Binary verdict: SUCCESS-ARTIFACT or BLOCKED-ARTIFACT or INSUFFICIENT`
5) `Single sharper missing theorem`
