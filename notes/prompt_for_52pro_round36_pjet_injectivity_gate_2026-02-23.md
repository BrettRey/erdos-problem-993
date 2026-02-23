# Round 36 (Web 5.2): P-Jet Injectivity Gate (No Multiset Uniqueness)

Context: Round35 correctly identified that child-factor multiset uniqueness is false.
This round targets the right object: injectivity at the polynomial `P` level.

## Authoritative packet
Canonical setup and definitions unchanged:
- canonical class: `d_leaf<=1`, canonical `(leaf,support,u)`, `deg(support)=2`, `B=T-{leaf,support}` rooted at `u`
- rooted DP: `P=dp0[u]=prod_c f_c`, `f_c=dp0[c]+dp1[c]`
- `m` leftmost mode of `I(T)`, `lambda=i_{m-1}/i_m`
- `K4=(deg(P), m, lambda, mu1..mu4)` at `lambda`
- `Threshold`, `B_max`, `R1_tail2`, `S` as previously defined

Moment/jet identity (fixed):
- `mu_r = lambda^r P^{(r)}(lambda)/P(lambda)`, `r=1..4`.
So at fixed `(m,lambda)`, K4 gives degree and 4-jet of `log P` at `lambda`.

## New hard facts
1) Full canonical frontier through `n<=23`: no `K4 -> P` split observed.
2) Multiset uniqueness is false (explicit canonical examples): same K4 and same P, different child-factor multisets.
   Therefore multiset uniqueness is not the theorem target.

## Required output (binary)
Return exactly one:

A) `SUCCESS-ARTIFACT`
- Provide a constructive canonical theorem route for `K4 -> P` (or `K4 -> Threshold`) at fixed `(m,lambda)`.
- Must explicitly avoid using child-factor multiset uniqueness.
- Give a lemma graph and exact dependency chain to `E_route1`.

B) `BLOCKED-ARTIFACT`
- Provide an explicit canonical pair `(T1,T2)` with same `K4` but different `P` (or different `Threshold`).
- Must include concrete g6 strings and exact values proving key equality and target inequality split.

If neither is achieved, output `INSUFFICIENT` and name the single missing theorem in one formal sentence.

## Hard constraints
- No unrestricted perturbation arguments detached from canonical DP.
- No geometric tails / no truncation assumptions.
- Must address why multiset-nonuniqueness does not block a possible `K4 -> P` theorem.

## Required output format
1) `Canonical inputs used`
2) `Artifact`
3) `Derivation`
4) `Binary verdict: SUCCESS-ARTIFACT or BLOCKED-ARTIFACT or INSUFFICIENT`
5) `Single missing theorem (if INSUFFICIENT)`
