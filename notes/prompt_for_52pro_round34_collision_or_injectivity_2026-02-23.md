# Round 34 (Web 5.2): Collision-Or-Injectivity Gate

Context: Round33 returned a plausible mechanism (cumulant-sum swap), but no explicit canonical instantiation.
This round must terminate that ambiguity.

## Authoritative packet (same definitions)
- Canonical class: `d_leaf<=1`, canonical `(leaf,support,u)` with `deg(support)=2`, `B=T-{leaf,support}` rooted at `u`.
- Rooted DP: `P=dp0[u]=prod_c f_c`, `f_c=dp0[c]+dp1[c]`.
- `m` leftmost mode of `I(T)`, `lambda=i_{m-1}/i_m`.
- `K4=(deg(P),m,lambda,mu1..mu4)` with factorial moments at `lambda`.
- `Threshold`, `B_max`, `R1_tail2`, `S` as before.

## New audited facts (treat as given)
1. Full canonical frontier through `n<=23`: no `K4->P` split observed.
   - `n<=21`: checked `175,722`, collisions `16,721`, `K4->P` splits `0`
   - `n=22`: checked `227,678`, collisions `16,271`, `K4->P` splits `0`
   - `n=23`: checked `528,196`, collisions `40,524`, `K4->P` splits `0`
2. Fixed `(m,lambda)=(7,1348/1357)` collision stress tests on rooted child-factor library:
   - unique factors from rooted trees up to size 10: `186`
   - exact pair-sum collision in vector `v(f)=(deg,kappa1..kappa4)`: none
   - exact triple-sum collision: none
   - size up to 12 (`874` factors): exact pair-sum collision: none
   - random 4-sum search (`800,000` trials): none found

## Required output (binary)
Return exactly one:

A) `SUCCESS-ARTIFACT`
- Provide a theorem-grade injectivity route for canonical child-factor sums:
  a precise statement (e.g., child-factor partition uniqueness at fixed `(m,lambda)`),
  and a lemma graph that would imply `K4 -> P` (or at least `K4 -> Threshold`).
- Must be formula-level, not analogy-level.

B) `BLOCKED-ARTIFACT`
- Provide one **explicit** collision witness in canonical symbols:
  two different multisets `{f_i}` and `{g_j}` of admissible child factors with
  `sum v(f_i)=sum v(g_j)` and `prod f_i != prod g_j`,
  including concrete polynomial tuples and exact equality checks.
- Generic “could occur at larger n” is not sufficient.

If you cannot produce A or B, output `INSUFFICIENT` and name the single theorem needed next.

## Hard constraints
- No fabricated computations.
- No unrestricted perturbation arguments detached from canonical child factors.
- No geometric tails or truncation assumptions.

## Output format
1) `Canonical inputs used`
2) `Artifact`
3) `Derivation`
4) `Binary verdict: SUCCESS-ARTIFACT or BLOCKED-ARTIFACT or INSUFFICIENT`
5) `Single next theorem target`
