# Round 35 (Web 5.2): Retarget After Multiset-Uniqueness Counterexample

Context: your previous INSUFFICIENT proposed this theorem target:
"child-factor multiset is uniquely determined by cumulative vector sum v=(deg,kappa1..kappa4)."
That target is now empirically false in canonical data.

## New hard facts (authoritative)
Within canonical class (`d_leaf<=1`, canonical bridge), there are explicit same-K4 same-P but different child-factor multisets:

Example 1 (n=12):
- g6_1: `K??CA?_CEOF_`
- g6_2: `K??CA?_cBOW_`
- shared: `deg(P)=5`, `m=4`, `lambda=124/147`,
  `P=(1,9,29,43,30,8)`
- different factor multisets:
  - `{(1,2), (1,7,15,13,4)}`
  - `{(1,1), (1,8,21,22,8)}`

Example 2 (n=14):
- g6_1: `M???C@?G?_@_q?K`?`
- g6_2: `M???C@?G?_@__`]??`
- shared: `deg(P)=6`, `m=4`, `lambda=45/73`,
  `P=(1,11,46,95,102,53,10)`
- different factor multisets:
  - `{(1,2), (1,9,28,39,24,5)}`
  - `{(1,1), (1,10,36,59,43,10)}`

Therefore, "multiset uniqueness from v-sum" is not the right theorem target.

## Keep these canonical definitions
- Canonical class, DP factorization, `K4=(deg(P),m,lambda,mu1..mu4)`,
- `Threshold`, `B_max`, `R1_tail2`, and accepted chain to `E_route1`.

## Required output (binary)
Return exactly one:

A) `SUCCESS-ARTIFACT`
- Provide a constructive symbolic route that targets `K4 -> P` (or at least `K4 -> Threshold`) **without** assuming child-factor multiset uniqueness.
- Must include a concrete lemma graph and explicit formulas, not analogy.

B) `BLOCKED-ARTIFACT`
- Provide an explicit canonical pair with same `K4` but different `P` (or different `Threshold`) if you can derive one.
- If no such pair can be produced, you must not claim BLOCKED.

Else output `INSUFFICIENT` with one single next theorem target.

## Hard constraints
- Do not reuse/refute by unrestricted perturbations.
- No geometric tails, no truncation assumptions.
- Must address the counterexamples above explicitly.

## Output format
1) `Canonical inputs used`
2) `Artifact`
3) `Derivation`
4) `Binary verdict: SUCCESS-ARTIFACT or BLOCKED-ARTIFACT or INSUFFICIENT`
5) `Single next theorem target`
