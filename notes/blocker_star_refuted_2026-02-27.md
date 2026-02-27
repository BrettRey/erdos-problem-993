# Refutation of Strong Blocker `(*)` (2026-02-27)

Context:
- In `notes/adjacent_separation_blocker_2026-02-27.md`, blocker `(*)` was proposed as a
  strictly stronger sufficient condition for adjacent canonical separation:
  no two root-admissible component multisets with same `W_lambda` and size difference `+1`.

Result:
- `(*)` is false.

Explicit counterexample at fixed `lambda=7/9` (all components from `d_leaf<=1` rooted-tree library):

Use four rooted component types:
- singleton: `g6='@'`, root `0`, size `1`, ratio `r=G/F=9/16`
- edge: `g6='A_'`, root `0`, size `2`, ratio `r=16/23`
- size-8 type: `g6='G?`@f?'`, root `2`, size `8`, ratio `r=6553/9129`
- size-9 type: `g6='H?AE@`g'`, root `1`, size `9`, ratio `r=6553/9129`

Root-admissible multisets (at most one singleton):
- A = {singleton, edge, size-8}
- B = {singleton, edge, size-9}

Then
- `N(A)=11`, `N(B)=12` (adjacent),
- product of ratios is equal (the differing component has identical ratio),
- so `rho=lambda*prod r_i` is equal for A and B.

Exact artifact:
- `results/root_admissible_product_adjacent_counterexample_lambda_7_9.json`

Important scope note:
- This refutes the strong local multiset rigidity `(*)` only.
- It does **not** by itself refute canonical adjacent separation
  `R_{m,lambda}(N) ∩ R_{m,lambda}(N+1)=empty`, because these multisets are not yet shown to
  be globally realizable with the same derived canonical `(m,lambda)` in full trees.
