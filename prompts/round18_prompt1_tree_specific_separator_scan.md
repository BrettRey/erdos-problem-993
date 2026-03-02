# Round 18, Prompt 1: Tree-Specific Separator Scan (Boundary-Correct)

**Target model:** Codex 5.3 (exhaustive computation)

---

## Why this round

Round 17 established (with boundary-correct CB indexing) that:

- `X_k >= 0` is false.
- Pairwise symmetric control (`S_{i,j}(k) >= 0`, `F(i,j) >= 0`) is false.
- `STP2(I,E)` holds at all rootings up to `n <= 18`.
- The real proof target is holistic: `D_k + X_k >= 0`.

A key ambiguity remains: which **tree-specific** constraints are doing the work (beyond abstract LC + STP2 axioms, which are insufficient).

---

## Conventions (must match exactly)

- Enumerate trees exhaustively with `networkx.nonisomorphic_trees(n)`.
- Support-step processing:
  - choose support root `r`,
  - absorb all leaf-neighbors first (initial `(1+x)^ell`),
  - then process non-leaf neighbors in increasing rooted-subtree size.
- Prefix at step `t`: `k < mode(I_new)` where mode is the **smallest** maximizing index.
- CB index set must include boundary column/index `-1` so decompositions are exact.

---

## Task 1: Reproduce boundary-correct Round 17 core numbers

Reproduce (or tightly confirm) the following on your own pipeline:

1. `X_k < 0` frequency for `n <= 18` (prefix regime).
2. `S_{i,j}(k) < 0` frequency for `n <= 15` with boundary-correct indexing.
3. Exactness check for split: `X_k = X_k^{pure} + X_k^{corr}` (0 mismatches).
4. `STP2(I,E)` all-rootings check through `n <= 18`.

If any mismatch with the posted Round 17 summary appears, report exact source (index range, prefix filter, step filter, or mode convention).

---

## Task 2: Test a candidate missing invariant on tree-realizable pairs

Define adjacent minor / P2 quantity for a pair `(P,Q)`:

`w_m(P,Q) := P(m)Q(m-1) - P(m-1)Q(m)`.

For rooted trees up to `n <= 22` (or max feasible exhaustive):

1. Scan all rootings and all valid `m` for `(I_v, E_v)`.
2. Report failures of `w_m(I_v,E_v) >= 0`.
3. If no failures, report global minimum margin and extremal witness.

Do the same for child factors `(I_c, E_c)` encountered in support-step updates.

Goal: determine whether `w >= 0` is genuinely tree-universal (empirically), and whether it separates tree-realizable data from known synthetic counterexamples.

---

## Task 3: Diagonal-Abel decomposition diagnostics for `X_k`

For each support-step prefix instance `(A,B,P,Q,k)` define diagonal index `s=i+j`:

- `u_i^{(s)} := A(i)B(s-i)`
- `Delta_{i,s-i}(A,B) = u_i^{(s)} - u_{i-1}^{(s)}`
- `W_i^{(s)} := P(k-i)Q(k-s+i)`
- `D_i^{(s)} := W_i^{(s)} - W_{i+1}^{(s)} = Delta_{k-i,k-s+i}(P,Q)`

Compute per diagonal:

- `S_k^{(s)} := sum_i Delta_{i,s-i}(A,B) W_i^{(s)} = sum_i u_i^{(s)} D_i^{(s)}`
- Split by midpoint `i <= floor(s/2)` vs `i >= ceil(s/2)`:
  - left/right mass of `u`,
  - left/right positive/negative mass of `D`,
  - sign-change profile of `D`.

Produce statistics specifically on cases with `X_k < 0` and on the worst 100 deficits.

Goal: isolate a numeric monotonicity/alignment pattern that is true for trees and could support a proof.

---

## Task 4: Holistic rescue metrics for `D_k + X_k`

For prefix support-step instances up to `n <= 18`:

1. Compute distribution of `D_k/|X_k|` on `X_k<0`.
2. Report min ratio and extremal witness.
3. Compute by-bucket profiles versus:
   - step `t`,
   - subtree-size gap,
   - mode distance `mode(I_new)-k`,
   - whether `X_k^{pure}<0` / `X_k^{corr}<0`.

Goal: identify strongest stable lower-bound shape (`D_k >= c*|X_k|` globally? by regime?).

---

## Output format

Provide:

1. A compact summary table of all scans.
2. Top-10 extremal witnesses (graph6, root/support, step, `k`, all relevant values).
3. A machine-readable JSON file with all aggregate statistics and extremals.
4. A short conclusion: which candidate invariant(s) look truly tree-specific and promising for proof.

