# Round 10, Prompt 1: W-Form Verification and s >= 2 Profiling

**Target model:** Codex 5.3 (exhaustive computation)

---

## Context

We are proving E >= J (ratio dominance) at support vertices of trees to establish unimodality of independent set sequences (Erdos Problem #993).

**CRITICAL UPDATE (n=28 counterexample):** Both SCC ((1+x)I >= E) and leaf-augmentation (I+xE >= E) are FALSE at n=28. The T_{3,4} broom (root v, 3 children w_i, each with 4 children x_{i,j}, each x with 1 leaf y) rooted at any x-vertex gives:
- SCC: d_14 = -2298
- Leaf-aug: g_14 = -334
- E >= J: HOLDS (all minors >= 0)
- W-form: HOLDS (s=1 case, both terms individually non-negative)

This means SCC is NOT a valid intermediate target. The proof goes directly via the W-form.

## The W-Form

At a support vertex r with ell leaf neighbors and non-leaf child subtrees T_1,...,T_s, process children incrementally:
- E^(0) = (1+x)^ell, J^(0) = [1]
- E^(t) = E^(t-1) * (E_t + x*J_t), J^(t) = J^(t-1) * E_t

At each step, setting f = E^(t-1), g = J^(t-1), q = E_t, r_t = J_t:
- A = f*q, B = f*r_t, C = g*q

The W-form says:
```
W_k = J_k * Delta_k(A, C) + J_{k+1} * d_{k-1}(B, C) >= 0
```

where J = J^(t) = C, Delta_k(A,C) = A_{k+1}*C_k - A_k*C_{k+1}, d_{k-1}(B,C) = B_k*C_{k-1} - B_{k-1}*C_k.

The s=1 case is PROVED (Karlin + transitivity). s >= 2 is OPEN.

## Your Tasks

### Task 1: Extend W-form verification to n = 22

The W-form has been verified 0 failures at n <= 20 (all trees, all support vertices, every incremental step). Extend this to n = 21 and n = 22.

Use `/opt/homebrew/bin/geng -c <n> <n-1>:<n-1> -q` to enumerate trees.

For each tree:
1. Find all support vertices (vertices adjacent to at least one leaf)
2. For each support vertex, do the incremental product build
3. At each step t, compute W_k for all k
4. Record any failures (W_k < 0)

Report: total trees, total support vertices, total incremental steps, total W_k checks, any failures.

Use `from indpoly import _polymul, _polyadd` for polynomial operations. Use the same DP code pattern as in `verify_gpt_counterexample.py`.

### Task 2: Profile s >= 2 extremals

For all trees n <= 20 (or n <= 22 if feasible), at each support vertex with s >= 2 non-leaf children:

1. Record the minimum W margin: min_k (W_term1 / |W_term2|) when W_term2 < 0
2. Record the absolute margin: min_k (W_term1 + W_term2)
3. Record the tree structure (degree sequence or Prufer code)
4. Group by s value

Find:
- The tree with the SMALLEST W ratio for each s = 2, 3, 4, ...
- The asymptotic behavior of the extremal ratio as s grows
- Whether the margin grows or shrinks with n

### Task 3: Factor-level checks at n = 28

For the T_{3,4} broom (n=28) and T_{3,3} broom (n=22, borderline):

At each support vertex, for each non-leaf child subtree (E_t, J_t):
1. Check: Is E_t LC? (c_k(E_t) >= 0 for all k)
2. Check: Does SCC hold for the factor? ((1+x)I_t >= E_t)
3. Check: Does E_t >= J_t?
4. Check: Does leaf-aug hold for the factor? (I_t + x*E_t >= E_t)

This tells us whether the FACTOR-level invariants hold even when the OVERALL tree invariants fail.

### Task 4: Check E >= J at ALL n = 28 broom rootings

We know E >= J holds at x-vertex rootings of T_{3,4}. Check it at ALL rootings (v, w, x, y vertices). If it fails at any non-support rooting, record which.

## Deliverables

1. W-form verification results for n = 21-22 (pass/fail, counts)
2. s >= 2 extremal families and margin profile (table by s and n)
3. Factor-level invariant check for T_{3,4} and T_{3,3}
4. E >= J status at all rootings of T_{3,4}
