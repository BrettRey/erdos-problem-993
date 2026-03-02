# Round 14, Prompt 1: W-Form Ratio Profiling for s ≥ 2

**Target model:** Codex 5.3 (exhaustive computation)

---

## Context

We are proving E ≽ J (prefix ratio dominance) at support vertices of trees, which implies unimodality of independent set sequences (Erdős Problem #993).

The W-form identity decomposes the increment at each step of processing non-leaf children:

```
w_k = d_k(E_acc·g, J_acc·g) + (B_k·C_k - B_{k-1}·C_{k+1})
     = Karlin(k) + correction(k)
```

where Karlin ≥ 0 (provably) and correction can be negative.

**Known results:**
- w_k ≥ 0 at ALL steps, ALL orderings, n ≤ 20 (2.44M steps, 0 failures)
- **s=1 case PROVED** (handles 63% of support vertices)
- The remaining open case is **s ≥ 2** (support vertex with ≥ 2 non-leaf children)
- Star+star (root + 1 leaf + star(a) + star(b)) is the extremal family for k=2 rescue ratio
  - k=2 ratio → 1+√3 ≈ 2.73 > 1
  - k=1 ratio → 2 > 1 (from GPT Round 13 analysis)
- STP2 (J(m+1)·E(n) ≤ J(m)·E(n+1) for m > n): 0 failures through n=20

**Critical question:** Does the minimum W ratio for s ≥ 2 stay bounded away from 1?

From earlier profiling (n ≤ 20):
```
s=1: min W ratio 1.195 (→ 1 as n → ∞, but PROVED)
s=2: min W ratio 1.664
s=3: min W ratio 1.783
s=4: min W ratio 1.742
```

If the s ≥ 2 ratio stays above some constant > 1, a ratio-based proof is viable. If it → 1, we need an absolute margin argument.

## Available Code

```python
# Core DP
from indpoly import dp, dp_all, _polymul

# Tree generation
# /opt/homebrew/bin/geng n n-1:n-1 -c  produces trees on n vertices
# res/mod partitioning: geng n n-1:n-1 -c res/mod
```

The W-form computation at each incremental step:
```python
# At support vertex root with ell leaves and non-leaf children c_1, ..., c_s:
# Start: E_acc = (1+x)^ell (binomial), J_acc = [1]
# For each child t:
#   g = E_t, h = J_t
#   A = polymul(E_acc, g)
#   B = polymul(E_acc, h)
#   C = polymul(J_acc, g)
#   Karlin_k = A[k+1]*C[k] - A[k]*C[k+1]
#   correction_k = B[k]*C[k] - B[k-1]*C[k+1]
#   w_k = Karlin_k + correction_k
#   E_acc = polymul(E_acc, I_t)  where I_t[k] = g[k] + h[k-1] (with h[-1]=0)
#   J_acc = polymul(J_acc, g)
```

Report the "W ratio" = Karlin(k)/|correction(k)| when correction(k) < 0.
Report the "absolute W margin" = min w_k (should be ≥ 0).

## Tasks

### Task 1: Extend w_k ≥ 0 verification to n=21-24

Using multiprocessing with geng res/mod partitioning, verify w_k ≥ 0 at all s ≥ 2 support vertices.

For each n:
- Total support vertices with s ≥ 2
- Total (step, k) positions checked
- Number of failures (should be 0)
- Min w_k value
- The tree/rooting/step/k achieving min w_k

Use ALL orderings of children (or prove that a specific ordering suffices — earlier scan showed ALL orderings give w_k ≥ 0, so any ordering is fine).

### Task 2: Profile min W ratio by (n, s)

For each n from 9 to 24 and each s from 2 to max_s(n):
- Find the minimum W ratio (Karlin/|correction| when correction < 0) across ALL trees at that n and s
- Report the extremal tree (graph6 string), the critical k, and the step

**Critical output:** A table:
```
n    s=2 ratio    s=3 ratio    s=4 ratio    ...
9    ???          ???          ???
10   ???          ???          ???
...
24   ???          ???          ???
```

Does the s=2 ratio decrease monotonically? Does it have a positive limit?

### Task 3: Profile min absolute W margin by (n, s)

Same scan, but report min w_k (absolute value, not ratio):
```
n    s=2 min_w    s=3 min_w    s=4 min_w    ...
9    ???          ???          ???
...
24   ???          ???          ???
```

Is min w_k ≥ 1 for all s ≥ 2? Does it grow with n?

### Task 4: Check chain STP2

The "chain STP2" condition is: STP2 holds for BOTH (E_t, J_t) AND (I_t, E_t) at every vertex t.

For each tree-realizable pair (I_t, E_t, J_t) at n=9-22:
- Check STP2(E_t, J_t): J_t(m+1)·E_t(n) ≤ J_t(m)·E_t(n+1) for all m > n
- Check STP2(I_t, E_t): E_t(m+1)·I_t(n) ≤ E_t(m)·I_t(n+1) for all m > n

Report any failures. If 0 failures, report the tightest slack (smallest positive value of RHS-LHS).

### Task 5: Identify the s=2 extremal family

For each n=9-24, among s=2 support vertices, identify the tree achieving minimum W ratio. Describe:
- Tree structure (graph6 + human-readable description)
- The two non-leaf children: their subtree types, sizes, E and J polynomials
- The specific (step, k) at which min ratio occurs
- Compare with the star+star family: is the extremal always a double-star?

## Deliverables

1. w_k ≥ 0 verification table (n=21-24)
2. Min W ratio table by (n, s) — the MOST IMPORTANT output
3. Min absolute w_k table by (n, s)
4. Chain STP2 verification results
5. s=2 extremal family identification for each n
6. Assessment: does s ≥ 2 ratio → 1 or stay bounded away from 1?
