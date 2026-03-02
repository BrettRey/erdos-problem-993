# Round 12, Prompt 1: STP2 Exhaustive Verification + Prefix E≽J at Larger n

**Target model:** Codex 5.3 (exhaustive computation)

---

## Context

We are proving prefix E ≽ J (ratio dominance for k < mode) at every support vertex of every tree. This implies unimodality of the IS sequence (Erdős Problem #993) when combined with the proved P3 (tail domination).

**Key Round 11 findings:**
- E≽J holds at ALL k through n≤22 (907M checks), but **FAILS globally at n=32** (d_15 = -3498, mode = 10). The failure is at k=15, well past the mode.
- **Prefix E≽J** (k < mode) still holds at n=32. This is the actual proof target.
- **STP2 condition**: r_{m+1}·q_n ≤ r_m·q_{n+1} for m > n ≥ 0, where q = E_t, r = J_t (factor polynomials from child subtree t). **0 failures** across 50,917 factor pairs (Round 11 Codex).
- The Karlin split: w_k = d_k(E_acc·g, J_acc·g) + d_k(x·E_acc·h, J_acc·g). First term provably ≥ 0.

## Available Code

Use the existing `indpoly.py` module for polynomial multiplication. Use `geng` at `/opt/homebrew/bin/geng` for tree enumeration:
```bash
geng -q <n> <n-1>:<n-1> -c   # generates all trees on n vertices
```

The DP at a rooted tree gives dp[v] = (E_v, J_v) where:
- E_v = exclude-root polynomial = product over children c of I_c
- J_v = include-root polynomial (before x-shift) = product over children c of E_c
- I_c = E_c + x·J_c (full IS polynomial of child subtree)

A support vertex is one adjacent to at least one leaf.

## Tasks

### Task 1: Exhaustive STP2 at All Incremental Steps (n ≤ 20)

At each support vertex with s ≥ 2 non-leaf children, process children one at a time. At each step t, extract the factor pair (q, r) = (E_t, J_t) from the child subtree rooted at child t.

Check STP2: for all 0 ≤ n < m ≤ max_degree, verify r_{m+1}·q_n ≤ r_m·q_{n+1}.

Count:
- Total factor pairs tested
- Total STP2 checks (individual (m,n) pairs)
- Any failures (report tree, rooting, child, (m,n) values)

This should be MORE thorough than the Round 11 scan (which tested 50,917 pairs). We want every factor pair at every incremental step at every support vertex.

### Task 2: Prefix E≽J at n = 23, 24, 25

For each tree on n vertices (via geng), for each support vertex r:
1. Compute E = dp[r][0], J = dp[r][1]
2. Compute mode = argmax of I = E + xJ
3. Check d_k = E_{k+1}·J_k - E_k·J_{k+1} ≥ 0 for k = 0, ..., mode-1 only (PREFIX)
4. Also check d_k for ALL k and report where failures occur (k vs mode)

Report:
- Total support vertices checked
- Prefix failures (should be 0)
- Full-range failures: how many, at what k relative to mode
- The smallest n where full-range E≽J first fails (if < 32)

Tree counts: n=23 has 1,903,890 trees, n=24 has 5,033,922, n=25 has 13,587,015. Use multiprocessing with geng's res/mod partitioning.

### Task 3: W-form Restricted to k < mode

At each incremental step (s ≥ 2 support vertices, n ≤ 18), compute:
- w_k for k < mode(I_new) — the PREFIX W-form
- w_k for k ≥ mode(I_new) — the TAIL W-form
- The Karlin term d_k(E_acc·g, J_acc·g) for each k
- The correction term d_k(x·E_acc·h, J_acc·g) for each k

Report:
- Min correction/Karlin ratio for k < mode vs k ≥ mode
- Does restricting to k < mode make the correction relatively smaller?
- Are there cases where w_k < 0 for k ≥ mode? (This would mean prefix-only is essential)

### Task 4: STP2 at the n=32 Failure Tree

Build the n=32 tree: support vertex r with 1 leaf, P_2 child (r--a--b), and T_{3,4} broom child (v with 3 children w_i, each w_i with 4 children x_{i,j}, each x with 1 leaf y).

At the support vertex r (s=2, ℓ=1):
1. Extract factor pairs (E_t, J_t) for both non-leaf children
2. Check STP2 for each factor
3. Compute the Karlin split at each incremental step
4. Identify which incremental step and which k produces d_15 = -3498
5. At that step: what is the Karlin term? What is the correction? Does the correction exceed the Karlin term?

## Deliverables

1. STP2 verification summary (total pairs, total checks, any failures)
2. Prefix E≽J verification at n=23-25 (failure count, failure locations relative to mode)
3. W-form prefix vs tail comparison
4. n=32 detailed Karlin split analysis
5. Any new extremal trees discovered (tightest prefix E≽J, tightest STP2)
