# Round 13, Prompt 1: Karlin Rescue Ratio Convergence

**Target model:** Codex 5.3 (exhaustive computation)

---

## Context

At each s≥2 support vertex, the incremental W-form gives:

```
C_k · w_k = Term1(Karlin, ≥0) + Term2(can be neg) + Term3(curvature, ≥0)
```

where A = E_acc·g, B = E_acc·h, C = J_acc·g, and:
- Term1 = C_k · d_k(A, C) (Karlin, provably ≥ 0)
- Term2 = C_{k+1} · d_{k-1}(B, C)
- Term3 = B_k · c_k(C) (curvature bonus, ≥ 0 since C is LC)

The shifted ratio condition ⋆ (B_k·C_k ≥ B_{k-1}·C_{k+1}, i.e. Term2+Term3 ≥ 0) FAILS at 3.4% of steps. But the full identity (Term1 + Term2 + Term3 ≥ 0) holds with 0 failures through n=18 (2.44M steps).

**Key finding**: the minimum Karlin rescue ratio (Term1/|Term2+Term3| when Term2+Term3 < 0) DECREASES with n:

```
n=9:  10.000
n=10:  5.625
n=11:  4.667
n=12:  4.684
n=13:  4.214
n=14:  3.947
n=15:  3.714
n=16:  3.556
n=17:  3.448
n=18:  3.333
n=19:  3.290
n=20:  3.196
```

The extremal is ALWAYS s=2, ell=1, k=2 (root with 1 leaf, 2 non-leaf children, failure at coefficient index 2).

**Questions**: Does this ratio converge to 3? Or does it keep decreasing past 3? If it converges to 3, the proof needs Term1 ≥ 3·|Term2+Term3|. If it goes to 1, the proof is possible but tight. If it goes below 1 at some n, the full identity fails.

## Available Code

Use `indpoly.py` for polynomial operations. Use `/opt/homebrew/bin/geng` for tree generation.

DP at rooted tree:
```python
# E_v = dp[root][0] = exclude-root poly
# J_v = dp[root][1] / x = include-root poly (no x-shift)
# I_v = E_v + x·J_v
```

Support vertex: root adjacent to ≥1 leaf. Children split into ell leaves and s≥2 non-leaf children.

Incremental processing: start with E_acc = (1+x)^ell, J_acc = [1]. For each non-leaf child t:
```
g = E_t, h = J_t
E_new = E_acc · (g + x·h) = E_acc · I_t
J_new = J_acc · g
```

At step ≥ 2, compute:
```
A = E_acc · g,  B = E_acc · h,  C = J_acc · g
Term1 = C_k · (A_{k+1}·C_k - A_k·C_{k+1})
Term2 = C_{k+1} · (B_k·C_{k-1} - B_{k-1}·C_k)
Term3 = B_k · (C_k² - C_{k-1}·C_{k+1})
```

Check for k = 1, ..., mode-1 only (prefix).

## Tasks

### Task 1: Extend min ratio to n=21-22

Continue the table above. n=21 has 2.1M trees, n=22 has 5.6M. Use multiprocessing with geng's res/mod partitioning for parallelism.

Report:
- Min rescue ratio at each n
- Exact fraction (numerator/denominator)
- The extremal tree (graph6 string)
- The child subtree sizes at the extremal

### Task 2: Identify the extremal family

The extremal is always s=2, ell=1, k=2. What are the two non-leaf children?

Hypothesis: one child is a path P_{n-k} and the other is also a path or caterpillar. The "worst" configuration might be a pendant star with a long path attached.

For each n=9..22, identify:
- Tree structure (adjacency list or description)
- The two child subtrees (type, size)
- E_acc, g, h polynomials at the critical step
- The exact values of Term1, Term2, Term3

### Task 3: Closed-form ratio for the extremal family

If the extremal family is identifiable (e.g., a specific parameterized construction), compute the ratio as a function of the parameter (path length, arm count, etc.).

Check whether the ratio is monotone decreasing in the parameter and whether it converges to a limit. If the limit is L:
- L > 1: the proof works with the 3-term identity
- L = 1: borderline (may need absolute margin argument)
- L < 1: the 3-term identity is insufficient

### Task 4: Run w_k ≥ 0 verification at n=19-22

The full check (all 3 terms) found 0 failures through n=18. Extend to n=22 using multiprocessing.

For each n, at each support vertex with s≥2:
- Process children in natural order (no need to try all orderings — earlier scan showed ALL orderings work)
- At each step ≥ 2, for each k < mode, compute total = Term1 + Term2 + Term3
- Report any case where total < 0

Also report:
- Min w_k value (should be ≥ 0)
- The tree/rooting achieving min w_k
- Whether min w_k = 0 or min w_k > 0 for n ≥ 13

## Deliverables

1. Min rescue ratio table extended to n=22 with exact fractions
2. Extremal family identification (tree structure at each n)
3. Closed-form ratio formula (if identifiable) and limiting value
4. w_k ≥ 0 verification at n=19-22 (pass/fail count)
