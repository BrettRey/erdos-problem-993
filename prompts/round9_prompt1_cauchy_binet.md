# Round 9, Prompt 1: Cauchy-Binet Expansion of the W Form

**Target model:** GPT 5.2 Pro (theory-heavy)

---

## Context

We are proving that the independent set sequence of every tree is unimodal (Erdos Problem #993). The proof reduces to showing E ≽ J (ratio dominance: E_{k+1}J_k ≥ E_k J_{k+1} for all k) at every support vertex of every tree, where E = dp[r][0] and J = dp[r][1]/x at a support vertex r.

At a support vertex r with ℓ leaf neighbours and non-leaf child subtrees T_1,...,T_s, the DP builds:
- E = (1+x)^ℓ · ∏ I_t  (where I_t = E_t + x·J_t is the IS poly of subtree t)
- J = ∏ E_t

We process non-leaf children incrementally:
- E^(0) = (1+x)^ℓ, J^(0) = [1]
- E^(t) = E^(t-1) · I_t = E^(t-1) · (E_t + x·J_t)
- J^(t) = J^(t-1) · E_t

Setting f = E^(t-1), g = J^(t-1), q = E_t, r = J_t:
- A = f·q
- B = f·r
- C = g·q = J^(t)
- E^(t) = A + xB

## What We Have Proved

1. **Karlin main term:** Δ_k(A, C) ≥ 0. (Because f ≽ g by inductive hypothesis, and TP2 closure under convolution with PF2 kernel q.)

2. **Id1 decomposition (algebraic identity):**
   ```
   J_k · Δ_k(E^(t), J^(t)) = J_k · Δ_k(A,J) + J_{k+1} · d_{k-1}(B,J) + B_k · c_k(J)
   ```
   where d_{k-1}(B,J) = B_k J_{k-1} - B_{k-1} J_k and c_k(J) = J_k² - J_{k-1}J_{k+1}.

3. **J^(t) is PF2** (proved: product of PF2 factors). So c_k(J) ≥ 0.

4. **s=1 case proved** (SCC + transitivity).

## What We Have Verified Computationally (0 failures)

- **W form (2-term):** J_k · Δ_k(A,J) + J_{k+1} · d_{k-1}(B,J) ≥ 0 (all trees n ≤ 20)
- **IF form (3-term):** same + B_k · c_k(J) ≥ 0 (all trees n ≤ 20)
- The W form holds even WITHOUT the curvature bonus. So the Karlin term alone suffices.
- Min W ratio (term1/|term2| when term2 < 0): 1.195 at n=20, s=1, pendant-star extremal.
- s=1 W ratio formula: ℓ(ℓ+1)/(ℓ-1)² → 1. Absolute margin: ℓ(3ℓ-1)/2 → ∞.

## Your Task: Cauchy-Binet Expansion

The W form says:
```
J_k · Δ_k(A, C) + J_{k+1} · d_{k-1}(B, C) ≥ 0
```
where A = f*q, B = f*r, C = g*q (convolutions).

**Step 1.** Write Δ_k(f*q, g*q) as a Cauchy-Binet sum:
```
Δ_k(A,C) = (f*q)_{k+1} · (g*q)_k - (f*q)_k · (g*q)_{k+1}
          = Σ_{i<j} (f_j g_i - f_i g_j) · K_{ij}(q, k)
```
where K_{ij}(q, k) is a 2×2 minor of the convolution kernel involving q.

Identify K_{ij}(q, k) explicitly. Since f ≽ g (inductive hypothesis), the weights w_{ij} = f_j g_i - f_i g_j ≥ 0 for i < j.

**Step 2.** Similarly write d_{k-1}(f*r, g*q) as:
```
d_{k-1}(B, C) = (f*r)_k · (g*q)_{k-1} - (f*r)_{k-1} · (g*q)_k
              = Σ_{i<j} (f_j g_i - f_i g_j) · L_{ij}(q, r, k)
```
where L_{ij} involves both kernels q and r.

Identify L_{ij}(q, r, k) explicitly.

**Step 3.** Substitute into W:
```
W = Σ_{i<j} w_{ij} · [J_k · K_{ij}(q,k) + J_{k+1} · L_{ij}(q,r,k)] ≥ 0
```

Since w_{ij} ≥ 0, it suffices to show the bracket is nonneg for each (i,j).

**Step 4.** Does the bracket simplify? Can you show that:
```
J_k · K_{ij}(q,k) + J_{k+1} · L_{ij}(q,r,k) ≥ 0
```
for all i < j, using:
- q is PF2 (log-concave, nonneg)
- r ≤ q coefficientwise
- SCC: (1+x)·(q + x·r) ≽ q, i.e., (1+x)I_t ≽ E_t at the child subtree
- J = g*q is PF2

If the pointwise bracket inequality holds, we have a proof of E ≽ J for all trees. If it fails, identify the failure cases — do they correspond to tree-realizable (q, r) pairs?

**Step 5.** If the pointwise inequality fails, consider the weighted version: can you show the sum is nonneg by pairing (i,j) contributions? (This would be a "Strassen coupling" or "transport" argument.)

## Key Constraints on (q, r) = (E_t, J_t)

These are DP values of a subtree rooted at the child:
- q = E_t is PF2 (LC, nonneg)
- r = J_t with r ≤ q coefficientwise (proved: J ≤ E)
- SCC holds: (1+x)I_t ≽ E_t where I_t = q + x·r
- E_t ≽ J_t (ratio dominance at child, by inductive hypothesis)
- Both q and r are products of PF2 factors (from deeper subtrees)

## Deliverables

1. Explicit formulas for K_{ij} and L_{ij}
2. Assessment: does the pointwise bracket inequality hold generically, or only for tree-realizable pairs?
3. If pointwise holds: proof. If not: quantify the failure and suggest coupling strategies.
4. Any simplification from the PF2 / SCC constraints on q, r.
