# Round 9, Prompt 3: Tree-Realizable Cone and SCC Product Closure

**Target model:** Codex 5.3 (computational + algebraic, exhaustive search strengths)

---

## Context

We are proving E ≽ J at support vertices of trees. The proof proceeds by induction: process non-leaf children one at a time, maintaining E^(t) ≽ J^(t).

At each step:
- E^(t) = E^(t-1) · (E_t + x·J_t)
- J^(t) = J^(t-1) · E_t

where (E_t, J_t) are the DP values of subtree t. The "factor" is the pair (E_t, J_t).

## The Problem

**Generic product closure is FALSE.** Example:
- E^(t-1) = J^(t-1) = [1,1], E_t = [1,1,1], J_t = [1]
- All 6 abstract constraints hold (PF2, J≤E, SCC, LC, E≽J at old, SCC at child)
- But Δ_2(E^(t), J^(t)) = -1 < 0. FAILURE.
- The issue: [1,1,1] is NOT a tree-realizable E_t.

**Tree-realizable pairs never fail.** Exhaustive test: 3.52M tree-realized cross-products at n ≤ 20, 0 violations.

## Your Task: Characterize the Tree-Realizable Cone

### Q1: What additional constraint distinguishes tree-realizable (E_t, J_t) from arbitrary PF2 pairs?

A tree-realizable pair (E_t, J_t) comes from rooting a subtree at vertex t. Then:
- E_t = ∏_{c ∈ children(t)} I_c  where I_c = E_c + x·J_c
- J_t = ∏_{c ∈ children(t)} E_c

Base case (leaf): E_t = [1], J_t = [1].
Single child c: E_t = E_c + x·J_c = I_c, J_t = E_c.
Two children: E_t = I_{c1} · I_{c2}, J_t = E_{c1} · E_{c2}.

The constraint is structural: E_t and J_t share factors (the E_c terms appear in both, via I_c = E_c + x·J_c in E_t and E_c directly in J_t).

### Q2: Enumerate tree-realizable (E_t, J_t) pairs for small n

For subtrees of size 1, 2, 3, ..., 8, enumerate ALL tree-realizable (E_t, J_t) pairs (up to graph isomorphism). For each pair, record:
- E_t (coefficient list)
- J_t (coefficient list)
- Which tree/rooting produced it
- Whether E_t ≽ J_t (should always hold)
- The "SCC ratio" min{Δ_k((1+x)I_t, E_t)} / max{|negative terms|}

### Q3: What property of the cone makes the correction term nonneg?

The synthetic counterexample has E_t = [1,1,1], J_t = [1]. This gives J_t/E_t ratios (1, 1/1, 1/1) that jump between indices — the ratio is not monotone.

Tree-realizable pairs always have E_t ≽ J_t (ratio dominance). Do they also satisfy:

(a) **Interlacing?** Are the roots of E_t and J_t interlaced?

(b) **Coefficient ratio monotonicity?** Is J_t,k / E_t,k nonincreasing?

(c) **Higher-order TP conditions?** E.g., is the 2×n matrix [E_t; J_t] totally positive (all 2×2 minors nonneg)? This is exactly E_t ≽ J_t.

(d) **SCC tightness?** How close is (1+x)I_t ≽ E_t to being tight? Does the margin predict the correction term sign?

### Q4: Can the W-form inequality be proved using only tree-realizable constraints?

The W form says:
```
J_k · Δ_k(A, C) + J_{k+1} · d_{k-1}(B, C) ≥ 0
```
where A = f·q, B = f·r, C = g·q (f = E_old, g = J_old, q = E_t, r = J_t).

Since f ≽ g (inductive hypothesis), we can write:
```
W = Σ_{i<j} (f_j g_i - f_i g_j) · [J_k · K_{ij} + J_{k+1} · L_{ij}]
```

The bracket depends on q, r, and k. Does the bracket being nonneg for all (i,j) reduce to a condition on (q, r) that tree-realizable pairs always satisfy?

### Q5: Test candidate invariants

Run the following checks for ALL tree-realizable (E_t, J_t) pairs at n ≤ 12:

1. **E_t ≽ J_t:** already known to hold.
2. **SCC:** (1+x)I_t ≽ E_t. Already known to hold.
3. **I_t ≽ (1+x)J_t:** Does I_t ratio-dominate (1+x)J_t? (From MEMORY: fails at ~12%)
4. **Leaf-augmentation at child:** I_t + x·E_t ≽ E_t. (Should hold by general conjecture.)
5. **g_k(child) = c_k(E_t) + d_k(I_t, E_t) ≥ 0?** (Should hold.)
6. **Stronger SCC:** Is Δ_k((1+x)I_t, E_t) ≥ Δ_k(x·I_t, E_t) for all k? (SCC margin vs contribution.)
7. **NEW: Is d_{k-1}(J_t, E_t) ≥ 0 for k ≤ mode?** (E_t ≽ J_t in prefix only vs everywhere.)
8. **NEW: Is the sequence (J_t,k / E_t,k) log-concave?** (Stronger than monotone.)

Report which invariants hold universally and which fail, with counterexamples.

## Deliverables

1. Full enumeration of tree-realizable (E_t, J_t) for subtrees of size 1-8.
2. Assessment of each candidate invariant (Q5 items 1-8).
3. Identification of the minimal additional constraint (beyond the 6 known) that makes the correction term nonneg.
4. If possible: a "tree-realizable cone lemma" statement that implies W.
