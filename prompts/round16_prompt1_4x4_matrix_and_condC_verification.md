# Round 16, Prompt 1: 4×4 N_c Matrix and Condition C Verification

**Target model:** Codex 5.3 (exhaustive computation)

---

## Context

We are proving E ≽ J (prefix ratio dominance) at support vertices of trees. The W-form induction processes non-leaf children one at a time. The key remaining gap is **chain STP2 multi-child closure**.

Two new frameworks emerged from Round 15 that need computational verification:

1. **4×4 update matrix N_c** with ladder-minor preservation
2. **Condition C**: a curvature-based bound on the negative term in the W-form identity

Both are promising but neither has been verified computationally beyond small examples.

## Known properties (background)

- At each step, the accumulated state is (E_acc, J_acc) with E_acc ≽ J_acc (prefix).
- Processing child c with triple (I_c, E_c, J_c):
  - E_new = E_acc · I_c
  - J_new = J_acc · E_c
- σ(F) = (F-1)/x for F(0)=1. σ(FG) = σ(F) + F·σ(G). σ(I) = σ(E) + J.
- M(F) = [F, 0; σ(F), 1] satisfies M(FG) = M(F)·M(G). Verified n≤12.
- Chain STP2: both STP2(E_t, J_t) and STP2(I_t, E_t) hold, 0 failures n≤22.
- STP2 diagonal form (Lemma A): under LC(f), STP2(f,g) ⟺ g(k+1)·f(k-1) ≤ g(k)·f(k) ∀k.

## Task 1: Construct and verify the 4×4 update matrix N_c

The augmented state vector tracks 4 polynomial quantities:
```
state = (σ(E_acc), E_acc, σ(J_acc), J_acc)
```

When processing child c with (I_c, E_c, J_c), the update rules are:
```
σ(E_new) = σ(E_acc · I_c) = σ(E_acc)·I_c + E_acc·σ(I_c)
                            = σ(E_acc)·I_c + E_acc·(σ(E_c) + J_c)
E_new = E_acc · I_c
σ(J_new) = σ(J_acc · E_c) = σ(J_acc)·E_c + J_acc·σ(E_c)
J_new = J_acc · E_c
```

Write this as a 4×4 matrix N_c (entries are polynomials in x) times the state vector:
```
[σ(E_new)]     [I_c    σ(E_c)+J_c   0       0      ] [σ(E_acc)]
[E_new   ]  =  [0      I_c          0       0      ] [E_acc   ]
[σ(J_new)]     [0      0            E_c     σ(E_c) ] [σ(J_acc)]
[J_new   ]     [0      0            0       E_c    ] [J_acc   ]
```

**Verify** this matrix equation holds for all tree-realizable triples (I_c, E_c, J_c) with n ≤ 12. For each tree T rooted at each vertex, compute the state after processing each child both directly (via the DP) and via repeated N_c multiplication. Confirm they agree.

## Task 2: Ladder-minor positivity — define and check

The STP2(E_acc, J_acc) diagonal condition is:
```
J_acc(k+1) · E_acc(k-1) ≤ J_acc(k) · E_acc(k)   ∀k
```

In terms of σ, this is:
```
σ(J_acc)(k) · E_acc(k-1) ≤ J_acc(k) · σ(E_acc)(k-1)   ∀k
```

which is equivalent to: for each k, the 2×2 determinant
```
det [σ(E_acc)(k-1)   σ(J_acc)(k) ]
    [E_acc(k-1)      J_acc(k)    ]  ≥ 0
```

This is the **ladder minor** of the state vector's coefficient matrix at index k.

Define:
```
Λ_k = σ(E_acc)(k-1) · J_acc(k) - E_acc(k-1) · σ(J_acc)(k)
    = E_acc(k) · J_acc(k) - E_acc(k-1) · J_acc(k+1)
```

**For each tree n ≤ 17 (exhaustive), at each support vertex, at each incremental step t:**
1. Compute the state vector (σ(E_acc), E_acc, σ(J_acc), J_acc)
2. Compute Λ_k for all k in the prefix range (k = 0, ..., mode-1)
3. Verify Λ_k ≥ 0 (this is just STP2 rephrased, so should always hold)
4. Compute Λ_k BEFORE and AFTER applying N_c at step t
5. Report: does Λ_k always increase or stay the same? Or can it decrease?
6. If it can decrease, what is the minimum ratio Λ_k(after) / Λ_k(before)?

The key question: **is the ladder-minor Λ_k monotonically nondecreasing through the N_c updates?** Or is there a more nuanced preservation mechanism?

## Task 3: Condition C verification

The curvature identity at each step t:
```
C_k · w_k = C_k · Karlin_k + C_{k+1} · d_{k-1}(B,C) + B_k · c_k(C)
```
where A = E_acc·g, B = E_acc·h, C = J_acc·g (g=E_t, h=J_t).

**Condition C** says: the negative term (Term 2) is controlled by the curvature bonus (Term 3):
```
-d_{k-1}(B,C) ≤ θ_k · (B_k / C_{k+1}) · c_k(C)
```
for some θ_k < 1. When this holds, w_k ≥ 0 follows.

**For each tree n ≤ 17 (exhaustive), at each support vertex, at each step t ≥ 2, at each k in the prefix:**
1. Compute Karlin_k, d_{k-1}(B,C), B_k, C_{k+1}, c_k(C)
2. When d_{k-1}(B,C) < 0, compute θ_k = -d_{k-1}(B,C) · C_{k+1} / (B_k · c_k(C))
3. Report: max θ_k across all trees/steps/k values
4. How does max θ_k vary by (n, step t, k)?
5. What tree/step/k achieves the maximum θ_k?
6. Does max θ_k approach 1 as n increases, or is it bounded away from 1?

**Critical:** if max θ_k < 1 with a finite gap, Condition C holds universally and the proof is done.

## Task 4: Condition C at star+star extremal

For the star+star family (root + ℓ leaves + star(a) + star(b)):
- E_1 = (1+x)^a, J_1 = [1], E_2 = (1+x)^b, J_2 = [1]
- At step 2: g = E_2 = (1+x)^b, h = J_2 = [1]
- B = E_acc · [1] = E_acc
- C = E_1 · E_2 = (1+x)^{a+b}

Compute θ_k for (a, b) ranging over 1 ≤ a ≤ b ≤ 100, ℓ = 1, at each k.

Report:
- max θ_k as a function of (a, b)
- Limiting behavior as a, b → ∞ (balanced a = b = n and unbalanced a = 1, b → ∞)
- Does θ_k → some limit < 1?
- The extremal (a, b, k) triple achieving the largest θ

## Task 5: N_c entries — positivity structure

For each tree-realizable triple (I_c, E_c, J_c) with n ≤ 15:
1. Construct N_c
2. Check: are all entries of N_c nonneg-coefficient polynomials?
3. Compute the Toeplitz matrix of each polynomial entry
4. For the 2×2 blocks of N_c, check: are the 2×2 Toeplitz blocks TN2?
5. Specifically: is the [σ(E_c)+J_c, I_c; 0, I_c] upper-left block's coefficient matrix TN2?

The goal: identify which positivity properties of N_c's entries can be proved from tree-realizability and might imply ladder-minor preservation.

## Deliverables

1. 4×4 N_c verification results (n ≤ 12)
2. Ladder-minor Λ_k preservation profile (n ≤ 17): always preserved, monotone, or more subtle?
3. Condition C θ_k profile by (n, step, k) — does θ stay bounded below 1?
4. Star+star θ_k for a, b ≤ 100 — limiting behavior
5. N_c positivity structure catalog
6. Assessment: which framework (ladder-minor vs Condition C) looks more tractable?
