# Round 17, Prompt 1: CB Pairwise Decomposition Scan

**Target model:** Codex 5.3 (exhaustive computation)

---

## Context

We are proving the independence set sequence of every tree is unimodal. The proof reduces to showing the **ladder minor** Λ_k = E(k)·J(k) - E(k-1)·J(k+1) stays nonneg through the support-vertex DP.

Round 16 established the **Cauchy-Binet expansion**:
```
Λ_k^{new} = Σ_{i,j} Δ_{i,j}(A,B) · P(k-i) · Q(k-j)
```
where:
- A = E_acc, B = J_acc (accumulated state)
- P = I_c = E_c + x·J_c, Q = E_c (child factor)
- Δ_{i,j}(A,B) = A(i)·B(j) - A(i-1)·B(j+1)

The **diagonal** (i=j) gives Δ_{i,i} = Λ_i^{old} (old ladder minor, ≥ 0 by IH), with weight P(k-i)·Q(k-i) ≥ 0.

The **cross terms** (i≠j) are the unknown. We need to verify that the cross sum X_k ≥ 0.

## Known properties

- STP2(A,B): Δ_{i,j}(A,B) ≥ 0 for j ≥ i (above-diagonal nonneg)
- STP2(P,Q) = STP2(I_c, E_c): holds for all tree-realizable triples
- LC(A), LC(B), LC(P), LC(Q): all hold
- J_c ≤ E_c coefficientwise (proved)
- I_c = E_c + x·J_c, so P(k) = Q(k) + R(k-1) where R = J_c

## Task 1: D_k / X_k decomposition for ladder minors

For each tree n ≤ 18 (exhaustive), at each support vertex, at each step t ≥ 2, at each k in the prefix:

1. Compute Λ_k^{new} via the DP (ground truth)
2. Compute D_k = Σ_i Δ_{i,i}(A,B) · P(k-i) · Q(k-i) = Σ_i Λ_i^{old} · P(k-i) · Q(k-i)
3. Compute X_k = Λ_k^{new} - D_k
4. Verify D_k ≥ 0 (should always hold since Λ_i^{old} ≥ 0 and P,Q nonneg)
5. **Report: is X_k ≥ 0 always?** Count failures if any.
6. Compute ratio D_k / X_k when both are positive. Profile min X_k / D_k.

## Task 2: Pairwise symmetric sums

For each (i,j) pair with i < j, define the **symmetric sum**:
```
S_{i,j}(k) = Δ_{i,j}(A,B) · P(k-i) · Q(k-j) + Δ_{j,i}(A,B) · P(k-j) · Q(k-i)
```

This pairs the above-diagonal term (where Δ_{i,j} ≥ 0 by STP2(A,B)) with the below-diagonal term (where Δ_{j,i} can be negative).

For each tree n ≤ 15 (exhaustive), at each support vertex, at each step t ≥ 2, at each k:

1. Compute S_{i,j}(k) for all i < j with nonzero terms
2. **Report: is S_{i,j}(k) ≥ 0 always?** Count failures.
3. If failures exist: which (i, j, k, tree, step) achieves the minimum S?
4. What is the gap distance j-i for the worst violations?
5. Profile: does S_{i,j}(k) ≥ 0 hold more often for smaller j-i?

## Task 3: Pure E / J correction split

Using P(k-i) = Q(k-i) + R(k-i-1) where R = J_c (with R(-1) = 0):

Split each term:
```
Δ_{i,j}(A,B) · P(k-i) · Q(k-j) = Δ_{i,j}(A,B) · Q(k-i) · Q(k-j)    [Pure E part]
                                  + Δ_{i,j}(A,B) · R(k-i-1) · Q(k-j)  [J correction]
```

The "pure E" part has factor Q(k-i)·Q(k-j) which is **symmetric in i↔j**.

Define:
```
X_k^{pure} = Σ_{i≠j} Δ_{i,j}(A,B) · Q(k-i) · Q(k-j)
X_k^{corr} = Σ_{i≠j} Δ_{i,j}(A,B) · R(k-i-1) · Q(k-j)
```
so X_k = X_k^{pure} + X_k^{corr}.

For each tree n ≤ 15 (exhaustive), at each support vertex, step ≥ 2, each k:

1. Compute X_k^{pure} and X_k^{corr} separately
2. **Report**: is X_k^{pure} ≥ 0? Is X_k^{corr} ≥ 0? Or only their sum?
3. For the pure E part, the symmetric pairwise sum is:
   ```
   S^{pure}_{i,j}(k) = [Δ_{i,j}(A,B) + Δ_{j,i}(A,B)] · Q(k-i) · Q(k-j)
   ```
   Is this always ≥ 0? (The Q·Q factor is nonneg, so this reduces to: is Δ_{i,j} + Δ_{j,i} ≥ 0 for i < j?)
4. Compute Δ_{i,j}(A,B) + Δ_{j,i}(A,B) = A(i)B(j) + A(j)B(i) - A(i-1)B(j+1) - A(j-1)B(i+1) for all i < j. Is this always ≥ 0?

## Task 4: STP2(I,E) derivation check

STP2(I,E) diagonal at index k says: E(k)·I(k) ≥ E(k+1)·I(k-1).

Using I(k) = E(k) + J(k-1), this becomes:
```
c_k(E) + E(k)·J(k-1) - E(k+1)·J(k-2) ≥ 0
```
where c_k(E) = E(k)² - E(k-1)·E(k+1) is the LC gap of E.

The term E(k)J(k-1) - E(k+1)J(k-2) can be negative (star with J=[1], k=2: gives -E(3)).

For each tree n ≤ 18 (exhaustive), at EVERY rooting (not just support vertices), at each k:
1. Compute c_k(E), E(k)J(k-1) - E(k+1)J(k-2), and their sum
2. Verify the sum ≥ 0 (confirming STP2(I,E))
3. When E(k)J(k-1) - E(k+1)J(k-2) < 0, compute the ratio c_k(E) / |E(k)J(k-1) - E(k+1)J(k-2)|
4. Report the minimum ratio across all trees/rootings/k
5. Which tree/rooting/k achieves the minimum ratio?

## Task 5: Monotonicity decomposition

Since Λ_k^{new} ≥ Λ_k^{old} (confirmed R16, 5.5M checks), and the diagonal term at i=j=k gives exactly Λ_k^{old} · P(0)Q(0) = Λ_k^{old}, the remaining terms sum to a nonneg quantity:
```
Λ_k^{new} - Λ_k^{old} = Σ_{i≠k or j≠k} [term_{i,j}]
```

Decompose this "increment" into:
```
Δ_k^{inc} = [D_k - Λ_k^{old}] + X_k
           = Σ_{i≠k} Λ_i^{old} · P(k-i)Q(k-i) + X_k
```

The first sum is ≥ 0 (nonneg terms). So if X_k ≥ 0, monotonicity is explained.

If X_k is sometimes negative: how much of the increment comes from D_k vs X_k?

For n ≤ 15: decompose the increment and report the X_k contribution as a fraction of total.

## Deliverables

1. D_k ≥ 0 confirmation and X_k sign profile (n ≤ 18)
2. Pairwise S_{i,j}(k) sign profile — are all symmetric pairs nonneg? (n ≤ 15)
3. Pure E / J correction split — which part carries the positivity? (n ≤ 15)
4. Symmetrized Δ_{i,j} + Δ_{j,i} sign profile (n ≤ 15)
5. STP2(I,E) derivation check at ALL rootings (n ≤ 18)
6. Monotonicity decomposition (n ≤ 15)
7. Extremal tree/step/k for each metric
