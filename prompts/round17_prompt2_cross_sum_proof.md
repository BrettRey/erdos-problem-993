# Round 17, Prompt 2: Cross Sum X_k ≥ 0 — Algebraic Proof

**Target model:** GPT 5.2 Pro (theory-heavy)

---

## Context

We are proving the independence set sequence of every tree is unimodal. The proof has been reduced to a single lemma: the **cross sum** X_k in the Cauchy-Binet expansion of the ladder minor is nonneg.

### The CB expansion

At each step of the support-vertex DP, processing child c:
```
E_new = E_acc · I_c,   J_new = J_acc · E_c
```

The ladder minor Λ_k = E(k)J(k) - E(k-1)J(k+1) has the CB expansion:
```
Λ_k^{new} = Σ_{i,j} Δ_{i,j}(A,B) · P(k-i) · Q(k-j)
```
where:
- A = E_acc, B = J_acc (accumulated state, polynomials with nonneg coefficients)
- P = I_c, Q = E_c (child factor)
- Δ_{i,j}(A,B) = A(i)B(j) - A(i-1)B(j+1)

Decomposition: Λ_k^{new} = D_k + X_k where
- D_k = Σ_i Λ_i^{old} · P(k-i)Q(k-i) ≥ 0 (diagonal, from IH)
- X_k = Σ_{i≠j} Δ_{i,j}(A,B) · P(k-i)Q(k-j) (cross terms)

**Target**: prove X_k ≥ 0.

### Verified computationally
- X_k ≥ 0: verified for w_k in Round 10 (410K support vertices, n≤18, 0 failures)
- Ladder-minor monotonicity (implies D_k + X_k ≥ Λ_k^{old} ≥ 0): R16, 5.5M checks, 0 violations

### Available properties
- **STP2(A,B)**: Δ_{i,j}(A,B) ≥ 0 for j ≥ i. This makes all above-diagonal cross terms nonneg.
- **STP2(P,Q) = STP2(I_c, E_c)**: holds for all tree-realizable triples.
- **LC(A), LC(B)**: A and B are log-concave (products of LC factors).
- **LC(P), LC(Q)**: P and Q are log-concave.
- **J_c ≤ E_c** coefficientwise (proved).
- **I_c = E_c + x·J_c**: the key structural identity.

## Task 1: The I = E + xJ split

Write P(k-i) = Q(k-i) + R(k-i-1) where R = J_c. Then:
```
X_k = Σ_{i≠j} Δ_{i,j}(A,B) · [Q(k-i)Q(k-j) + R(k-i-1)Q(k-j)]
    = X_k^{pure} + X_k^{corr}
```

### Pure E part

```
X_k^{pure} = Σ_{i≠j} Δ_{i,j}(A,B) · Q(k-i)·Q(k-j)
```

The factor Q(k-i)·Q(k-j) is **symmetric in i↔j**. Grouping (i,j) with (j,i) for i < j:
```
X_k^{pure} = Σ_{i<j} [Δ_{i,j} + Δ_{j,i}] · Q(k-i)·Q(k-j)
```

where Q(k-i)Q(k-j) > 0 (when both in range). So X_k^{pure} ≥ 0 iff each symmetrized bracket is nonneg:
```
Δ_{i,j}(A,B) + Δ_{j,i}(A,B) = [A(i)B(j) + A(j)B(i)] - [A(i-1)B(j+1) + A(j-1)B(i+1)]
```

**Key question**: for A, B with nonneg coefficients and STP2(A,B), is this symmetrized bracket always ≥ 0?

Note: A(i)B(j) + A(j)B(i) is the permanent of the 2×2 matrix [[A(i), A(j)], [B(i), B(j)]]. The bracket says: the permanent at indices (i,j) ≥ permanent at indices (i-1, j+1). This is a "Schur-concavity" type condition on the permanent.

**Investigate**: does LC(A) + LC(B) + STP2(A,B) + A ≥ B coefficientwise imply this?

### J correction part

```
X_k^{corr} = Σ_{i≠j} Δ_{i,j}(A,B) · R(k-i-1) · Q(k-j)
```

This is NOT symmetric in i↔j. The factor R(k-i-1)·Q(k-j) depends on i,j differently.

For i < j: Δ_{i,j} ≥ 0, R(k-i-1) is at a higher index than R(k-j-1), and Q(k-j) is at a lower index than Q(k-i). Since R = J_c ≤ E_c = Q coefficientwise, R terms are "smaller."

**Investigate**: can X_k^{corr} be bounded in terms of X_k^{pure}? I.e., is |X_k^{corr}| ≤ X_k^{pure}?

## Task 2: The symmetrized bracket

Focus on proving:
```
F(i,j) := A(i)B(j) + A(j)B(i) - A(i-1)B(j+1) - A(j-1)B(i+1) ≥ 0   for i < j
```

where A, B are nonneg LC sequences with STP2(A,B) and A ≥ B coefficientwise.

### Approach 1: Direct manipulation

Write F(i,j) = [A(i)B(j) - A(i-1)B(j+1)] + [A(j)B(i) - A(j-1)B(i+1)]
             = Δ_{i,j}(A,B) + Δ_{j,i}(A,B)

The first term Δ_{i,j} ≥ 0 (STP2, since j > i). The second Δ_{j,i} can be negative.

Can we bound |Δ_{j,i}| ≤ Δ_{i,j}? STP2(A,B) gives Δ_{i,j} ≥ 0 but says nothing directly about Δ_{j,i}.

### Approach 2: LC of A and B

By LC(A): A(i)A(j) ≥ A(i-1)A(j+1) for i ≤ j (the "interlacing" property).
By LC(B): B(i)B(j) ≥ B(i-1)B(j+1) for i ≤ j.

Multiply: A(i)A(j)B(i)B(j) ≥ A(i-1)A(j+1)B(i-1)B(j+1).

This gives: [A(i)B(i)][A(j)B(j)] ≥ [A(i-1)B(i-1)][A(j+1)B(j+1)], i.e., AB is LC.

But F involves A(i)B(j) + A(j)B(i) (mixed products), not A(i)B(i).

### Approach 3: Schur-convexity / majorization

The function g(i,j) = A(i)B(j) + A(j)B(i) for fixed A,B is symmetric. It's a weighted permanent.

For LC nonneg sequences, g might be Schur-concave on the integer lattice: more "equal" (i,j) gives larger g.

The condition F(i,j) ≥ 0 says g(i,j) ≥ g(i-1,j+1) for i < j. Since (i-1,j+1) is more spread than (i,j) (same sum, larger gap), this IS Schur-concavity.

**Prove or disprove**: A(i)B(j) + A(j)B(i) is Schur-concave for LC nonneg A, B.

Note: for A = B (both LC), this is 2A(i)A(j) ≥ 2A(i-1)A(j+1) for i ≤ j, which is exactly LC! So for A = B, it holds.

For A ≠ B, we need the mixed version. This is related to the Hadamard product of LC sequences.

### Approach 4: STP2 as extra leverage

STP2(A,B) says B(j+1)/B(j) ≤ A(i)/A(i-1) for j > i-1, i.e., j ≥ i.

For the first bracket: A(i)B(j) - A(i-1)B(j+1) = A(i-1)B(j)[A(i)/A(i-1) - B(j+1)/B(j)] ≥ 0.

For the second bracket: A(j)B(i) - A(j-1)B(i+1) = A(j-1)B(i)[A(j)/A(j-1) - B(i+1)/B(i)].

The sign depends on A(j)/A(j-1) vs B(i+1)/B(i). Since j > i, and LC makes ratios decrease:
- A(j)/A(j-1) ≤ A(i)/A(i-1) (LC of A)
- B(i+1)/B(i) ≥ B(j+1)/B(j) (LC of B, i < j)

STP2 gives B(j+1)/B(j) ≤ A(i)/A(i-1). Combined: B(i+1)/B(i) ≥ B(j+1)/B(j) and A(j)/A(j-1) ≤ A(i)/A(i-1). So:

We need A(j)/A(j-1) ≥ B(i+1)/B(i), which is STP2 at (m=i, n=j-1) — but this requires m > n, i.e., i > j-1, i.e., i ≥ j. Since i < j, this is EXACTLY the wrong direction.

So STP2 does NOT directly control the second bracket. The below-diagonal term resists pointwise control.

### Approach 5: Aggregation

Maybe individual pairwise sums can be negative, but the full sum over (i,j) is nonneg.

Group the cross terms by "gap" g = j - i > 0:
```
X_k = Σ_{g=1}^{∞} Σ_i Δ_{i,i+g}(A,B) · P(k-i)Q(k-i-g) + Δ_{i+g,i}(A,B) · P(k-i-g)Q(k-i)
```

For each gap g, the inner sum involves A,B at pairs (i, i+g). Does the sum over i have a sign?

**Key insight**: summing Δ_{i,i+g}(A,B) over i with weights gives:
```
Σ_i Δ_{i,i+g}(A,B) · w_i = Σ_i [A(i)B(i+g) - A(i-1)B(i+g+1)] · w_i
```
which telescopes if w_i has special structure.

**Investigate** whether the factor weights P(k-i)Q(k-j) or Q(k-i)Q(k-j) enable a useful telescoping or Abel summation.

## Task 3: The "complementary STP2" argument

The CB expansion has TWO STP2 conditions:
1. STP2(A,B) controls the accumulator minors Δ_{i,j}
2. STP2(P,Q) = STP2(I_c, E_c) controls the factor contribution

For the full CB sum to be nonneg, we need the negative contributions (from the "wrong" quadrant of each STP2) to be dominated by the positive contributions (from the "right" quadrant).

The key insight: the two STP2 conditions cover COMPLEMENTARY regions.
- When i < j: Δ_{i,j}(A,B) ≥ 0 (STP2(A,B)), but factor minors uncertain
- When i > j: Δ_{i,j}(A,B) uncertain, but factor has P at lower index and Q at higher index

**Can you formalize this complementarity?** Is there a matrix identity or convolution argument that combines the two STP2 conditions to give nonnegativity of the full sum?

One approach: rewrite X_k as a quadratic form in the coefficient vector of some auxiliary sequence, and show the associated matrix is positive semidefinite using both STP2 conditions.

## Task 4: Connection to w_k

The ladder minor Λ_k = E(k)J(k) - E(k-1)J(k+1) and the P2 minor w_k = E(k)J(k-1) - E(k-1)J(k) are related but not identical (different index shift).

Round 16 P2 noted: Λ_k(E_new,J_new) = Λ_k(A*Q, B*Q) + Λ_{k-1}(A*R, B*Q).

The first term is the "Karlin" part (same denominator Q for both), and the second is the "correction" from R = J_c.

**Show**: the "Karlin" part Λ_k(A*Q, B*Q) has the same CB structure as X_k^{pure}, and the "correction" Λ_{k-1}(A*R, B*Q) relates to X_k^{corr}.

This would unify the two decompositions and potentially allow the Karlin ≥ 0 proof to be lifted to the full setting.

## Deliverables

1. Analysis of X_k^{pure} and X_k^{corr} — which is harder to control?
2. The symmetrized bracket F(i,j): proof or counterexample of Schur-concavity for LC sequences
3. If F(i,j) ≥ 0 fails: which pairwise terms can be negative? Can they be aggregated away?
4. The "complementary STP2" argument: formalization attempt
5. Aggregation approaches: Abel summation, telescoping, or quadratic form
6. Connection between X_k decomposition and Karlin/correction decomposition
7. Assessment: which approach is most promising? What's the minimal additional property needed?
