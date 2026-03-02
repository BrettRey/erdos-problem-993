# Round 13, Prompt 3: STP2 Algebraic Proof

**Target model:** GPT 5.2 Pro (theory-heavy)

---

## Context

STP2 is a condition on tree-realizable factor pairs (E_t, J_t) that has 0 failures across 921K+ factor instances (exhaustive through n=18). It's the most promising algebraic condition for controlling the mixed brackets in the Cauchy-Binet expansion.

### The condition

**STP2**: J_t(m+1) · E_t(n) ≤ J_t(m) · E_t(n+1) for all m > n ≥ 0.

Equivalently: the ratio J_t(k)/E_t(k) is nonincreasing in k (when both positive).

### What we know about (E_t, J_t)

At a vertex t with children c_1, ..., c_d:
```
E_t = ∏ I_{c_i} = ∏ (E_{c_i} + x·J_{c_i})
J_t = ∏ E_{c_i}
```

Properties of tree-realizable (E_t, J_t) pairs:
1. **E_t ≽ J_t** at support vertices (prefix ratio dominance, 0 failures n≤22)
2. **J_t ≤ E_t** coefficientwise (PROVED: J^{(t)} ≤ E^{(t)} by induction)
3. **E_t is LC** (product of LC polynomials)
4. **J_t is LC** (product of E_{c_i}, each LC)
5. **I_t = E_t + x·J_t** is the IS polynomial of the subtree
6. **E_t(0) = J_t(0) = 1** always (empty set is independent, root excluded)
7. **deg(J_t) = deg(E_t) - 1** (J_t has one less degree since root is included)
8. **E_t - J_t has nonneg coefficients** (from property 2)

### STP2 in terms of E and J

STP2 says: for k < l, J_t(l)/E_t(l) ≤ J_t(k)/E_t(k).

Since both J/E ratios are between 0 and 1 (property 2), STP2 says the "inclusion fraction" J/E decreases with k. At small k, J ≈ E (both close to 1 at k=0). At large k, J << E (J has shorter degree).

**Stronger than E ≽ J**: E ≽ J says J(k+1)/E(k+1) ≤ J(k)/E(k) (consecutive indices). STP2 says this for ALL pairs m > n, not just consecutive.

### What STP2 controls

In the Cauchy-Binet expansion of w_k, the mixed brackets are:
```
L_{i,j} = h_{k-i} · g_{k-1-j} - h_{k-1-i} · g_{k-j}
```
where g = E_t, h = J_t.

For i < j (so k-i > k-j), setting a = k-i, b = k-j:
```
L_{i,j} = h_a · g_{b-1} - h_{a-1} · g_b
```

STP2 with m = a-1 > n = b-1 gives h_a · g_{b-1} ≤ h_{a-1} · g_b, so **L_{i,j} ≤ 0**.

This means: under STP2, all mixed brackets have definite sign. The correction term is systematically negative, bounded, and controlled.

## Your Tasks

### Part 1: Prove STP2 for Single-Child Vertices

When t has one child c (so t is degree 2 in the rooted tree):
```
E_t = I_c = E_c + x·J_c
J_t = E_c
```

STP2 becomes: E_c(m+1) · (E_c + x·J_c)(n) ≤ E_c(m) · (E_c + x·J_c)(n+1) for m > n.

I.e.: E_c(m+1) · (E_c(n) + J_c(n-1)) ≤ E_c(m) · (E_c(n+1) + J_c(n))

This is: E_c(m+1)·E_c(n) - E_c(m)·E_c(n+1) ≤ E_c(m)·J_c(n) - E_c(m+1)·J_c(n-1)

LHS ≤ 0 (since E_c is LC, with m > n this is a Toeplitz minor ≤ 0).

RHS: E_c(m)·J_c(n) - E_c(m+1)·J_c(n-1). For this to be ≥ 0, need STP2 of (E_c, J_c) at indices (m, n-1).

So STP2 for (I_c, E_c) follows from LC of E_c + STP2 of (E_c, J_c). This is an inductive reduction.

Can you complete this induction? The base case: at a leaf, E = [1], J = [1], STP2 is vacuous.

### Part 2: Prove STP2 for Multi-Child Vertices

When t has children c_1, ..., c_d:
```
E_t = ∏ I_{c_i},  J_t = ∏ E_{c_i}
```

If STP2 holds for each (I_{c_i}, E_{c_i}), does it hold for (∏I_{c_i}, ∏E_{c_i})?

The product structure: if two pairs (f_1, g_1) and (f_2, g_2) satisfy STP2, does (f_1·f_2, g_1·g_2)?

This is a product closure question. For individual factors, STP2 says the Toeplitz-like matrix T(f, g) with T_{m,n} = f(m)·g(n) - f(n)·g(m) is non-negative for m > n.

The convolution product: (f_1·f_2)(k) = Σ_i f_1(i)·f_2(k-i). The STP2 condition for the convolution product involves sums of products of the original T matrices.

This is related to the **Schur product theorem** (Hadamard product of positive-semidefinite matrices is PSD) or **Karlin's total positivity** results (products of TP matrices).

**Key question**: Is the STP2 matrix positive-semidefinite? If so, the Schur product theorem gives product closure.

Actually, STP2 says the matrix M with M_{m,n} = J(m)/E(m) (for n < m, with some modification) is "monotone." This might be related to totally positive matrices.

### Part 3: Does STP2 Follow from Known Conditions?

Consider whether STP2 follows from:
(a) E ≽ J (ratio dominance) + LC of both E and J
(b) J ≤ E coefficientwise + LC of both
(c) I = E + xJ has nonneg coefficients + LC of I

Construct test: find a pair (E, J) satisfying (a)-(c) but violating STP2. If such pairs exist, STP2 is a genuinely new condition requiring tree-specific structure.

If no such pairs exist (at least for "reasonable" coefficient sizes), STP2 may follow from the standard conditions.

### Part 4: Alternative — Direct Proof Without STP2

If STP2 is hard to prove, consider proving the weaker statement directly:

```
d_k(E_acc·g, J_acc·g) ≥ |(E_acc·h)_k · (J_acc·g)_k - (E_acc·h)_{k-1} · (J_acc·g)_{k+1}|
```

for k < mode, at support vertices.

This is: Karlin term ≥ |correction|. The correction involves h = J_t, and h ≤ g = E_t. The Karlin term involves g alone.

Since h ≤ g coefficientwise, (E_acc·h)_k ≤ (E_acc·g)_k = A_k. So the correction is bounded:

```
|B_k·C_k - B_{k-1}·C_{k+1}| ≤ max(B_k·C_k, B_{k-1}·C_{k+1}) ≤ A_k·C_k + A_{k-1}·C_{k+1}
```

And the Karlin term:
```
d_k(A, C) = A_{k+1}·C_k - A_k·C_{k+1}
```

So we need: A_{k+1}·C_k - A_k·C_{k+1} ≥ something involving A_k·C_k and A_{k-1}·C_{k+1}.

For k < mode in the growth phase, A_{k+1}/A_k is large (coefficients are growing). This makes the Karlin term large relative to the correction. Can you formalize this growth-phase argument?

## Deliverables

1. STP2 proof for single-child vertices (inductive step)
2. STP2 product closure analysis
3. Whether STP2 follows from standard conditions (E≽J + LC + coeff-wise)
4. If STP2 is blocked: direct proof of Karlin ≥ |correction| in the growth phase
5. Assessment: what's the most promising path to close the proof?
