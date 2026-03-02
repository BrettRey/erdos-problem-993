# Round 14, Prompt 2: Cauchy-Binet Pairwise Domination

**Target model:** GPT 5.2 Pro (theory-heavy)

---

## Context

We are proving E ≽ J at support vertices of trees (Erdős Problem #993). The W-form at each incremental step gives:

```
w_k = Karlin(k) + correction(k) ≥ 0
```

where:
- Karlin(k) = d_k(E_acc·g, J_acc·g) = A_{k+1}C_k - A_kC_{k+1} ≥ 0 (provably, from E_acc ≽ J_acc + LC of g)
- correction(k) = B_kC_k - B_{k-1}C_{k+1} (can be negative)
- A = E_acc·g, B = E_acc·h, C = J_acc·g, where g = E_child, h = J_child

**Computationally verified:** 0 failures through n=24, 10M+ steps. The Karlin term ALWAYS dominates |correction| when correction < 0.

**The proof reduces to:** show Karlin(k) ≥ |correction(k)| for all prefix k, at all s ≥ 2 support vertices.

## Key Structural Properties

1. **E_acc ≽ J_acc** (inductive hypothesis): d_k(E_acc, J_acc) ≥ 0 for k in prefix
2. **h ≤ g coefficientwise** (proved: J_child ≤ E_child)
3. **g is LC** (E_child is LC at all rootings)
4. **h is LC** (J_child is a product of LC polynomials)
5. **STP2(g, h)**: h(m+1)·g(n) ≤ h(m)·g(n+1) for m > n (verified 0 failures, 921K+ instances)
6. **g(0) = h(0) = 1** (empty set is independent)

## The Cauchy-Binet Expansion

### Karlin term

By the CB formula for convolution LR minors:

```
d_k(E_acc·g, J_acc·g) = Σ_{i<j} W_{ij} · K_{ij}(k)
```

where:
- W_{ij} = E_acc(j)·J_acc(i) - E_acc(i)·J_acc(j) ≥ 0 (from E_acc ≽ J_acc, since i < j)
- K_{ij}(k) = g(k+1-i)·g(k-j) - g(k-i)·g(k+1-j) ≥ 0 (from LC of g, since k-j < k+1-i)

Each pair (i,j) contributes a non-negative amount W_{ij}·K_{ij} to the Karlin term.

### Correction term

The correction is B_kC_k - B_{k-1}C_{k+1} = d_k(x·B, C) where (x·B)(m) = B(m-1).

Since B = E_acc·h and C = J_acc·g:

```
correction = Σ_i Σ_j [E_acc(i)·h(k-i)] · [J_acc(j)·g(k-j)]
           - Σ_i Σ_j [E_acc(i)·h(k-1-i)] · [J_acc(j)·g(k+1-j)]
```

Rearranging with the same (i,j) pairing trick:

```
correction = Σ_{i,j} E_acc(i)·J_acc(j) · [h(k-i)·g(k-j) - h(k-1-i)·g(k+1-j)]
```

Split into symmetric/antisymmetric parts over i < j:

```
correction = Σ_{i<j} W_{ij} · L_{ij}(k) + Σ_{i<j} S_{ij} · M_{ij}(k) + Σ_i diag_i
```

where:
- W_{ij} = E_acc(j)J_acc(i) - E_acc(i)J_acc(j) (same as in Karlin: ≥ 0)
- S_{ij} = E_acc(j)J_acc(i) + E_acc(i)J_acc(j) (symmetric part, always ≥ 0)
- L_{ij}(k) = [h(k-j)g(k-i) - h(k-1-j)g(k+1-i)] - [h(k-i)g(k-j) - h(k-1-i)g(k+1-j)] (antisymmetric mixed bracket)
- M_{ij}(k) = [h(k-j)g(k-i) - h(k-1-j)g(k+1-i)] + [h(k-i)g(k-j) - h(k-1-i)g(k+1-j)] (symmetric mixed bracket)

Note: the exact index arithmetic needs care. Please derive the correct CB decomposition of the correction term and identify the bracket types.

## Your Tasks

### Part 1: Derive the CB expansion of the correction

Write the correction B_kC_k - B_{k-1}C_{k+1} as a sum over pairs (i,j) with i < j, identifying:
- Which terms share the SAME W_{ij} weights as the Karlin term
- Which terms have different (symmetric) weights
- The exact bracket expressions for each

### Part 2: STP2 sign control

Under STP2(g, h) (i.e., h(m+1)g(n) ≤ h(m)g(n+1) for m > n), determine the sign of the mixed brackets.

Specifically, for the bracket h(a)g(b) - h(a-1)g(b+1) with a > b:
- STP2 gives h(a)g(b) ≤ h(a-1)g(b+1) (setting m=a-1, n=b in STP2: h(a)g(b) ≤ h(a-1)g(b+1)). So this bracket is ≤ 0.

Work out the sign of each bracket type in the correction's CB expansion under STP2.

### Part 3: Pairwise domination attempt

For each pair (i,j) with i < j, we have:
- Karlin contribution: W_{ij} · K_{ij}(k)  (≥ 0, since both factors ≥ 0)
- Correction contribution: W_{ij} · L_{ij}(k) + ... (sign depends on brackets)

**Key question:** For each (i,j), is |L_{ij}(k)| ≤ K_{ij}(k)?

If yes, then pairwise: Karlin_pair ≥ |correction_pair|, and summing gives w_k ≥ 0.

This would mean: the pure (g,g) bracket at each pair dominates the mixed (h,g) bracket.

Since h ≤ g coefficientwise, each h factor is "smaller" than the corresponding g factor. Intuitively, replacing one g by h in a bracket should make it smaller. Can you formalize this?

Specifically:
```
K_{ij}(k) = g(k+1-i)g(k-j) - g(k-i)g(k+1-j)
L_{ij}(k) ≈ h(a)g(b) - h(c)g(d)  (for appropriate indices)
```

If h(a) ≤ g(a) and h(c) ≤ g(c), then:
```
|h(a)g(b) - h(c)g(d)| ≤ max(h(a)g(b), h(c)g(d)) ≤ max(g(a)g(b), g(c)g(d))
```

But K_{ij}(k) = g(a')g(b') - g(c')g(d') is a DIFFERENCE, not a max. So we need K_{ij} ≥ max(...), which is stronger than K_{ij} ≥ 0.

**Does this hold?** For the pure bracket, can we show g(a)g(b) - g(c)g(d) ≥ g(a)g(b) (i.e., g(c)g(d) ≤ 0)? No, since g has positive coefficients. So pairwise domination in this crude form fails.

**Better approach:** Use the specific index structure. The pure bracket K_{ij}(k) uses consecutive index shifts (k+1-i, k-j vs k-i, k+1-j), while the mixed bracket L_{ij}(k) uses shifted indices (k-i, k-j vs k-1-i, k+1-j). The LC structure of g might make K_{ij} inherently larger.

### Part 4: Alternative — aggregate bound via Schur complement

Instead of pairwise, try an aggregate bound. The W-form is:

```
w_k = Σ_{i<j} W_{ij} · [K_{ij}(k) + L_{ij}(k)] + (symmetric terms) + (diagonal terms)
```

The W_{ij} weights are the same. If we can show K_{ij} + L_{ij} ≥ 0 for each pair, we're done.

K_{ij} + L_{ij} ≥ 0 means: the combined (Karlin + correction) contribution of each (i,j) pair is non-negative. This is a LOCAL condition on 4 values of g and 2 values of h.

Check computationally (using tree data) whether K_{ij} + L_{ij} ≥ 0 for each pair. If this FAILS, the pairwise approach is dead and we need a global argument.

### Part 5: The star+star extremal insights

At the star+star extremal (g = (1+x)^a, h = [1]):
- K_{ij}(k) = binom(a,k+1-i)·binom(a,k-j) - binom(a,k-i)·binom(a,k+1-j) (Newton's inequality form)
- L_{ij}(k) = δ_{k-i,0}·binom(a,k-j) - δ_{k-1-i,0}·binom(a,k+1-j) (very sparse since h = [1])

Can you exploit this sparsity to prove K + L ≥ 0 for the extremal family, then argue that general trees are "further from the boundary"?

## Deliverables

1. Correct CB expansion of the correction term (exact index formulas)
2. Sign analysis of mixed brackets under STP2
3. Assessment of pairwise domination: does |L_{ij}| ≤ K_{ij} hold?
4. If not: assessment of K_{ij} + L_{ij} ≥ 0 (combined pair positivity)
5. Star+star case analysis (where h = [1] simplifies everything)
6. Recommendation: is the CB approach viable, or should we try something else?
