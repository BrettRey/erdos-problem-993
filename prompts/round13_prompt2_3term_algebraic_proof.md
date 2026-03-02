# Round 13, Prompt 2: Three-Term Identity Algebraic Proof

**Target model:** GPT 5.2 Pro (theory-heavy)

---

## Context

We are proving prefix E ≽ J at support vertices of trees (Erdős Problem #993, unimodality of IS sequences). The W-form induction processes non-leaf children one at a time.

### The curvature-augmented identity (VERIFIED, 0 failures n≤20)

At each incremental step, with E_acc ≽ J_acc inductively:

```
C_k · w_k = Term1 + Term2 + Term3
```

where A = E_acc·g, B = E_acc·h, C = J_acc·g (g=E_t, h=J_t from child subtree t):

- **Term1** = C_k · d_k(A, C) = C_k · (A_{k+1}·C_k - A_k·C_{k+1}) ≥ 0
  - Proof: E_acc ≽ J_acc, g is LC → Karlin (TP2 convolution)

- **Term2** = C_{k+1} · d_{k-1}(B, C) = C_{k+1} · (B_k·C_{k-1} - B_{k-1}·C_k)
  - Can be negative: B = E_acc·h mixes E_acc with h=J_t

- **Term3** = B_k · c_k(C) = B_k · (C_k² - C_{k-1}·C_{k+1}) ≥ 0
  - Proof: C = J_acc·g, both LC → C is LC → c_k(C) ≥ 0, B_k ≥ 0

### Computational findings

- **Term1 + Term2 + Term3 ≥ 0**: 0 failures, 2.44M steps, n≤18
- **Term2 + Term3 < 0** (⋆ fails): 3.4% of steps. Karlin rescues all.
- **Min Karlin rescue ratio** (Term1/|Term2+Term3| when Term2+Term3 < 0):
  - Decreasing with n: 10 → 5.6 → 4.7 → ... → 3.33 → 3.29 → 3.20
  - Always at s=2, ell=1, k=2
  - Appears to converge to ~3

### What's proved

- **s=1 (one non-leaf child)**: PROVED. Handles 62.8% of support vertices.
- **P3 (tail domination)**: PROVED via leaf-swap injection.
- **Term1 ≥ 0**: PROVED (Karlin/TP2).
- **Term3 ≥ 0**: PROVED (LC of C).
- **J ≤ E coefficientwise**: PROVED.
- **E is LC at support vertices**: PROVED (product of LC polynomials).

### What's open

- **Term1 + Term2 + Term3 ≥ 0 for s≥2**: the identity holds computationally, but we need an algebraic proof.

## Your Tasks

### Part 1: Structural Analysis of Term2

Term2 = C_{k+1} · (B_k · C_{k-1} - B_{k-1} · C_k)

where B = E_acc · h and C = J_acc · g.

Since E_acc ≽ J_acc and h ≤ g coefficientwise (J_t ≤ E_t), we have B ≤ A coefficientwise but B and C are not directly comparable by ≽.

The sign of d_{k-1}(B, C) = B_k·C_{k-1} - B_{k-1}·C_k depends on whether B ≽ C. Since B = E_acc · h and C = J_acc · g, this is related to whether (E_acc/J_acc) · (h/g) is ratio-decreasing.

E_acc/J_acc is ratio-decreasing (by E_acc ≽ J_acc). But h/g = J_t/E_t is also ratio-decreasing (by E_t ≽ J_t inductively). The product of two ratio-decreasing sequences is NOT necessarily ratio-decreasing.

**Question**: Under what conditions does d_{k-1}(B, C) ≥ 0? This would make all three terms non-negative.

### Part 2: Why the Sum Term1 + Term2 + Term3 Works

Even when Term2 < 0 (which drags Term2+Term3 < 0 sometimes), Term1 always rescues. The min ratio is ~3 at large n.

Can you show Term1 ≥ |Term2| (ignoring Term3)? This would be:

```
C_k · d_k(A, C) ≥ |C_{k+1} · d_{k-1}(B, C)|
```

Expanding: A = E_acc·g, B = E_acc·h, C = J_acc·g.

Use the Cauchy-Binet expansion for d_k(A, C):

```
d_k(A, C) = Σ_{i<j} (E_acc_i · J_acc_j - E_acc_j · J_acc_i) · (g_{k+1-i}·g_{k-j} - g_{k-i}·g_{k+1-j})
```

And d_{k-1}(B, C):

```
d_{k-1}(B, C) = Σ_{i<j} (E_acc_i · J_acc_j - E_acc_j · J_acc_i) · (h_{k-i}·g_{k-1-j} - h_{k-1-i}·g_{k-j})
```

The (E_acc, J_acc) brackets are the SAME in both. The difference is in the (g,g) vs (h,g) brackets. Since h ≤ g, the (h,g) brackets should be "smaller" than (g,g) brackets.

Can you make this precise? Each pair (i,j) contributes:
- To Term1: C_k · (E_acc bracket) · (g,g bracket at k)
- To Term2: C_{k+1} · (E_acc bracket) · (h,g bracket at k-1)

The E_acc brackets have definite sign (from E_acc ≽ J_acc). So the question reduces to comparing:
```
C_k · (g,g bracket at k) vs C_{k+1} · (h,g bracket at k-1)
```
for each pair (i,j).

### Part 3: The Role of STP2

The STP2 condition (J_t(m+1)·E_t(n) ≤ J_t(m)·E_t(n+1) for m > n) controls the sign of the mixed brackets. Under STP2, the (h,g) brackets have controlled sign.

Does STP2 + LC of g + h ≤ g imply the pointwise bound from Part 2? Work out the algebra for a specific pair (i,j) and see if STP2 gives the needed inequality.

### Part 4: The Ratio ≥ 3 Phenomenon

The rescue ratio converges to ~3. Why 3?

In the curvature-augmented identity, the three terms have specific polynomial degrees and coefficient structures. The ratio 3 might come from:
- The number of terms in the expansion (three terms, one dominant)
- A specific Cauchy-Binet cancellation pattern
- An extremal tree family where the contributions balance in a 3:1 ratio

Can you identify the algebraic source of "3"?

### Part 5: Alternative Identity

Instead of the curvature-augmented identity, try the ORIGINAL w_k directly:

```
w_k = d_k(E_acc · f, J_acc · g) where f = g + x·h
```

Use the linearity of d_k in the first argument:
```
w_k = d_k(E_acc · g, J_acc · g) + d_k(x · E_acc · h, J_acc · g)
     = Karlin + correction
```

The Karlin term ≥ 0. The correction:
```
d_k(x·E_acc·h, J_acc·g) = (x·E_acc·h)_{k+1} · (J_acc·g)_k - (x·E_acc·h)_k · (J_acc·g)_{k+1}
                         = (E_acc·h)_k · C_k - (E_acc·h)_{k-1} · C_{k+1}
                         = B_k · C_k - B_{k-1} · C_{k+1}
```

So the correction is exactly the ⋆ quantity! And:

```
w_k = Karlin + (B_k · C_k - B_{k-1} · C_{k+1})
```

But this is just 2 terms (Karlin + ⋆), without the curvature. The curvature identity splits Term1 differently... or does it?

Wait, there's an inconsistency. Let me reconcile:
- 2-term: w_k = d_k(E_acc·g, J_acc·g) + (B_k·C_k - B_{k-1}·C_{k+1})
- 3-term: C_k·w_k = C_k·d_k(A,C) + C_{k+1}·d_{k-1}(B,C) + B_k·c_k(C)

Multiply the 2-term by C_k:
```
C_k·w_k = C_k·Karlin + C_k·(B_k·C_k - B_{k-1}·C_{k+1})
         = C_k·Karlin + B_k·C_k² - B_{k-1}·C_k·C_{k+1}
```

The 3-term version:
```
C_k·w_k = C_k·Karlin + C_{k+1}·(B_k·C_{k-1} - B_{k-1}·C_k) + B_k·(C_k² - C_{k-1}·C_{k+1})
         = C_k·Karlin + B_k·C_{k-1}·C_{k+1} - B_{k-1}·C_k·C_{k+1} + B_k·C_k² - B_k·C_{k-1}·C_{k+1}
         = C_k·Karlin + B_k·C_k² - B_{k-1}·C_k·C_{k+1}
```

They're the same! So the 3-term identity is just the 2-term identity multiplied by C_k, with the correction split into two pieces (Term2 and Term3).

This means the actual content is just:
```
w_k = Karlin + (B_k·C_k - B_{k-1}·C_{k+1})
```

The "curvature rescue" in the 3-term form is an artifact of the split — it's the SAME information rearranged.

**So the real proof target is**: Karlin(k) + B_k·C_k - B_{k-1}·C_{k+1} ≥ 0, where Karlin = d_k(E_acc·g, J_acc·g) ≥ 0 and the ⋆ term = B_k·C_k - B_{k-1}·C_{k+1} can be negative.

Equivalently: d_k(E_acc·g, J_acc·g) ≥ B_{k-1}·C_{k+1} - B_k·C_k.

Can you prove this 2-term version directly?

## Deliverables

1. Analysis of when d_{k-1}(B, C) ≥ 0 (making all 3 terms non-negative)
2. Attempt to show Term1 ≥ |Term2| pair-by-pair in the CB expansion
3. Role of STP2 in controlling mixed brackets
4. Algebraic source of the ratio ≥ 3 phenomenon
5. The 2-term equivalence (confirm or correct my reconciliation)
6. Direct proof attempt for: d_k(E_acc·g, J_acc·g) ≥ |B_k·C_k - B_{k-1}·C_{k+1}|
