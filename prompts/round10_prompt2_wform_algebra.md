# Round 10, Prompt 2: Algebraic Structure of the W-Form for s >= 2

**Target model:** GPT 5.2 Pro (theory-heavy)

---

## Context

We are proving E >= J (ratio dominance) at every support vertex of every tree. This implies unimodality of the independent set sequence (Erdos Problem #993).

**The proof reduces to the W-form.** At each incremental step when processing non-leaf children of a support vertex:
```
W_k = J_k * Delta_k(A, C) + J_{k+1} * d_{k-1}(B, C) >= 0
```

where A = f*q, B = f*r, C = g*q, with:
- f = E^(t-1), g = J^(t-1) (accumulated from previous steps; f >= g by induction)
- q = E_t, r = J_t (DP pair of current child subtree)
- J = J^(t) = C = g*q

**Term 1** (Karlin): J_k * Delta_k(A,C) >= 0 always. PROVED via Karlin's theorem (f >= g implies f*q >= g*q in ratio dominance, since q is PF2).

**Term 2** (correction): J_{k+1} * d_{k-1}(B,C) can be negative. d_{k-1}(B,C) = B_k*C_{k-1} - B_{k-1}*C_k = (f*r)_k*(g*q)_{k-1} - (f*r)_{k-1}*(g*q)_k.

**The s=1 case is PROVED:** When there's only one non-leaf child, both terms are individually non-negative.

**The s >= 2 case is OPEN:** The Cauchy-Binet pointwise bracket approach FAILS on tree-realizable data (tested in Round 9). Generic product closure is FALSE.

## CRITICAL NEW FINDING: SCC is FALSE at n=28

The T_{3,4} broom (n=28) rooted at degree-2 support vertices gives:
- SCC ((1+x)I >= E): d_14 = -2298. FALSE.
- Leaf-augmentation (I+xE >= E): g_14 = -334. FALSE.
- E >= J: HOLDS (all minors >= 0)
- W-form: HOLDS (s=1, both terms individually non-negative)

This kills SCC as an intermediate target. The W-form bypasses SCC entirely.

## Computational Profile of the W-form (n <= 20)

- 0 failures across all trees, all support vertices, all incremental steps
- Min ratio (term1/|term2| when term2 < 0): **1.195** at n=20, s=1
- For s=1: ratio formula = ell*(ell+1)/(ell-1)^2 -> 1 as ell -> infinity
- For s >= 2: min ratio >= 1.664 (n <= 20)
- Absolute margin (term1 + term2) grows with n
- Curvature term (B_k * c_k(J)) is always non-negative but never needed

## Your Task: Find an Algebraic Proof of W >= 0 for s >= 2

### Approach A: Factor-Pair Structure

For each step t, the pair (q, r) = (E_t, J_t) comes from a child subtree. Known properties:
- q is PF2 (LC, nonneg)
- r <= q coefficientwise (J <= E is PROVED)
- q >= r in ratio dominance (E >= J at child, by inductive hypothesis at support vertices OR by recursive assumption)
- SCC at factor: (1+x)(q + xr) >= q (verified 0 fails n <= 22; may fail for large n — unknown)
- Leaf-aug at factor: (q + xr) + xq >= q (verified 0 fails n <= 12; may fail for large n — unknown)

The key question: WHICH properties of (q, r) are needed to make W_k >= 0?

### Approach B: The Correction as a Cross-Convolution Minor

Write:
```
d_{k-1}(f*r, g*q) = sum_{i,j} (f_i * r_{k-j} * g_j * q_{k-1-i} - f_i * r_{k-1-j} * g_j * q_{k-i})
```

This can be reorganized as:
```
d_{k-1}(B, C) = sum_{a < b} (f_a g_b - f_b g_a) * [r_{k-b} * q_{k-1-a} - r_{k-1-b} * q_{k-a}]
                + sum_{a < b} (f_b g_a - f_a g_b) * [... symmetric ...]
```

Wait — since f >= g (ratio dominance), the weights (f_a g_b - f_b g_a) have a definite sign pattern (>= 0 for a > b? or a < b? — depends on convention). The bracket involves a 2x2 minor of the matrix formed by shifting r and q.

Can you:
1. Write out the Cauchy-Binet expansion of d_{k-1}(B,C) explicitly
2. Combine with the Cauchy-Binet expansion of Delta_k(A,C) (which is the Karlin term)
3. Show that the combined sum J_k * [Karlin CB] + J_{k+1} * [correction CB] has all non-negative summands, OR
4. Show a pairing/coupling of terms that makes the sum non-negative

### Approach C: Recursion on s

For s = 2: the accumulated pair after step 1 is (E^(1), J^(1)) where E^(1) = (1+x)^ell * I_1, J^(1) = E_1. At step 2, we process (E_2, J_2).

The W-form at step 2 is:
```
W_k = J^(2)_k * Delta_k(E^(1)*E_2, J^(2)) + J^(2)_{k+1} * d_{k-1}(E^(1)*J_2, J^(2))
```

Can you factor this using the known structure of E^(1) = (1+x)^ell * I_1 and J^(1) = E_1?

Is there a "two-child identity" that reduces the s=2 W-form to quantities involving only (E_1, J_1), (E_2, J_2), and ell?

### Approach D: Monotone Transport

The W-form can be written as:
```
W_k = sum_j J_j * [Delta if j=k, d if j=k+1, 0 otherwise] applied to (A, B, C)
```

This is a weighted average of LR minors at consecutive positions. Can you interpret this as a transport problem where:
- "Supply" comes from the Karlin term (non-negative)
- "Demand" comes from the correction term (possibly negative)
- The J_k / J_{k+1} ratio controls the coupling

If J is LC (which is PROVED), the ratios J_k/J_{k+1} are non-increasing. Does this monotonicity help?

## Deliverables

1. Explicit Cauchy-Binet expansion of W_k as a double sum over (i,j)
2. Assessment: which approach (A-D) is most promising for s >= 2?
3. If possible: a proof of W >= 0 for s = 2 (simplest open case)
4. If not: identification of the precise obstruction and what additional constraint on (q,r) would resolve it
