# Round 12, Prompt 3: Completing the s ≥ 2 Proof

**Target model:** GPT 5.2 Pro (theory-heavy)

---

## Context

We are proving unimodality of the IS sequence for all trees (Erdős Problem #993). The proof reduces to showing prefix E ≽ J at support vertices, via a W-form induction processing non-leaf children one at a time.

**What's proved:**
- s = 1 (one non-leaf child): PROVED via Karlin + transitivity. Handles 62.8% of support vertices.
- P3 (tail domination): PROVED via leaf-swap injection.
- Karlin term ≥ 0: PROVED (TP2 convolution theorem).

**What's open:**
- s ≥ 2: the correction term d_k(x·E_acc·h, J_acc·g) needs bounding by the Karlin term.

**Extremal data (s ≥ 2):**

| s | min ratio (Karlin/|corr|) | absolute margin | extremal family |
|---|--------------------------|----------------|-----------------|
| 1 | 1.195 (n=20) | 2m → ∞ | pendant-star |
| 2 | 1.264 (n=18) | 1 | various |
| 3 | 1.320 (n=18) | 1 | various |
| 4 | 1.406 (n=18) | 1 | various |
| ≥5 | > 1.5 | 1 | rare |

Key: the ratio INCREASES with s (s=1 is hardest, PROVED). Absolute margin = 1 for all s.

**s = 2 closed forms (balanced double-star):**
- Ratio → 3/2 as arm length L → ∞
- Formula: ((2L-1)(3L²+5L+4)) / (4(L-1)(L²+2L-4))

**n = 32 failure tree:**
- Support vertex with ℓ=1, s=2 (P_2 child + T_{3,4} broom child)
- Full E≽J fails at k=15, mode=10, d_15 = -3498
- Prefix E≽J (k < 10) still holds

## Your Tasks

### Part 1: The s=1 Proof — What Makes It Work?

The s=1 proof: when there's exactly one non-leaf child t, E_acc = (1+x)^ℓ and J_acc = [1].
Then w_k = d_k((1+x)^ℓ · I_t, E_t).

By Karlin: if I_t ≽ E_t (SCC), then (1+x)^ℓ · I_t ≽ (1+x)^ℓ · E_t ≽ E_t.

But SCC fails at n=28! The s=1 proof actually uses a different argument: at the single-child case, E_acc = (1+x)^ℓ and J_acc = [1], so E_acc ≽ J_acc is trivial. Then the Karlin term d_k(g, g) = 0 (degenerate!) and the full w_k reduces to d_k((1+x)^ℓ · f, g) which is shown non-negative by the binomial smoothing lemma.

Wait — that doesn't sound right either. Let me state what's actually proved:

For s=1, E_new = (1+x)^ℓ · I_t and J_new = E_t. Then:
d_k(E_new, J_new) = d_k((1+x)^ℓ · I_t, E_t)

The proof uses: (1+x)^ℓ · I_t = (1+x)^ℓ · E_t + x(1+x)^ℓ · J_t.
Binomial smoothing: (1+x)^ℓ · E_t ≽ E_t (since E_t is LC).
The xJ_t part contributes a non-negative correction: x(1+x)^ℓ · J_t has non-negative coefficients, and adding it to (1+x)^ℓ · E_t only increases the ratios.

**Clarify**: write out the exact s=1 proof so we can see which step fails for s ≥ 2.

### Part 2: Why s = 2 is Harder

For s = 2 with children t_1, t_2:

After processing child 1:
```
E_acc = (1+x)^ℓ · I_1
J_acc = E_1
```

Processing child 2: f = I_2, g = E_2, h = J_2.
```
E_new = E_acc · (E_2 + xJ_2) = (1+x)^ℓ · I_1 · I_2
J_new = J_acc · E_2 = E_1 · E_2
```

So w_k = d_k((1+x)^ℓ · I_1 · I_2, E_1 · E_2).

The Karlin term: d_k(E_acc · g, J_acc · g) = d_k((1+x)^ℓ · I_1 · E_2, E_1 · E_2).
Since E_acc = (1+x)^ℓ · I_1 and E_acc ≽ J_acc = E_1 (this is exactly the s=1 result!), and g = E_2 is LC, Karlin gives Karlin term ≥ 0.

The correction: d_k(x · E_acc · h, J_acc · g) = d_k(x · (1+x)^ℓ · I_1 · J_2, E_1 · E_2).

This is a cross-term mixing I_1 and J_2 with E_1 and E_2. The difficulty: the x-shift and the mixing of child-1 and child-2 polynomials.

**Task**: Bound this correction. Use:
- E_1 ≽ J_1 and E_2 ≽ J_2 (by inductive hypothesis at support vertices, or just E ≽ J from s=1 proof)
- J_i ≤ E_i coefficientwise
- All polynomials LC
- (1+x)^ℓ provides binomial smoothing

### Part 3: Child Ordering

The incremental proof processes children in some order. For s = 2, the two orderings (t_1 first vs t_2 first) give the same w_k (since E_new and J_new don't depend on ordering), but the Karlin/correction split DOES depend on ordering.

**Question**: Is there always an ordering where the correction is "small enough"? The correction at step 2 depends on E_acc from step 1:
- If the "easier" child (smaller subtree? more leaves?) is processed second, E_acc from step 1 is larger/smoother, making the Karlin term dominate more.
- Conversely, processing the "harder" child second might give a worse correction.

Can we exploit ordering freedom? If one ordering always makes the correction bounded, that suffices.

### Part 4: Product Closure of the Ratio Bound

For s ≥ 3, we process children iteratively. At each step, the ratio (Karlin/|correction|) needs to stay above 1. Empirically it INCREASES with s (1.26, 1.32, 1.41, ...).

**Conjecture**: if E_acc ≽ J_acc with ratio gap δ (meaning min_k d_k(E_acc, J_acc)/(E_acc_k · J_acc_k) ≥ δ), then after processing one more child, the new ratio gap δ' ≥ δ·(1 + ε) for some ε > 0 depending on the child.

This would give: each child processed INCREASES the ratio gap, making the correction relatively smaller. The bottleneck is s=1 (smallest ratio), which is PROVED.

Can you formalize this? The ratio gap at s=1 is the binomial smoothing gap from (1+x)^ℓ, which is at least ℓ/(mode+1) or similar. Each subsequent child adds to this gap via the Karlin mechanism.

### Part 5: The Absolute Margin = 1 Phenomenon

For ALL s ≥ 2, the tightest w_k = 1 (not 0, but exactly 1). This happens at small trees where the polynomials have small coefficients.

If w_k ≥ 1 always (for s ≥ 2), and the minimum is achieved at specific small trees, perhaps the proof can be structured as:
- Base cases: verify w_k ≥ 1 for all trees up to some threshold n_0
- Induction: for n > n_0, show w_k ≥ 1 using the growth of coefficients

The growth argument: as the tree gets larger, E and J have larger coefficients, so the Karlin term grows faster than the correction. The minimum w_k = 1 is "pinned" at small trees and cannot decrease.

Is this viable? The pendant-star (s=1) has w_k = 2m → ∞, which is monotone in tree size. For s ≥ 2, the absolute margin stays at 1 but occurs at different (small) trees for each n.

## Deliverables

1. Clean statement of the s=1 proof, identifying which steps fail for s ≥ 2
2. Explicit bound on the s=2 correction term
3. Analysis of child ordering: does ordering freedom help?
4. Product closure conjecture: does the ratio gap increase with each child?
5. Absolute margin argument: can w_k ≥ 1 be proved by growth + base cases?
6. If you can close the s=2 case: a complete proof sketch
