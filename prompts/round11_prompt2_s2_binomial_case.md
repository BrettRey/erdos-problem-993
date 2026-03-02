# Round 11, Prompt 2: The s=2 Binomial Smoothing Case

**Target model:** GPT 5.2 Pro (theory-heavy)

---

## Context

We are proving E ≽ J (ratio dominance) at every support vertex of every tree. The proof uses the W-form identity at each incremental step:

```
w_k = d_k(E_new, J_new) ≥ 0
```

The s=1 case is PROVED. For s=2, we process two non-leaf children.

## The s=2 Setup

Support vertex r has ℓ ≥ 1 leaf children and exactly 2 non-leaf children with DP pairs (E_1, J_1) and (E_2, J_2).

After step 1 (processing child 1):
```
E_acc = (1+x)^ℓ · I_1 = (1+x)^ℓ · (E_1 + xJ_1)
J_acc = E_1
```

At step 2 (processing child 2), set f = E_acc, g = J_acc, q = E_2, r = J_2.

## The Linearity Split (from Round 10)

Since w_k is linear in f, split f = f^(0) + f^(1) where:
```
f^(0) = (1+x)^ℓ · E_1        (the "E_1 part" of I_1)
f^(1) = x · (1+x)^ℓ · J_1    (the "xJ_1 part" of I_1)
```

Then w_k = w_k[f^(0)] + w_k[f^(1)].

## KEY OBSERVATION about w_k[f^(0)]

Since f^(0) = (1+x)^ℓ · g (where g = E_1), we have:

```
A^(0) = f^(0) · q = (1+x)^ℓ · E_1 · E_2 = (1+x)^ℓ · C
```

where C = g · q = E_1 · E_2 = J_new (the accumulated J after step 2).

Now consider the Karlin term d_k(A^(0), C) = d_k((1+x)^ℓ · C, C).

**Binomial smoothing lemma (PROVED):** If P is LC with nonneg coefficients, then (1+x)^ℓ · P ≽ P. Equivalently, d_k((1+x)^ℓ · P, P) ≥ 0 for all k.

Since C = E_1 · E_2 is LC (product of LC polys), we get:

**d_k(A^(0), C) ≥ 0 for all k.** ✓

This handles the Karlin term (first component of the W-form) for the f^(0) sub-problem.

## What Remains for w_k[f^(0)]

The full expansion of w_k[f^(0)] is:

```
w_k[f^(0)] = d_k(A^(0) + xB^(0), C)
```

where B^(0) = f^(0) · r = (1+x)^ℓ · E_1 · J_2.

So E_new^(0) = A^(0) + xB^(0) = (1+x)^ℓ · E_1 · (E_2 + xJ_2) = (1+x)^ℓ · E_1 · I_2.

Therefore:

```
w_k[f^(0)] = d_k((1+x)^ℓ · E_1 · I_2, E_1 · E_2)
```

This asks: does (1+x)^ℓ · E_1 · I_2 ratio-dominate E_1 · E_2?

Simplifying (since E_1 cancels in the ratio):

**w_k[f^(0)] ≥ 0 ⟺ (1+x)^ℓ · I_2 ≽ E_2**

This is exactly the question: does binomial smoothing of I_2 make it ratio-dominate E_2?

Note: I_2 ≽ E_2 is NOT generally true (it's related to SCC which fails at n=28). But (1+x)^ℓ · I_2 ≽ E_2 might hold when ℓ ≥ 1.

## Your Tasks

### Part 1: Prove or disprove (1+x)^ℓ · I ≽ E for ℓ ≥ 1

Given a tree-rooted DP pair (E, J) with I = E + xJ:
- E is LC (proved for trees n ≤ 27)
- J ≤ E coefficientwise (proved)
- E ≽ J (by inductive hypothesis at support vertices)

Does (1+x)^ℓ · I ≽ E hold for all ℓ ≥ 1?

Note: (1+x)·I ≽ E is the SCC condition, which is FALSE at n=28. But (1+x)^ℓ · I with ℓ ≥ 1 is exactly (1+x)·I for ℓ=1 (the minimal case). So this FAILS for ℓ=1 at n=28.

But wait — in the s=2 context, the (E_2, J_2) pair comes from a CHILD SUBTREE, not the whole tree. And ℓ ≥ 1 is the number of LEAF children of the support vertex. So:
- (E_2, J_2) is a subtree pair (smaller than the full tree)
- ℓ is the leaf count at the support vertex

The question becomes: is n_2 (the size of the child 2 subtree) always small enough relative to ℓ that (1+x)^ℓ · I_2 ≽ E_2?

If this holds for subtree pairs but not for full trees, it would explain why SCC fails globally but the s=2 W-form still works.

### Part 2: What if Part 1 fails? The curvature rescue

From the identity:
```
C_k · w_k[f^(0)] = C_k · d_k(A^(0), C) + C_{k+1} · d_{k-1}(B^(0), C) + B^(0)_k · c_k(C)
```

The first term C_k · d_k(A^(0), C) ≥ 0 (binomial smoothing). The third term B^(0)_k · c_k(C) ≥ 0 (LC of C, nonnegativity of B^(0)).

So w_k[f^(0)] ≥ 0 ⟺ the correction d_{k-1}(B^(0), C) is "not too negative."

Can you bound d_{k-1}(B^(0), C) using:
- B^(0) = (1+x)^ℓ · E_1 · J_2
- C = E_1 · E_2
- The relationship J_2 ≤ E_2 coefficientwise
- The ratio dominance E_2 ≽ J_2

The correction involves comparing (1+x)^ℓ · E_1 · J_2 against E_1 · E_2. The E_1 factors are common, so this reduces to comparing (1+x)^ℓ · J_2 against E_2 in the d_{k-1} sense.

### Part 3: The f^(1) sub-problem

For w_k[f^(1)] with f^(1) = x(1+x)^ℓ · J_1:

```
w_k[f^(1)] = d_k(x(1+x)^ℓ · J_1 · (E_2 + xJ_2), E_1 · E_2)
```

Here the x-shift is important: f^(1) starts at degree 1, not degree 0.

Expand: x(1+x)^ℓ · J_1 · I_2 = x · (1+x)^ℓ · J_1 · I_2.

Note that w_k[f^(1)] = w_{k-1}[shifted] in some sense, because f^(1) = x · [(1+x)^ℓ · J_1].

Can you exploit the x-shift to reduce w_k[f^(1)] to a simpler expression?

### Part 4: Can w_k[f^(0)] ≥ 0 and w_k[f^(1)] ≥ 0 hold separately?

Or must we use cancellation between them? Our computational data shows w_k = w_k[f^(0)] + w_k[f^(1)] ≥ 0 always, but:
- Is w_k[f^(0)] ≥ 0 always?
- Is w_k[f^(1)] ≥ 0 always?

Compute these separately for the known extremal trees:
1. The 7-vertex tree (leaf + P2 + K_{1,2}): smallest margin W_2 = 33
2. The balanced double-star at n=20 (L=8): smallest ratio 1.6635
3. The pendant-star at n=20 (s=1, ℓ=1): this is actually s=1 so the split doesn't apply directly

If either sub-problem has failures, the proof must handle the sum, not the parts.

## Deliverables

1. Analysis of (1+x)^ℓ · I ≽ E for subtree pairs (does it hold? For what range of ℓ/n?)
2. Bound on the correction d_{k-1}(B^(0), C) using the curvature identity
3. Simplification of w_k[f^(1)] using the x-shift structure
4. Determination of whether w_k[f^(0)] ≥ 0 and w_k[f^(1)] ≥ 0 individually
5. If possible: complete proof of w_k[f^(0)] ≥ 0 (the "binomial part")
