# Round 15, Prompt 3: Direct Proof of w_k ≥ 0 at Step 2

**Target model:** GPT 5.2 Pro (theory-heavy)

---

## Context

We are proving E ≽ J (prefix ratio dominance) at support vertices of trees via the W-form induction. The induction processes non-leaf children one at a time.

**Step 1** (one non-leaf child processed): **PROVED.** Handles 63% of support vertices.

**Step 2** (two non-leaf children processed): **THE HARDEST STEP.** Exhaustive computation through n=17 and sampling to n=24 confirms that the minimum W ratio (Karlin/|correction|) across ALL trees, ALL steps, ALL k values is always achieved at step 2.

The minimum ratio at step 2 converges to **1+√3 ≈ 2.732**, achieved by the star+star extremal family. This is bounded away from 1, so a ratio-based proof is viable.

**Steps 3+** have higher ratios (the accumulated pair becomes "smoother" with more factors).

So: **proving w_k ≥ 0 at step 2 is the key remaining challenge.**

## Step-2 Structure

At step 2, the accumulated pair is:
```
E_acc = (1+x)^ℓ · I_1 = (1+x)^ℓ · (E_1 + x·J_1)
J_acc = E_1
```

where (E_1, J_1) is the tree-realizable pair from the first non-leaf child, and ℓ ≥ 1 is the number of leaves at the support vertex.

Now we process the second non-leaf child with pair (g, h) = (E_2, J_2):
```
A = E_acc · g = (1+x)^ℓ · I_1 · E_2
B = E_acc · h = (1+x)^ℓ · I_1 · J_2
C = J_acc · g = E_1 · E_2
```

The W-form:
```
w_k = d_k(A, C) + (B_k·C_k - B_{k-1}·C_{k+1})
    = Karlin_k + correction_k
```

## Known Properties

1. **E_acc ≽ J_acc** (from step 1 proof): d_k((1+x)^ℓ·I_1, E_1) ≥ 0 for prefix k
2. **g is LC**: E_2 is LC (product of LC factors)
3. **h ≤ g coefficientwise**: J_2 ≤ E_2 (proved)
4. **h is LC**: J_2 is LC (product of LC factors)
5. **g(0) = h(0) = 1**
6. **STP2(g, h)**: J_2(m+1)/J_2(m) ≤ E_2(n+1)/E_2(n) for m > n (verified, 0 failures)
7. **I_1 = E_1 + x·J_1**: the first child's IS polynomial
8. **C = E_1·E_2** is a product of two LC polynomials, hence LC
9. **(1+x)^ℓ has binomial coefficients**: highly structured

## The Star+Star Extremal

When child 1 and child 2 are both stars (centers with only leaf children):
```
E_1 = (1+x)^a,  J_1 = [1],  I_1 = (1+x)^a + x
E_2 = (1+x)^b,  J_2 = [1],  I_2 = (1+x)^b + x
```

At this extremal:
- h = J_2 = [1] (trivial include-root polynomial)
- B = E_acc · [1] = E_acc
- C = (1+x)^a · (1+x)^b = (1+x)^{a+b}
- The correction B_k·C_k - B_{k-1}·C_{k+1} involves only E_acc and binomial C

**Result (PROVED):** At the star+star, for k=2:
```
w_2 = s(s-1)·Q(a,b)/6 where Q = (a-(b-1)/2)² + (3b²+6b+35)/4 > 0
```
Rescue ratio = 1+√3 at the limit (infimum, never achieved).

## Your Tasks

### Part 1: Simplify the step-2 Karlin term

The Karlin term is:
```
Karlin_k = d_k(E_acc·g, E_1·g)
         = d_k((1+x)^ℓ·I_1·g, E_1·g)
```

Since I_1 = E_1 + x·J_1:
```
E_acc·g = (1+x)^ℓ·(E_1 + x·J_1)·g = (1+x)^ℓ·E_1·g + x·(1+x)^ℓ·J_1·g
```

So:
```
Karlin_k = d_k((1+x)^ℓ·E_1·g + x·(1+x)^ℓ·J_1·g, E_1·g)
```

By linearity of d_k in the first argument:
```
Karlin_k = d_k((1+x)^ℓ·E_1·g, E_1·g) + d_k(x·(1+x)^ℓ·J_1·g, E_1·g)
```

The first part: d_k((1+x)^ℓ·F, F) where F = E_1·g = C.
This equals: [(1+x)^ℓ·F]_{k+1}·F_k - [(1+x)^ℓ·F]_k·F_{k+1}.

Since (1+x)^ℓ·F ≽ F iff F is LC (proved: the binomial smoothing lemma), this first part ≥ 0.

The second part: d_k(x·(1+x)^ℓ·J_1·g, E_1·g).
This involves J_1 (the include-root poly of child 1), which is ≤ E_1 coefficientwise.

**Can you bound the second part?** Is it always ≥ 0? Is it bounded below by some explicit function of the first part?

### Part 2: Simplify the step-2 correction

The correction is:
```
correction_k = B_k·C_k - B_{k-1}·C_{k+1}
```

where B = (1+x)^ℓ·I_1·J_2 and C = E_1·E_2.

Since I_1 = E_1 + x·J_1:
```
B = (1+x)^ℓ·(E_1 + x·J_1)·J_2 = (1+x)^ℓ·E_1·J_2 + x·(1+x)^ℓ·J_1·J_2
```

Split the correction similarly:
```
correction_k = [(1+x)^ℓ·E_1·J_2]_k·C_k - [(1+x)^ℓ·E_1·J_2]_{k-1}·C_{k+1}
             + [x·(1+x)^ℓ·J_1·J_2]_k·C_k - [x·(1+x)^ℓ·J_1·J_2]_{k-1}·C_{k+1}
```

The first part: correction from (1+x)^ℓ·E_1·J_2 against C = E_1·E_2.
Since J_2 ≤ E_2, we have (1+x)^ℓ·E_1·J_2 ≤ (1+x)^ℓ·E_1·E_2 = (1+x)^ℓ·C coefficientwise.

**Key question:** Can you show that the correction from the first part is dominated by the corresponding Karlin from the first part? I.e., a "factor-by-factor" domination?

### Part 3: The w_k identity in terms of (E_1, J_1, E_2, J_2)

Combine Parts 1 and 2 to write w_k entirely in terms of the two child pairs:

```
w_k = [binomial Karlin from (1+x)^ℓ]
    + [cross Karlin from J_1]
    + [correction from E_1·J_2]
    + [correction from J_1·J_2]
```

Each piece involves products of child polynomials with binomial weights.

At the star+star extremal (J_1 = J_2 = [1]):
- Cross Karlin = d_k(x·(1+x)^ℓ, E_1·E_2) (involves just the binomial and the product C)
- Correction from E_1·J_2 = [(1+x)^ℓ·E_1]_k·C_k - [(1+x)^ℓ·E_1]_{k-1}·C_{k+1}
- Correction from J_1·J_2 = [x·(1+x)^ℓ]_k·C_k - [x·(1+x)^ℓ]_{k-1}·C_{k+1}

This is the simplest possible version (J = [1] eliminates most complexity).

**Can you prove w_k ≥ 0 for the star+star at general k?** Not just k=2 (already done) but all k in the prefix. The computational scan shows min w_k = 1 at (a,b)=(2,2), k=4.

### Part 4: Why does nonzero J make things easier?

The star+star extremal has J_1 = J_2 = [1]. For general trees, J_i has more nonzero coefficients.

Computationally: the W ratio is HIGHER for trees with nontrivial J_i (the star+star is the minimum). Why?

Intuition: when J_i has more nonzero coefficients:
- h = J_2 is "closer" to g = E_2 (more overlap in coefficient patterns)
- This makes the correction more "aligned" with the Karlin term
- The mismatch between B and C is reduced

Can you formalize: "the correction is bounded by Karlin · J_2(1)/E_2(1)" or some similar monotone function of the J/E ratio?

At the extremal: J_2 = [1], so J_2(1)/E_2(1) = 0/b = 0 (minimal "alignment").
For a path: J/E ratio is higher → ratio is higher. Makes sense.

### Part 5: Induction on J-complexity

Consider the following proof strategy:

**Base case:** J_1 = J_2 = [1] (star+star). PROVED (w_2 ≥ 0, all-k verified).

**Inductive step:** If w_k ≥ 0 holds when J_i has coefficients up to degree d-1, then it holds when J_i has degree d.

Alternatively: parametrize by t ∈ [0,1] and consider J_i(t) = (1-t)·[1] + t·J_i. Show w_k(t) is nondecreasing in t. Then w_k(1) ≥ w_k(0) ≥ 0.

Is this monotonicity true? At t=0, J_i = [1] (star+star, proved). At t=1, J_i = actual J_i.

This would reduce the general case to the star+star base case.

### Part 6: Factor-level sufficient condition

Instead of proving w_k ≥ 0 directly, find a sufficient condition on individual factors that implies w_k ≥ 0 and is preserved under tree DP.

Candidate: "the correction contributed by factor i is bounded by a fixed fraction of the Karlin contributed by factor i."

Since the W ratio is ≥ 1+√3 at the extremal, there is room for a factor-level bound. If each factor contributes correction ≤ (1/(1+√3)) × Karlin, summing gives total correction ≤ (1/(1+√3)) × total Karlin < total Karlin.

But this assumes the contributions decompose additively by factor, which is only approximately true.

## Deliverables

1. Simplified Karlin term at step 2 (binomial part + cross part)
2. Simplified correction at step 2 (E_1·J_2 part + J_1·J_2 part)
3. Star+star all-k proof (extend beyond k=2)
4. Why nontrivial J increases the ratio (formal or semi-formal argument)
5. Assessment of the J-complexity induction approach
6. Factor-level sufficient condition analysis
7. Overall assessment: which approach is most likely to close step 2?
