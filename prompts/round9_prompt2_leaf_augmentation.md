# Round 9, Prompt 2: Prove Leaf-Augmentation LR Dominance

**Target model:** GPT 5.2 Pro or Codex 5.3 (both could contribute)

---

## The Conjecture

**Leaf-augmentation LR dominance:** For every tree T rooted at any vertex r, define:
- E = dp[r][0] (exclude-root IS polynomial)
- J = dp[r][1]/x (include-root IS polynomial divided by x)
- I = E + xJ (full IS polynomial)
- I_tilde = I + xE

Then I_tilde ≽ E (ratio dominance), i.e.:
```
(a_{k+1} + b_k) · b_k  ≥  (a_k + b_{k-1}) · b_{k+1}   for all k
```

where a_k = coeff of x^k in I, b_k = coeff of x^k in E.

Equivalently: g_k := c_k(E) + d_k(I, E) ≥ 0 for all k.

## Computational Status

**Verified:** 0 failures through n = 22 (ALL 9.1M trees, ALL rootings, not just support vertices). This is approximately 100M+ coefficient checks.

## Tree Interpretation

I + xE is the IS polynomial of the tree T⁺ obtained by attaching one new leaf to the root r. Specifically:
- At root r of T⁺, the exclude polynomial is (1+x)E (the new leaf can be in or out)
- The include polynomial is xJ (same as before, since the new leaf is adjacent to r)
- I(T⁺) = (1+x)E + xJ = E + xE + xJ = I + xE = I_tilde

So the conjecture says: for every tree T and every vertex r, the IS polynomial of T⁺ (with one extra pendant leaf at r) ratio-dominates E(T, r).

## Why This Matters

1. **Curvature rescue:** g_k ≥ 0 implies T3 ≥ |T1| in the Strong Condition C decomposition (curvature alone rescues the negative current minor). This is verified with 0 failures.

2. **Connection to E ≽ J:** The identity
   ```
   T1 + T3 = b_{k-1}(c_k + d_k) + j_{k-2} c_k = b_{k-1} g_k + j_{k-2} c_k
   ```
   shows that g_k ≥ 0 plus c_k ≥ 0 (LC of E) gives T1+T3 ≥ 0 immediately.

3. **Self-contained target:** This is a property of a single rooted tree, not of a product. No incremental induction is needed.

## Structural Properties Available

At a rooted tree (T, r):

1. **E is PF2** (LC, nonneg): proved (product of PF2 factors from subtrees)
2. **J ≤ E coefficientwise:** proved (J^(t) ≤ E^(t) by induction)
3. **I = E + xJ:** definition
4. **I_tilde = I + xE = (1+x)E + xJ**
5. **If r is a leaf:** E = [1], J = dp of parent's subtree excluding r. Trivial.
6. **If r has children c_1,...,c_d:**
   ```
   E = ∏ I_c  (product of subtree IS polys)
   J = ∏ E_c  (product of subtree exclude polys)
   ```
7. **E/J structure:** E_k/J_k = ∏ (I_c)_k / ∏ (E_c)_k doesn't factor nicely

## Suggested Approaches

### Approach A: Induction on the Number of Vertices

Base cases (n=1,2) are trivial. For the induction step, remove a leaf ℓ.

If ℓ is not adjacent to r: T − ℓ has the same root, and I_tilde and E change by removing ℓ's contribution from the product. Can we show g_k is preserved?

If ℓ is adjacent to r: E gains a (1+x) factor, J gains an E_ℓ = [1] factor (trivial). The relationship between (I_tilde, E) for T and for T − ℓ involves removing the (1+x) from E.

### Approach B: Product Expansion

Since E = ∏ I_c and J = ∏ E_c:
```
I_tilde = (1+x) ∏ I_c + x ∏ E_c
E = ∏ I_c
```

Can we show d_k(I_tilde, E) ≥ 0 by expanding the product and using properties of the factors?

For a single child (d=1): I_tilde = (1+x)I_1 + xE_1 = (1+x)(E_1+xJ_1) + xE_1 = (1+2x)E_1 + x(1+x)J_1, and E = E_1+xJ_1 = I_1. Then g_k = c_k(I_1) + d_k((1+x)I_1+xE_1, I_1). This involves the LC of I_1 and the LR minor.

### Approach C: Direct Coefficient Analysis

Write a_k = e_k + j_{k-1} (where I = E + xJ). Then:
```
I_tilde_k = a_k + b_{k-1} = e_k + j_{k-1} + e_{k-1}
```

and g_k = (e_{k+1}+j_k+e_k)e_k - (e_k+j_{k-1}+e_{k-1})e_{k+1}
= e_k² - e_{k-1}e_{k+1} + j_k e_k - j_{k-1}e_{k+1}
= c_k(E) + d_k(xJ, E)
= c_k(E) + (j_k e_k - j_{k-1} e_{k+1})

So g_k ≥ 0 iff the LC surplus of E at k exceeds the "J grows faster than E" deficit at k.

This says: either E is "LC enough" at k, or J doesn't outrun E from k-1 to k.

### Approach D: TP2 / Variation Diminishing

Note that I_tilde = (1+x)E + xJ = E + x(E+J) = E + xI. The matrix [I_tilde; E] = [E+xI; E] has the structure of adding a "shifted I" to the first row. Can the theory of variation-diminishing operators help?

## Deliverables

1. A proof of g_k ≥ 0 for all tree rootings, or identification of the precise obstruction.
2. If the proof uses induction on n: the induction step and which leaf removal works.
3. If the proof uses product structure: the key inequality at the factor level.
4. Assessment of which approach (A-D) is most promising and why.
