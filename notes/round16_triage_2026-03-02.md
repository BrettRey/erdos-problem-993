# Round 16 Triage (2026-03-02)

**Models:** Codex 5.3 (P1), GPT 5.2 Pro (P2, P3)

## P1: 4×4 N_c Matrix and Condition C Verification (Codex)

### N_c correction (CRITICAL)
- **Original N_c was WRONG**: σ(FG) = σ(F) + F·σ(G), NOT σ(F)·G + F·σ(G).
- Corrected N_c has 1's on (1,1) and (3,3) diagonal (not I_c, E_c).
- Counterexample: 3-vertex star, σ((1+x)²) = 2+x but wrong formula gives 2+2x.
- **Verified**: 20,038 child steps n≤12, 0 failures with corrected matrix.

### Ladder-minor monotonicity (KEY FINDING)
- Λ_k = E_acc(k)·J_acc(k) - E_acc(k-1)·J_acc(k+1)
- **ALWAYS ≥ 0**: confirmed at all steps.
- **MONOTONICALLY NONDECREASING**: 5,502,481 full-range (k,step) pairs n≤17, 0 decreases.
- Also 2,590,788 prefix-only checks, 0 decreases.
- This is strictly STRONGER than preservation (preservation = stays ≥ 0; monotonicity = grows).

### Condition C: DEAD
- θ_k EXCEEDS 1 starting at n=8.
- Grows to ~37.85 at n=17 (max across all trees/steps/k).
- Star+star extremal: θ diverges quadratically (~m²/4 balanced, ~3m/2 unbalanced).
- Global max θ for a,b≤100: ~4875.87 at (98,100,99).
- **Condition C cannot work.** The negative term is NOT controlled by curvature rescue.

### Other findings
- Block-Toeplitz TN2 of N_c entries fails structurally.
- **Assessment**: ladder-minor approach >> Condition C.

## P2: Ladder-Minor Cauchy-Binet (GPT)

### Exact CB formula
```
Λ_k^{new} = Σ_{i,j} Δ_{i,j}(A,B) · P(k-i)·Q(k-j)
```
where Δ_{i,j}(A,B) = A(i)B(j) - A(i-1)B(j+1), A=E_acc, B=J_acc, P=I_c, Q=E_c.

### Sign structure
- STP2(A,B) → j ≥ i ⟹ Δ_{i,j} ≥ 0 ("above-diagonal" nonneg)
- STP2(P,Q) → factor minors nonneg for i ≥ j region
- On diagonal (i=j): Δ_{i,i} = Λ_i^{old} ≥ 0 by IH
- "Complementary" sign control: the two STP2 conditions cover overlapping regions

### Factor minors
- μ_k(I_c, E_c) = c_k(E_c) + [J_c(k-1)E_c(k) - J_c(k)E_c(k-1)] + ...
- ν_k(I_c, E_c) for off-diagonal

### Step-2 specialization
- A = (1+x)^ℓ · I_1, B = E_1
- Δ_{i,j}(A,B) = binomial-weighted sum of (I_1,E_1) cross-minors

### W-form connection
- Λ_k(E_new,J_new) = Λ_k(A*Q, B*Q) + Λ_{k-1}(A*R, B*Q) — "Karlin + correction"
- Unifies ladder-minor and W-form frameworks

### Missing piece
- Below-diagonal control lemma: show Σ_{i>j} terms dominated by Σ_{j>i} terms
- Planar network interpretation feasible for width-2 but needs explicit construction

## P3: Condition C Closure (GPT)

### Expansion matches P2
- d_{k-1}(B,C) expansion as double sum agrees with P2's structure.

### Cauchy-Binet for LC gaps
- c_k(F*G) = Σ_{p<q} [Toeplitz minor of F] · [Toeplitz minor of G]
- Clean sum of products of nonneg minors (when F, G are LC).

### Induction loop analysis
- **Step 6 BROKEN**: E≽J + LC + J≤E does NOT imply STP2(E,J).
- STP2 must be proved DIRECTLY as a tree lemma, not derived from P2.
- The R14 counterexample fails E≽J prefix, so it doesn't apply when all 4 conditions hold.
- But the implication still needs a proof, not just absence of counterexample.

### STP2+LC+J≤E does NOT suffice for Condition C abstractly
- Counterexample: g=[1,2,3], h=[1,1].

### Assessment
- Condition C needs tree-realizability; cannot be proved from abstract axioms.
- Extremal principle (star+star = worst) is key missing piece but θ diverges anyway.
- **Condition C is DEAD. Focus entirely on ladder-minor CB.**

## Strategic Summary

The proof picture is now very clean:

1. **Target**: STP2(E,J) at every rooting of every tree (via ladder-minor Λ_k ≥ 0)
2. **Tool**: Cauchy-Binet expansion of Λ_k^{new}
3. **Known**: diagonal sum D_k ≥ 0 (from IH + nonneg weights)
4. **Known**: above-diagonal (j>i) cross terms have Δ_{i,j} ≥ 0 (from STP2(A,B))
5. **Remaining**: show cross sum X_k ≥ 0, or show D_k dominates -X_k
6. **Key structure**: I = E + xJ gives "pure E" (symmetric) + "J correction"
7. **Bonus**: ladder-minor is MONOTONICALLY nondecreasing (5.5M checks)

The single missing lemma is: the off-diagonal cross sum X_k in the CB expansion is nonneg.
