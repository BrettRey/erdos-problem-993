# THE OPEN PROBLEM: STP2 Closure Under Convolution

## What this is

This is the core unsolved step in proving Erdős Problem #993: that the independent set sequence of every tree is unimodal. Everything else is proved. This one theorem would close it.

## The problem

Given two pairs of nonneg integer-coefficient polynomials (I₁, E₁) and (I₂, E₂), each satisfying the STP2 condition (ladder minors nonneg), prove that their convolution products also satisfy STP2:

```
STP2(I₁, E₁)  ∧  STP2(I₂, E₂)  ⟹  STP2(I₁*I₂, E₁*E₂)
```

where * is polynomial multiplication (convolution of coefficient sequences), and STP2(f,g) means f(k)·g(k) ≥ f(k-1)·g(k+1) for all k.

## Available hypotheses (all verified computationally)

The sequences aren't arbitrary — they come from the tree DP recurrence. Known properties:

1. **STP2(Iₘ, Eₘ)** for each child m — the ratio I(k)/E(k) is nonincreasing
2. **Log-concavity** of all four sequences I₁, E₁, I₂, E₂
3. **Nonneg coefficients** throughout
4. **Coefficientwise domination**: E ≤ I (because I = E + xJ with J ≥ 0)
5. **I(0) = 1** for all tree factors
6. **E is PF2** (Pólya frequency of order 2) — proved
7. **J ≤ E coefficientwise** — proved
8. **STP2(E, J)** also holds at each child — confirmed

## The Cauchy-Binet approach (most promising direction)

Define the ladder minor: Λ(f,g)(k) = f(k)·g(k) - f(k-1)·g(k+1).

The CB expansion gives an identity:

```
Λ(I₁*I₂, E₁*E₂)(k) = Σᵢ Σⱼ [I₁(i)E₁(j) - I₁(j)E₁(i)] · I₂(k-i) · E₂(k+1-j)
```

Split into:
- **D_k** (diagonal, j ≥ i): each term nonneg by STP2(I₁, E₁). Always positive.
- **X_k** (cross, j < i): each term nonpositive. Can be negative.

**The gap: prove D_k + X_k ≥ 0.**

Empirically (n ≤ 17, exhaustive): min D_k/|X_k| = 4.577, median 11.25. D_k ALWAYS dominates by a large margin. But no proof.

**Dual form (Form 2):** Same sum, different grouping:
```
Λ(I₁*I₂, E₁*E₂)(k) = Σᵢ Σⱼ I₁(i)E₁(j) · [I₂(k-i)E₂(k-j) - I₂(k-j)E₂(k-i+1)]
```
Nonneg for i ≥ j (by STP2(I₂, E₂)). Covers the complementary region.

**Key structural fact:** Gap-1 terms (|i-j| = 1) are ALWAYS nonneg (0 failures). Only gaps ≥ 2 can be negative.

## What has been tried and FAILED

**Do not attempt these — they are proven dead ends:**

1. ❌ **Each CB term nonneg**: FALSE. Cross terms X_k CAN be negative (2015 failures at n≤15).
2. ❌ **Symmetrized bracket F(i,j) ≥ 0**: FALSE (47.5% failure rate!).
3. ❌ **Pairwise symmetric sum S(i,j,k) ≥ 0**: FALSE (3.1% failures).
4. ❌ **Averaged CB (Form1+Form2)/2 for all (i,j)**: Fails for j < i, gap ≥ 2.
5. ❌ **Prefix-cumulative gap rescue**: FAILS with correct indexing.
6. ❌ **STP2 from abstract axioms alone**: FALSE (explicit counterexample exists: E=(1,19,77,72,50,22), J=(1,17,21,24,25)).
7. ❌ **STP2 generic product closure**: FALSE (not all STP2 pairs compose).
8. ❌ **Compound matrix Λ²**: Scalar sparsity doesn't survive coefficient-level unfolding (100% mismatch).
9. ❌ **Condition C / SCC**: Fails at n=28.
10. ❌ **Step-ratio monotonicity**: FALSE (22 finite-finite violations n≤17).

## What might work (unexplored or partially explored)

1. **Total positivity theory**: The sequences are PF2. Convolution of PF2 sequences is PF2. But STP2 is about PAIRS, not individual sequences. Is there a variation of the Cauchy-Binet theorem for TP matrices that handles this directly?

2. **Generating function approach**: Work with the generating functions directly rather than coefficient-by-coefficient. The W(x) = (1+x)E'(x) - DE(x) approach worked for star+star.

3. **Induction on polynomial degree**: For degree-1 factors (leaves), STP2 closure is trivial. Can a degree-induction work?

4. **Real-rootedness**: If I and E are real-rooted, their ratio I/E has interlacing roots, which might force STP2. Are tree-realizable I,E always real-rooted? (Stars are: (1+x)^n is real-rooted.)

5. **Schur functions / symmetric function theory**: The convolution product has representation-theoretic interpretations.

## What I need

Fill the `sorry` in `stp2_conv_closure` in `Formal/STP2Closure.lean`. The auxiliary sorries (conv_comm, conv_nonneg, lc_conv, cb_expansion_form1) can be filled or used as axioms — the main target is the closure theorem.

If you can't prove the full theorem, partial progress is valuable:
- Prove it for degree-1 factors (I = 1+x, E = x)
- Prove it for degree-2 factors
- Prove the gap-1 terms are nonneg
- Identify which hypothesis is doing the heavy lifting
- Suggest a proof strategy we haven't tried

This is a genuine open problem in combinatorics. A proof would resolve a 1987 Erdős conjecture.
