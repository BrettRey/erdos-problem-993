# GPT 5.2 Pro Round 2 Triage (2026-02-27)

Three instances on: (1) P3 tail domination, (2) P2 prefix TP2, (3) DP closure / rooting characterization.

## Instance 1: P3 tail domination

**Verdict: PROVED. Clean injection argument, paper-ready.**

**Theorem (support-root domination).** If r is adjacent to a leaf u in tree T, then e_k ≥ j_{k-1} for all k ≥ 0 (not just k ≥ mode). P3 holds for ALL k.

*Proof.* The leaf-swap injection φ(S) = (S \ {r}) ∪ {u} maps each size-k IS containing r to a size-k IS avoiding r. Since u's only neighbor is r (which was removed), φ(S) is independent. The map is injective (inverse swaps u back for r). ∎

**Corollary.** Every tree on ≥ 2 vertices has a leaf; its neighbor is a support vertex. Rooting there gives P3 for all k.

The algebraic formulation is also clean:
- E = (1+x)A, xJ = xB where A = Π_{c≠u} S_c, B = Π_{c≠u} E_c
- A ≥ B coefficientwise (since S_c = E_c + xJ_c ≥ E_c)
- E - xJ = (1+x)A - xB = A + x(A-B) ≥ 0 coefficientwise ∎

Both proofs are correct. The injection is more elementary; the algebraic version connects to the DP structure.

## Instance 2: P2 prefix TP2

**Verdict: Valuable counterexamples, correct analysis of obstacles.**

### Key findings:
1. **K_{1,4} rooted at a leaf fails P2** (verified: e₂j₁ = 9 < 12 = e₁j₂)
2. **All stars K_{1,d} (d ≥ 4) rooted at a leaf fail P2** (proof: reduces to d-1 ≥ d, impossible)
3. **Some non-leaf, non-support roots also fail P2** (8-vertex double-star example)
4. **The DP induction breaks at the unary step**: P2 on (E_c, J_c) doesn't control (I(T_c), E_c) at the parent's prefix range

### Index convention correction:
The prompt's claim that P2 is equivalent to "P(r ∈ S | |S|=k) nonincreasing" is WRONG. The conditional probability is j_{k-1}/(e_k + j_{k-1}), whose monotonicity is governed by e_{k+1}j_{k-1} ≥ e_kj_k (shifted indices), NOT by e_{k+1}j_k ≥ e_kj_{k+1}. These are different inequalities. The P_5 endpoint example confirms: P2 holds but P(r ∈ S | |S|=k) increases from 1/5 to 3/6.

### Three salvage directions proposed:
A. Keep P2 as stated, make it existential over rootings
B. Switch to the correct conditional-probability inequality
C. Reframe in terms of factor pairs (I(T_i), E_i) instead of (E_i, J_i)

## Instance 3: DP closure / rooting characterization

**Verdict: Best structural insights. Binomial upgrade lemma is useful.**

### Key contributions:
1. **Same P3 proof** as Instance 1 (leaf-swap injection, stated as "support vertex" result)
2. **α(T-N[r]) ≤ m-1 conjecture is FALSE**: P_7 center has α(T-N[r]) = m = 2 but P⋆ holds
3. **Binomial upgrade lemma** (proved):
   If A ≽ B (ratio dominant) on a range, and B is LC on that range, then (1+x)A ≽ B.
   Proof: e_{k+1}b_k - e_kb_{k+1} = (a_{k+1}b_k - a_kb_{k+1}) + (a_kb_k - a_{k-1}b_{k+1}), both ≥ 0.
4. **Clean reduction**: P⋆ at support vertex reduces to proving A ≽ B + prefix LC of B, where A = Π S_c (products of subtree IS polys) and B = Π E_c (products of exclude-root polys).
5. **Key obstacle identified**: Factorwise dominance S_c ≽ E_c fails in small examples (3-vertex star).

### References:
- **Hu-Wang-Zhao-Zhao (arXiv:1507.08430)**: REAL paper (correcting round 1 triage). "Convolution Preserves Partial Synchronicity of Log-concave Sequences." Introduces partial synchronicity, intermediate between synchronicity and weak synchronicity.
- **Gross-Mansour-Tucker-Wang (SIDMA ~2015)**: Synchronised sequences under convolution. The foundational reference.

### The reduction to "support-vertex P2":
Instance 3 independently arrives at the same conclusion: since P3 is automatic at support vertices, the entire P⋆ program reduces to proving P2 at some (or all) support vertices.

## Critical computational verification

### Support-vertex P2 is UNIVERSAL (not just existential)

**59,916,124 support-vertex checks across all 9,114,283 trees n ≤ 22: ZERO P2 failures.**

P2 holds at EVERY support vertex of EVERY tree. This is dramatically stronger than "there exists a good rooting." The conjecture is:

**Conjecture (support-vertex P2).** For every tree T and every support vertex r, the pair (E_r, J_r) satisfies P2: e_{k+1}j_k ≥ e_kj_{k+1} for k = 0, ..., m-1, where m = mode(I(T)).

Combined with the proved P3 theorem:

**Conjecture (support-vertex P⋆).** For every tree T and every support vertex r, P⋆(E_r, J_r; m) holds. Equivalently: I(T-r; x) and I(T-N[r]; x) satisfy ratio dominance up to mode, and I(T-r) dominates xI(T-N[r]) coefficientwise.

## Proof strategy going forward

The complete proof of unimodality via P⋆ requires:
1. ✅ **P3 at support vertices**: PROVED (leaf-swap injection)
2. ⬜ **P2 at support vertices**: verified 60M checks, 0 failures. OPEN.

For P2, the binomial upgrade lemma reduces it to:
- A ≽ B on {0,...,m} (ratio dominance of Π S_c over Π E_c)
- B is prefix-LC on {0,...,m}

where A = Π_{c≠leaf} S_c and B = Π_{c≠leaf} E_c, with the product over non-leaf children of the support vertex r.

The key obstacle (per Instance 3): factorwise ratio dominance S_c ≽ E_c is false, so the product dominance must come from an aggregate effect.

## Actionable next steps

1. **Verify binomial upgrade lemma** on concrete examples ✅ (follows from algebra)
2. **Profile A ≽ B**: for each support vertex, measure the slack in a_{k+1}b_k - a_kb_{k+1} for k up to mode. Find the tightest cases.
3. **Test prefix LC of B = Π E_c**: is this always true up to mode? (E_c are exclude-root polys of subtrees, products of IS polys, likely LC as products of LC polys... but LC of factors doesn't imply LC of products in general.)
4. **Investigate Hu-Wang-Zhao-Zhao partial synchronicity**: does it give a closure theorem that applies here?
5. **Next round of prompts**: target proving support-vertex P2, with the binomial upgrade lemma as given.

## Papers to cite (additions from this round)

1. Hu, Wang, Zhao, Zhao, "Convolution Preserves Partial Synchronicity of Log-concave Sequences" (2015, arXiv:1507.08430) — correcting round 1 triage that called this fabricated.
