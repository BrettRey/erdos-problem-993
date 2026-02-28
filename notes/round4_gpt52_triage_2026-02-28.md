# Round 4 GPT 5.2 Triage (2026-02-28)

## Outputs received

Two outputs (SV1 prompt cancelled after our scan showed it false).

### Output 1: SV1 + Canonical max-ℓ (Instance 2)

**Confirmed findings:**
- SV1 is false (gives explicit 12-vertex counterexample, matching our scan)
- "A ≥ B + B LC ⟹ P2 for ℓ≥2" is FALSE as a purely algebraic claim (counterexample: B=1+2x+3x²+4x³, A=B+x, verified)
- Max-ℓ=1 trees are exactly "no strong support vertex" (leaf edges form a matching)

**New algebra (potentially valuable):**
- Binomial-minor expansion: Δ_k^(ℓ) = Σ_{t=0}^ℓ C(ℓ,t) · M_k^(t)
- Shifted minor recursion: b_{k-1} · M_k^(t) = b_k · M_{k-1}^(t-1) + a_{k-t} · c_k
- Explicit ℓ=2 decomposition (equation 3) with two curvature terms c_k and c_{k-1}
- "2-step compensated deficit" inequality (4) as sufficient condition for ℓ≥2
- Degree-2 support vertex tie-break for ℓ=1 regime

**Assessment:** Solid algebra, useful for the paper. The binomial-minor expansion generalizes our single-ℓ identity to arbitrary ℓ. The ℓ=2 counterexample correctly identifies that tree-specific structure is needed.

### Output 2: Condition C / Delta identity (Instance 3)

**Confirmed findings:**
- Clean algebraic proof of identity (†) -- matches our sympy verification
- Determinantal packaging: (C_k) = det(M^(1)_k) + det(M^(2)_k) -- new, useful

**CRITICAL ERROR in "counterexample to product closure":**
- GPT 5.2 tested the WEAK version of Condition C (d_k + (b_k/b_{k-1})·d_{k-1} + c_k ≥ 0)
- Our scans test the STRONG version (b_{k-1}·d_k + b_k·d_{k-1} + a_{k-1}·c_k ≥ 0)
- GPT's counterexample (A1=[1,2,1,6,1], B1=[1,1,1,1,1], etc.) PASSES the strong version
- The weak version fails at k=6 (value = -6), but the strong version gives +26 (a_{k-1}/b_{k-1} = 20 at that point)

**However, GPT's conclusion is directionally correct:**
- The STRONG version of Condition C is also NOT generically product-closed
- Random synthetic pairs (I=E+xJ, E LC, Cond C holds): 166 failures / 125K pairs
- BUT: adding the tree constraint J ≤ E coefficientwise eliminates ALL failures (0 / 125K)

**The J ≤ E constraint is the key.**

## The J ≤ E constraint

In the tree DP, for a non-leaf child c of root r:
- E_c = dp0[c] = ∏_{gc} I_{gc}  (product of grandchild IS polys)
- J_c = dp1s[c] = ∏_{gc} E_{gc}  (product of grandchild exclude polys)
- Since I_{gc} ≥ E_{gc} coefficientwise, the products satisfy E_c ≥ J_c coefficientwise

This means I_c - E_c = x·J_c where J_c ≤ E_c. The "include bonus" (from including c) is bounded by the "exclude" polynomial. This is a natural constraint: including the root adds at most as many sets as excluding it.

**Full constraint set for product closure (all verified computationally):**
1. E[0] = 1
2. I = E + x·J (DP structure)
3. J ≤ E coefficientwise (tree constraint)
4. E is log-concave
5. Strong Condition C holds for (I, E)

With 1-5: product closure holds (0 failures, 125K random pairs + 701K tree pairs).
Without constraint 3: fails (166 / 125K).

## Updated proof structure

**Theorem (to prove):** If (I_1, E_1) and (I_2, E_2) both satisfy constraints 1-5, then (I_1·I_2, E_1·E_2) satisfies constraints 1, 2, 4, and 5 (NOT 3).

**CORRECTION (verified):** Constraint 3 (J ≤ E) is NOT product-closed. Two leaf factors (J=E=1) give J'=[2,1] > E'=[1]. But Condition C holds for the product anyway. So the induction is: each FACTOR satisfies 1-5 (with J_c ≤ E_c from tree DP), and the PRODUCT satisfies 1, 2, 4, 5 (Condition C holds, but J' may exceed E').

**Note:** Constraints 1, 2, 4 are obviously multiplicative:
- (E_1·E_2)[0] = 1·1 = 1
- I_1·I_2 = (E_1+xJ_1)(E_2+xJ_2) = E_1E_2 + x(J_1E_2 + E_1J_2 + xJ_1J_2)
  So I = E + x·J' where J' = J_1E_2 + E_1J_2 + xJ_1J_2
- J' ≤ E_1E_2? **NO.** Even two leaf factors (J=E=1) give J'=[2,1] > E'=[1].
  The x·J_1·J_2 cross-term destroys J≤E at the product level.
  **BUT**: Condition C still holds for the product! So J≤E is needed only at
  the FACTOR level, not as a product-closed invariant.
- Constraint 4 (LC of product) is preserved by Cauchy products of LC sequences.
- Constraint 5 (strong Condition C) is the key step.

## Summary of GPT 5.2 Round 4

| Claim | Status |
|-------|--------|
| SV1 is false | **CONFIRMED** (our scan + GPT's explicit example) |
| A≥B + B LC ⟹ P2 for ℓ≥2 | **FALSE** (GPT's algebraic counterexample, verified) |
| Identity (†) | **CONFIRMED** (trivial algebra) |
| Determinantal packaging | **NEW**, useful for TP approach |
| Weak Cond C not product-closed | **CONFIRMED** (GPT's example) |
| Strong Cond C not product-closed | **TRUE for arbitrary pairs**, FALSE for tree-like pairs with J≤E |
| J ≤ E is the missing constraint | **NEW DISCOVERY** (our verification) |
| Binomial-minor expansion | **NEW**, useful generalization |
| ℓ=1 degree-2 tie-break | **NEW**, useful structural insight |
