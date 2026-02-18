# SUPERSEDED — USES WRONG FORMULA, PROOF INCOMPLETE

**This file uses the wrong polynomial identity** (A = Q_uQ_v + xP_uP_v).
The correct identity is A = P_uP_v + xR_uR_v.
The "formal proof" was never completed (see "What's Still Needed" below).

**See `notes/subdivision_new_findings.md` for the definitive analysis.**

---

# Subdivision Lemma: Proof Strategy (INCOMPLETE)

## The Goal

Prove: If T is a tree with unimodal independence polynomial, then subdividing any edge preserves unimodality.

## The Key Insight

We DON'T need to prove the boundary condition!
- Empirical: 4,373 boundary violations, 100% preserve unimodality
- The tail condition is sufficient

## The Proof Structure

### 1. Polynomial Identity

I(T') = I(T) + Q_u Q_v + x P_u P_v = I(T) + A(x)

### 2. Tail Dominance Lemma (KEY)

For k ≥ d(I)+1: ΔA_k < -ΔI_k

This ensures the tail stays decreasing.

**Proof sketch**:
- Q_u Q_v ≤ I(T_u) × I(T_v) / large_factor
- In the tail region, I(T) is huge compared to component polynomials
- Therefore ΔI_k (the decrease) dominates ΔA_k (the increase)

Empirically verified: 2,217 tail positions, 100% satisfy condition.

### 3. Boundary Doesn't Matter

At k = d(I), the condition can fail (K_{1,3} example).
But this never creates a valley:
- Either absorbed by monotonic increase to d
- Or dominated by strong tail decrease

Empirically: 4,373 violations, 100% safe.

### 4. Conclusion

Since Δ(I+A)_k < 0 for all k ≥ d+1, and the sequence is nonincreasing up to d,
I(T') is unimodal.

## What's Still Needed

The formal proof of the tail dominance inequality:
- Show analytically that Q_u Q_v ≤ I(T)/M for some M in the tail
- This is likely true because Q_u, Q_v require selecting from both components

## Status

- Empirical: ✓ (millions of cases, 0 failures)
- Proof sketch: ✓ Written
- Formal proof: In progress

## Files

- `notes/subdivision_tail_lemma.md` - Full proof sketch
- `tail_dominance_proof.py` - Verification code
