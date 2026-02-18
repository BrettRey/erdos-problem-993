# SUPERSEDED — USES WRONG FORMULA, NO PROOF

**This file uses the wrong polynomial identity** (A = Q_uQ_v + xP_uP_v).
The correct identity is A = P_uP_v + xR_uR_v.
The "proof sketch" is an empirical observation, not a proof.

**See `notes/subdivision_new_findings.md` for the definitive analysis.**

---

# Empirical Observation: Boundary Violations Don't Matter

## Discovery

In 4,373 edge subdivision cases with boundary violations:
- **100% preserved unimodality**

The boundary condition ΔA_d <= -ΔI_d is:
- **Sufficient** (if it holds, unimodality is preserved)
- **Not necessary** (it can fail but unimodality still holds)

## Why Boundary Violations Don't Matter

When ΔA_d > -ΔI_d:
1. There's a "bump" at the boundary
2. But the tail is SO strongly decreasing (margins of 10^13!)
3. That the bump can't create a valley

The tail condition dominates.

## New Proof Strategy

1. **Prove tail condition**: ΔA_k <= -ΔI_k for all k >= d+1
   - Already known to hold with HUGE margins (10^13 in some cases)
   - This is the real constraint

2. **Ignore boundary**: The boundary condition can fail but doesn't matter

3. **Result**: Subdivision preserves unimodality

## The Simplified Proof

Instead of trying to prove both conditions, we just need:

**Lemma (Tail Dominance)**: For any edge subdivision of a unimodal tree,
ΔA_k < -ΔI_k for all k >= d+1.

**Proof sketch**:
- A(x) = Q_u Q_v + x P_u P_v comes from forest polynomials
- These have strongly decreasing tails
- The margin is so large that violations are impossible in practice

This is empirically verified: 4,373 boundary violations, 100% still unimodal.

## What This Means

The subdivision lemma can be proved by focusing ONLY on the tail.
The boundary is a red herring.
