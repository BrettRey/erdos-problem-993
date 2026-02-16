# SUPERSEDED — USES WRONG FORMULA

**This file uses the wrong polynomial identity** (A = Q_uQ_v + xP_uP_v).
The correct identity is A = P_uP_v + xR_uR_v.

**See `notes/subdivision_new_findings.md` for the definitive analysis.**

---

# Proof Attack: Subdivision Lemma (SUPERSEDED)

## Current Status

The subdivision lemma states: if T is unimodal, then subdividing any edge preserves unimodality.

Empirically verified through n=19 (5.7M edge subdivisions) with no counterexamples.

## The Key Identity

For edge uv subdivided with new vertex w:
- I(T') = I(T) + Q_u Q_v + x P_u P_v = I(T) + A(x)

Where:
- P_u = polynomial for sets excluding u
- Q_u = polynomial for sets including u (with x factor)
- Same for v

## Sufficient Conditions for Unimodality of I + A

Let d = first descent index of I(T).

We need:
1. **Mode condition**: d(A) >= d (A doesn't introduce earlier descent)
2. **Tail dominance**: ΔA_k <= -ΔI_k for all k >= d+1

If both hold, then I+A is unimodal.

## What We Know Empirically

For all tested cases:
- Condition 1 holds: A is nondecreasing up to d(I)
- Condition 2 holds in the tail region
- The boundary at k=d can fail (K_{1,3} counterexample)

## The Path to Proof

### Key Insight: Forest Polynomial Structure

The polynomials p = P_u P_v and j = Q_u Q_v + P_u Q_v + Q_u P_v come from FORESTS (disconnected components after removing edges).

**Conjecture**: For forest independence polynomials, the coefficient ratios are monotonically decreasing in the tail region.

### Why This Might Be True

1. **Forest = collection of trees**: Each tree component has unimodal polynomial
2. **Products of unimodal sequences**: When you multiply two unimodal sequences, the result tends to be unimodal in the tail
3. **Levit-Mandrescu theorem**: The tail of tree independence polynomials is monotone decreasing

### Formal Approach

Let F be a forest with components C_1, ..., C_m. Then:
I(F) = ∏ I(C_i)

For each component, we know:
- The coefficients are unimodal
- In the tail (after the mode), they're strictly decreasing
- The ratio r_k = i_{k+1}/i_k is < 1 and decreasing

For the product: i_k = ∑_{k1+...+km=k} ∏ i_{kj}

In the tail region, the dominant contributions come from tuples where each component is in its decreasing regime. This might allow a proof by induction on the number of components.

## Attempted Proof Strategy

### Lemma: Product of Two Tail-Decreasing Sequences

Let (a_k) and (b_k) be two sequences that are:
- Positive for all k
- Unimodal with mode at m_a and m_b
- Strictly decreasing for k >= max(m_a, m_b) + 1

Then the product sequence (c_k) = (a * b)_k is also strictly decreasing for k >= max(m_a, m_b) + 1.

**Proof sketch**: 
For k >= max(m_a, m_b) + 1, both a_k and b_k are decreasing.
Consider c_{k+1}/c_k = (a_{k+1}b_{k+1})/(a_k b_k).
We need to show this <= 1.

Since both sequences are decreasing:
a_{k+1} <= a_k and b_{k+1} <= b_k

Therefore:
c_{k+1}/c_k = (a_{k+1}/a_k) * (b_{k+1}/b_k) <= 1 * 1 = 1

The inequality is strict unless both ratios equal 1, which happens only at the mode.

### Application to Forest Polynomials

If we can prove that P_u and P_v (the "exclude-root" polynomials) are tail-decreasing, then their product p = P_u P_v is tail-decreasing.

Similarly for the j polynomial if we can establish the right properties.

## Next Steps

1. **Prove tail-decreasing for P polynomials**: Show that for any rooted tree, the P polynomial (sets excluding root) has decreasing ratios in the tail.

2. **Handle the boundary**: The boundary condition at k=d can fail. We need to show this doesn't create a valley - maybe by analyzing the specific structure of A.

3. **Finite kernel fallback**: If the general proof is too hard, show it for the finite class of "leaf-light" forests, which covers all potential counterexamples.

## Empirical Support

- n<=19: 5.7M edge subdivisions, 0 failures
- Finite core (b0=8, λ0=4): 19.7M edge checks, 0 failures

The pattern is extremely consistent, suggesting a clean proof exists.
