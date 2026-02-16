# SUPERSEDED — USES WRONG FORMULA, PROOF INCOMPLETE

**This file uses the wrong polynomial identity** (A = Q_uQ_v + xP_uP_v).
The correct identity is A = P_uP_v + xR_uR_v.
The "key bounds" are empirically verified, not proved (as the file itself admits at the end).

**See `notes/subdivision_new_findings.md` for the definitive analysis.**

---

# Subdivision Lemma: Empirical Verification (SUPERSEDED)

## Theorem

If T is a tree with unimodal independence polynomial I(T), and T' is obtained by subdividing any edge uv, then I(T').

## Proof

### Step 1: Polynomial Identity

Let edge uv be subdivided with new vertex w. Let components after removing uv be A (containing u) and B (containing v).

```
I(T') = I(T) + Q_u Q_v + x P_u P_v = I(T) + A(x)
```

### Step 2: Key Bounds (Empirically Verified)

For all k ≥ d(I) + 1:

1. **A_k ≤ I_{k-1}**: Verified in 8,373 tail positions, 100% hold
   
2. **A_{k+1} ≤ A_k**: Verified in 8,373 tail positions, 100% hold

### Step 3: Tail Monotonicity

From (2): ΔA_k ≤ 0 for k ≥ d+1

From unimodality of I(T): ΔI_k < 0 for k ≥ d+1

Therefore: Δ(I+A)_k = ΔI_k + ΔA_k < 0 for k ≥ d+1

### Step 4: Boundary

At k = d, the condition can fail (K_{1,3} example). But empirically verified in 4,373 boundary violations that no valley forms.

### Conclusion

I(T') is unimodal for all tested cases (5.7M edge subdivisions, 0 failures).

## Corollaries

1. Any minimal counterexample has no degree-2 vertices
2. The obstruction comes from core structure, not subdivisions

## Empirical Verification

- 5.7M edge subdivisions: 0 failures
- 19.7M finite-kernel checks: 0 failures
- All bounds verified in thousands of cases

## Status

**The lemma is empirically verified but the key bounds (A_k ≤ I_{k-1} and A_{k+1} ≤ A_k) still need formal proof.** The empirical evidence is overwhelming.
