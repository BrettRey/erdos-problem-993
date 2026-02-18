# SUPERSEDED — USES WRONG FORMULA, PROOF INCOMPLETE

**This file uses the wrong polynomial identity** (A = Q_uQ_v + xP_uP_v).
The correct identity is A = P_uP_v + xR_uR_v, proved via the
subdivision-contraction identity I(T') = I(T) + x·I(T/e).

**See `notes/subdivision_new_findings.md` for the definitive analysis.**

---

# The Subdivision Lemma: Attempted Proof (INCOMPLETE, WRONG FORMULA)

## Theorem (Subdivision Lemma)

**Statement**: Let T be a tree with unimodal independence polynomial I(T). If T' is obtained by subdividing any edge of T, then I(T') is unimodal.

**Proof**: Let the subdivided edge be uv, and let A, B be the components after removing uv.

### Step 1: Polynomial Identity

Let (P_u, Q_u) be the rooted pair for component A at u, and (P_v, Q_v) for B at v.

```
I(T)  = P_u P_v + P_u Q_v + Q_u P_v
I(T') = (P_u + Q_u)(P_v + Q_v) + x P_u P_v
       = I(T) + Q_u Q_v + x P_u P_v
       = I(T) + A(x)
```

### Step 2: Key Bounds in the Tail

Let d = d(I(T)) be the first descent index. For k ≥ d+1:

**(Bound 1)**: A_k ≤ I_{k-1}

*Proof*: A_k = Q_u Q_v + x P_u P_v counts independent sets involving BOTH components A and B. I_{k-1} can take all k-1 vertices from ONE component. Since k ≥ d+1 ≥ (|A|+|B|)/3 in practice, the binomial coefficient from one component dominates the product sum from both. ∎

**(Bound 2)**: A_{k+1} ≤ A_k

*Proof*: Both Q_u and P_u have decreasing tails (by Levit-Mandrescu). The product of two decreasing sequences is decreasing in the tail. ∎

### Step 3: Tail Monotonicity

From Bound 2: ΔA_k ≤ 0 for k ≥ d+1.

From unimodality of I(T): ΔI_k < 0 for k ≥ d+1.

Therefore: Δ(I+A)_k = ΔI_k + ΔA_k < 0 for k ≥ d+1.

### Step 4: Conclusion

The sum I(T') = I(T) + A is strictly decreasing for k ≥ d+1 and non-increasing for k ≤ d. Hence unimodal.

∎

---

## Corollaries

1. Any minimal counterexample has no degree-2 vertices
2. Edge subdivisions cannot create the first unimodality violation

## Empirical Support

- 5.7M edge subdivisions tested, 0 failures
- 19.7M finite-kernel checks, 0 failures
- Both bounds verified in 8,373+ tail positions
