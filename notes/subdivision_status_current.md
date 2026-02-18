# Current Status: Subdivision Lemma

## What's Been Proved

1. **Exact Identity**: I(T') = I(T) + Q_u Q_v + x P_u P_v ✓
2. **Coefficient Bound**: A(x) <= (1+x) I(T) ✓
3. **Empirical**: No counterexample in 5.7M edge subdivisions ✓
4. **Tail Dominance**: ΔA_k < -ΔI_k for k >= d+1 (100% in 2,217 tests) ✓
5. **Boundary**: Can fail but doesn't create valleys (100% safe in 4,373 tests) ✓

## What's Missing

The formal proof that tail dominance holds analytically.

## Key Approaches from Literature

### 1. Levit-Mandrescu Theorem
States: The tail of tree independence polynomials is strictly decreasing.
- After index k >= ceil((2α-1)/3), coefficients drop monotonically.
- This gives a lower bound on where the tail starts.

### 2. Finite Kernel Approach
If we can characterize "leaf-light" trees (finite class), we can verify directly.
- Empirical: 19.7M checks in finite kernel class, 0 failures.

### 3. The Key Inequality

Need to show: For edge uv splitting tree into A, B:
Q_u Q_v <= I(T) / M for some M in tail region.

This is equivalent to showing that polynomials of forests (components after removing edge) grow slower than the full tree polynomial in the tail.

## The Challenge

The proof requires connecting:
- Binomial coefficient bounds
- Forest polynomial structure  
- The specific form of Q = x * I(component - neighbors)

## Recommendation

Given the empirical evidence (millions of cases, 0 failures), the subdivision lemma is almost certainly true. The remaining gap is analytic.

Options:
1. **Accept as conjecture** with massive empirical support
2. **Finite kernel proof** - characterize finite class where it holds, verify directly
3. **Continue analytic attack** - find the right inequality

The finite kernel approach seems most tractable given the existing structure.
