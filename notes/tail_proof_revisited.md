# The Correct Proof Approach

## The Key Observation

The empirical data shows:
- Smaller component → larger ratio |ΔA|/|ΔI| (1.73 avg)
- Larger component → still only ~2.0 avg

But the DOMINANCE still holds 100%. This means ΔI_k is negative enough that even with |ΔA| being ~2x |ΔI| in ratio, the actual values make ΔI_k + ΔA_k < 0.

## The Real Proof Strategy

Instead of bounding ratios, we should argue about the ACTUAL VALUES:

### Lemma: Tail Decreases Faster Than A Can Increase

For k >= d+1:
- I(T) is in the DESCENT region, dropping rapidly
- A(x) = Q_u Q_v + x P_u P_v is in the ASCENT region (or flat)

The key: In the tail, I(T) drops by MANY units per step, while A can only increase by a few units.

## Formal Argument

1. **I(T) in tail**: After the mode, binomial coefficients drop exponentially.
   For k >= d+1 >= alpha/2:
   I(T)_k ≈ C(n, k) in worst case, decreasing exponentially with k.

2. **A(x) structure**: 
   - Q_u Q_v = x^2 * I(A-N[u]) * I(B-N[v])
   - P_u P_v = I(A-u) * I(B-v)
   
   These are products of smaller polynomials, so their coefficients grow more slowly.

3. **The inequality**:
   
   In the tail region, consider the ratio:
   
   r_k = A_k / I(T)_k
   
   We need to show r_{k+1} < r_k (A grows slower than I in tail)
   OR show ΔI_k + ΔA_k < 0 directly.

4. **Direct approach**:
   
   ΔI_k = I_{k+1} - I_k is negative and large in magnitude.
   ΔA_k = A_{k+1} - A_k is bounded.
   
   We need |ΔI_k| > |ΔA_k|.
   
   Empirically this always holds. Theoretically:
   - |ΔI_k| is at least the drop in binomial coefficients
   - |ΔA_k| is bounded by the size of the smaller component

## The Simple Bounding Argument

Let s = min(|A|, |B|) = size of smaller component.

Claim: |ΔA_k| <= s^2 for all k.

Reason: A_k counts independent sets involving BOTH components.
The maximum contribution comes from choosing s-2 vertices from the smaller component and s-2 from the larger, which is bounded by s^2.

Meanwhile, |ΔI_k| >= something much larger in the tail.

For the path case (worst for I):
- I(T) ≈ Fibonacci-like, still decreasing by 1 per step
- For s >= 3, |ΔI_k| >= 3 > s^2 = 9? No, this doesn't work.

Actually we need a different approach.

## The Actual Empirical Finding

The key is that the FIRST difference ΔA_d can be positive, but subsequent ΔA_k become negative quickly.

The PROPER approach:
1. Show the tail of A(x) is eventually decreasing
2. Show the transition at d is not steep enough to create valley

This is what the empirical 100% success shows.

## Next Step

We need to characterize exactly WHY A(x) has a decreasing tail in all cases. This might relate to the forest polynomial structure.
