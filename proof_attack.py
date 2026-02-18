#!/usr/bin/env python3
"""Analytical attack on the correction condition.

KEY INSIGHT I WANT TO TEST:

At a gap=+1 merge, mode(total) = m_A + 1. This means:
  (A+B)[m_A+1] >= (A+B)[k] for all k
where A = dp0, B = dp1 = x*P, P = prod(dp0[ci]).

So A[m_A+1] + P[m_A] >= A[m_A] + P[m_A-1]  ... (*)

Since P[m_A] >= P[m_A-1] (from gap=+1), P is still increasing at m_A.
Let m_P = mode(P). Then m_A < m_P, so m_A <= m_P - 1.

Now B = x*P, so mode(B) = 1 + m_P.

ATTEMPT: Rewrite the excess formula in terms of P and A.

excess(total) = (m_A + 1) - mean(total)
             = (m_A + 1) - [|A|*mean(A) + |B|*mean(B)] / (|A|+|B|)
             = (m_A + 1) - [|A|*mean(A) + |P|*(1+mean(P))] / (|A|+|P|)

Now |A| = prod |total[ci]|, |P| = prod |dp0[ci]|.

Define: w = |P|/(|A|+|P|), so w_A = 1-w, w_B = w.

excess = (m_A+1) - (1-w)*mean(A) - w*(1+mean(P))
       = (m_A+1) - mean(A) + w*mean(A) - w - w*mean(P)
       = (m_A+1-mean(A)) - w*(1 + mean(P) - mean(A))
       = (1 + excess(A)) - w*(1 + mean(P) - mean(A))

Now: mean(P) = sum_i mean(dp0[ci])
     mean(A) = sum_i mean(total[ci])
     mean(A) - mean(P) = sum_i [mean(total[ci]) - mean(dp0[ci])] = sum_i delta_i > 0

So: 1 + mean(P) - mean(A) = 1 - sum_i delta_i

And: excess = 1 + excess(A) - w*(1 - sum_i delta_i)

For excess < 1: w*(1 - sum_i delta_i) > excess(A)

APPROACH 1: Express w in terms of the DP quantities and find a lower bound.

w = prod |dp0[ci]| / (prod |total[ci]| + prod |dp0[ci]|)

For each child: |total[ci]| = |dp0[ci]| + |dp1[ci]|, so |total[ci]| > |dp0[ci]|.
ratio_i = |dp0[ci]| / |total[ci]| = 1 - |dp1[ci]|/|total[ci]| < 1.

w = prod ratio_i / (1 + prod ratio_i)

If there's one child: w = r/(1+r) where r = |dp0[c]|/|total[c]|.

APPROACH 2: Direct coefficient-level argument.

For the gap=+1 case, mode(total) = m_A+1 means (A+B)[m_A+1] is the largest coefficient.

Since (A+B)[m_A+1] >= (A+B)[m_A]:
This is exactly A[m_A+1]+P[m_A] >= A[m_A]+P[m_A-1]
i.e., P[m_A]-P[m_A-1] >= A[m_A]-A[m_A+1]

RHS >= 0 (since m_A is mode of A).
LHS >= 0 (since m_A < mode(P)).

KEY INEQUALITY: P[m_A] - P[m_A-1] >= A[m_A] - A[m_A+1]

Now I want to use this to bound the excess. Let me define:
Delta_A(j) = A[j] - A[j+1] (the "drop" in A going right)
Delta_P(j) = P[j] - P[j-1] (the "rise" in P going right)

At position m_A: Delta_P(m_A) >= Delta_A(m_A) >= 0.

This means: the rise in P at position m_A is at least as large as 
the drop in A at position m_A.

Can this constrain the excess?

Let me just compute everything carefully for the single-child case
and see if I can find a closed-form proof.
"""

import math

def polyadd(a, b):
    n = max(len(a), len(b))
    result = [0] * n
    for i in range(len(a)): result[i] += a[i]
    for i in range(len(b)): result[i] += b[i]
    return result

def polymul(a, b):
    if not a or not b: return [0]
    n = len(a) + len(b) - 1
    result = [0] * n
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            result[i + j] += ai * bj
    return result

def poly_mode(poly):
    return max(range(len(poly)), key=lambda k: poly[k])

def poly_mean(poly):
    total = sum(poly)
    return sum(k * poly[k] for k in range(len(poly))) / total

def poly_excess(poly):
    return poly_mode(poly) - poly_mean(poly)

def main():
    print("=" * 80)
    print("  SINGLE-CHILD ATTACK: Path recurrence")
    print("=" * 80)
    print()
    
    # For a path P_n, the DP is:
    # dp0[v_n] = total[v_{n-1}]
    # dp1[v_n] = x * dp0[v_{n-1}]
    # total[v_n] = dp0[v_n] + dp1[v_n] = total[v_{n-1}] + x * dp0[v_{n-1}]
    
    # Start: dp0[v_0] = [1], dp1[v_0] = [0,1], total[v_0] = [1,1]
    
    dp0 = [1]     # dp0 of current vertex's child
    dp1 = [0, 1]  # dp1 of current vertex's child
    total = [1, 1] # total of current vertex's child
    
    print(f"{'step':>4} {'gap':>3} {'excess':>8} {'excess_A':>8} {'w':>6} "
          f"{'1-Σδ':>6} {'correction':>10} {'slack':>8} {'formula':>8}")
    print("-" * 75)
    
    for step in range(2, 35):
        # At the new vertex:
        A = total[:]  # dp0 of new vertex = total of child
        B = [0] + dp0[:]  # dp1 of new vertex = x * dp0 of child
        new_total = polyadd(A, B)
        
        m_A = poly_mode(A)
        m_total = poly_mode(new_total)
        gap = m_total - m_A
        
        if gap == 1:
            excess_total = poly_excess(new_total)
            excess_A = poly_excess(A)
            
            S_A = sum(A)
            S_B = sum(B)
            w = S_B / (S_A + S_B)
            
            mu_B = poly_mean(B)
            mu_A = poly_mean(A)
            mu_diff = mu_B - mu_A
            
            correction = w * mu_diff
            slack = correction - excess_A
            formula = 1 + excess_A - correction
            
            print(f"{step:4d} {gap:+3d} {excess_total:8.4f} {excess_A:+8.4f} {w:6.4f} "
                  f"{mu_diff:6.4f} {correction:10.4f} {slack:8.4f} {formula:8.4f}")
            
            # PROOF ATTEMPT: For single child, express everything in terms of 
            # dp0 and dp1 of the child
            # A = f + g, B = x*f where f = dp0[child], g = dp1[child]
            
            f = dp0[:]
            g = dp1[:]
            
            S_f = sum(f)
            S_g = sum(g)
            mu_f = poly_mean(f) if sum(f) > 0 else 0
            mu_fg = poly_mean(polyadd(f, g))
            
            # w = S_f / (S_f + S_g + S_f) = S_f / (2*S_f + S_g)
            w_check = S_f / (2 * S_f + S_g)
            
            # mu_B = 1 + mu_f
            # mu_A = mu(f+g) = (S_f*mu_f + S_g*mu_g) / (S_f + S_g)
            # mu_diff = 1 + mu_f - mu(f+g) = 1 - [mu(f+g) - mu_f]
            #         = 1 - S_g*(mu_g - mu_f)/(S_f + S_g)
            
            # Let's define delta = mu(f+g) - mu_f:
            delta = mu_fg - mu_f
            mu_diff_check = 1 - delta
            
        elif gap == 0:
            excess_total = poly_excess(new_total)
            print(f"{step:4d} {gap:+3d} {excess_total:8.4f}  (gap=0, automatic)")
        else:
            excess_total = poly_excess(new_total)
            print(f"{step:4d} {gap:+3d} {excess_total:8.4f}")
        
        # Update for next step
        dp0 = total[:]  # dp0[new] = total[child]
        dp1 = B[:]       # dp1[new] = x * dp0[child]
        total = new_total[:]
    
    print()
    
    # Now try to find the pattern in the single-child case
    # For a path, everything is determined by the Fibonacci-like recurrence
    print("PATH EXCESS: Does excess(P_n) oscillate and stay bounded?")
    print()
    
    dp0 = [1]
    total = [1, 1]
    
    excesses = []
    for step in range(2, 201):
        A = total[:]
        B = [0] + dp0[:]
        new_total = polyadd(A, B)
        
        ex = poly_excess(new_total)
        excesses.append((step, ex))
        
        dp0 = total[:]
        total = new_total[:]
    
    max_ex = max(ex for _, ex in excesses)
    min_ex = min(ex for _, ex in excesses)
    
    print(f"  Path excess range: [{min_ex:.6f}, {max_ex:.6f}]")
    print(f"  Max excess < 1: {max_ex < 1}")
    print()
    
    # Show the oscillation
    print("  Excess every 10 steps:")
    for step, ex in excesses:
        if step % 10 == 0:
            print(f"    P_{step}: excess = {ex:+.6f}")
    
    print()
    
    # THEORETICAL OBSERVATION for paths:
    # The IS poly of P_n satisfies I(P_n) = I(P_{n-1}) + x*I(P_{n-2})
    # If I(P_{n-1}) has excess e_{n-1} and we're in a gap=+1 merge:
    #   excess(P_n) = 1 + e_{n-1} - w*(1-delta)
    # 
    # Since the path recurrence gives w and delta in terms of Fibonacci-like quantities,
    # we should be able to show this is bounded.
    
    print("=" * 80)

if __name__ == "__main__":
    main()
