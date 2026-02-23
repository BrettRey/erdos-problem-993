#!/usr/bin/env python3
"""Verify TP_2 closure for the cross-Turan determinant under convolution.

Question: if D(f1, g1, j) >= 0 and D(f2, g2, j) >= 0 for all j,
where D(f, g, j) = f_j * g_{j-1} - f_{j+1} * g_{j-2},
does D(f1*f2, g1*g2, j) >= 0 for all j?

Equivalently: if the 2-row matrices
  [f1_0 f1_1 ...]    and    [f2_0 f2_1 ...]
  [g1_0 g1_1 ...]           [g2_0 g2_1 ...]
(with g shifted by 1) are each TP_2, is the convolution-product
  [F_0  F_1  ...]
  [G_0  G_1  ...]
also TP_2?

This would follow from a variant of the Cauchy-Binet formula for Toeplitz matrices,
but let me verify it computationally first.

NOTE: The standard Cauchy-Binet for TP_2 applies to matrix products. Convolution
IS matrix multiplication for Toeplitz matrices. So:

If A = Toeplitz(f1) and B = Toeplitz(f2) are TP_2, then A*B = Toeplitz(f1*f2) is TP_2.
Similarly for g.

But the cross-Turan involves TWO different sequences (f and g). The question is
whether the *pair* (f, g) being cross-TP_2 is preserved under convolution.

This is equivalent to: the matrix
  [Toeplitz(f)]
  [Toeplitz(g) shifted by 1 row]
being TP_2.

Actually, for our specific application, we have:
  f_c = dp0[c] + dp1[c]
  g_c = dp0[c]

So f_c = g_c + dp1[c]. The cross-Turan D(f_c, g_c, j) = D(g_c + dp1[c], g_c, j)
= D(g_c, g_c, j) + D(dp1[c], g_c, j)
= 0 + D(dp1[c], g_c, j)  [since D(g,g,j) = g_j*g_{j-1} - g_{j+1}*g_{j-2} = LC surplus]

Wait no: D(g,g,j) = g_j * g_{j-1} - g_{j+1} * g_{j-2}, which is a cross-index LC-like
quantity, NOT the standard LC surplus g_j^2 - g_{j+1}*g_{j-1}. These are different!

D(g,g,j) = g_j * g_{j-1} - g_{j+1} * g_{j-2}: this is >= 0 iff g_j/g_{j+1} >= g_{j-2}/g_{j-1},
i.e., the lambda-ratios are non-increasing. For a LC sequence, the lambda-ratios
lambda_k = g_{k-1}/g_k form a non-decreasing sequence (by LC: g_k^2 >= g_{k-1}*g_{k+1}
means g_{k-1}/g_k <= g_k/g_{k+1}, so lambda_k <= lambda_{k+1}).

So g_j/g_{j+1} <= g_{j+1}/g_{j+2} (ratios non-decreasing).
D(g,g,j) >= 0 iff g_j/g_{j+1} >= g_{j-2}/g_{j-1}, i.e., lambda_{j+1} >= lambda_{j-1}.
Since lambdas are non-decreasing (by LC), lambda_{j+1} >= lambda_{j-1}. YES!
So D(g,g,j) >= 0 for any LC sequence g.

Actually wait, let me recheck. g_j*g_{j-1} - g_{j+1}*g_{j-2}:
= g_{j-1}*(g_j - g_{j+1}*g_{j-2}/g_{j-1})
= g_{j-1}*g_j*(1 - (g_{j+1}/g_j)*(g_{j-2}/g_{j-1}))
= g_{j-1}*g_j*(1 - (lambda_{j+1})^{-1} * lambda_{j-1})    [where lambda_k = g_{k-1}/g_k]

For LC: lambda_{k+1} >= lambda_k >= ... >= lambda_{j-1}. So lambda_{j+1} >= lambda_{j-1}.
Thus (lambda_{j+1})^{-1} <= (lambda_{j-1})^{-1}, so
(lambda_{j+1})^{-1} * lambda_{j-1} <= 1.
So D(g,g,j) >= 0. Confirmed.

OK so even D(g,g,j) >= 0 for LC sequences. And the full cross-Turan:
D(f,g,j) = D(g+dp1, g, j) = D(g,g,j) + D(dp1, g, j)

Both parts are non-negative (D(g,g,j) >= 0 by LC of g; D(dp1,g,j) >= 0 by the
inductive cross-Turan). So D(f,g,j) >= 0.

This is already a proof that the SINGLE-FACTOR cross-Turan holds!
And since f_c = g_c + dp1[c], we have:
D(f_c, g_c, j) = D(g_c, g_c, j) + D(dp1[c], g_c, j)

where D(g_c, g_c, j) >= 0 (by LC of g_c) and D(dp1[c], g_c, j) is the recursive
cross-Turan one level down.

For the PRODUCT cross-Turan, we need D(prod f_c, prod g_c, j) >= 0.

Approach: Since each f_c = g_c + dp1[c] with dp1[c] having non-negative coefficients,
and each g_c is LC, perhaps there's a direct induction on the number of factors.

Base case (1 factor): D(f_1, g_1, j) >= 0 (proved above).

Step: Given D(F, G, j) >= 0 where F = prod_{i<k} f_i, G = prod_{i<k} g_i,
show D(F*f_k, G*g_k, j) >= 0.

(F*f_k)_j = sum_r F_r * f_k(j-r)
(G*g_k)_j = sum_r G_r * g_k(j-r)

D(F*f_k, G*g_k, j)
= sum_{r,s} [F_r * f_k(j-r) * G_s * g_k(j-1-s) - F_r * f_k(j+1-r) * G_s * g_k(j-2-s)]
= sum_{r,s} F_r * G_s * [f_k(j-r) * g_k(j-1-s) - f_k(j+1-r) * g_k(j-2-s)]

For fixed (r, s) with r <= s:
  The bracketed term is f_k(a)*g_k(b) - f_k(a+1)*g_k(b-1) where a=j-r, b=j-1-s.
  This is D(f_k, g_k, a) when b = a - 1 + (r-s) = a - 1 - (s-r).

Hmm, this isn't quite D(f_k, g_k, a) unless s = r.

When s = r: f_k(j-r)*g_k(j-1-r) - f_k(j+1-r)*g_k(j-2-r) = D(f_k, g_k, j-r) >= 0.

When s != r: the cross terms don't factor as simple D(f_k, g_k, ...) terms.

So the Cauchy-Binet approach doesn't directly apply to the "shifted-pair" TP_2 condition.

Let me think about this differently.

Actually, let me just test it computationally for small cases. Take two pairs
(f_1, g_1) and (f_2, g_2) from actual trees and check.
"""

from __future__ import annotations

import subprocess
from indpoly import _polyadd, _polymul


def D_cross(f, g, j):
    """D(f, g, j) = f_j * g_{j-1} - f_{j+1} * g_{j-2}."""
    fj = f[j] if 0 <= j < len(f) else 0
    fj1 = f[j+1] if 0 <= j+1 < len(f) else 0
    gj1 = g[j-1] if 0 <= j-1 < len(g) else 0
    gj2 = g[j-2] if 0 <= j-2 < len(g) else 0
    return fj * gj1 - fj1 * gj2


def check_pair_tp2(f, g, label=""):
    """Check D(f, g, j) >= 0 for all j."""
    max_j = max(len(f), len(g)) + 2
    for j in range(2, max_j):
        d = D_cross(f, g, j)
        if d < 0:
            print(f"  FAIL {label}: j={j}, D={d}, f={f}, g={g}")
            return False
    return True


def main():
    # Test with specific tree pairs
    print("=== Testing TP_2 closure under convolution ===")
    print()

    # Pair 1: single edge (1 vertex subtree)
    # dp0 = [1], dp1 = [0,1], I = [1,1]
    f1 = [1, 1]
    g1 = [1]

    # Pair 2: path of 2 (2-vertex subtree)
    # dp0 = [1,1], dp1 = [0,1], I = [1,2,1]
    # Wait: for a path c-w, rooted at c:
    # dp0[w] = [1], dp1[w] = [0,1]
    # dp0[c] = dp0[w]+dp1[w] = [1,1]
    # dp1[c] = x*dp0[w] = [0,1]
    # I(T_c) = dp0[c]+dp1[c] = [1,2]? No:
    # Wait, the subtree rooted at c with child w (edge c-w):
    # dp0[c] = product over children of (dp0+dp1) = dp0[w]+dp1[w] = [1]+[0,1] = [1,1]
    # dp1[c] = x * product over children of dp0 = x*[1] = [0,1]
    # I(subtree_c) = dp0[c] + dp1[c] = [1,1]+[0,1] = [1,2]
    # This is correct: 2 vertices with 1 edge, IS: {}, {c}, {w} = 1+2x.
    f2 = [1, 2]
    g2 = [1, 1]

    # Pair 3: star K_{1,3} subtree (3-leaf star rooted at center)
    # Children are 3 leaves, each has dp0=[1], dp1=[0,1].
    # dp0[center] = product of ([1]+[0,1])^3 = [1,1]^3 = [1,3,3,1]
    # dp1[center] = x * [1]^3 = [0,1]
    # I(star) = [1,3,3,1] + [0,1] = [1,4,3,1]
    f3 = [1, 4, 3, 1]
    g3 = [1, 3, 3, 1]

    # Pair 4: path of 3 (3-vertex path rooted at one end)
    # c - w1 - w2
    # dp0[w2]=[1], dp1[w2]=[0,1]
    # dp0[w1] = [1]+[0,1] = [1,1]; dp1[w1] = x*[1] = [0,1]
    # dp0[c] = [1,1]+[0,1] = [1,2]; dp1[c] = x*[1,1] = [0,1,1]
    # I(subtree_c) = [1,2]+[0,1,1] = [1,3,1]
    f4 = [1, 3, 1]
    g4 = [1, 2]

    pairs = [(f1, g1, "edge"), (f2, g2, "P2"), (f3, g3, "K13"), (f4, g4, "P3")]

    # Check each pair
    for f, g, label in pairs:
        ok = check_pair_tp2(f, g, label)
        print(f"Single pair ({label}): {'OK' if ok else 'FAIL'}")

    print()

    # Check all pairwise products
    for i, (f1, g1, l1) in enumerate(pairs):
        for j, (f2, g2, l2) in enumerate(pairs):
            F = _polymul(f1, f2)
            G = _polymul(g1, g2)
            label = f"{l1} x {l2}"
            ok = check_pair_tp2(F, G, label)
            if not ok:
                print(f"  Product ({label}): FAIL! F={F}, G={G}")
            else:
                print(f"  Product ({label}): OK")

    print()

    # Check triple products
    for i, (f1, g1, l1) in enumerate(pairs):
        for j, (f2, g2, l2) in enumerate(pairs):
            for k, (f3, g3, l3) in enumerate(pairs):
                F = _polymul(_polymul(f1, f2), f3)
                G = _polymul(_polymul(g1, g2), g3)
                label = f"{l1} x {l2} x {l3}"
                ok = check_pair_tp2(F, G, label)
                if not ok:
                    print(f"  Triple product ({label}): FAIL! F={F}, G={G}")

    print("Triple products: all done (only failures shown above)")

    # Now check: is cross = 2*p1*q1 - pm*q0 - p0*qm >= 0 when (P, P') satisfies TP_2?
    # And is R = cross + mismatch >= 0?
    # These are DIFFERENT from the TP_2 condition on (P, P').
    # The TP_2 condition gives D(P, P', j) >= 0, which is Term A.
    # R also includes Terms B and C.

    # Let me check: for the product pairs, what are D, Term A, B, C, R?
    print("\n=== Detailed decomposition for products ===")
    for i, (f1, g1, l1) in enumerate(pairs):
        for j, (f2, g2, l2) in enumerate(pairs):
            P = _polymul(f1, f2)  # = prod I(T_c)
            Pprime = _polymul(g1, g2)  # = prod dp0[c]
            # Q = x * P', so q_k = P'_{k-1}
            Q = [0] + Pprime

            for m in range(2, len(P)):
                p0 = P[m-2] if m-2 < len(P) else 0
                p1 = P[m-1] if m-1 < len(P) else 0
                pm = P[m] if m < len(P) else 0
                q0 = Q[m-2] if m-2 < len(Q) else 0
                q1 = Q[m-1] if m-1 < len(Q) else 0
                qm = Q[m] if m < len(Q) else 0

                termA = p1*q1 - pm*q0
                termB = p1*(q1 - q0)
                termC = p0*(q1 - qm)
                R = termA + termB + termC

                b0 = p0 + q0
                b1 = p1 + q1
                b2 = pm + qm
                lc_P = p1*p1 - pm*p0
                lc_Q = q1*q1 - qm*q0
                combined = lc_P + lc_Q + R

                if R < 0:
                    print(f"  R < 0! {l1}x{l2}, m={m}: R={R}, A={termA}, B={termB}, C={termC}")
                    print(f"    P={P}, Q={Q}")
                    print(f"    p0={p0} p1={p1} pm={pm} q0={q0} q1={q1} qm={qm}")

    print("\nDone. (Only R < 0 cases shown above)")


if __name__ == "__main__":
    main()
