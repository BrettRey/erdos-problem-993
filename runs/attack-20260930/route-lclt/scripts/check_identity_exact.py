#!/usr/bin/env python3
"""Exact check (Fractions) of the smoothing identities on all trees n <= NMAX:
 (I1) P(X=j) = sum_{J subset R} w_J b_{B_J,p}(j - |J|),  w_J = lam^|J|(1+lam)^{B_J}/Z
 (I2) 2p_k - p_{k-1} - p_{k+1} = sum_J w_J * (-Delta^2 b_{B_J,p})(k - |J|)
 (I3) fixed kernel: for B0 >= 1, sum over J with B_J >= B0 of w_J (b_{B_J} (j-|J|))
      = sum_m P(M'=m, E0) b_{B0}(j-m),  M' = |J| + Bin(B_J - B0, p)
at lam = i_{k-1}/i_k for every window k, both bipartition sides."""
import sys, subprocess
from fractions import Fraction
from math import comb
from lclt_lib import parent_line_to_adj, tree_polys, window

def binpmf(B, p, j):
    if j < 0 or j > B: return Fraction(0)
    return comb(B, j) * p**j * (1 - p)**(B - j)

def check(adj):
    n = len(adj)
    Z, ZL, col = tree_polys(adj)
    alpha = len(Z) - 1
    lo, hi, q = window(n, alpha)
    cnt = 0
    for k in range(lo, hi + 1):
        lam = Fraction(Z[k - 1], Z[k]); p = lam / (1 + lam)
        Zl = sum(c * lam**j for j, c in enumerate(Z))
        pk = [Z[j] * lam**j / Zl for j in range(alpha + 1)] + [Fraction(0)] * 3
        for side in (0, 1):
            R = [v for v in range(n) if col[v] != side]; L = [v for v in range(n) if col[v] == side]
            terms = []
            for mask in range(1 << len(R)):
                Js = {R[i] for i in range(len(R)) if mask >> i & 1}
                B = sum(1 for v in L if not any(w in Js for w in adj[v]))
                terms.append((len(Js), B, lam**len(Js) * (1 + lam)**B))
            W = sum(t[2] for t in terms)
            assert W == Zl
            for j in range(alpha + 1):
                assert pk[j] == sum(w * binpmf(B, p, j - s) for s, B, w in terms) / W
            lhs = 2 * pk[k] - pk[k - 1] - pk[k + 1]
            rhs = sum(w * (2 * binpmf(B, p, k - s) - binpmf(B, p, k - 1 - s) - binpmf(B, p, k + 1 - s)) for s, B, w in terms) / W
            assert lhs == rhs
            B0 = max(1, min(t[1] for t in terms) + 1)
            for j in range(alpha + 1):
                a = sum(w * binpmf(B, p, j - s) for s, B, w in terms if B >= B0)
                bsum = Fraction(0)
                for s, B, w in terms:
                    if B < B0: continue
                    for m in range(s, s + B - B0 + 1):
                        bsum += w * binpmf(B - B0, p, m - s) * binpmf(B0, p, j - m)
                assert a == bsum
            cnt += 1
    return cnt

if __name__ == "__main__":
    nmax = int(sys.argv[1])
    tot = 0; trees = 0
    for n in range(3, nmax + 1):
        out = subprocess.run(["gentreeg", "-p", "-q", str(n)], capture_output=True, text=True).stdout
        for line in out.splitlines():
            par = list(map(int, line.split()))
            if len(par) != n: continue
            tot += check(parent_line_to_adj(par)); trees += 1
    print(f"identities I1-I3 hold exactly: {trees} trees, {tot} (tree,k,side) cases, n<=%d" % nmax)
