#!/usr/bin/env python3
"""KT4 (DIAGNOSTIC, float64): how Gaussian is the smoothing-shift variable?

Condition on J = I cap R.  Then X = |J| + Bin(B_J, p), p = lam/(1+lam), and the
exact identity
    2 p_k - p_{k-1} - p_{k+1} = E_J[ -Delta^2 b_{B_J,p}(k - |J|) ]
reduces local curvature to the law of the shift M = |J| + p B_J (and of B_J).
The Gaussian-comparison sufficient condition is
    d_K(M, Normal)  <  c(rho) := (1/sqrt(2 pi)) / TV(phi''') * rho^{3/2}
                     = 0.2640... * rho^{3/2},   TV(phi''') = int |phi'''| = 1.5107...
(ignoring O(1/v) binomial-vs-normal and B-fluctuation corrections).
This script computes the joint law of (|J|, B_J) under the R-marginal
(weights lam^|J| (1+lam)^B_J) by a tree DP with 2-D arrays (float64,
nonnegative entries, per-step renormalised), then d_K(M, N(EM, Var M)),
the relative spread sd(B)/E(B), and the ratio  r = d_K / c(rho).
The side L (binomial side) is chosen as the one with larger E|I cap side|.
lam is the symmetric fugacity sqrt(i_{k-1}/i_{k+1}) (float).
"""
import sys, json, math
import numpy as np
from scipy.signal import fftconvolve
from scipy.stats import norm
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from lclt_lib import tree_polys, window

def conv(a, b):
    if a.size * b.size < 4096:
        out = np.zeros((a.shape[0] + b.shape[0] - 1, a.shape[1] + b.shape[1] - 1))
        for i in range(a.shape[0]):
            for j in range(a.shape[1]):
                if a[i, j] != 0.0:
                    out[i:i + b.shape[0], j:j + b.shape[1]] += a[i, j] * b
        return out
    out = fftconvolve(a, b)
    out[out < 0] = 0.0
    return out

def norm_(a):
    s = a.sum()
    return a / s, s

def joint_JB(adj, Lside_is_even, lam):
    """Return 2-D array P[j, B] (normalised) for J = I cap R, B = #available L."""
    n = len(adj)
    parent = [-1] * n; depth = [0] * n; order = [0]; seen = [False] * n; seen[0] = True
    i = 0
    while i < len(order):
        u = order[i]; i += 1
        for w in adj[u]:
            if not seen[w]:
                seen[w] = True; parent[w] = u; depth[w] = depth[u] + 1; order.append(w)
    isL = [((d % 2 == 0) == Lside_is_even) for d in depth]
    one = np.ones((1, 1))
    # For R vertex: (Gout, Gin) arrays; for L vertex: (F0, F1) = parent not in J / parent in J
    A = [None] * n; Bv = [None] * n
    y = 1.0 + lam
    for u in reversed(order):
        ch = [c for c in adj[u] if c != parent[u]]
        if isL[u]:
            allp = one; outp = one
            for c in ch:
                Gout, Gin = A[c], Bv[c]
                tot = np.zeros((max(Gout.shape[0], Gin.shape[0]), max(Gout.shape[1], Gin.shape[1])))
                tot[:Gout.shape[0], :Gout.shape[1]] += Gout; tot[:Gin.shape[0], :Gin.shape[1]] += Gin
                allp = conv(allp, tot); outp = conv(outp, Gout)
                s = max(allp.max(), 1e-300); allp /= s; outp /= s
            # F1: parent in J -> u unavailable: all children free
            F1 = allp.copy()
            # F0: available iff no child in J: outp shifted in B by one with weight y; plus (allp - outp)
            F0 = np.zeros((allp.shape[0], allp.shape[1] + 1))
            F0[:allp.shape[0], :allp.shape[1]] += allp
            F0[:outp.shape[0], :outp.shape[1]] -= outp
            F0[F0 < 0] = 0.0
            F0[:outp.shape[0], 1:outp.shape[1] + 1] += y * outp
            s = max(F0.max(), F1.max()); A[u], Bv[u] = F0 / s, F1 / s
        else:
            p0 = one; p1 = one
            for c in ch:
                p0 = conv(p0, A[c]); p1 = conv(p1, Bv[c])
                s = max(p0.max(), p1.max()); p0 /= s; p1 /= s
            Gin = np.zeros((p1.shape[0] + 1, p1.shape[1])); Gin[1:, :] = lam * p1
            s = max(p0.max(), Gin.max()); A[u], Bv[u] = p0 / s, Gin / s
    r = order[0]
    if isL[r]:
        P = A[r]
    else:
        Gout, Gin = A[r], Bv[r]
        P = np.zeros((max(Gout.shape[0], Gin.shape[0]), max(Gout.shape[1], Gin.shape[1])))
        P[:Gout.shape[0], :Gout.shape[1]] += Gout; P[:Gin.shape[0], :Gin.shape[1]] += Gin
    return P / P.sum()

def dK_shift(P, p):
    js, bs = np.nonzero(P > 0)
    w = P[js, bs]
    m = js + p * bs
    EM = (w * m).sum(); VM = (w * (m - EM) ** 2).sum()
    EB = (w * bs).sum(); VB = (w * (bs - EB) ** 2).sum()
    order = np.argsort(m, kind="stable")
    ms = m[order]; ws = w[order]
    # merge equal atoms
    uniq, idx = np.unique(np.round(ms, 12), return_index=True)
    wsum = np.add.reduceat(ws, idx)
    F = np.cumsum(wsum); Fm = F - wsum
    if VM <= 0:
        return 0.0, EM, VM, EB, VB
    z = norm.cdf((uniq - EM) / math.sqrt(VM))
    d = max(np.max(np.abs(F - z)), np.max(np.abs(Fm - z)))
    return float(d), EM, VM, EB, VB

TVphi3 = 2 * (2 * norm.pdf(math.sqrt(3)) + norm.pdf(0) + 2 * norm.pdf(math.sqrt(3)))
C0 = (1 / math.sqrt(2 * math.pi)) / TVphi3

def analyse_W(adj, ks=None, npts=5):
    n = len(adj)
    Z, ZL, col = tree_polys(adj)
    alpha = len(Z) - 1
    lo, hi, q = window(n, alpha)
    if hi < lo: return None
    if ks is None:
        ks = sorted(set(lo + round(i * (hi - lo) / (npts - 1)) for i in range(npts)))
    rows = []
    for k in ks:
        lam = math.sqrt(Z[k - 1] / Z[k + 1]) if Z[k+1] < 10**300 else math.exp(0.5*(math.log(Z[k-1]) - math.log(Z[k+1])))
        p = lam / (1 + lam)
        best = None
        for even in (True, False):
            P = joint_JB(adj, even, lam)
            d, EM, VM, EB, VB = dK_shift(P, p)
            EV = p * (1 - p) * EB
            sig2 = VM + EV    # total variance = Var E[X|J] + E Var[X|J]
            rho = EV / sig2
            rec = dict(k=k, lam=lam, side_even=even, dK=d, rho=rho, c_rho=C0 * rho ** 1.5,
                       ratio=d / (C0 * rho ** 1.5), sdB_over_EB=math.sqrt(VB) / EB if EB > 0 else None,
                       sigma2=sig2)
            if best is None or rec["ratio"] < best["ratio"]:
                best = rec
        rows.append(best)
    return dict(n=n, alpha=alpha, rows=rows, max_ratio=max(r["ratio"] for r in rows),
                max_dK=max(r["dK"] for r in rows), min_rho=min(r["rho"] for r in rows))

if __name__ == "__main__":
    print(json.dumps(dict(TV_phi3=TVphi3, C0=C0)))


# ---------------------------------------------------------------------------
# KT5 (DIAGNOSTIC, float64): the fixed-kernel smoothing certificate.
# On E0 = {B_J >= B0}: X = M' + Y, Y ~ Bin(B0, p) independent of M',
# M' = |J| + Bin(B_J - B0, p).  With h(j) = -Delta^2 b_{B0,p}(j),
#   D_k = E[h(k - M'); E0] + (-Delta^2 P(X = ., E0^c))(k),
#   |second| <= 2 P(E0^c).
# Gaussian comparison (Abel summation):
#   E[h(k - M') | E0] >= E h(k - Gd) - TV(h) * d_K(M' | E0, Gd),
# Gd the discretised normal with the conditional mean/variance of M'.
# Cert := P(E0) [E h(k-Gd) - TV(h) dK] - 2 P(E0^c).  Cert > 0 certifies
# log-concavity at k by this mechanism alone (given the true d_K).
# Also reports Exact := D_k computed from the true law (sanity: must be > 0).
# ---------------------------------------------------------------------------
from scipy.stats import binom as _binom

def cert_one(P, p, k, B0):
    J, Bn = P.shape
    Bs = np.arange(Bn)
    PE0 = P[:, B0:].sum() if B0 < Bn else 0.0
    if PE0 <= 0:
        return None
    maxm = J + Bn
    law = np.zeros(maxm + 1)
    for B in range(B0, Bn):
        col = P[:, B]
        if col.sum() == 0: continue
        kern = _binom.pmf(np.arange(B - B0 + 1), B - B0, p)
        law[:J + B - B0] += np.convolve(col, kern)
    law /= law.sum()
    ms = np.arange(law.size)
    mbar = (law * ms).sum(); s2 = (law * (ms - mbar) ** 2).sum(); s = math.sqrt(max(s2, 1e-300))
    Gd = norm.cdf((ms + 0.5 - mbar) / s) - norm.cdf((ms - 0.5 - mbar) / s)
    Gd[0] = norm.cdf((0.5 - mbar) / s); Gd[-1] = 1 - norm.cdf((ms[-1] - 0.5 - mbar) / s)
    dK = np.max(np.abs(np.cumsum(law) - np.cumsum(Gd)))
    b = _binom.pmf(np.arange(B0 + 1), B0, p)
    bb = np.concatenate([np.zeros(3), b, np.zeros(3)])
    d2 = bb[2:] - 2 * bb[1:-1] + bb[:-2]          # Delta^2 b at j = -2..B0+2 (index shift 2)
    h = -d2                                        # h(j) for j = -2 .. B0+2
    TV = np.abs(np.diff(np.concatenate([[0.0], h, [0.0]]))).sum()
    def Eh(lawvec):
        tot = 0.0
        for m, w in enumerate(lawvec):
            if w == 0: continue
            j = k - m
            if -2 <= j <= B0 + 2:
                tot += w * h[j + 2]
        return tot
    EhG = Eh(Gd); EhM = Eh(law)
    cert = PE0 * (EhG - TV * dK) - 2 * (1 - PE0)
    return dict(B0=B0, PE0=PE0, dK_Mprime=float(dK), TV_h=float(TV), EhG=float(EhG), EhM=float(EhM),
                cert=float(cert), cert_rel=float(cert / EhG) if EhG > 0 else None, s2_Mprime=float(s2))

def certificate(adj, npts=5, etas=(0.1, 0.2, 0.35, 0.5)):
    n = len(adj)
    Z, ZL, col = tree_polys(adj)
    alpha = len(Z) - 1
    lo, hi, q = window(n, alpha)
    if hi < lo: return None
    ks = sorted(set(lo + round(i * (hi - lo) / (npts - 1)) for i in range(npts)))
    rows = []
    for k in ks:
        lam = math.exp(0.5 * (math.log(Z[k - 1]) - math.log(Z[k + 1])))
        p = lam / (1 + lam)
        # exact curvature of p_k at this lam (float of exact-integer ratio): D_k / p_k
        best = None
        for even in (True, False):
            P = joint_JB(adj, even, lam)
            EB = (P.sum(axis=0) * np.arange(P.shape[1])).sum()
            for eta in etas:
                B0 = max(1, int(math.floor((1 - eta) * EB)))
                c = cert_one(P, p, k, B0)
                if c is None: continue
                c.update(k=k, eta=eta, side_even=even, lam=lam, EB=float(EB))
                if best is None or c["cert_rel"] is not None and (best["cert_rel"] is None or c["cert_rel"] > best["cert_rel"]):
                    best = c
        rows.append(best)
    return dict(n=n, alpha=alpha, window=[lo, hi, q], rows=rows,
                min_cert_rel=min(r["cert_rel"] for r in rows if r and r["cert_rel"] is not None),
                all_cert=all(r and r["cert"] > 0 for r in rows))
