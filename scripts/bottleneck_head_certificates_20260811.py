#!/usr/bin/env python3
"""Bottleneck family (1,1) x (2h+1,2h+1): symbolic head-window analysis.

Reflected coefficients of B and C at fixed depth are polynomials in h.
This script derives them exactly from binomial formulas, validates them
against exact numeric instances, and computes the Turan and
partial-synchronization margins at head depths as polynomials in h.

Conventions:
  n = 2h+1.  u = P_n, v = P_{n-1} (path independence polynomials).
  Q_a=(1+x)^2, R_a=1, Q_b=u^2, R_b=v^2.
  B = (1+x)^2 u^2 + x u^2 + x (1+x)^2 v^2      (deg N = 2h+4)
  C = (1+x)^2 u^2 + x v^2                       (deg N)
Reflected coefficient at depth d means [x^{N-d}].

Exact facts used (verified numerically below as well):
  [x^{h+1-j}] P_{2h+1} = binom(h+1+j, 2j)        (depth j in u, deg u = h+1)
  [x^{h-j}]   P_{2h}   = binom(h+1+j, 2j+1)      (depth j in v, deg v = h)
"""
import sympy as sp

h = sp.symbols('h', nonnegative=True)
DMAX = 8          # symbolic head depths 0..DMAX
HCHECK = range(0, 31)

# reflected coefficient generators, as polynomials in h
def ru(j):        # depth j of u = P_{2h+1}
    return sp.expand_func(sp.binomial(h + 1 + j, 2 * j))

def rv(j):        # depth j of v = P_{2h}
    return sp.expand_func(sp.binomial(h + 1 + j, 2 * j + 1))

def refl_square(r, d):
    """depth-d reflected coefficient of S^2 given S's reflected gen r."""
    return sp.expand(sum(r(i) * r(d - i) for i in range(d + 1)))

def refl_shift121(seq_fn, d):
    """reflected of (1+x)^2 * S at depth d: S_d + 2 S_{d-1} + S_{d-2}."""
    out = seq_fn(d)
    if d >= 1:
        out += 2 * seq_fn(d - 1)
    if d >= 2:
        out += seq_fn(d - 2)
    return sp.expand(out)

u2 = lambda d: refl_square(ru, d) if d >= 0 else 0
v2 = lambda d: refl_square(rv, d) if d >= 0 else 0

# component degrees: (1+x)^2 u^2 -> N;  x u^2 -> N-1;  x (1+x)^2 v^2 -> N-1;
# x v^2 -> N-3.   (N = 2h+4, deg u^2 = 2h+2, deg v^2 = 2h.)
def Bd(d):
    out = refl_shift121(u2, d)                      # (1+x)^2 u^2, offset 0
    if d >= 1:
        out += u2(d - 1)                            # x u^2, offset 1
        out += refl_shift121(v2, d - 1)             # x(1+x)^2 v^2, offset 1
    return sp.expand(out)

def Cd(d):
    out = refl_shift121(u2, d)
    if d >= 3:
        out += v2(d - 3)                            # x v^2, offset 3
    return sp.expand(out)

Bpoly = [sp.Poly(Bd(d), h) for d in range(DMAX + 1)]
Cpoly = [sp.Poly(Cd(d), h) for d in range(DMAX + 1)]

# ---- numeric validation against exact polynomial arithmetic ----
def pmul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    return out

def padd(a, b):
    m = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
            for i in range(m)]

def P(nn):
    a, b = [1], [1, 1]      # P_0, P_1
    if nn == 0:
        return a
    for _ in range(nn - 1):
        a, b = b, padd(b, [0] + a)
    return b

def BC_num(hh):
    n = 2 * hh + 1
    u, v = P(n), P(n - 1)
    one2 = [1, 2, 1]
    u2n, v2n = pmul(u, u), pmul(v, v)
    Bn = padd(padd(pmul(one2, u2n), [0] + u2n), [0] + pmul(one2, v2n))
    Cn = padd(pmul(one2, u2n), [0] + v2n)
    return Bn, Cn

ok = True
for hh in HCHECK:
    Bn, Cn = BC_num(hh)
    N = len(Bn) - 1
    assert N == 2 * hh + 4
    for d in range(min(DMAX, N) + 1):
        sb = Bpoly[d].eval(hh)
        sc = Cpoly[d].eval(hh)
        if Bn[N - d] != sb or Cn[N - d] != sc:
            ok = False
            print(f"MISMATCH h={hh} d={d}: num B {Bn[N-d]} sym {sb}; "
                  f"num C {Cn[N-d]} sym {sc}")
print("symbolic reflected coefficients validated on h=0..30:", ok)

print("\nReflected coefficients as polynomials in h:")
for d in range(DMAX + 1):
    print(f"  B_{{N-{d}}} = {Bpoly[d].as_expr()}")
for d in range(DMAX + 1):
    print(f"  C_{{N-{d}}} = {Cpoly[d].as_expr()}")

# ---- head-window margins as polynomials in h ----
def at(polys, d):
    return polys[d].as_expr() if 0 <= d <= DMAX else sp.Integer(0)

print("\nTuran gaps of B at reflected depth d (positive <=> LC there):")
for d in range(1, DMAX):
    # reflected Turan at depth d: B_{N-d}^2 - B_{N-d+1} B_{N-d-1}
    g = sp.expand(at(Bpoly, d) ** 2 - at(Bpoly, d - 1) * at(Bpoly, d + 1))
    gp = sp.Poly(g, h)
    roots = sp.real_roots(gp) if gp.degree() > 0 else []
    print(f"  d={d}: {sp.factor(g)}")
    print(f"        real roots: {[sp.nsimplify(r, rational=False) for r in roots]}")

# partial sync xC ~_p B and C ~_p B at head pairs.
# Work in reflected depths: absolute index m = N - dm etc.
# For polynomials A (=xC or C, deg N+1 or N) and B (deg N), the sync
# inequality at (m, n), m >= n:
#   A_m B_n + A_n B_m >= A_{m+1} B_{n-1} + A_{n-1} B_{m+1}
# xC has reflected depth d meaning [x^{N+1-d}] xC = C_{N-d+... }:
# [x^k] xC = C_{k-1}, so reflected depth d of xC (deg N+1) = C at depth d.
def sync_margin(Adepth, Bdepth, dm, dn, degA_minus_degB):
    """margin at (m,n) with m = degB - dm + degA_minus_degB ... use depths.
    We parametrize by reflected depths (dm, dn) OF EACH polynomial in its
    own reflected indexing, requiring the absolute indices to satisfy the
    sync pattern. For A of degree N+e and B of degree N, absolute m for A
    means depth dA = N+e-m; absolute n for B depth dB = N-n.
    Inequality terms: A_m B_n + A_n B_m - A_{m+1} B_{n-1} - A_{n-1} B_{m+1}.
    """
    e = degA_minus_degB
    Am = lambda m_abs_depth: Adepth(m_abs_depth)         # depth in A
    Bm = lambda d: Bdepth(d)
    # absolute: m has A-depth dm, B-depth dm - e... caller passes consistent.
    return None  # placeholder; explicit cases below

# Explicit: xC vs B. deg xC = N+1, deg B = N. Absolute index m ~ A-depth
# a(m) = N+1-m, B-depth b(m) = N-m = a(m)-1.
XC = lambda d: at(Cpoly, d)      # depth d of xC equals depth d of C? check:
# [x^{N+1-d}] xC = [x^{N-d}] C  -> yes, depth d of xC = depth d of C.
Br = lambda d: at(Bpoly, d)
Cr = lambda d: at(Cpoly, d)

def sync_margin_xCB(dm, dn):
    """Sync inequality for A=xC (deg N+1), B (deg N) at absolute (m,n) with
    m at A-depth dm i.e. m = N+1-dm, n at B-depth dn i.e. n = N-dn, m>=n.
    A_m = XC(dm); A_n = XC(N+1-n)=XC(dn+1); B_m = Br(N-m)=Br(dm-1);
    B_n = Br(dn); A_{m+1}=XC(dm-1); B_{n-1}=Br(dn+1); A_{n-1}=XC(dn+2);
    B_{m+1}=Br(dm-2)."""
    lhs = XC(dm) * Br(dn) + XC(dn + 1) * Br(dm - 1)
    rhs = XC(dm - 1) * Br(dn + 1) + XC(dn + 2) * Br(dm - 2)
    return sp.expand(lhs - rhs)

def sync_margin_CB(dm, dn):
    """A=C (deg N), B (deg N), absolute m = N-dm, n = N-dn, m >= n <=> dm <= dn.
    A_m=Cr(dm); B_n=Br(dn); A_n=Cr(dn); B_m=Br(dm);
    A_{m+1}=Cr(dm-1); B_{n-1}=Br(dn+1); A_{n-1}=Cr(dn+1); B_{m+1}=Br(dm-1)."""
    lhs = Cr(dm) * Br(dn) + Cr(dn) * Br(dm)
    rhs = Cr(dm - 1) * Br(dn + 1) + Cr(dn + 1) * Br(dm - 1)
    return sp.expand(lhs - rhs)

print("\nxC ~_p B margins at head (A-depth dm, B-depth dn), m>=n:")
for dm in range(0, 5):
    for dn in range(max(dm - 1, 0), 5):   # m>=n  <=> N+1-dm >= N-dn <=> dn >= dm-1
        g = sync_margin_xCB(dm, dn)
        gp = sp.Poly(g, h)
        if gp.is_zero:
            continue
        roots = [sp.N(r, 6) for r in sp.real_roots(gp)] if gp.degree() > 0 else []
        lead = gp.LC()
        print(f"  (dm={dm},dn={dn}): deg {gp.degree()}, lead {lead}, "
              f"real roots {roots}")

print("\nC ~_p B margins at head, m>=n (dm <= dn):")
for dm in range(0, 5):
    for dn in range(dm, 5):
        g = sync_margin_CB(dm, dn)
        gp = sp.Poly(g, h)
        if gp.is_zero:
            continue
        roots = [sp.N(r, 6) for r in sp.real_roots(gp)] if gp.degree() > 0 else []
        print(f"  (dm={dm},dn={dn}): deg {gp.degree()}, lead {gp.LC()}, "
              f"real roots {roots}")
