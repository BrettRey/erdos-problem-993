"""Order-of-magnitude N0 for arXiv:2609.20961 under several instantiations of its constants.

Two routes are modelled.
  PAPER route (Sections 5-6 as written): centroid decomposition with b = ceil(n^{1/4}),
    Berry-Esseen with Raic's constant 58 (arXiv:1802.06475v4, Thm 1.1, d = 1: 42 d^{1/4} + 16),
    normal comparison (5.7), complementary event, then Kolmogorov -> characteristic function by
    clipping at +-R' (|chi - e^{-t^2/2}| <= 2 R' |t| D + 4/R'^2, Chebyshev), then the curvature
    integral (6.4) quantified exactly as the Lean `curvature_limit` does it.
  LEAN route (the formalisation's `stdCharFn_tendsto_gaussian`): fixed threshold b with
    b^{a-1} <= eta, direct characteristic-function comparison, N = ceil(Lambda^2/c) + 1.

Constants for Prop. 4.1 are produced by the paper's own proof (p. 5-7): given (p, rho, L=lambda_max),
  C0 = 1/(1-rho), u = p/(2-p), q* = L/(1+L), kappa = q*^{u-p} rho, L_ = q*^{u-1} + G kappa < G,
  H = L^{u-1} sup_y e^{-(u-1)y}(1+C0 y)^{2u/p},  C1 = H^{1/u}/(G^{1/u} - L_^{1/u}),
  C = L * max( sup e^{-y}(1+C0 y)^{2/p},  C1^2 sup e^{-y}(1+G y)^{2/u} ),  G optimised.
Numerics: python-flint arb; optimisation over p uses float evaluations of F (diagnostic); the
chosen p is then certified (F_upper < 1 by ball covering) and its certified upper bound is used
as rho.  All outputs are log10 N0; they are orders of magnitude, not sharp thresholds.
"""
from flint import arb, ctx, fmpq
import json, math

ctx.prec = 256
LN10 = arb(10).log()
pi = arb.pi()
def lg(x): return x.log() / LN10
def f(x, d=6):
    try:
        return float(x.mid())
    except Exception:
        return x.mid().str(d, radius=False)

# ---------------- scalar sup F_L(p) --------------------------------------------------------
def g(y, p, L):
    r = L * (-y).exp(); q = r / (1 + r)
    return y * q ** p / r.log1p()

def F_float(p, L):
    best = 0.0
    for i in range(1, 3001):
        y = 40.0 * i / 3000
        r = L * math.exp(-y); q = r / (1 + r)
        best = max(best, y * q ** p / math.log1p(r))
    return best

def F_upper(p, L, Y=60, cells=6000, thresh=1, maxdepth=14):
    p = arb(p); L = arb(L)
    assert Y > 1 / (float(p.mid()) - 1)
    ub = arb(0); w = arb(Y) / cells
    stack = [((w * i), (w * (i + 1)), 0) for i in range(cells)]
    while stack:
        y0, y1, d = stack.pop()
        v = g(y0.union(y1), p, L)
        if v.upper() >= thresh and d < maxdepth:
            m = (y0 + y1) / 2; stack += [(y0, m, d + 1), (m, y1, d + 1)]; continue
        if v.upper() > ub.upper(): ub = arb(v.upper())
    tail = 2 * arb(Y) * L ** (p - 1) * (-(p - 1) * arb(Y)).exp()
    return ub if ub.upper() > tail.upper() else arb(tail.upper())

# ---------------- sup_{y>=0} e^{-alpha y} (1+beta y)^k ------------------------------------------
def supexp(alpha, beta, k):
    alpha, beta, k = arb(alpha), arb(beta), arb(k)
    ys = k / alpha - 1 / beta
    if ys <= 0: return arb(1)
    return (-alpha * ys).exp() * (k * beta / alpha) ** k

def C_paper(p, rho, L):
    p, rho, L = arb(p), arb(rho), arb(L)
    C0 = 1 / (1 - rho); u = p / (2 - p); qs = L / (1 + L)
    kappa = qs ** (u - p) * rho
    assert kappa < 1
    H = L ** (u - 1) * supexp(u - 1, C0, 2 * u / p)
    Cd = L * supexp(1, C0, 2 / p)
    best = None
    Gmin = qs ** (u - 1) / (1 - kappa)
    for j in range(-100, 400):
        G = arb(10) ** (arb(j) / 10)
        if not G > Gmin * arb(1.001): continue
        Lg = qs ** (u - 1) + G * kappa
        C1 = H ** (1 / u) / (G ** (1 / u) - Lg ** (1 / u))
        Cg = L * C1 ** 2 * supexp(1, G, 2 / u)
        Cm = Cd if Cd > Cg else Cg
        if best is None or Cm < best: best = Cm
    return best

def consts_K(L):
    """paper's c (5.3) and c (Lemma 6.1) for K = [1/4, L]"""
    L = arb(L); q = arb(1) / 4
    c_lin = min(float((q / (2 * (1 + q) ** 4)).mid()), float((L / (2 * (1 + L) ** 4)).mid()))
    beta = min(0.2, float((2 * L / (1 + L) ** 2).mid()), float((L / (1 + L)).mid()))
    c_ch = beta / float((2 * (1 + L) ** 2).mid())
    return arb(c_lin), arb(c_ch)

# ---------------- shared curvature step ----------------------------------------------------------
def curvature_targets(Cup, c_ch):
    cprime = c_ch / (pi ** 2 * Cup)
    c0 = cprime if cprime < arb(1) / 2 else arb(1) / 2
    eps = 1 / (2 * (2 * pi).sqrt())
    tgt = pi * eps / 3
    def tail(R):
        R = arb(R)
        return R * (-c0 * R * R).exp() / c0 + arb.const_sqrt_pi() / (2 * c0 ** 1.5) * (c0.sqrt() * R).erfc()
    lo, hi = arb(0), arb(300)
    for _ in range(300):
        m = (lo + hi) / 2
        if tail(arb(10) ** m) <= tgt: hi = m
        else: lo = m
    R = arb((arb(10) ** hi).mid()).ceil()
    eps1 = pi * eps / (4 * R ** 3)
    return R, eps1, cprime

# ---------------- LEAN route -----------------------------------------------------------------------
def lean_route(a, C, c_lin, c_ch):
    a, C = arb(a), arb(C)
    potL = 1 / (arb(2) ** (1 - a) - 1)
    Cup = arb(1) / 4 + C * (1 + potL)
    R, eps1, cprime = curvature_targets(Cup, c_ch)
    K1 = C * (1 + potL); K2 = (C + C.sqrt()) * (1 + potL)
    e1 = (eps1 / (4 * R)) ** 2 * c_lin / K1; e2 = eps1 * c_lin / (2 * R ** 2 * K2)
    eta = e1 if e1 < e2 else e2
    log10_b = -lg(eta) / (1 - a)
    log10_Lam = log10_b + lg(R + R ** 3 / eps1 + R ** 2 / eps1.sqrt())
    N = 2 * log10_Lam - lg(c_lin)
    return dict(log10_N0=f(N), log10_C=f(lg(C)), one_minus_a=f(1 - a), log10_Cup=f(lg(Cup)),
                log10_cprime=f(lg(cprime)), log10_R=f(lg(R)), log10_eta=f(lg(eta)))

# ---------------- PAPER route ----------------------------------------------------------------------
def paper_route(a, C, c_lin, c_ch, BE=58):
    a, C = arb(a), arb(C)
    Cup = arb(1) / 4 + C / (1 - arb(2) ** (a - 1))
    R, eps1, cprime = curvature_targets(Cup, c_ch)
    Rp = (8 / eps1).sqrt()                         # 4/R'^2 = eps1/2
    D_tgt = eps1 / (4 * Rp * R)                    # 2 R' R D <= eps1/2
    k1 = 1 / (1 - arb(2) ** (a - 1))
    k2 = (1 / (1 - arb(2) ** (1 - 2 * a))).sqrt() if a > arb(1) / 2 else None
    def D(log10n):
        n = arb(10) ** log10n; b = n ** (arb(1) / 4)
        be = BE * b * arb(2).sqrt() / (c_lin * n).sqrt()
        Eh = (C * b ** (a - 1) * k1 / c_lin).sqrt()
        ES = ((C.sqrt() * k2 * n ** a) + C * n * b ** (a - 1) * k1) / (c_lin * n)
        return be + Eh / pi.sqrt() + ES / (2 * pi * arb(1).exp()).sqrt() + 2 * ES
    lo, hi = arb(1), arb(10) ** 9
    assert D(hi) < D_tgt
    for _ in range(200):
        m = (lo + hi) / 2
        if D(m) <= D_tgt: hi = m
        else: lo = m
    return dict(log10_N0=f(hi), log10_C=f(lg(C)), one_minus_a=f(1 - a), log10_Cup=f(lg(Cup)),
                log10_cprime=f(lg(cprime)), log10_R=f(lg(R)), log10_D_target=f(lg(D_tgt)))

out = {}
L12 = 12
c_lin12, c_ch12 = consts_K(L12)
out['paper_c_lin(K=[1/4,12])'] = f(c_lin12); out['paper_c_char(K=[1/4,12])'] = f(c_ch12)

# S1: Lean as formalised (cross-check of lean_chain_n0.py)
b_sc = arb(fmpq(999, 1000)); th = 1 - (1 - b_sc) / 8; pL = arb(3) / 2 + th / 2
rhoL = b_sc ** th * arb(5) ** (1 - th); uL = pL / (2 - pL); aL = 1 - 1 / uL
C0 = 1 / (1 - rhoL); Gc = 2 / (1 - rhoL)
logHc = 32001 * (1 + C0 * 32001).log()
C1 = (((arb(12).log() + logHc) / uL).exp() + 1) / (Gc ** (1 / uL) - (1 + Gc * rhoL) ** (1 / uL))
C_Lean = 12 * ((1 + 3 * C0) ** 3 + C1 ** 2 * (1 + 2 * Gc) ** 2)
out['S1_lean_as_formalised'] = lean_route(aL, C_Lean, 1 / arb(8 * 13 ** 4), 1 / arb(114244))
# S2: paper route, the Lean's Sec.-4 constants, paper's c's
out['S2_paper_route_LeanSec4'] = paper_route(aL, C_Lean, c_lin12, c_ch12)
# S2b: paper route, Lean's p but the paper's own (sharper) C formula
out['S2b_paper_route_Leanp_paperC'] = paper_route(aL, C_paper(pL, rhoL, 12), c_lin12, c_ch12)

# S3: paper route, best p allowed by the paper's own relaxation (b just above 0.99397)
bstar = arb(0.9940)
K32 = 2 * arb(12).sqrt() / arb(1).exp()
best = None
for j in range(1, 60):
    p = arb(2) - arb(0.0032) * arb(j) / 60      # p in (1.9968, 2)
    rho = K32 ** (4 - 2 * p) * bstar ** (2 * p - 3)
    if not rho < 1: continue
    a = 2 - 2 / p
    r = paper_route(a, C_paper(p, rho, 12), c_lin12, c_ch12)
    if best is None or r['log10_N0'] < best[1]['log10_N0']: best = (float(p.mid()), r)
out['S3_paper_route_best_p_relaxation'] = dict(p=best[0], **best[1])

# S4/S5: best p from the true sup at L = 12 (optimise with float F, then certify)
def optimise(L, route, c_lin, c_ch, pgrid):
    best = None
    for p in pgrid:
        rho = F_float(p, L) * 1.0001
        if rho >= 1: continue
        a = 2 - 2 / p
        try:
            r = route(a, C_paper(p, rho, L), c_lin, c_ch)
        except AssertionError:
            continue
        if best is None or r['log10_N0'] < best[2]['log10_N0']: best = (p, rho, r)
    p, rho, r = best
    rho_cert = F_upper(p, L, thresh=rho)        # certify F_L(p) <= rho (the value used)
    assert rho_cert <= arb(rho), (p, rho, rho_cert)
    return dict(p=p, rho_used=rho, rho_certified_upper=float(rho_cert.upper()), **r)

pg12 = [1.8915 + 0.1085 * j / 60 for j in range(1, 60)]
out['S4_paper_route_true_sup_L12'] = optimise(12, paper_route, c_lin12, c_ch12, pg12)
out['S5_lean_route_true_sup_L12'] = optimise(12, lean_route, c_lin12, c_ch12, pg12)
# S6: hypothetical lambda_max = 3 (Section 7 sharpened; stars need lambda > 2)
c_lin3, c_ch3 = consts_K(3)
pg3 = [1.5068 + 0.49 * j / 60 for j in range(1, 60)]
out['S6_lean_route_true_sup_L3'] = optimise(3, lean_route, c_lin3, c_ch3, pg3)
out['S6b_paper_route_true_sup_L3'] = optimise(3, paper_route, c_lin3, c_ch3, pg3)
# S6c/S6d: lambda_max = 4 and 6 (lambda_max = 3 is refuted for forests by mean_alpha_ratio.py:
# a 17-vertex tree has E_3 X / alpha = 0.6542 < 2/3, and disjoint copies keep the ratio)
for Lh, p_lo in [(4, 1.5657), (6, 1.6657)]:
    cl, cc = consts_K(Lh)
    pg = [p_lo + (2 - p_lo) * j / 60 for j in range(1, 60)]
    out[f'S6_lean_route_true_sup_L{Lh}'] = optimise(Lh, lean_route, cl, cc, pg)
    out[f'S6_paper_route_true_sup_L{Lh}'] = optimise(Lh, paper_route, cl, cc, pg)
# S7: floor of the Fourier/curvature machinery: hypothetical exponent a and constant C
for a_h, C_h in [(0.5, 10), (0.0, 1)]:
    out[f'S7_lean_route_hypothetical_a={a_h}_C={C_h}_L12'] = lean_route(a_h, C_h, c_lin12, c_ch12)
    out[f'S7_lean_route_hypothetical_a={a_h}_C={C_h}_L3'] = lean_route(a_h, C_h, c_lin3, c_ch3)
print(json.dumps(out, indent=1, default=str))
