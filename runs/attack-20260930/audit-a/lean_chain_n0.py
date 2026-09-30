"""Instantiate every existential in the Lean chain of arXiv:2609.20961's formalisation
(github.com/junwei-lu/Erdos_993_Tree_Independent_Set_Unimodality, commit b2a1d3e) with the
witness the Lean proof constructs, and report log10 of the resulting N0.

Every step is taken from the Lean source (file:line given in comments; lean_src/ is a clone
of b2a1d3e).  Ball arithmetic (python-flint arb, 400 bits) throughout; huge quantities are
carried as logarithms.  The only existential without an explicit witness in the Lean is the
Gaussian-tail cutoff R (Curvature.lean:304 `exists_tail_small`, obtained from a limit); we
use the least natural R that satisfies the stated inequality, via the closed form of
  int_{|t|>R} t^2 exp(-c0 t^2) dt = R exp(-c0 R^2)/c0 + sqrt(pi)/(2 c0^{3/2}) erfc(sqrt(c0) R).
Any valid witness is >= this one, so the N0 reported is the smallest this proof can deliver.

Run: venv/bin/python lean_chain_n0.py   (from the erdos-problem-993 project root)
"""
from flint import arb, ctx, fmpq
import json, sys

ctx.prec = 400
LN10 = arb(10).log()
def lg(x):            # log10 of a positive arb
    return x.log() / LN10
def s(x, d=6):        # short printable midpoint
    return x.mid().str(d, radius=False)

out = {}
pi = arb.pi()

# --- RootMoments.lean -------------------------------------------------------------
b_sc  = arb(fmpq(999, 1000))                     # numeric_b, :446 ; exists_scalar_two :456
theta = 1 - (1 - b_sc) / 8                       # RootMoments.lean:535
p     = arb(3) / 2 + theta / 2                   # RootMoments.lean:538  -> 31999/16000
rho   = b_sc ** theta * arb(5) ** (1 - theta)    # RootMoments.lean:538
u     = p / (2 - p)                              # RootMoments.lean:1110 -> 31999
a     = 1 - 1 / u                                # RootMoments.lean:1133 -> 31998/31999
C0    = 1 / (1 - rho)
Gc    = 2 / (1 - rho)                            # RootMoments.lean:1115
# exists_poly_mul_exp_le (:648): M = (1 + c n)^n with n = ceil(k)+1
k_H   = 2 * u / p                                # exactly 32000
n_H   = 32001                                    # ceil(32000)+1 (k_H is an exact integer)
assert (k_H - 32000).contains(0)
logHc = n_H * (1 + C0 * n_H).log()               # Hc = (1 + C0*n_H)^n_H   RootMoments.lean:1121 via :648
# C1 = ((12 Hc)^{1/u} + 1) / (Gc^{1/u} - (1 + Gc rho)^{1/u})     RootMoments.lean:1125
num   = ((arb(12).log() + logHc) / u).exp() + 1
den   = Gc ** (1 / u) - (1 + Gc * rho) ** (1 / u)
C1    = num / den
# M1: c = C0, k = 2/p (=32000/31999, ceil = 2) -> n = 3 ; M2: c = Gc, k = 2/u (ceil 1) -> n = 2
M1    = (1 + C0 * 3) ** 3
M2    = (1 + Gc * 2) ** 2
Crm   = 12 * (M1 + C1 ** 2 * M2)                 # root_moments constant C, RootMoments.lean:1133
out['sec4'] = dict(p=s(p, 10), one_minus_a=s(1 - a, 8), u=s(u, 8), rho=s(rho, 8),
                   C0=s(C0), log10_Hc=s(lg(arb(1)) + logHc / LN10), C1=s(C1),
                   M1=s(M1), M2=s(M2), log10_C=s(lg(Crm)))

# --- CLT.lean ----------------------------------------------------------------------
potL  = 1 / (arb(2) ** (1 - a) - 1)             # CLT.lean:622
Cv    = arb(1) / 4 + Crm * (1 + potL)            # CLT.lean:1067
c_lin = 1 / arb(8 * 13 ** 4)                     # CLT.lean:1080
K1    = Crm * (1 + potL)                         # CLT.lean:1870
K2    = (Crm + Crm.sqrt()) * (1 + potL)          # CLT.lean:1871
# --- Fourier.lean / Curvature.lean ------------------------------------------------
c_ch  = 1 / arb(114244)                          # Fourier.lean:700 (=1/338^2)
cprime = c_ch / (pi ** 2 * Cv)                   # Curvature.lean:80
c0    = cprime if cprime < arb(1) / 2 else arb(1) / 2   # Curvature.lean:341
eps   = 1 / (2 * (2 * pi).sqrt())                # Curvature.lean:634
target_tail = pi * eps / 3                       # Curvature.lean:348

def tail(R):
    R = arb(R)
    return R * (-c0 * R * R).exp() / c0 + arb.const_sqrt_pi() / (2 * c0 ** 1.5) * (c0.sqrt() * R).erfc()

# least natural R >= 1 with tail(R) <= target (tail is decreasing in R); bisection on log R
lo, hi = arb(0), arb(200)                        # log10 R bounds
for _ in range(400):
    mid = (lo + hi) / 2
    if tail(arb(10) ** mid) <= target_tail:
        hi = mid
    else:
        lo = mid
R = (arb(10) ** hi).mid()
R = arb(R).ceil()                                # natural number, tail(R) <= target
assert tail(R) <= target_tail
eps1  = pi * eps / (4 * R ** 3)                  # Curvature.lean:350
e1 = (eps1 / (4 * R)) ** 2 * c_lin / K1          # CLT.lean:1874
e2 = eps1 * c_lin / (2 * R ** 2 * K2)
eta = e1 if e1 < e2 else e2
log10_b = -lg(eta) * u                           # b = ceil(eta^{1/(a-1)}) + 1  CLT.lean:1887  (1/(1-a) = u)
log10_Lam = log10_b + lg(R + R ** 3 / eps1 + R ** 2 / eps1.sqrt())   # CLT.lean:1899
log10_N1 = 2 * log10_Lam - lg(c_lin)             # N1 = ceil(Lam^2/c)+1  CLT.lean:1910
X = R ** 2 + 5 * R ** 5 / (12 * pi * eps)        # Curvature.lean:353
log10_Ncurv_aux = lg(X / c_lin)
out['clt_curv'] = dict(log10_potL=s(lg(potL)), log10_Cv=s(lg(Cv)), c_lin=s(c_lin),
                       c_char=s(c_ch), log10_cprime=s(lg(cprime)), log10_R=s(lg(R)),
                       log10_eps1=s(lg(eps1)), log10_eta=s(lg(eta)), log10_b_threshold=s(log10_b),
                       log10_Lambda=s(log10_Lam), log10_N1=s(log10_N1),
                       log10_X_over_c=s(log10_Ncurv_aux))
out['log10_N0_lean'] = s(log10_N1, 8)
out['decomposition'] = dict(
    note="log10 N0 ~= 2*u*(-log10 eta) + 2*log10(R^3/eps1) - log10 c_lin",
    two_u=s(2 * u, 8), minus_log10_eta=s(-lg(eta)))
print(json.dumps(out, indent=1))
