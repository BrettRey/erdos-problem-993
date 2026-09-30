"""Certified (arb ball arithmetic) checks on the two scalar inputs that pin the exponent a of
arXiv:2609.20961, Prop. 4.1, and on the Section 7 mean-range constant.

(1) Lemma 4.2 (p. 5): the paper needs b<1 with 12 e^{-1-3b/2} < b "possible because 12 < e^{5/2}".
    We certify 12 < e^{5/2}, bracket the least admissible b, and the least p for which the
    paper's own interpolation constant rho = (2 sqrt12/e)^{4-2p} b^{2p-3} is < 1.
(2) The true one-dimensional supremum the lemma needs,
        F_L(p) = sup_{0<r<=L} log(L/r) * (r/(1+r))^p / log(1+r),     L = lambda_max,
    (with y = log(lambda/r); the sup over lambda in [1/4, L] is attained at lambda = L because
    the expression increases with lambda at fixed r).  The lemma holds with exponent p iff
    F_L(p) < 1.  Certified: lower bounds by point evaluation, upper bounds by covering
    y = log(L/r) in [0, Y] with ball intervals plus an analytic tail bound for y > Y.
(3) Section 7 (p. 13): E_12 X >= log(13) alpha / Q*(12) with
        Q*(lambda) = sqrt((1+lambda)/lambda) * (3 log2/sqrt2 + 2 log((sqrt(lambda)+sqrt(1+lambda))/(1+sqrt2))).
    Certify log13/Q*(12) > 17/25 and locate the least lambda with log(1+lambda)/Q*(lambda) > 2/3.

Run from project root: venv/bin/python runs/attack-20260930/audit-a/scalar_and_meanrange.py
"""
from flint import arb, ctx, fmpq
import json

ctx.prec = 200
e = arb(1).exp()
out = {}

# ---- (1) paper's relaxation ------------------------------------------------------------
out['e^{5/2}'] = arb(fmpq(5, 2)).exp().str(12)
out['12 < e^{5/2}'] = bool(arb(12) < arb(fmpq(5, 2)).exp())

def h(b):  # 12 e^{-1-3b/2} - b ; decreasing in b
    return 12 * (-1 - 3 * b / 2).exp() - b
lo, hi = arb(fmpq(9, 10)), arb(1)
assert h(lo) > 0 and h(hi) < 0
for _ in range(120):
    m = (lo + hi) / 2
    if h(m) > 0: lo = m
    else: hi = m
b_min = hi
out['paper_b_min (12e^{-1-3b/2}=b)'] = b_min.str(15)

K32 = 2 * arb(12).sqrt() / e     # the paper's exponent-3/2 constant 2 sqrt12 / e
out['2sqrt12/e'] = K32.str(12)
def rho_paper(p, b):
    return K32 ** (4 - 2 * p) * b ** (2 * p - 3)
# least p with rho_paper < 1 at b = b_min (rho decreasing in p near 2)
lo, hi = arb(fmpq(3, 2)), arb(2)
for _ in range(200):
    m = (lo + hi) / 2
    if rho_paper(m, b_min) < 1: hi = m
    else: lo = m
p_relax = hi
a_relax = 2 - 2 / p_relax
out['paper_relaxation_p_min'] = p_relax.str(12)
out['paper_relaxation_1-a_max'] = (1 - a_relax).str(8)
out['paper_relaxation_u_min=1/(1-a)'] = (1 / (1 - a_relax)).str(8)

# ---- (2) the true sup F_L(p) --------------------------------------------------------------
def g(y, p, L):
    r = L * (-y).exp()
    q = r / (1 + r)
    return y * q ** p / r.log1p()

def F_lower(p, L, grid=4000, Y=40.0):
    """certified lower bound: max over grid points (each evaluation is a valid lower bound
    for the sup once we take the lower endpoint of the ball)."""
    best, arg = arb(0), 0.0
    for i in range(1, grid + 1):
        y = arb(Y * i / grid)
        v = g(y, p, L)
        if v.lower() > best.lower():
            best, arg = arb(v.lower()), Y * i / grid
    # golden refine around arg (still certified lower bounds)
    a0, b0 = max(arg - Y / grid, 1e-9), arg + Y / grid
    for _ in range(80):
        m1 = a0 + (b0 - a0) * 0.382; m2 = a0 + (b0 - a0) * 0.618
        v1 = g(arb(m1), p, L); v2 = g(arb(m2), p, L)
        if v1.mid() > v2.mid(): b0 = m2
        else: a0 = m1
        for v in (v1, v2):
            if v.lower() > best.lower(): best = arb(v.lower())
    return best, (a0 + b0) / 2

def F_upper(p, L, Y=60, cells=6000, thresh=1, maxdepth=14):
    """certified upper bound: cover y in [0, Y] by ball intervals, subdividing adaptively any
    cell whose enclosure reaches `thresh`; the tail y > Y is bounded analytically:
    r <= L e^{-Y} <= 1, q <= r, log(1+r) >= r/2  =>  g <= 2 y L^{p-1} e^{-(p-1)y},
    decreasing for y > 1/(p-1); we require Y > 1/(p-1)."""
    assert Y > 1 / (float(p.mid()) - 1)
    ub = arb(0)
    w = arb(Y) / cells
    stack = [((w * i), (w * (i + 1)), 0) for i in range(cells)]
    while stack:
        y0, y1, d = stack.pop()
        v = g(y0.union(y1), p, L)
        if v.upper() >= thresh and d < maxdepth:
            m = (y0 + y1) / 2
            stack.append((y0, m, d + 1)); stack.append((m, y1, d + 1))
            continue
        if v.upper() > ub.upper(): ub = arb(v.upper())
    tail = 2 * arb(Y) * arb(L) ** (p - 1) * (-(p - 1) * arb(Y)).exp()
    if tail.upper() > ub.upper(): ub = arb(tail.upper())
    return ub

res = {}
for L in [12, 11, 8, 6, 4, 3, fmpq(5, 2)]:
    La = arb(L)
    rec = {}
    lo_b, _ = F_lower(arb(2), La)
    rec['F(2)_lower'] = lo_b.str(8)
    # bisection for p* using midpoints of lower bound (diagnostic), then certify a bracket
    plo, phi = 1.02, 2.0
    for _ in range(40):
        pm = (plo + phi) / 2
        v, _ = F_lower(arb(pm), La, grid=1500)
        if float(v.mid()) >= 1: plo = pm
        else: phi = pm
    # certify: F(plo) >= 1 (so p* >= plo) and F(phi + 0.002) < 1 (so p* <= phi + 0.002)
    vlo, _ = F_lower(arb(plo), La)
    p_up = phi + 0.002
    vup = F_upper(arb(p_up), La)
    rec['p*_certified_lower'] = plo if vlo >= 1 else None
    rec['F(p_lower)_lower'] = vlo.str(8)
    rec['p_certified_admissible'] = p_up if vup < 1 else None
    rec['F(p_admissible)_upper'] = vup.str(8)
    rec['1-a at p_admissible'] = (2.0 / p_up - 1)
    rec['u=1/(1-a) at p_admissible'] = p_up / (2 - p_up)
    res[str(L)] = rec
out['true_sup'] = res

# ---- (3) Section 7 ----------------------------------------------------------------------
def Qstar(l):
    l = arb(l)
    s2 = arb(2).sqrt()
    return ((1 + l) / l).sqrt() * (3 * arb(2).log() / s2 + 2 * ((l.sqrt() + (1 + l).sqrt()) / (1 + s2)).log())
ratio12 = arb(13).log() / Qstar(12)
out['log13/Q*(12)'] = ratio12.str(10)
out['log13/Q*(12) > 17/25'] = bool(ratio12 > arb(fmpq(17, 25)))
out['margin log13/Q*(12) - 17/25'] = (ratio12 - arb(fmpq(17, 25))).str(6)
out['log13/Q*(12) > 64/95 (Lean)'] = bool(ratio12 > arb(fmpq(64, 95)))
def ratio(l): return (1 + arb(l)).log() / Qstar(l)
lo, hi = arb(2), arb(12)
assert ratio(lo) < arb(2) / 3 and ratio(hi) > arb(2) / 3
for _ in range(100):
    m = (lo + hi) / 2
    if ratio(m) > arb(2) / 3: hi = m
    else: lo = m
out['least lambda with log(1+l)/Q*(l) > 2/3 (Sec. 7 method)'] = hi.str(10)
out['ratio at lambda=11'] = ratio(11).str(8)
print(json.dumps(out, indent=1, default=str))
