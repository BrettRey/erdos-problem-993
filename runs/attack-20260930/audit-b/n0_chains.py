#!/usr/bin/env python3
"""Effective-constant chains for arXiv:2609.20961 (Fang-Lu-Nevo-Yao-Zheng).

DIAGNOSTIC, ORDER OF MAGNITUDE. Rational inputs fixed by the sources are kept
as Fractions; everything transcendental is evaluated with mpmath at 60 digits.
Nothing here is a certified inequality; the only certified claim of this run is
in certify_lemma42.py (python-flint / Arb).

Scenarios
  L  : the Lean development at commit b2a1d3e, every existential replaced by the
       witness its proof constructs (exists_tail_small replaced by the minimal
       integer tail cut).  Route: constant decomposition size b, direct
       characteristic-function CLT, n >= Lambda^2 / c.
  P  : the PDF as written: Lemma 4.2 by the interpolation in its proof
       (parameters b_L, s = 2 - p optimised over a grid), b = ceil(n^{1/4}),
       Kolmogorov CLT (5.2) with Raic's constant 58, Kolmogorov -> ch.f. on
       compacts by clipping, then the Fourier curvature step (6.4).
  S  : Lemma 4.2 sharpened to the true supremum (certified pair (p, rho) from
       certify_lemma42.py), otherwise the paper's constants, Lean-style route.
  SP : same sharpened lemma, paper route (b = n^{1/4}).
  F  : route floor: a = 2/3 and all root-moment constants = 1 (C_up = 1/4 + 3.4),
       paper's c_sigma and c_F, Lean-style route.  Not attainable by any known
       proof; it shows what this architecture costs even with ideal inputs.
"""
import json, sys
from fractions import Fraction as Fr
import mpmath as mp

mp.mp.dps = 60
PI = mp.pi
LOG10 = lambda x: mp.log10(x)

# ---------------------------------------------------------------- common pieces
EPS_CURV = 1 / (2 * mp.sqrt(2 * PI))       # Lean mean_lc: curvature error budget, = half of 1/sqrt(2 pi)


def tail(R, c0):
    """2 * int_R^inf t^2 exp(-c0 t^2) dt, closed form."""
    R = mp.mpf(R); c0 = mp.mpf(c0)
    return R / c0 * mp.e ** (-c0 * R * R) + mp.sqrt(PI) * mp.erfc(R * mp.sqrt(c0)) / (2 * c0 ** mp.mpf(1.5))


def min_tail_cut(c0, target):
    """minimal integer R >= 1 with tail(R, c0) <= target (bisection on integers)."""
    lo, hi = mp.mpf(1), mp.mpf(1)
    if tail(lo, c0) <= target:
        return lo
    while tail(hi, c0) > target:
        hi *= 2
    while hi - lo > 1:
        mid = mp.floor((lo + hi) / 2)
        if tail(mid, c0) <= target:
            hi = mid
        else:
            lo = mid
    return hi


def curvature_budget(c_prime):
    """Lean curvature_limit: c0 = min(c', 1/2), R = tail cut at pi*eps/3, eps1 = pi*eps/(4R^3),
    X = R^2 + 5R^5/(12 pi eps) (needs sigma^2 >= X)."""
    c0 = min(c_prime, mp.mpf(1) / 2)
    R = min_tail_cut(c0, PI * EPS_CURV / 3)
    eps1 = PI * EPS_CURV / (4 * R ** 3)
    X = R ** 2 + 5 * R ** 5 / (12 * PI * EPS_CURV)
    return dict(c0=c0, R=R, eps1=eps1, X=X)


def lean_route_threshold(a, K1, K2, c_sigma, R, eps1):
    """Lean stdCharFn_tendsto_gaussian(R, eps1): eta, b = ceil(eta^{1/(a-1)})+1,
    Lambda = b (R + R^3/eps1 + R^2/sqrt(eps1)), N = ceil(Lambda^2/c)+1.  Returns log10 values."""
    eta1 = (eps1 / (4 * R)) ** 2 * c_sigma / K1
    eta2 = eps1 * c_sigma / (2 * R ** 2 * K2)
    eta = min(eta1, eta2)
    log10_b = LOG10(eta) / (a - 1)            # b ~ eta^{1/(a-1)} (the +1 and ceiling are negligible)
    log10_Lam = log10_b + LOG10(R + R ** 3 / eps1 + R ** 2 / mp.sqrt(eps1))
    log10_N = 2 * log10_Lam - LOG10(c_sigma)
    return dict(eta1=eta1, eta2=eta2, eta=eta, log10_b=log10_b, log10_Lambda=log10_Lam, log10_N=log10_N,
                binding=("mean-shift (eta1)" if eta1 <= eta2 else "variance-shift (eta2)"))


# ------------------------------------------------ paper constants for Prop 4.1
def sup_exp_poly(C, k, w=mp.mpf(1)):
    """sup_{y>=0} exp(-w y) (1 + C y)^k, exact maximiser."""
    ystar = k / w - 1 / C
    if ystar <= 0:
        return mp.mpf(1)
    return mp.e ** (-w * ystar) * (1 + C * ystar) ** k


def paper_prop41_constants(p, rho, G=None):
    """Explicit versions of the constants in the proof of Prop 4.1 (PDF pp. 6-7).
    C_delta: q(1-q) delta^2 <= C_delta n^a;  C_gamma: q(1-q) gamma^2 <= C_gamma n^{2a}."""
    p = mp.mpf(p); rho = mp.mpf(rho)
    a = 2 - 2 / p
    u = p / (2 - p)
    C0 = 1 / (1 - rho)
    C_delta = 12 * sup_exp_poly(C0, 2 / p)
    qs = mp.mpf(12) / 13
    # G: the proof needs q*^{u-1} + G q*^{u-p} rho < G; pick G minimising C1 over a small grid
    best = None
    for Gc in ([G] if G else [mp.mpf(10) ** e for e in range(-3, 7)]):
        L = qs ** (u - 1) + Gc * qs ** (u - p) * rho
        if not L < Gc:
            continue
        # H^{1/u} = 12^{(u-1)/u} * sup_y exp(-(u-1) y/u) (1 + C0 y)^{2/p}
        Hu = mp.mpf(12) ** ((u - 1) / u) * sup_exp_poly(C0, 2 / p, (u - 1) / u)
        C1 = Hu / (Gc ** (1 / u) - L ** (1 / u))
        C_gamma = 12 * C1 ** 2 * sup_exp_poly(Gc, 2 / u)
        if best is None or C_gamma < best[2]:
            best = (Gc, C1, C_gamma)
    Gc, C1, C_gamma = best
    return dict(p=p, rho=rho, a=a, u=u, C0=C0, C_delta=C_delta, G=Gc, C1=C1, C_gamma=C_gamma)


C_SIGMA_PAPER = mp.mpf(12) / (2 * 13 ** 4)                    # (5.3) at lambda = 12
C_F_PAPER = (mp.mpf(24) / 169) / 338                          # Lemma 6.1: inf beta / 338
RAIC = 58                                                     # Raic 2019 Thm 1.1, 42 d^{1/4} + 16, d = 1
PHI_SCALE = 1 / mp.sqrt(2 * PI * mp.e)                        # sup |z| phi(z); bounds |f'(s)| for s >= 1/2


def lean_style_chain(consts, c_sigma, c_F):
    a = consts["a"]
    one_plus_potL = 1 / (1 - mp.mpf(2) ** (a - 1))
    K1 = consts["C_delta"] * one_plus_potL
    K2 = (consts["C_delta"] + mp.sqrt(consts["C_gamma"])) * one_plus_potL
    C_up = mp.mpf(1) / 4 + K1
    c_prime = c_F / (PI ** 2 * C_up)
    cb = curvature_budget(c_prime)
    th = lean_route_threshold(a, K1, K2, c_sigma, cb["R"], cb["eps1"])
    log10_N0 = max(th["log10_N"], LOG10(cb["X"] / c_sigma), LOG10(1000))
    return dict(K1=K1, K2=K2, C_up=C_up, c_prime=c_prime, **cb, **th, log10_N0=log10_N0)


def paper_route_chain(consts, c_sigma, c_F, detail=False):
    """PDF route: b = n^{1/4}; d_K bound from end of Sec. 5; ch.f. on |t|<=R from d_K by clipping."""
    a = consts["a"]
    K_M = consts["C_delta"] / (1 - mp.mpf(2) ** (a - 1))      # E(M-mu)^2 <= K_M n b^{a-1}
    K_S = mp.sqrt(consts["C_gamma"] / (1 - mp.mpf(2) ** (1 - 2 * a)))  # E|sum gamma(xi-q)| <= K_S n^a
    C_up = mp.mpf(1) / 4 + K_M
    c_prime = c_F / (PI ** 2 * C_up)
    cb = curvature_budget(c_prime)
    R, eps1 = cb["R"], cb["eps1"]
    Rc = mp.sqrt(8 / eps1)                                     # clipping radius: 4/Rc^2 = eps1/2
    dK = eps1 / (4 * R * Rc)                                   # need 2 Rc R dK <= eps1/2
    # four terms of the Kolmogorov bound, each required <= dK/4; solve each for log10 n
    t = dK / 4
    # (i) Berry-Esseen: RAIC * b / sqrt(S), S >= c_sigma n / 2, b = n^{1/4}  ->  RAIC sqrt(2/c) n^{-1/4}
    l_BE = 4 * LOG10(RAIC * mp.sqrt(2 / c_sigma) / t)
    # (ii) mean shift: E|h|/sqrt(pi) <= sqrt(K_M/c) b^{(a-1)/2}/sqrt(pi) ; b = n^{1/4}
    l_mean = 4 * (2 / (1 - a)) * LOG10(mp.sqrt(K_M / (PI * c_sigma)) / t)
    # (iii) variance ratio + complementary event: (2 + PHI_SCALE) E|S - sigma^2| / sigma^2
    #       <= (2+PHI)(K_S n^{a-1} + K_M b^{a-1}) / c  -> two pieces, each <= t/2
    l_var_b = 4 * (1 / (1 - a)) * LOG10(2 * (2 + PHI_SCALE) * K_M / (c_sigma * t))
    l_var_n = (1 / (1 - a)) * LOG10(2 * (2 + PHI_SCALE) * K_S / (c_sigma * t))
    l_curv = LOG10(cb["X"] / c_sigma)
    parts = dict(berry_esseen=l_BE, mean_shift=l_mean, var_shift_b=l_var_b, var_shift_n=l_var_n,
                 curvature_sigma=l_curv, patching=LOG10(1000))
    log10_N0 = max(parts.values())
    out = dict(K_M=K_M, K_S=K_S, C_up=C_up, c_prime=c_prime, Rc=Rc, dK_target=dK, **cb,
               per_step_log10_n=parts, log10_N0=log10_N0,
               dominant=max(parts, key=lambda k: parts[k]))
    return out


# ------------------------------------------------------------------- scenarios
def scenario_L():
    bL = Fr(999, 1000)
    theta = 1 - (1 - bL) / 8                    # Lean exists_scalar_p
    p = Fr(3, 2) + theta / 2                    # = 31999/16000
    u = p / (2 - p)                             # = 31999
    a = 1 - 1 / u                               # = 31998/31999
    k_H = 2 * u / p                             # = 32000 exactly
    assert p == Fr(31999, 16000) and u == 31999 and k_H == 32000
    rho = mp.mpf(bL.numerator) / bL.denominator
    rho = rho ** (mp.mpf(theta.numerator) / theta.denominator) * mp.mpf(5) ** (1 - mp.mpf(theta.numerator) / theta.denominator)
    C0 = 1 / (1 - rho)
    Gc = 2 / (1 - rho)
    nH = int(k_H) + 1                           # ceil(k)+1 with k integral: ceil(32000)+1 = 32001
    # exists_poly_mul_exp_le: M = (1 + c n)^n,  n = ceil(k)+1
    log10_Hc = nH * LOG10(1 + C0 * nH)
    n1 = 3                                      # k = 2/p = 1.0000312..., ceil = 2, +1 = 3
    n2 = 2                                      # k = 2/u = 6.25e-5, ceil = 1, +1 = 2
    M1 = (1 + C0 * n1) ** n1
    M2 = (1 + Gc * n2) ** n2
    uf = mp.mpf(31999)
    twelveHc_pow = mp.power(10, (LOG10(12) + log10_Hc) / uf)       # (12 Hc)^{1/u}
    den = Gc ** (1 / uf) - (1 + Gc * rho) ** (1 / uf)
    C1 = (twelveHc_pow + 1) / den
    C = 12 * (M1 + C1 ** 2 * M2)
    af = mp.mpf(31998) / 31999
    potL = 1 / (mp.mpf(2) ** (1 - af) - 1)
    K1 = C * (1 + potL)
    K2 = (C + mp.sqrt(C)) * (1 + potL)
    C_v = mp.mpf(1) / 4 + C * (1 + potL)
    c_sigma = mp.mpf(1) / (8 * 13 ** 4)          # Lean linear_variance
    c_F = mp.mpf(1) / 114244                     # Lean charFn_bound, = 1/338^2
    c_prime = c_F / (PI ** 2 * C_v)               # Lean norm_charFn_std_le
    cb = curvature_budget(c_prime)
    th = lean_route_threshold(af, K1, K2, c_sigma, cb["R"], cb["eps1"])
    log10_N0 = max(th["log10_N"], LOG10(cb["X"] / c_sigma), LOG10(1000))
    return dict(p=str(p), u=str(u), a=str(a), one_minus_a=str(1 - a), rho=rho, C0=C0, Gc=Gc,
                log10_Hc=log10_Hc, M1=M1, M2=M2, C1=C1, C=C, potL=potL, K1=K1, K2=K2, C_v=C_v,
                c_sigma=c_sigma, c_F=c_F, c_prime=c_prime, **cb, **th, log10_N0=log10_N0)


def scenario_P(grid=True):
    A = 2 * mp.sqrt(12) / mp.e
    bmin = mp.findroot(lambda b: 12 * mp.e ** (-1 - 3 * b / 2) - b, 0.994)
    best = None
    rows = []
    for fb in [mp.mpf(x) for x in ("0.05", "0.1", "0.2", "0.35", "0.5", "0.65", "0.8", "0.95")]:      # position of b_L in (bmin, 1)
        bL = bmin + fb * (1 - bmin) if fb != 0 else bmin
        smax = -mp.log(bL) / (2 * (mp.log(A) - mp.log(bL)))
        for fs in [mp.mpf(x) for x in ("0.2", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "0.95", "0.98", "0.99", "0.995", "0.999")]:
            s = fs * smax
            p = 2 - s
            rho = A ** (2 * s) * bL ** (1 - 2 * s)
            if not rho < 1:
                continue
            cst = paper_prop41_constants(p, rho)
            ch = paper_route_chain(cst, C_SIGMA_PAPER, C_F_PAPER)
            rows.append(dict(bL=bL, s=s, p=p, one_minus_a=1 - cst["a"], rho=rho, C0=cst["C0"],
                             log10_N0=ch["log10_N0"], dominant=ch["dominant"]))
            if best is None or ch["log10_N0"] < best[1]["log10_N0"]:
                best = (cst, ch, dict(bL=bL, s=s))
    return dict(bmin=bmin, A=A, grid=rows, best_params=best[2], best_consts=best[0], best_chain=best[1])


def scenario_S(p, rho):
    cst = paper_prop41_constants(p, rho)
    return dict(consts=cst, lean_style=lean_style_chain(cst, C_SIGMA_PAPER, C_F_PAPER),
                paper_route=paper_route_chain(cst, C_SIGMA_PAPER, C_F_PAPER))


def scenario_F():
    cst = dict(a=mp.mpf(2) / 3, C_delta=mp.mpf(1), C_gamma=mp.mpf(1))
    return dict(consts=cst, lean_style=lean_style_chain(cst, C_SIGMA_PAPER, C_F_PAPER),
                paper_route=paper_route_chain(cst, C_SIGMA_PAPER, C_F_PAPER))


def fmt(x):
    if isinstance(x, mp.mpf):
        return mp.nstr(x, 8)
    if isinstance(x, dict):
        return {k: fmt(v) for k, v in x.items()}
    if isinstance(x, list):
        return [fmt(v) for v in x]
    return x


if __name__ == "__main__":
    cert = json.load(open("data/lemma42_certificates.json"))
    out = {}
    out["L_lean_as_formalised"] = scenario_L()
    out["P_paper_as_written"] = scenario_P()
    out["S_sharp_lemma42"] = {}
    for c in cert["certified"]:
        if Fr(c["p"]) >= 2:          # p = 2 certificate is for comparison with (4.5) only; Prop 4.1 needs p < 2
            continue
        key = f"p={c['p']}"
        out["S_sharp_lemma42"][key] = scenario_S(mp.mpf(Fr(c["p"]).numerator) / Fr(c["p"]).denominator,
                                                 mp.mpf(Fr(c["rho_upper"]).numerator) / Fr(c["rho_upper"]).denominator)
    out["F_route_floor"] = scenario_F()
    json.dump(fmt(out), open("data/n0_chains.json", "w"), indent=1)
    L = out["L_lean_as_formalised"]; P = out["P_paper_as_written"]
    print("L (Lean): 1-a =", L["one_minus_a"], " C =", mp.nstr(L["C"], 5), " C_v =", mp.nstr(L["C_v"], 5),
          " c' =", mp.nstr(L["c_prime"], 5), " R =", mp.nstr(L["R"], 5), " eps1 =", mp.nstr(L["eps1"], 5),
          " eta =", mp.nstr(L["eta"], 5), L["binding"], " log10 b =", mp.nstr(L["log10_b"], 6),
          " log10 N0 =", mp.nstr(L["log10_N0"], 6))
    bc = P["best_chain"]; bp = P["best_params"]; cs = P["best_consts"]
    print("P (paper): best bL =", mp.nstr(bp["bL"], 8), " s =", mp.nstr(bp["s"], 5), " 1-a =", mp.nstr(1 - cs["a"], 5),
          " C0 =", mp.nstr(cs["C0"], 5), " C_delta =", mp.nstr(cs["C_delta"], 5), " C_gamma =", mp.nstr(cs["C_gamma"], 5))
    print("   C_up =", mp.nstr(bc["C_up"], 5), " c' =", mp.nstr(bc["c_prime"], 5), " R =", mp.nstr(bc["R"], 5),
          " dK target =", mp.nstr(bc["dK_target"], 5))
    print("   per-step log10 n:", {k: mp.nstr(v, 6) for k, v in bc["per_step_log10_n"].items()})
    print("   log10 N0 =", mp.nstr(bc["log10_N0"], 6), " dominant:", bc["dominant"])
    for k, v in out["S_sharp_lemma42"].items():
        cst = v["consts"]
        print(f"S {k}: 1-a = {mp.nstr(1-cst['a'],5)} C0 = {mp.nstr(cst['C0'],5)} C_delta = {mp.nstr(cst['C_delta'],5)}"
              f" C_gamma = {mp.nstr(cst['C_gamma'],5)} | Lean-style log10 N0 = {mp.nstr(v['lean_style']['log10_N0'],6)}"
              f" ({v['lean_style']['binding']}), paper-route log10 N0 = {mp.nstr(v['paper_route']['log10_N0'],6)}"
              f" ({v['paper_route']['dominant']})")
    Fl = out["F_route_floor"]
    print("F floor: Lean-style log10 N0 =", mp.nstr(Fl["lean_style"]["log10_N0"], 6), " R =", mp.nstr(Fl["lean_style"]["R"], 5),
          "| paper-route log10 N0 =", mp.nstr(Fl["paper_route"]["log10_N0"], 6), Fl["paper_route"]["dominant"],
          {k: mp.nstr(v, 5) for k, v in Fl["paper_route"]["per_step_log10_n"].items()})
