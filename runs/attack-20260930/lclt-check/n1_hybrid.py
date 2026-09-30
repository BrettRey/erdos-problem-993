#!/usr/bin/env python3
"""Best explicit N1 for the hybrid: route-lclt binomial-smoothing certificate (C)
+ Fang et al. Section-5 centroid Kolmogorov chain, with the wave-1 audits' constants.

DIAGNOSTIC ORDER-OF-MAGNITUDE (mpmath 60 digits), like audit-b/n0_chains.py, whose
constant functions are imported unchanged (paper_prop41_constants, C_SIGMA_PAPER,
C_F_PAPER, RAIC, PHI_SCALE) so the comparison with Fang's N0 is like for like.

Certificate requirement (claim C, re-derived in lclt-check):
   d_K(M' | E0, G_d) <= theta * eps_K,   eps_K = c_* ((1-eta) rho_0)^{3/2},
   c_* = phi(0)/int|phi'''| = 0.2641979111 (arb-certified here), theta = 1/2
   (half the budget left for the tail term and the O(1/v), O(x^2) corrections).
Kolmogorov chain for M' ASSUMED to have the structure of Fang (5.2)-(5.7):
   (i)   Berry-Esseen   RAIC_d * b / sqrt(S'),  S' >= eta * c_sigma * n / 2
   (ii)  mean shift     sqrt(K_M / (pi eta c_sigma)) b^{(a-1)/2}
   (iii) variance shift (2 + PHI)(K_S n^{a-1} + K_M b^{a-1}) / (eta c_sigma)
   Each of (i), (ii) <= t = theta eps_K / 4 and each piece of (iii) <= t/2 (audit-b's split).
   Why these transfer: given the centroid reveal, E[M'|rev] = E[X|rev] - p B0 and
   Var(M'|rev) = Var(X|rev) - v0, so (ii),(iii) are the SAME random quantities as for X
   (up to the E0 conditioning); only the variance normaliser drops from sigma^2 to
   Var M' >= eta sigma^2 (approx).  Step (i) does NOT transfer: given the reveal, M' is
   a binomial mixture over a sum of independent 2-vectors (|J_j|, B_j), not a sum of
   independent summands; we price it at Raic's d=2 constant 42*2^{1/4}+16 = 65.94 as a
   hypothesis.
Routes:  paper route b = n^{1/4};  fixed-b route (b chosen from (ii),(iii-b), then n from (i),(iii-n)).
rho_0:   (a) empirical 0.4448 (route-lclt KT2, census n<=20; unproved)
         (b) provable from the audits: p(1-p)E B_J = E|I cap L|/(1+lam) >= (k/2)/13 >= n/104
             (larger side, k >= n/4, lam <= 12 by Fang Prop 7.1) and sigma^2 <= C_up n with
             C_up = 1/4 + K_M (audit-b, Fang Prop 5.1 upper bound), so rho_0 >= 1/(104 C_up).
"""
import sys, json
sys.path.insert(0, "/Users/brettreynolds/projects/LLM-CLI-projects/papers/queue/erdos-problem-993/runs/attack-20260930/audit-b")
import mpmath as mp
from n0_chains import paper_prop41_constants, C_SIGMA_PAPER, C_F_PAPER, RAIC, PHI_SCALE
mp.mp.dps = 60
PI = mp.pi; L10 = mp.log10
C_STAR = mp.mpf("0.2641979111219313291")
RAIC2 = 42 * mp.mpf(2) ** mp.mpf(0.25) + 16

def chain(cst, rho0, eta, theta=mp.mpf(1) / 2, c_sigma=C_SIGMA_PAPER, be=RAIC2):
    a = cst["a"]
    K_M = cst["C_delta"] / (1 - mp.mpf(2) ** (a - 1))
    K_S = mp.sqrt(cst["C_gamma"] / (1 - mp.mpf(2) ** (1 - 2 * a)))
    C_up = mp.mpf(1) / 4 + K_M
    if rho0 == "provable":
        rho = 1 / (104 * C_up)
    else:
        rho = mp.mpf(rho0)
    epsK = C_STAR * ((1 - eta) * rho) ** mp.mpf(1.5)
    t = theta * epsK / 4
    cs = eta * c_sigma
    # paper route, b = n^{1/4}
    l_BE = 4 * L10(be * mp.sqrt(2 / cs) / t)
    l_mean = 4 * (2 / (1 - a)) * L10(mp.sqrt(K_M / (PI * cs)) / t)
    l_var_b = 4 * (1 / (1 - a)) * L10(2 * (2 + PHI_SCALE) * K_M / (cs * t))
    l_var_n = (1 / (1 - a)) * L10(2 * (2 + PHI_SCALE) * K_S / (cs * t))
    paper = dict(berry_esseen=l_BE, mean_shift=l_mean, var_shift_b=l_var_b, var_shift_n=l_var_n)
    # fixed-b route
    lb_mean = (2 / (1 - a)) * L10(mp.sqrt(K_M / (PI * cs)) / t)
    lb_var = (1 / (1 - a)) * L10(2 * (2 + PHI_SCALE) * K_M / (cs * t))
    lb = max(lb_mean, lb_var)
    ln_BE = 2 * (L10(be) + lb - L10(t)) + L10(2 / cs)
    fixed = dict(log10_b=lb, b_binding=("mean_shift" if lb_mean >= lb_var else "var_shift_b"),
                 berry_esseen_given_b=ln_BE, var_shift_n=l_var_n)
    return dict(one_minus_a=1 - a, K_M=K_M, K_S=K_S, C_up=C_up, rho0=rho, epsK=epsK, t=t,
                paper_route=paper, paper_log10_N1=max(paper.values()), paper_dominant=max(paper, key=lambda k: paper[k]),
                fixed_route=fixed, fixed_log10_N1=max(ln_BE, l_var_n),
                fixed_dominant=("berry_esseen with b from " + fixed["b_binding"]) if ln_BE >= l_var_n else "var_shift_n")

def per_reveal_floor(cst, frac=mp.mpf("0.1"), c_sigma=C_SIGMA_PAPER):
    """Not the hybrid: if the curvature were averaged per reveal, the mean shift would only
    need E(M-mu)^2/sigma^2 <= frac (a constant).  b from K_M b^{a-1}/c_sigma <= frac."""
    a = cst["a"]; K_M = cst["C_delta"] / (1 - mp.mpf(2) ** (a - 1))
    return (1 / (1 - a)) * L10(K_M / (c_sigma * frac))

if __name__ == "__main__":
    scen = {
        "S_p=19/10_rho<=199/200 (certified sharpened Lemma 4.2, audit-b)": paper_prop41_constants(mp.mpf(19) / 10, mp.mpf(199) / 200),
        "S_p=193/100_rho<=49/50": paper_prop41_constants(mp.mpf(193) / 100, mp.mpf(49) / 50),
        "P_paper_as_written_best (s=0.0029905,b_L=0.99427146)": None,
        "F_floor_a=2/3_C=1 (hypothetical)": dict(a=mp.mpf(2) / 3, C_delta=mp.mpf(1), C_gamma=mp.mpf(1)),
    }
    # P scenario: rebuild best params reported by audit-b (n0_chains_stdout.txt)
    A = 2 * mp.sqrt(12) / mp.e
    s = mp.mpf("0.0029905"); bL = mp.mpf("0.99427146")
    scen["P_paper_as_written_best (s=0.0029905,b_L=0.99427146)"] = paper_prop41_constants(2 - s, A ** (2 * s) * bL ** (1 - 2 * s))
    fang = {"S_p=19/10_rho<=199/200 (certified sharpened Lemma 4.2, audit-b)": dict(lean_style=2069.69, paper_route=5375.03),
            "S_p=193/100_rho<=49/50": dict(lean_style=2827.88, paper_route=7430.42),
            "P_paper_as_written_best (s=0.0029905,b_L=0.99427146)": dict(paper_route=243772.0),
            "F_floor_a=2/3_C=1 (hypothetical)": dict(lean_style=240.155, paper_route=537.077)}
    out = {}
    for name, cst in scen.items():
        out[name] = dict(fang_log10_N0_audit_b=fang[name], per_reveal_mean_shift_log10_b=float(per_reveal_floor(cst)))
        for rho0 in ("0.4448", "provable"):
            for eta in (mp.mpf(3) / 10, mp.mpf(1) / 2):
                r = chain(cst, rho0, eta)
                key = f"rho0={rho0},eta={float(eta)}"
                out[name][key] = {k: (float(v) if isinstance(v, mp.mpf) else ({kk: (float(vv) if isinstance(vv, mp.mpf) else vv) for kk, vv in v.items()} if isinstance(v, dict) else v)) for k, v in r.items()}
    # route-lclt's own formula N_K = (K / (c_* ((1-eta) rho0)^{3/2}))^2 : what K would the Fang chain correspond to?
    out["route_lclt_NK_formula"] = {f"eta={e}": float((1 / (C_STAR * ((1 - mp.mpf(e)) * mp.mpf("0.4448")) ** mp.mpf(1.5))) ** 2) for e in ("0.3", "0.5")}
    json.dump(out, open("n1_hybrid.json", "w"), indent=1)
    for name, d in out.items():
        if name == "route_lclt_NK_formula":
            print(name, d); continue
        print("==", name, "| Fang N0 (audit-b):", d["fang_log10_N0_audit_b"], "| per-reveal mean-shift log10 b:", round(d["per_reveal_mean_shift_log10_b"], 1))
        for key, r in d.items():
            if not key.startswith("rho0"): continue
            print(f"  {key}: 1-a={r['one_minus_a']:.4g} K_M={r['K_M']:.4g} C_up={r['C_up']:.4g} rho0={r['rho0']:.4g} epsK={r['epsK']:.4g}"
                  f" | paper route log10 N1={r['paper_log10_N1']:.1f} ({r['paper_dominant']}) {{" +
                  ", ".join(f"{k}:{v:.1f}" for k, v in r['paper_route'].items()) + "}"
                  f" | fixed-b log10 N1={r['fixed_log10_N1']:.1f} (log10 b={r['fixed_route']['log10_b']:.1f}, {r['fixed_dominant']})")
