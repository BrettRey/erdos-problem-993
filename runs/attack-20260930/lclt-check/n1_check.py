#!/usr/bin/env python3
"""Independent recomputation of the hybrid threshold N1 (lclt-check, resumed run).
DIAGNOSTIC ORDER-OF-MAGNITUDE ARITHMETIC (Python floats, log10 scale). Does NOT import
audit-b/n0_chains.py; it re-implements the chain from the formulas stated in
audit-b/RETURN.json (steps 1-6) and checks the result against the partial n1_hybrid.json.

Hybrid = route-lclt certificate (C) + a Fang-Section-5-type centroid Kolmogorov chain
applied to M' (ASSUMED transfer; see report).

Certificate tolerance (claim C):  d_K(M'|E0, G) < eps_K := c_* ((1-eta) rho0)^{3/2},
   c_* = 1/(2 + 8 e^{-3/2}) = 0.264197911... (Arb-certified in kernel_check_exact.py).
Budget: theta = 1/2 of eps_K for the Kolmogorov chain (rest for tail, discretisation and
O(1/v) kernel corrections), split into 4 terms, each <= t := theta*eps_K/4.
Chain terms (normaliser Var M' >= eta*sigma^2 >= eta*c_sigma*n):
 (i)   Berry-Esseen      RAIC * b * sqrt(2/(eta c_sigma n))                  <= t
 (ii)  mean shift        sqrt(K_M b^{a-1} / (pi eta c_sigma))                 <= t
 (iii) variance shift b  (2+PHI) K_M b^{a-1} /(eta c_sigma)                   <= t/2
 (iv)  variance shift n  (2+PHI) K_S n^{a-1} /(eta c_sigma)                   <= t/2
 K_M = C_delta/(1-2^{a-1}),  C_delta = 12 sup_y e^{-y}(1+C0 y)^{2/p},  C0 = 1/(1-rho),
 a = 2 - 2/p;  K_S = sqrt(C_gamma/(1-2^{1-2a})).  Sources: audit-b RETURN.json steps 1-6
 (Fang et al. arXiv:2609.20961 Sec. 4-5, eqs. (4.3), (5.3), (5.5)-(5.7)).
Routes: fixed-b (Lean-style: b constant, then n) and paper (b = n^{1/4}).
"""
import json, math
L10 = math.log10
C_STAR = 1 / (2 + 8 * math.exp(-1.5))
C_SIGMA = 12 / (2 * 13 ** 4)          # Fang (5.3) at lam = 12  (audit-b step 3)
RAIC1 = 58.0                          # Raic d=1: 42+16 (audit-a/b read arXiv:1802.06475v4 eq. 1.2)
RAIC2 = 42 * 2 ** 0.25 + 16           # Raic d=2 (hypothetical 2-D version for (|J|, B_J))
PHI = 1 / math.sqrt(2 * math.pi * math.e)

def c_delta(p, rho):
    C0 = 1 / (1 - rho); ystar = max(0.0, 2 / p - 1 / C0)
    return 12 * math.exp(-ystar) * (1 + C0 * ystar) ** (2 / p)

def chain(a, Cd, Cg, rho0, eta, theta=0.5, raic=RAIC1, c_sigma=C_SIGMA):
    K_M = Cd / (1 - 2 ** (a - 1))
    a_n = max(a, 0.55)   # (iv) formula needs 2a > 1; for a <= 1/2 use a_n = 0.55 (conservative; term is subdominant)
    K_S = math.sqrt(Cg / (1 - 2 ** (1 - 2 * a_n)))
    C_up = 0.25 + K_M
    rho = 1 / (104 * C_up) if rho0 == "provable" else rho0
    epsK = C_STAR * ((1 - eta) * rho) ** 1.5
    t = theta * epsK / 4
    cs = eta * c_sigma
    # fixed-b route
    lb_mean = L10(K_M / (math.pi * cs * t * t)) / (1 - a)
    lb_var = L10(2 * (2 + PHI) * K_M / (cs * t)) / (1 - a)
    lb = max(lb_mean, lb_var)
    ln_BE = 2 * (L10(raic) + lb - L10(t)) + L10(2 / cs)
    ln_var_n = L10(2 * (2 + PHI) * K_S / (cs * t)) / (1 - a_n)
    fixed = max(ln_BE, ln_var_n)
    # paper route b = n^{1/4}
    paper = dict(BE=4 * L10(raic * math.sqrt(2 / cs) / t), mean=4 * lb_mean, var_b=4 * lb_var, var_n=ln_var_n)
    # decomposition of the mean-shift log factor
    dec = dict(log10_K_M=L10(K_M), log10_1_over_c_sigma=-L10(c_sigma), log10_1_over_pi_eta=-L10(math.pi * eta),
               log10_1_over_t2=-2 * L10(t), total=L10(K_M / (math.pi * cs * t * t)))
    return dict(a=a, one_minus_a=1 - a, C_delta=Cd, K_M=K_M, K_S=K_S, rho0=rho, epsK=epsK, t=t,
                fixed_log10_b=lb, fixed_log10_N1=fixed,
                fixed_dominant="BE with b set by " + ("mean shift" if lb_mean >= lb_var else "variance shift (b)") if ln_BE >= ln_var_n else "variance shift (n)",
                paper_terms=paper, paper_log10_N1=max(paper.values()), paper_dominant=max(paper, key=paper.get),
                mean_shift_factor_decomposition=dec)

if __name__ == "__main__":
    scen = {}
    p1, r1 = 19 / 10, 199 / 200
    scen["S1: a=18/19 certified sharpened Lemma 4.2 (p=19/10, rho<=199/200), C_gamma=2.08e9 (audit-b)"] = (2 - 2 / p1, c_delta(p1, r1), 2.08e9)
    p2, r2 = 193 / 100, 49 / 50
    scen["S2: a=0.9637 (p=193/100, rho<=49/50), C_gamma=audit-b-scale 1e9 (subdominant)"] = (2 - 2 / p2, c_delta(p2, r2), 1e9)
    A = 2 * math.sqrt(12) / math.e; sP = 0.0029905; bL = 0.99427146
    pP = 2 - sP; rP = A ** (2 * sP) * bL ** (1 - 2 * sP)
    scen["P: PDF Lemma 4.2 as written, best params (audit-b s=0.0029905, b=0.99427146), C_gamma=3.0e12 (audit-b)"] = (2 - 2 / pP, c_delta(pP, rP), 3.0e12)
    scen["H88: a=0.88 (audit-a growth slope at lam=12), C_delta=C_gamma=1 (idealised)"] = (0.88, 1.0, 1.0)
    scen["H80: a=0.80 (audit-b growth slope), C_delta=C_gamma=1"] = (0.80, 1.0, 1.0)
    scen["H67: a=2/3, C=1"] = (2 / 3, 1.0, 1.0)
    scen["H00: a=0, C=1 (ideal root moments)"] = (0.0, 1.0, 1.0)
    out = dict(constants=dict(c_star=C_STAR, c_sigma=C_SIGMA, RAIC1=RAIC1, RAIC2=RAIC2, PHI=PHI,
                              C_delta_S1=c_delta(p1, r1), C_delta_S2=c_delta(p2, r2)))
    for name, (a, Cd, Cg) in scen.items():
        d = {}
        for rho0 in (0.4448, "provable"):
            for eta in (0.3, 0.5):
                d[f"rho0={rho0},eta={eta}"] = chain(a, Cd, Cg, rho0, eta)
        d["rho0=0.4448,eta=0.3,RAIC_d2"] = chain(a, Cd, Cg, 0.4448, 0.3, raic=RAIC2)
        out[name] = d
    # hypothetical lam_max = 4 (unproved; audit-a: lam=3 refuted for forests, lam=4 open), c_sigma = 4/(2*5^4),
    # a = 0.69 (audit-a growth slope at lam = 4), C = 1
    out["L4: lam_max=4 (unproved), a=0.69, C=1, c_sigma=4/(2*5^4)"] = {
        "rho0=0.4448,eta=0.3": chain(0.69, 1.0, 1.0, 0.4448, 0.3, c_sigma=4 / (2 * 5 ** 4))}
    # route-lclt's hypothetical Lemma L(K): d_K(M') <= K n^{-1/2}  ->  n >= (K/tol)^2
    hyp = {}
    for eta in (0.3, 0.5):
        epsK = C_STAR * ((1 - eta) * 0.4448) ** 1.5
        for K in (0.267, 1.0, 10.0):
            hyp[f"eta={eta},K={K}"] = dict(N_full_budget=(K / epsK) ** 2, N_half_budget=(K / (epsK / 2)) ** 2)
    out["hypothetical_LemmaL_K_over_sqrt_n"] = hyp
    # Fang N0 reference values (audit-a / audit-b, log10)
    out["fang_N0_reference_log10"] = {"Lean witnesses (both audits)": 1.68e7, "PDF best as written (audit-b)": 2.44e5,
        "sharpened Lemma 4.2 p=19/10, fixed-b (audit-b)": 2069.69, "sharpened, paper route (audit-b)": 5375.03,
        "architecture floor a=0,C=1, lam_max=12, fixed-b (audit-a)": 103}
    json.dump(out, open("rerun/n1_check.json", "w"), indent=1)
    for name, d in out.items():
        if not isinstance(d, dict) or name in ("constants", "fang_N0_reference_log10", "hypothetical_LemmaL_K_over_sqrt_n"):
            print(name, json.dumps(d)); continue
        print("==", name)
        for k, r in d.items():
            print(f"  {k}: 1-a={r['one_minus_a']:.4f} K_M={r['K_M']:.4g} rho0={r['rho0']:.3g} epsK={r['epsK']:.3g} t={r['t']:.3g}"
                  f" | fixed-b: log10 b={r['fixed_log10_b']:.1f}, log10 N1={r['fixed_log10_N1']:.1f} ({r['fixed_dominant']})"
                  f" | paper: {r['paper_log10_N1']:.1f} ({r['paper_dominant']})"
                  f" | meanshift factor {r['mean_shift_factor_decomposition']['total']:.2f} = KM {r['mean_shift_factor_decomposition']['log10_K_M']:.2f}"
                  f" + 1/cs {r['mean_shift_factor_decomposition']['log10_1_over_c_sigma']:.2f} + 1/(pi eta) {r['mean_shift_factor_decomposition']['log10_1_over_pi_eta']:.2f}"
                  f" + 1/t^2 {r['mean_shift_factor_decomposition']['log10_1_over_t2']:.2f}")
