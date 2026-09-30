"""Heterogeneous star-of-hubs: centre z joined to c_j hubs of type (m_j, s_j).
Exact closed forms; used for a capped search for R1 counterexamples smaller than n = 249.
Only the top 16 window levels are scanned in the search (all SH failures seen sit there);
any hit is then rescanned on the full window and rechecked with the generic DP."""
import sys, json, itertools, time
from flint import fmpz_poly
import fam
X = fam.X

def mixed(types):
    """types: list of (count, m, s)"""
    A = {}; B = {}; U = {}; Up = {}
    for c, m, s in types:
        a, b = fam.AB(s); A[(m, s)] = a; B[(m, s)] = b
        U[(m, s)] = a ** m + X * b ** m; Up[(m, s)] = a ** (m - 1) + X * b ** (m - 1)
    def prodU(skip=None):
        p = fmpz_poly([1])
        for c, m, s in types:
            cc = c - (1 if (m, s) == skip else 0)
            p *= U[(m, s)] ** cc
        return p
    def prodA(skip=None):
        p = fmpz_poly([1])
        for c, m, s in types:
            cc = c - (1 if (m, s) == skip else 0)
            p *= A[(m, s)] ** (m * cc)
        return p
    P, Q = prodU(), prodA()
    I = P + X * Q
    h = sum(c for c, _, _ in types)
    n = 1 + sum(c * (1 + m * (s + 1)) for c, m, s in types)
    cl = {"z": (1, h, Q)}; nb = {"z": {"z": 1}}
    for c, m, s in types:
        t = f"{m}_{s}"
        Pj, Qj = prodU((m, s)), prodA((m, s))
        cl["hub" + t] = (c, m + 1, B[(m, s)] ** m * Pj)
        cl["mid" + t] = (c * m, s + 1, A[(m, s)] ** (m - 1) * (Pj + X * Qj))
        cl["leaf" + t] = (c * m * s, 1, (1 + X) ** (s - 1) * (Up[(m, s)] * Pj + X * A[(m, s)] ** (m - 1) * Qj))
        nb["z"]["hub" + t] = c
        nb["hub" + t] = {"hub" + t: 1, "z": 1, "mid" + t: m}
        nb["mid" + t] = {"mid" + t: 1, "hub" + t: 1, "leaf" + t: s}
        nb["leaf" + t] = {"leaf" + t: 1, "mid" + t: 1}
    return fam.Fam("MSH" + str(types), n, I, cl, nb)

if __name__ == "__main__":
    NMAX = int(sys.argv[1]); TLIM = float(sys.argv[2])
    t0 = time.time(); cands = []
    ms = range(7, 17)
    for h in range(3, 10):
        for ma, mb in itertools.combinations_with_replacement(ms, 2):
            for ca in range(1, h + 1):
                cb = h - ca
                types = [(ca, ma, 2)] + ([(cb, mb, 2)] if cb and mb != ma else [])
                if cb and mb == ma: continue
                n = 1 + sum(c * (3 * m + 1) for c, m, s in types)
                if n < NMAX: cands.append(types)
    for c, m, s in [(h, m, 3) for h in range(3, 10) for m in range(7, 17)]:
        if 1 + c * (1 + m * 4) < NMAX: cands.append([(c, m, s)])
    cands.sort(key=lambda T: 1 + sum(c * (1 + m * (s + 1)) for c, m, s in T))
    print(json.dumps({"n_candidates": len(cands)}), flush=True)
    done = 0
    for T in cands:
        if time.time() - t0 > TLIM: break
        F = mixed(T)
        rows = F.scan(range(max(F.lo, F.q - 15), F.q + 1))
        done += 1
        fails = [(r["k"], r["r1_fail"]) for r in rows if r["r1_fail"]]
        if fails:
            print(json.dumps({"hit": str(T), "n": F.n, "fails": fails}), flush=True)
    print(json.dumps({"evaluated": done, "of": len(cands), "seconds": round(time.time() - t0)}), flush=True)
