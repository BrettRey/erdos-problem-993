#!/usr/bin/env python3
"""B = one bipartition class (T - B independent => binomial components with EXACT
curvature; this is the decomposition behind Fang et al. Lemma 6.1).  Exact rho."""
import json, random
import pbmix_killtest as P
import pbmix_families2 as F

def classes(adj, n):
    col = [-1]*n; col[0] = 0; st = [0]
    while st:
        v = st.pop()
        for w in adj[v]:
            if col[w] < 0:
                col[w] = 1-col[v]; st.append(w)
    return [{v for v in range(n) if col[v]==c} for c in (0,1)]

out = []
def run(name, adj, n):
    for ci, B in enumerate(classes(adj, n)):
        r = P.analyse_binomial_family(f"{name} B=class{ci}", adj, n, B, cert=False)
        rh = [x["rho_float"] for x in r["rows"]]; bf = [x["between_frac_float"] for x in r["rows"]]
        print(r["family"], "n", n, "|B|", len(B), "window", r["window"], "allLC", all(x["LC_exact"] for x in r["rows"]),
              "rho[min,max]", (round(min(rh),4), round(max(rh),4)), "bf[min,max]", (round(min(bf),4), round(max(bf),4)), flush=True)
        out.append(r)

for d in range(3, 7):
    adj, n = P.complete_binary(d); run(f"CBT d={d}", adj, n)
rng = random.Random(4242)
for n in (40, 80, 120):
    adj = P.random_tree(n, rng); run(f"rand n={n}", adj, n)
for L in (40, 80):
    for e in (1, 2):
        adj, n = F.comb_tree(L, e); run(f"comb L={L} every={e}", adj, n)
json.dump(out, open("results_bipartition.json", "w"), indent=1)
