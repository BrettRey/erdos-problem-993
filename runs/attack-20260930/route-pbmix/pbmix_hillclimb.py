#!/usr/bin/env python3
"""Adversarial hill-climb: maximise max_{k in [ceil(n/4), q]} rho_k(T, B=branch set).
rho is computed exactly (Fraction) and compared exactly; the objective is recorded
as float only for logging.  Mutation: move a leaf to a new parent."""
import json
import random
import sys
import time
from fractions import Fraction as Fr

import pbmix_killtest as P


def score(adj, n):
    B = {v for v in range(n) if len(adj[v]) >= 3}
    if len(B) > 10:
        return None, None
    comps = P.components_branch_set(adj, n, B)
    f = [0] * (n + 1)
    for mult, s, Q in comps:
        for i, c in enumerate(Q):
            f[s + i] += mult * c
    while len(f) > 1 and f[-1] == 0:
        f.pop()
    alpha = P.alpha_of(f)
    lo, q, _, _ = P.windows(n, alpha)
    best, bestk, info = Fr(-10), None, None
    for k in range(lo, q + 1):
        st = P.mixture_stats(comps, k, want_cert=True)
        if st["rho"] is not None and st["rho"] > best:
            best, bestk = st["rho"], k
            info = {"cert_pb_pass": st["cert_pb"] > 0, "cert_best_pass": st["cert_best"] > 0,
                    "LC": st["delta_f"] > 0}
    return best, (bestk, len(B), info)


def mutate(adj, n, rng):
    adj = [list(a) for a in adj]
    leaves = [v for v in range(n) if len(adj[v]) == 1]
    v = rng.choice(leaves)
    p = adj[v][0]
    adj[p].remove(v)
    adj[v] = []
    cands = [w for w in range(n) if w != v and (len(adj[w]) > 0 or n == 1)]
    w = rng.choice(cands)
    adj[w].append(v)
    adj[v].append(w)
    return adj


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 22
    budget = float(sys.argv[2]) if len(sys.argv) > 2 else 240
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    rng = random.Random(seed)
    t0 = time.time()
    global_best = (Fr(-10), None, None)
    restarts = 0
    while time.time() - t0 < budget:
        restarts += 1
        while True:
            adj = P.random_tree(n, rng)
            s, meta = score(adj, n)
            if s is not None:
                break
        stall = 0
        while stall < 150 and time.time() - t0 < budget:
            adj2 = mutate(adj, n, rng)
            s2, meta2 = score(adj2, n)
            if s2 is None:
                stall += 1
                continue
            if s2 >= s:
                if s2 > s:
                    stall = 0
                else:
                    stall += 1
                adj, s, meta = adj2, s2, meta2
            else:
                stall += 1
        if s > global_best[0]:
            edges = sorted({(min(a, b), max(a, b)) for a in range(n) for b in adj[a]})
            global_best = (s, meta, edges)
            print(json.dumps({"n": n, "restart": restarts, "rho_max_float": float(s),
                              "k": meta[0], "B": meta[1], "info": meta[2],
                              "degseq": sorted((len(a) for a in adj), reverse=True)[:8]}),
                  flush=True)
    s, meta, edges = global_best
    json.dump({"n": n, "seed": seed, "restarts": restarts, "rho_max": str(s),
               "rho_max_float": float(s), "k": meta[0], "B": meta[1], "info": meta[2],
               "edges": edges}, open(f"results_hillclimb_n{n}_s{seed}.json", "w"), indent=1)


if __name__ == "__main__":
    main()
