"""R1 family sweep: structured adversarial families, 30 <= n <= NMAX.
Exact R1 verdict per tree (integer test), float normalised loads A,B,C,E,F.
Usage: python families_r1.py PART NPARTS NMAX > data/families_partX.jsonl
"""
import json
import sys
import time
from r1lib import B, evaluate


def hubstar_into(b, root, m, t):
    for _ in range(m):
        c = b.add(root)
        for _ in range(t):
            b.add(c)


def H(m, t, extra_leaves=0, tail=0, tail_at_leaf=0):
    b = B()
    hubstar_into(b, 0, m, t)
    for _ in range(extra_leaves):
        b.add(0)
    if tail:
        b.path(0, tail)
    if tail_at_leaf:
        b.path(b.n - 1 if not tail else 2, tail_at_leaf)
    return b.adj()


def hub_chain(h, m, t, bridge):
    """h hubs, consecutive hubs joined by a path with `bridge` edges (1 = adjacent)"""
    b = B()
    hub = 0
    for i in range(h):
        hubstar_into(b, hub, m, t)
        if i < h - 1:
            hub = b.path(hub, bridge)
    return b.adj()


def star_of_hubs(h, m, t, sub=1, centre_leaves=0):
    """centre joined by paths of `sub` edges to h hubs, each carrying m t-cherries"""
    b = B()
    for _ in range(h):
        hub = b.path(0, sub)
        hubstar_into(b, hub, m, t)
    for _ in range(centre_leaves):
        b.add(0)
    return b.adj()


def hub_of_hubs(m, s, t):
    b = B()
    for _ in range(m):
        c = b.add(0)
        for _ in range(s):
            g = b.add(c)
            for _ in range(t):
                b.add(g)
    return b.adj()


def hub_with_hub_neighbours(h, m, t, own_cherries):
    """centre hub adjacent to h hubs (each with m t-cherries) and own_cherries t-cherries"""
    b = B()
    for _ in range(h):
        w = b.add(0)
        hubstar_into(b, w, m, t)
    hubstar_into(b, 0, own_cherries, t)
    return b.adj()


def caterpillar(L, a):
    b = B()
    spine = [0] + [None] * (L - 1)
    for i in range(1, L):
        spine[i] = b.add(spine[i - 1])
    for s in spine:
        for _ in range(a):
            b.add(s)
    return b.adj()


def cherry_caterpillar(L, m, t):
    """spine path of L hubs, each carrying m t-cherries"""
    return hub_chain(L, m, t, 1)


def double_broom(L, a, c):
    b = B()
    end = b.path(0, L)
    for _ in range(a):
        b.add(0)
    for _ in range(c):
        b.add(end)
    return b.adj()


def star_of_stars(m, s, centre_leaves=0):
    b = B()
    for _ in range(m):
        c = b.add(0)
        for _ in range(s):
            b.add(c)
    for _ in range(centre_leaves):
        b.add(0)
    return b.adj()


def spider(legs):
    b = B()
    for L in legs:
        b.path(0, L)
    return b.adj()


def bundle_child(b, parent, arms):
    c = b.add(parent)
    for _ in range(arms):
        b.path(c, 2)
    return c


def galvin_T(m, tt):
    b = B()
    for _ in range(m):
        bundle_child(b, 0, tt)
    return b.adj()


def kl_T3(m, n2):
    b = B()
    for arms in (3, m, n2):
        bundle_child(b, 0, arms)
    return b.adj()


def TG(m, tt):
    b = B()
    for _ in range(m):
        v = b.add(0)
        for _ in range(3):
            bundle_child(b, v, tt)
    b.add(0)
    return b.adj()


def hubs_via_deg2(h, m, t):
    """hub u with h pendant paths u-p-w, each w a hub with m t-cherries; u also carries m t-cherries"""
    b = B()
    hubstar_into(b, 0, m, t)
    for _ in range(h):
        p = b.add(0)
        w = b.add(p)
        hubstar_into(b, w, m, t)
    return b.adj()


def specs(nmax):
    S = []
    for t in range(1, 13):
        for m in range(3, 80):
            n = 1 + m * (t + 1)
            if 30 <= n <= nmax:
                S.append((f"H({m},{t})", ("H", (m, t))))
    for m in range(6, 40, 3):
        for t in (2, 3, 4, 5):
            for e in (1, 3, 8, 15):
                S.append((f"H({m},{t})+{e}leaves", ("H", (m, t, e))))
            for tail in (1, 2, 4, 8, 16, 30):
                S.append((f"H({m},{t})+tail{tail}", ("H", (m, t, 0, tail))))
                S.append((f"H({m},{t})+leaftail{tail}", ("H", (m, t, 0, 0, tail))))
    for h in (2, 3, 4, 5, 6):
        for m in (3, 5, 7, 9, 12, 16, 20):
            for t in (1, 2, 3, 5):
                for br in (1, 2, 3, 4):
                    S.append((f"hubchain(h{h},m{m},t{t},br{br})", ("hub_chain", (h, m, t, br))))
                for sub in (1, 2, 3):
                    for cl in (0, 5):
                        S.append((f"starofhubs(h{h},m{m},t{t},sub{sub},cl{cl})", ("star_of_hubs", (h, m, t, sub, cl))))
                for own in (0, 3, 9):
                    S.append((f"hubhubnbrs(h{h},m{m},t{t},own{own})", ("hub_with_hub_neighbours", (h, m, t, own))))
                S.append((f"hubsviadeg2(h{h},m{m},t{t})", ("hubs_via_deg2", (h, m, t))))
    for m in range(3, 16):
        for s in range(1, 6):
            for t in (1, 2, 3):
                S.append((f"hubofhubs({m},{s},{t})", ("hub_of_hubs", (m, s, t))))
    for L in range(4, 80, 3):
        for a in (1, 2, 3, 4, 6):
            S.append((f"caterpillar({L},{a})", ("caterpillar", (L, a))))
    for L in range(2, 60, 4):
        for a in (2, 5, 10, 20, 40):
            for c in (a, 1, 3):
                S.append((f"doublebroom({L},{a},{c})", ("double_broom", (L, a, c))))
    for m in range(3, 30, 2):
        for s in (2, 3, 4, 6, 9):
            for cl in (0, 5, 15, 30):
                S.append((f"starofstars({m},{s},cl{cl})", ("star_of_stars", (m, s, cl))))
    for m in range(3, 40, 2):
        for legs in ([2] * m, [3] * m, [4] * m, [1] * m + [2] * m, [2] * m + [5], [1] * (2 * m) + [2] * m):
            S.append((f"spider{legs[:3]}..x{len(legs)}", ("spider", (legs,))))
    for m in range(2, 30, 2):
        for tt in (2, 3, 4, 5, 7, 10):
            S.append((f"galvinT({m},{tt})", ("galvin_T", (m, tt))))
            S.append((f"TG({m},{tt})", ("TG", (m, tt))))
    for m in range(3, 70, 3):
        S.append((f"klT3({m},{m})", ("kl_T3", (m, m))))
    return S


FUN = {"H": H, "hub_chain": hub_chain, "star_of_hubs": star_of_hubs, "hub_of_hubs": hub_of_hubs,
       "hub_with_hub_neighbours": hub_with_hub_neighbours, "caterpillar": caterpillar,
       "double_broom": double_broom, "star_of_stars": star_of_stars, "spider": spider,
       "galvin_T": galvin_T, "kl_T3": kl_T3, "TG": TG, "hubs_via_deg2": hubs_via_deg2}


def main():
    part, nparts, nmax = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    S = specs(nmax)
    for i, (name, (fn, args)) in enumerate(S):
        if i % nparts != part:
            continue
        adj = FUN[fn](*args)
        n = len(adj)
        if n < 30 or n > nmax:
            continue
        t0 = time.time()
        r = evaluate(adj)
        r["family"] = name
        r["ctor"] = [fn, list(args) if fn != "spider" else [list(args[0])]]
        r["n_viol"] = len(r["viol"])
        r["viol"] = r["viol"][:5]
        r["secs"] = round(time.time() - t0, 3)
        print(json.dumps(r), flush=True)


if __name__ == "__main__":
    main()
