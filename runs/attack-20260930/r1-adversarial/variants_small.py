"""Structured variants of the star of hubs to push the smallest R1 violation below n=249.
Centre (vertex 0) joined to hubs; hub i carries m2 2-cherries, m3 3-cherries, m1 pendant P2s and e leaves.
Nothing is added to the centre except in the 'cl' runs (centre leaves), recorded separately.
Usage: python variants_small.py NMAX > data/variants_small.jsonl"""
import itertools, json, sys
from r1lib import B, evaluate

def build(hubs, centre_leaves=0):
    b = B()
    for (m2, m3, m1, e) in hubs:
        w = b.add(0)
        for t, cnt in ((2, m2), (3, m3), (1, m1)):
            for _ in range(cnt):
                c = b.add(w)
                for _ in range(t):
                    b.add(c)
        for _ in range(e):
            b.add(w)
    for _ in range(centre_leaves):
        b.add(0)
    return b.adj()

def size(hubs, cl=0):
    return 1 + cl + sum(1 + 3 * m2 + 4 * m3 + 2 * m1 + e for (m2, m3, m1, e) in hubs)

def main():
    nmax = int(sys.argv[1])
    specs = []
    for h in (5, 6, 7, 8, 9):
        for m2 in range(6, 13):
            for m3 in range(0, 4):
                for m1 in range(0, 4):
                    for e in range(0, 4):
                        specs.append((f"uni h{h} m2{m2} m3{m3} m1{m1} e{e}", [(m2, m3, m1, e)] * h, 0))
        for m in range(7, 12):
            for a in range(1, h):
                specs.append((f"mix h{h} {a}x{m}+{h-a}x{m+1}", [(m, 0, 0, 0)] * a + [(m + 1, 0, 0, 0)] * (h - a), 0))
                specs.append((f"mix h{h} {a}x{m}+{h-a}x{m+2}", [(m, 0, 0, 0)] * a + [(m + 2, 0, 0, 0)] * (h - a), 0))
        for m in range(8, 11):
            for cl in (1, 2, 4):
                specs.append((f"uni h{h} m2{m} cl{cl}", [(m, 0, 0, 0)] * h, cl))
    for name, hubs, cl in specs:
        n = size(hubs, cl)
        if n > nmax or n < 150:
            continue
        r = evaluate(build(hubs, cl))
        assert r["n"] == n
        print(json.dumps({"name": name, "hubs": hubs, "cl": cl, "n": n, "alpha": r["alpha"], "window": r["window"],
                          "n_viol": len(r["viol"]), "viol": r["viol"][:3], "viol_k": sorted({v["k"] for v in r["viol"]}),
                          "lc_fail": r["lc_fail"], "F": r["F"], "F_at": r["F_at"], "C": r["C"]}), flush=True)


if __name__ == '__main__':
    main()
