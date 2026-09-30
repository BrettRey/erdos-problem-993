"""Try 1-3 extra pendant leaves at representative positions of near-violating stars of hubs (n <= 236),
exact R1 verdict via r1lib. Usage: python leaf_tweaks.py > data/leaf_tweaks.jsonl"""
import itertools, json
from r1lib import adj_from_edges, evaluate
rows = [json.loads(l) for l in open('data/multiset_centre_scan.jsonl')]
near = sorted([r for r in rows if not r['viol_k'] and r['n'] <= 235], key=lambda r: -r['F_centre'])[:12]
def build(ms):
    edges = []; nxt = 1; pos = {"hub": [], "cc": [], "leaf": []}
    for m in ms:
        w = nxt; nxt += 1; edges.append((0, w)); pos["hub"].append(w)
        for _ in range(m):
            c = nxt; nxt += 1; edges.append((w, c)); pos["cc"].append(c)
            for _ in range(2):
                edges.append((c, nxt)); pos["leaf"].append(nxt); nxt += 1
    return edges, nxt, pos
for r in near:
    ms = r['ms']
    base, n0, pos = build(ms)
    # representative positions: first cherry centre / leaf of the first hub of each distinct size; hubs of each size
    reps = []
    seen = set()
    hub_of = {}
    idx = 0
    for hi, m in enumerate(ms):
        if m in seen: 
            idx += m; continue
        seen.add(m)
        reps += [("hub", pos["hub"][hi]), ("cc", pos["cc"][idx]), ("leaf", pos["leaf"][2 * idx])]
        idx += m
    reps.append(("centre", 0))
    for extra in (1, 2, 3):
        for combo in itertools.combinations_with_replacement(range(len(reps)), extra):
            if n0 + extra > 236: continue
            edges = list(base); nxt = n0
            for c in combo:
                edges.append((reps[c][1], nxt)); nxt += 1
            ev = evaluate(adj_from_edges(nxt, edges))
            print(json.dumps({"ms": ms, "added_at": [reps[c][0] + str(reps[c][1]) for c in combo], "n": nxt,
                              "window": ev["window"], "n_viol": len(ev["viol"]), "viol": ev["viol"][:2], "F": ev["F"], "F_at": ev["F_at"],
                              "edges": edges if ev["viol"] else None}), flush=True)
