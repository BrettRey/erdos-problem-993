"""Generic forest independence polynomial (copied from ../r1_third_check_20260930.py)."""
def pmul(p, q):
    r = [0] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        if x:
            for j, y in enumerate(q): r[i + j] += x * y
    return r
def padd(p, q):
    L = max(len(p), len(q)); return [(p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0) for i in range(L)]
def indpoly(adj, alive):
    alive = set(alive); seen = set(); total = [1]
    for r in alive:
        if r in seen: continue
        order, parent, stack = [], {r: None}, [r]; seen.add(r)
        while stack:
            v = stack.pop(); order.append(v)
            for w in adj[v]:
                if w in alive and w not in seen: seen.add(w); parent[w] = v; stack.append(w)
        out, inn = {}, {}
        for v in reversed(order):
            ex, ic = [1], [0, 1]
            for w in adj[v]:
                if w in alive and parent.get(w) == v: ex = pmul(ex, padd(out[w], inn[w])); ic = pmul(ic, out[w])
            out[v], inn[v] = ex, ic
        total = pmul(total, padd(out[r], inn[r]))
    while len(total) > 1 and total[-1] == 0: total.pop()
    return total
