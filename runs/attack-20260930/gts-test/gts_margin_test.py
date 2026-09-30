"""Does Csikvari's generalized tree shift (GTS) move the minimum central
log-concavity margin monotonically toward the star?

GTS (Csikvari, 'On a poset of trees', Def. 2.1, p. 2): x, y vertices of tree T
with every interior vertex of the x-y path of degree 2; z = y's neighbour on the
path; move every edge y-w (w in N(y) - {z}) to x-w. Proper iff neither x nor y
is a leaf; then T' has one more leaf. Star = unique maximum, path = unique
minimum of the induced poset (Thm 2.4).

m(T) = min over k in W(T) of delta_k = 1 - i_{k-1} i_{k+1} / i_k^2 (exact),
W(T) = [ceil(n/4), q], q = ceil((2 alpha - 1)/3) (also reported for ceil(n/5)).
Strong form S: m(T') <= m(T) for every proper GTS T -> T'.
Weak form W: every non-star T has SOME proper GTS image T' with m(T') <= m(T).
W implies m(T) >= m(star) > 0 for all T by chaining up the poset, hence
log-concavity on W(T) and unimodality.
"""
import sys, json, networkx as nx
from fractions import Fraction
from math import ceil

def indpoly(T):
    root = next(iter(T)); order = list(nx.dfs_preorder_nodes(T, root)); par = {root: None}
    for u in order:
        for w in T[u]:
            if w != par[u]: par[w] = u
    ex, ic = {}, {}
    def mul(p, q):
        r = [0] * (len(p) + len(q) - 1)
        for i, a in enumerate(p):
            if a:
                for j, b in enumerate(q): r[i + j] += a * b
        return r
    for v in reversed(order):
        e, c = [1], [0, 1]
        for w in T[v]:
            if w != par[v]:
                s = [(ex[w][i] if i < len(ex[w]) else 0) + (ic[w][i] if i < len(ic[w]) else 0) for i in range(max(len(ex[w]), len(ic[w])))]
                e = mul(e, s); c = mul(c, ex[w])
        ex[v], ic[v] = e, c
    s = [(ex[root][i] if i < len(ex[root]) else 0) + (ic[root][i] if i < len(ic[root]) else 0) for i in range(max(len(ex[root]), len(ic[root])))]
    while s[-1] == 0: s.pop()
    return s

def margin(T, lowdiv):
    I = indpoly(T); n = T.number_of_nodes(); a = len(I) - 1
    q = ceil((2 * a - 1) / 3); lo = max(1, ceil(n / lowdiv))
    best = None
    for k in range(lo, min(q, a - 1) + 1):
        d = 1 - Fraction(I[k - 1] * I[k + 1], I[k] ** 2)
        if best is None or d < best: best = d
    return best

def proper_gts_images(T):
    leaves = {v for v in T if T.degree(v) == 1}
    for x in T:
        if x in leaves: continue
        for y in T:
            if y == x or y in leaves: continue
            path = nx.shortest_path(T, x, y)
            if any(T.degree(v) != 2 for v in path[1:-1]): continue
            z = path[-2]
            U = T.copy()
            for w in list(T[y]):
                if w != z:
                    U.remove_edge(y, w); U.add_edge(x, w)
            yield (x, y), U

N = int(sys.argv[1]); out = {}
for lowdiv in (4, 5):
    res = {}
    for n in range(4, N + 1):
        trees = list(nx.nonisomorphic_trees(n))
        star_m = None; s_viol = 0; w_fail = []; shifts = 0; worst = None
        for T in trees:
            if max(d for _, d in T.degree()) == n - 1:
                star_m = margin(T, lowdiv); continue
            mT = margin(T, lowdiv); ok_w = False
            for xy, U in proper_gts_images(T):
                assert nx.is_tree(U) and sum(1 for v in U if U.degree(v) == 1) == sum(1 for v in T if T.degree(v) == 1) + 1
                shifts += 1; mU = margin(U, lowdiv)
                if mT is None or mU is None: continue
                if mU <= mT: ok_w = True
                else:
                    s_viol += 1
                    r = mU / mT
                    if worst is None or r > worst[0]: worst = (r, nx.to_graph6_bytes(T, header=False).decode().strip(), xy)
            if not ok_w and mT is not None: w_fail.append(nx.to_graph6_bytes(T, header=False).decode().strip())
        res[n] = dict(trees=len(trees), proper_shifts=shifts, strong_violations=s_viol, weak_failures=len(w_fail),
                      weak_failure_examples=w_fail[:5], star_margin=str(star_m),
                      worst_strong_ratio=(str(worst[0]), worst[1], str(worst[2])) if worst else None)
        print(f'W=[n/{lowdiv},q] n={n} trees={len(trees)} shifts={shifts} S-viol={s_viol} W-fail={len(w_fail)} star={star_m} worstS={float(worst[0]) if worst else None}', flush=True)
    out[f'n/{lowdiv}'] = res
json.dump(out, open(f'gts_margin_n4_{N}.json', 'w'), indent=1)
