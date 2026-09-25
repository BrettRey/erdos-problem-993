"""Exhaustive Mason-(b) check driven by gentreeg -p parent arrays.

usage:  gentreeg -p N res/mod | python mason_b_census_gentreeg_20260902.py N res

Mason (b) at index k:  (k+2) i_{k+2} i_k <= (k+1) i_{k+1}^2,
equivalently mu_{k+1} <= mu_k for mu_k = (k+1) i_{k+1}/i_k.
Prints every failure (parent array, polynomial, failing k) and a DONE line
with the tree count.  Exact integer arithmetic throughout.
"""
import sys


def ipoly(par):
    """par: gentreeg -p line, entry v-1 is the parent of vertex v (0 = root)."""
    n = len(par)
    ch = [[] for _ in range(n + 1)]
    root = None
    for v in range(1, n + 1):
        if par[v - 1] == 0:
            root = v
        else:
            ch[par[v - 1]].append(v)
    order = []
    st = [root]
    while st:
        u = st.pop()
        order.append(u)
        st.extend(ch[u])
    P = [None] * (n + 1)
    R = [None] * (n + 1)
    for u in reversed(order):
        p = [0, 1]
        r = [1]
        for c in ch[u]:
            Rc = R[c]
            Pc = P[c]
            out = [0] * (len(p) + len(Rc) - 1)
            for i, x in enumerate(p):
                if x:
                    for j, y in enumerate(Rc):
                        out[i + j] += x * y
            p = out
            s = [(Pc[i] if i < len(Pc) else 0) + (Rc[i] if i < len(Rc) else 0)
                 for i in range(max(len(Pc), len(Rc)))]
            out = [0] * (len(r) + len(s) - 1)
            for i, x in enumerate(r):
                if x:
                    for j, y in enumerate(s):
                        out[i + j] += x * y
            r = out
        P[u] = p
        R[u] = r
    p, r = P[root], R[root]
    return [(p[i] if i < len(p) else 0) + (r[i] if i < len(r) else 0)
            for i in range(max(len(p), len(r)))]


def main():
    n = int(sys.argv[1])
    res = sys.argv[2]
    cnt = 0
    fails = 0
    for line in sys.stdin:
        par = [int(x) for x in line.split()]
        if len(par) != n:
            continue
        cnt += 1
        p = ipoly(par)
        a = len(p) - 1
        bad = [k for k in range(a - 1)
               if (k + 2) * p[k + 2] * p[k] > (k + 1) * p[k + 1] * p[k + 1]]
        if bad:
            fails += 1
            print("FAIL", n, res, par, p, bad, flush=True)
    print("DONE", n, res, "trees", cnt, "mason_b_failures", fails, flush=True)


if __name__ == "__main__":
    main()
