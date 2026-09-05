"""Exploratory audit (NOT a verification): brute-force check of the fivefold
claim on small forests.  Results here are heuristic evidence only.
"""
import sys

def canon(adj, n):
    def bfs(s):
        dist = [-1]*n; dist[s]=0; order=[s]
        for u in order:
            for v in adj[u]:
                if dist[v] < 0:
                    dist[v] = dist[u]+1; order.append(v)
        return dist
    d0 = bfs(0)
    a = max(range(n), key=lambda i: d0[i])
    da = bfs(a)
    b = max(range(n), key=lambda i: da[i])
    db = bfs(b)
    path_len = da[b]
    centers = [v for v in range(n) if da[v]+db[v]==path_len and min(da[v],db[v])==path_len//2]
    def enc(root):
        parent = [-1]*n; order=[root]; parent[root]=root
        for u in order:
            for v in adj[u]:
                if parent[v]==-1:
                    parent[v]=u; order.append(v)
        lab = ['']*n
        for u in reversed(order):
            ch = sorted(lab[v] for v in adj[u] if parent[v]==u and v!=root)
            lab[u] = '('+''.join(ch)+')'
        return lab[root]
    return min(enc(c) for c in centers)

def gen_trees(maxn):
    all_t = {1: [tuple()]}
    for n in range(2, maxn+1):
        seen = {}
        for edges in all_t[n-1]:
            adj = [[] for _ in range(n)]
            for u,v in edges:
                adj[u].append(v); adj[v].append(u)
            for u in range(n-1):
                e2 = edges + ((u, n-1),)
                adj[u].append(n-1); adj[n-1].append(u)
                c = canon(adj, n)
                adj[u].pop(); adj[n-1].pop()
                if c not in seen:
                    seen[c] = e2
        all_t[n] = list(seen.values())
    return all_t

def profile(n, edges):
    nbr = [0]*n
    for u,v in edges:
        nbr[u] |= 1<<v; nbr[v] |= 1<<u
    indep = []
    for m in range(1<<n):
        ok = True; mm = m
        while mm:
            bb = mm & -mm; i = bb.bit_length()-1; mm ^= bb
            if nbr[i] & m: ok=False; break
        if ok: indep.append(m)
    alpha = max(bin(m).count('1') for m in indep)
    maxsets = [m for m in indep if bin(m).count('1')==alpha]
    e = [0]*(alpha+2); b_ = [0]*(alpha+2)
    for m in indep:
        d = alpha - bin(m).count('1')
        if any(m & M == m for M in maxsets): e[d]+=1
        else: b_[d]+=1
    return alpha, e, b_

def combine(f):
    off = 0; E=[]
    for (k, edges) in f:
        for u,v in edges: E.append((u+off, v+off))
        off += k
    return off, E

if __name__ == '__main__':
    maxn = int(sys.argv[1]) if len(sys.argv)>1 else 11
    T = gen_trees(maxn)
    for n in sorted(T): print('trees order', n, len(T[n]))
    items = [(k,e) for k in range(1,maxn+1) for e in T[k]]
    def rec2(start, remaining, acc):
        if acc: yield list(acc)
        for idx in range(start, len(items)):
            k,e = items[idx]
            if k <= remaining:
                acc.append(items[idx])
                yield from rec2(idx, remaining-k, acc)
                acc.pop()
    bad = 0; cnt=0
    for f in rec2(0, maxn, []):
        n, E = combine(f)
        alpha, e, b_ = profile(n, E)
        if alpha < 5 or b_[2] != 0: continue
        cnt += 1
        if not (5*b_[3] <= e[3] and 5*b_[4] <= e[4]):
            bad += 1
            print('FAIL', n, E, alpha, e[:6], b_[:6])
    print('eligible forests checked:', cnt, 'failures:', bad)
