"""Structured family sweep for the constant-sigma rule S(sigma).
Usage: python families.py PART NPARTS NMAX OUT.jsonl
Each line: family, params, n, exact interval ends (num/den), binding (u,k,deg),
float ranking fields, S(2/3) exact violation count."""
import sys, json, time
from slib import analyse, B, edges_of

def hubstar(b, p, m, s):
    h = b.add(p)
    for _ in range(m):
        x = b.add(h)
        for _ in range(s):
            b.add(x)
    return h

def t_msh(hubs, svec, cl=0, cp=0, sub=1):
    b = B()
    for m, s in zip(hubs, svec):
        p = b.path(0, sub - 1) if sub > 1 else 0
        hubstar(b, p, m, s)
    for _ in range(cl): b.add(0)
    if cp: b.path(0, cp)
    return b.adj()

def t_nest(h, g, m, s, own=0):
    """centre -> h hubs -> each hub g sub-hubs H(m,s) (+ own supports with s leaves)"""
    b = B()
    for _ in range(h):
        hub = b.add(0)
        for _ in range(g): hubstar(b, hub, m, s)
        for _ in range(own):
            x = b.add(hub)
            for _ in range(s): b.add(x)
    return b.adj()

def t_twocentre(h1, h2, m, s, L):
    b = B()
    for _ in range(h1): hubstar(b, 0, m, s)
    c2 = b.path(0, L)
    for _ in range(h2): hubstar(b, c2, m, s)
    return b.adj()

def t_catspine(L, c, m, s):
    b = B(); spine = [0]
    for i in range(1, L): spine.append(b.add(spine[-1]))
    for v in spine:
        for _ in range(c): hubstar(b, v, m, s)
    return b.adj()

def t_H(m, s, cl=0):
    b = B()
    for _ in range(m):
        x = b.add(0)
        for _ in range(s): b.add(x)
    for _ in range(cl): b.add(0)
    return b.adj()

def specs(nmax):
    S = []
    def n_msh(hubs, svec, cl=0, cp=0, sub=1): return 1 + sum(1 + m * (s + 1) + (sub - 1) for m, s in zip(hubs, svec)) + cl + cp
    # A: uniform MSH
    for s in (1, 2, 3):
        for m in range(3, 26):
            for h in range(1, 41):
                if 30 <= n_msh([m] * h, [s] * h) <= nmax:
                    S.append(('MSH', dict(h=h, m=m, s=s), ('msh', ([m] * h, [s] * h))))
    # B: two hub sizes
    for s in (1, 2, 3):
        for m1 in range(3, 26, 2):
            for m2 in range(m1 + 1, 26, 3):
                for c1 in range(1, 9):
                    for c2 in range(1, 9):
                        hubs = [m1] * c1 + [m2] * c2
                        if 30 <= n_msh(hubs, [s] * len(hubs)) <= nmax:
                            S.append(('MSHmix', dict(m1=m1, c1=c1, m2=m2, c2=c2, s=s), ('msh', (hubs, [s] * len(hubs)))))
    # C: mixed leaf counts per hub
    for m in (6, 8, 10, 12, 15):
        for h in (4, 6, 8, 10):
            for pat in ((1, 2), (1, 3), (2, 3), (1, 2, 3)):
                sv = [pat[i % len(pat)] for i in range(h)]
                if 30 <= n_msh([m] * h, sv) <= nmax:
                    S.append(('MSHsvar', dict(m=m, h=h, pat=list(pat)), ('msh', ([m] * h, sv))))
    # D: centre decorations and subdivided centre-hub edges
    for s in (1, 2, 3):
        for m in (5, 8, 10, 12, 15, 20):
            for h in (2, 4, 6, 8, 10):
                for cl in (1, 3, 8, 20):
                    if n_msh([m] * h, [s] * h, cl=cl) <= nmax:
                        S.append(('MSH+cl', dict(h=h, m=m, s=s, cl=cl), ('msh', ([m] * h, [s] * h, cl))))
                for cp in (1, 2, 3, 5):
                    if n_msh([m] * h, [s] * h, cp=cp) <= nmax:
                        S.append(('MSH+cp', dict(h=h, m=m, s=s, cp=cp), ('msh', ([m] * h, [s] * h, 0, cp))))
                for sub in (2, 3):
                    if n_msh([m] * h, [s] * h, sub=sub) <= nmax:
                        S.append(('MSH+sub', dict(h=h, m=m, s=s, sub=sub), ('msh', ([m] * h, [s] * h, 0, 0, sub))))
    # E: nesting
    for s in (1, 2):
        for m in (3, 4, 5, 6, 8, 10):
            for g in (2, 3, 4, 6):
                for h in (2, 3, 4, 6, 8):
                    for own in (0, 3):
                        n = 1 + h * (1 + g * (1 + m * (s + 1)) + own * (s + 1))
                        if 30 <= n <= nmax:
                            S.append(('nest', dict(h=h, g=g, m=m, s=s, own=own), ('nest', (h, g, m, s, own))))
    # F: two centres
    for s in (1, 2):
        for m in (5, 8, 10, 12):
            for h1 in (2, 4, 6, 8):
                for h2 in (1, 2, 4, 8):
                    for L in (1, 2, 3, 5):
                        n = 1 + L + (h1 + h2) * (1 + m * (s + 1))
                        if 30 <= n <= nmax:
                            S.append(('twocentre', dict(h1=h1, h2=h2, m=m, s=s, L=L), ('two', (h1, h2, m, s, L))))
    # G: caterpillar spines carrying hub-stars
    for s in (1, 2):
        for m in (3, 5, 8, 10, 12):
            for c in (1, 2, 3, 4):
                for L in (2, 3, 4, 6, 8, 12):
                    n = L * (1 + c * (1 + m * (s + 1)))
                    if 30 <= n <= nmax:
                        S.append(('catspine', dict(L=L, c=c, m=m, s=s), ('cat', (L, c, m, s))))
    # H: single hub-stars with centre leaves
    for s in range(1, 8):
        for m in range(3, 80):
            for cl in (0, 2, 5, 10):
                n = 1 + m * (s + 1) + cl
                if 30 <= n <= nmax:
                    S.append(('H', dict(m=m, s=s, cl=cl), ('H', (m, s, cl))))
    return S

CT = {'msh': t_msh, 'nest': t_nest, 'two': t_twocentre, 'cat': t_catspine, 'H': t_H}

if __name__ == '__main__':
    part, nparts, nmax, out = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    S = specs(nmax)
    with open(out, 'w') as f:
        for i, (fam, par, (ct, args)) in enumerate(S):
            if i % nparts != part: continue
            adj = CT[ct](*args); t0 = time.time(); r = analyse(adj)
            r.update(family=fam, params=par, ctor=[ct, list(args)], secs=round(time.time() - t0, 3))
            f.write(json.dumps(r) + '\n'); f.flush()
