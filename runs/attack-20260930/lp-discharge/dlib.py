"""Shared helpers: trees, window defect rows, families. Defects via the wave-1
pv_lib (validated exactly against C for n <= 13)."""
import sys, json, subprocess
from fractions import Fraction
sys.path.insert(0, '../route-freecount')
from pv_lib import tree_data, co, window, parse_parent_array

def rows_for(adj):
    """For each k in the window (k <= alpha-1): exact defects D_v(k) and their
    floats normalized by k*i_{k-1}*i_k. D_v(k) = k i_{k-1} j^v_k - (k+1) i_k j^v_{k-1}."""
    n = len(adj); I, J, E, a = tree_data(adj); lo, q = window(n, a)
    out = []
    for k in range(lo, min(q, a - 1) + 1):
        D = [k * co(I, k - 1) * co(J[v], k) - (k + 1) * co(I, k) * co(J[v], k - 1) for v in range(n)]
        assert sum(D) == k * (k + 1) * (co(I, k - 1) * co(I, k + 1) - co(I, k) ** 2)
        nrm = k * co(I, k - 1) * co(I, k)
        out.append((k, D, [float(Fraction(d, nrm)) for d in D]))
    return out

def gentreeg(n):
    p = subprocess.run(['/opt/homebrew/bin/gentreeg', '-p', '-q', str(n)], capture_output=True, text=True, check=True)
    for line in p.stdout.split('\n'):
        if line.strip(): yield parse_parent_array(line.split())

def from_edges(n, edges):
    adj = [[] for _ in range(n)]
    for a, b in edges: adj[a].append(b); adj[b].append(a)
    return adj

def build(spec):
    """Tree from nested spec: list of child specs; [] is a leaf."""
    edges = []; cnt = [1]
    def rec(node, children):
        for ch in children:
            c = cnt[0]; cnt[0] += 1; edges.append((node, c)); rec(c, ch)
    rec(0, spec)
    return from_edges(cnt[0], edges)

leaf = []
def star(m): return [leaf] * m
def hubstar(m, s): return [star(s)] * m                      # H(m,s) rooted at hub
def sos(m, s): return [star(s)] * m                          # star of stars: centre -> m hubs -> s leaves
def msh(hubs, s=2): return [hubstar(m, s) for m in hubs]      # centre -> hub-stars H(m,s)

def families():
    F = {}
    for m in range(3, 16):
        for s in (1, 2, 3):
            F[f'H({m},{s})'] = build(hubstar(m, s))
    for m in range(3, 8):
        for s in (2, 3, 4):
            F[f'SoS({m},{s})'] = build(sos(m, s))
    for hubs in ([9, 9], [9] * 4, [9] * 6, [9] * 8, [9] * 4 + [10] * 4, [8] * 8, [10] * 8, [5] * 8, [12] * 6):
        F['MSH(' + ','.join(map(str, hubs)) + ')'] = build(msh(hubs))
    for a in range(3, 20, 2):
        F[f'bistar_mid({a},{a})'] = build([[star(a)[0]] * 0 + star(a), star(a)]) if False else build([star(a), star(a)])
    return F

def load_r1_killer():
    d = json.load(open('../r1-analytic/data/recheck_MSH_4x9_4x10_generic.json'))
    return from_edges(d['n'], d['edge_list_0indexed'])
