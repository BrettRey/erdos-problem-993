"""R1 adversarial lane: exact rerooting evaluator for the closed-neighbourhood
load L_u(k) = sum_{v in N[u]} D_v(k)/(deg v + 1), with
D_v(k) = k i_{k-1}(T) j^v_k - (k+1) i_k(T) j^v_{k-1},  j^v = I(T - N[v]).

Algorithm (independent of route-freecount/pv_lib.py, which runs a fresh DP on
each forest T - N[v]): one rooted DP plus one rerooting pass, polynomials as
python-flint fmpz_poly (exact integers).
  P[b] = I(subtree b), Q[b] = I(subtree b - b), R[b] = prod_children Q[c]
  Pu[c] = I(T - subtree c), Qu[c] = I(T - subtree c - parent c)
  j^v = R[v] * Qu[v].
All inequality verdicts are exact integer tests (L_u * M with M = lcm(deg+1)).
Floats appear only in the labelled normalised diagnostics:
  normA = L_u / (k i_{k-1} i_k)                         (task normalisation)
  normB = L_u / (sum_v |D_v| / n)                       (wave-1 normalisation)
  normC = pos_u / neg_u, pos/neg = weighted positive / negative D mass in N[u]
          (absorption ratio; R1 fails at (u,k) iff normC > 1)
  normE = L_u / sum_{v in N[u]} (T1_v + T2_v)/(deg v+1), T1 = k i_{k-1} j_k,
          T2 = (k+1) i_k j_{k-1}, D = T1 - T2 (gross-term normalisation, in
          (-1,1), R1 fails iff normE > 0; continuous even when all D_v < 0)
  normF = (k+1) L_u / sum_{v in N[u]} T2_v/(deg v+1)   (load in units of the
          built-in 1/k slack: a flat free-probability profile on N[u] gives
          exactly -1; R1 fails iff normF > 0.  Primary steering objective.)
"""
from math import lcm
import flint

X = flint.fmpz_poly([0, 1])
ONE = flint.fmpz_poly([1])


def adj_from_edges(n, edges):
    adj = [[] for _ in range(n)]
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    return adj


def adj_from_parent(par):
    """0-indexed parent array, par[root] = -1."""
    n = len(par)
    return adj_from_edges(n, [(i, p) for i, p in enumerate(par) if p >= 0])


def coeffs(p):
    return [int(c) for c in p.coeffs()]


def all_seqs(adj):
    """Return I (list of ints) and J (list over v of list of ints)."""
    n = len(adj)
    parent = [-1] * n
    order = [0]
    seen = [False] * n
    seen[0] = True
    for u in order:
        for w in adj[u]:
            if not seen[w]:
                seen[w] = True
                parent[w] = u
                order.append(w)
    assert len(order) == n, "not connected"
    children = [[w for w in adj[u] if w != parent[u]] for u in range(n)]
    P = [None] * n
    Q = [None] * n
    R = [None] * n
    for b in reversed(order):
        q = ONE
        r = ONE
        for c in children[b]:
            q = q * P[c]
            r = r * Q[c]
        Q[b] = q
        R[b] = r
        P[b] = q + X * r
    Pu = [None] * n
    Qu = [None] * n
    Pu[0] = ONE
    Qu[0] = ONE
    for a in order:
        ch = children[a]
        d = len(ch)
        if d == 0:
            continue
        preP = [ONE] * (d + 1)
        preQ = [ONE] * (d + 1)
        for i, c in enumerate(ch):
            preP[i + 1] = preP[i] * P[c]
            preQ[i + 1] = preQ[i] * Q[c]
        sufP = ONE
        sufQ = ONE
        for i in range(d - 1, -1, -1):
            c = ch[i]
            exP = preP[i] * sufP
            exQ = preQ[i] * sufQ
            Qu[c] = exP * Pu[a]
            Pu[c] = Qu[c] + X * exQ * Qu[a]
            sufP = sufP * P[c]
            sufQ = sufQ * Q[c]
    I = coeffs(P[0])
    J = [coeffs(R[v] * Qu[v]) for v in range(n)]
    # exact double-count identity: k i_k = sum_v j^v_{k-1}
    alpha = len(I) - 1
    for k in range(1, alpha + 2):
        s = sum(Jv[k - 1] if k - 1 < len(Jv) else 0 for Jv in J)
        assert k * (I[k] if k < len(I) else 0) == s, ("double count", k)
    return I, J


def co(s, k):
    return s[k] if 0 <= k < len(s) else 0


def window(n, alpha):
    return max(1, -((-n) // 4)), -((-(2 * alpha - 1)) // 3)


def evaluate(adj, full=False):
    """Exact R1 check on the window + float diagnostics.
    Returns dict with exact violation list and the maxima of normA/B/C."""
    n = len(adj)
    I, J = all_seqs(adj)
    alpha = len(I) - 1
    lo, q = window(n, alpha)
    deg = [len(a) for a in adj]
    M = 1
    for d in set(deg):
        M = lcm(M, d + 1)
    W = [M // (d + 1) for d in deg]
    best = {"A": (-1e300, None), "B": (-1e300, None), "C": (-1e300, None), "E": (-1e300, None), "F": (-1e300, None)}
    viol = []
    lc_fail = []
    for k in range(lo, q + 1):
        a, b = k * co(I, k - 1), (k + 1) * co(I, k)
        T1 = [a * co(Jv, k) for Jv in J]
        T2 = [b * co(Jv, k - 1) for Jv in J]
        D = [T1[v] - T2[v] for v in range(n)]
        WG = [W[v] * (T1[v] + T2[v]) for v in range(n)]
        W2 = [W[v] * T2[v] for v in range(n)]
        sD = sum(D)
        if sD > 0:
            lc_fail.append(k)
        denA = k * co(I, k - 1) * co(I, k)
        sabs = sum(abs(x) for x in D)
        WD = [W[v] * D[v] for v in range(n)]
        for u in range(n):
            LM = WD[u]
            pos = WD[u] if WD[u] > 0 else 0
            neg = -WD[u] if WD[u] < 0 else 0
            G = WG[u]
            G2 = W2[u]
            for v in adj[u]:
                G += WG[v]
                G2 += W2[v]
                x = WD[v]
                LM += x
                if x > 0:
                    pos += x
                else:
                    neg -= x
            if LM > 0:
                viol.append({"u": u, "k": k, "deg_u": deg[u], "L_u": f"{LM}/{M}"})
            fA = LM / (M * denA)
            fB = (LM * n) / (M * sabs) if sabs else 0.0
            fC = pos / neg if neg else (float("inf") if pos else 0.0)
            fE = LM / G if G else 0.0
            fF = (k + 1) * LM / G2 if G2 else (float("inf") if LM > 0 else -1e300)
            if fF > best["F"][0]:
                best["F"] = (fF, (u, k))
            if fE > best["E"][0]:
                best["E"] = (fE, (u, k))
            if fA > best["A"][0]:
                best["A"] = (fA, (u, k))
            if fB > best["B"][0]:
                best["B"] = (fB, (u, k))
            if fC > best["C"][0]:
                best["C"] = (fC, (u, k))
    out = {"n": n, "alpha": alpha, "window": [lo, q], "viol": viol, "lc_fail": lc_fail,
           "A": best["A"][0], "A_at": best["A"][1], "B": best["B"][0], "B_at": best["B"][1],
           "C": best["C"][0], "C_at": best["C"][1], "E": best["E"][0], "E_at": best["E"][1],
           "F": best["F"][0], "F_at": best["F"][1]}
    for key in ("A_at", "B_at", "C_at", "E_at", "F_at"):
        if out[key] is not None:
            u, k = out[key]
            out[key] = {"u": u, "deg_u": deg[u], "k": k, "k_rel": (k - lo) / max(1, q - lo)}
    return out


class B:
    """tree builder: vertex 0 is the root"""
    def __init__(self):
        self.n = 1
        self.edges = []

    def add(self, p):
        v = self.n
        self.n += 1
        self.edges.append((p, v))
        return v

    def path(self, p, L):
        cur = p
        for _ in range(L):
            cur = self.add(cur)
        return cur

    def adj(self):
        return adj_from_edges(self.n, self.edges)
