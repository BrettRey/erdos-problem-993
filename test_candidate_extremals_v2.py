"""
Test more subtree families for the Karlin rescue ratio.

From test_candidate_extremals.py, path+path gives NO ⋆ failures and
path+star gives ratio ~5.8-6 at best. The actual extremal (ratio 3.33
at n=18) must involve a different subtree type.

New candidates:
- Caterpillar rooted at one end
- Double star (two adjacent star centers)
- Path with pendant edges ("combs")
- Binary trees
"""

from fractions import Fraction

try:
    from numpy import convolve as _np_convolve
    _HAS_NUMPY = True
except ImportError:
    _HAS_NUMPY = False

_INT64_SAFE = 2**62


def _polymul_python(a, b):
    la, lb = len(a), len(b)
    out = [0] * (la + lb - 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    return out


def _polymul(a, b):
    if not a or not b:
        return []
    if _HAS_NUMPY:
        max_a = max(abs(x) for x in a)
        max_b = max(abs(x) for x in b)
        max_terms = min(len(a), len(b))
        if max_a > 0 and max_b > 0 and max_terms * max_a * max_b < _INT64_SAFE:
            return _np_convolve(a, b).astype(int).tolist()
    return _polymul_python(a, b)


def _polyadd(a, b):
    la, lb = len(a), len(b)
    out = [0] * max(la, lb)
    for i in range(la):
        out[i] += a[i]
    for i in range(lb):
        out[i] += b[i]
    return out


def _coeff(p, k):
    return p[k] if 0 <= k < len(p) else 0


def subtree_EJ(children_list):
    """Compute (E, J) for a rooted tree given as adjacency via children_list.

    children_list[v] = list of children of vertex v.
    Root is vertex 0.
    """
    n = len(children_list)
    # Post-order traversal
    order = []
    stack = [(0, False)]
    while stack:
        v, done = stack.pop()
        if done:
            order.append(v)
            continue
        stack.append((v, True))
        for c in children_list[v]:
            stack.append((c, False))

    dp0 = [None] * n  # E
    dp1s = [None] * n  # J
    for v in order:
        if not children_list[v]:
            dp0[v] = [1]
            dp1s[v] = [1]
        else:
            pS = [1]
            pE = [1]
            for c in children_list[v]:
                Ic = _polyadd(dp0[c], [0] + dp1s[c])
                pS = _polymul(pS, Ic)
                pE = _polymul(pE, dp0[c])
            dp0[v] = pS
            dp1s[v] = pE
    return dp0[0], dp1s[0]


def make_path(m):
    """Path P_m as children list, rooted at one end."""
    children = [[] for _ in range(m)]
    for i in range(m-1):
        children[i].append(i+1)
    return children


def make_star(k):
    """Star K_{1,k} as children list, rooted at center."""
    children = [[] for _ in range(k+1)]
    for i in range(1, k+1):
        children[0].append(i)
    return children


def make_pendant_star(p, k):
    """Pendant star: path of length p, then star K_{1,k} at the end.

    v0 - v1 - ... - v_{p-1} - center - (k leaves)
    Rooted at v0.
    Total vertices: p + 1 + k
    """
    total = p + 1 + k
    children = [[] for _ in range(total)]
    # Path
    for i in range(p-1):
        children[i].append(i+1)
    # Connect path end to center
    center = p
    if p > 0:
        children[p-1].append(center)
    # Leaves
    for i in range(k):
        children[center].append(p + 1 + i)
    return children


def make_caterpillar(spine_len, pendants_per_vertex):
    """Caterpillar: spine of length spine_len, each internal spine vertex has
    pendants_per_vertex pendant leaves.

    Rooted at one end of spine.
    """
    total = spine_len + (spine_len - 2) * pendants_per_vertex if spine_len > 2 else spine_len
    if spine_len <= 2:
        return make_path(spine_len)
    children = [[] for _ in range(total)]
    # Spine
    for i in range(spine_len - 1):
        children[i].append(i + 1)
    # Pendants on internal vertices (1 to spine_len-2)
    next_v = spine_len
    for i in range(1, spine_len - 1):
        for _ in range(pendants_per_vertex):
            children[i].append(next_v)
            next_v += 1
    return children


def make_double_star(k1, k2):
    """Double star: two adjacent centers, one with k1 leaves, other with k2 leaves.

    Rooted at center1.
    v0 (center1) has children: center2, plus k1 leaves.
    center2 has k2 leaf children.
    Total: 2 + k1 + k2
    """
    total = 2 + k1 + k2
    children = [[] for _ in range(total)]
    center2 = 1
    children[0].append(center2)
    for i in range(k1):
        children[0].append(2 + i)
    for i in range(k2):
        children[center2].append(2 + k1 + i)
    return children


def compute_rescue_ratio(E1, J1, E2, J2, ell=1, k=2):
    """Compute the Karlin rescue ratio at given k."""
    E_acc = [1]
    for _ in range(ell):
        E_acc = _polymul(E_acc, [1, 1])
    J_acc = [1]

    g0, h0 = E1, J1
    I_c0 = _polyadd(g0, [0] + h0)
    E_acc = _polymul(E_acc, I_c0)
    J_acc = _polymul(J_acc, g0)

    g1, h1 = E2, J2
    A = _polymul(E_acc, g1)
    B = _polymul(E_acc, h1)
    C = _polymul(J_acc, g1)

    Ck = _coeff(C, k)
    Ckm = _coeff(C, k-1)
    Ckp = _coeff(C, k+1)
    Bk = _coeff(B, k)
    Bkm = _coeff(B, k-1)

    star_margin = Bk * Ck - Bkm * Ckp
    if star_margin >= 0:
        return None

    Ak = _coeff(A, k+1)
    Akm = _coeff(A, k)
    dk_AC = Ak * Ck - Akm * Ckp
    Term1 = Ck * dk_AC
    t23 = Ckp * (Bk * Ckm - Bkm * Ck) + Bk * (Ck * Ck - Ckm * Ckp)

    if t23 >= 0 or Term1 <= 0:
        return None

    ratio = Fraction(Term1, abs(t23))
    w_k = (Term1 + t23) // Ck if Ck > 0 else 0
    return ratio, Term1, t23, w_k


def tree_size(children_list):
    return len(children_list)


def main():
    best_overall = {}

    # 1. Pendant star + pendant star
    print("=== Pendant star + pendant star ===")
    print(f"{'n':>4s} {'p1':>3s} {'k1':>3s} {'p2':>3s} {'k2':>3s} {'ratio':>12s} {'float':>10s}")
    print("-" * 55)
    for p1 in range(0, 10):
        for k1 in range(1, 15):
            sz1 = p1 + 1 + k1 if p1 > 0 else 1 + k1
            for p2 in range(0, 10):
                for k2 in range(1, 15):
                    sz2 = p2 + 1 + k2 if p2 > 0 else 1 + k2
                    n = 2 + sz1 + sz2
                    if n > 35:
                        continue
                    if p1 == 0:
                        E1, J1 = subtree_EJ(make_star(k1))
                    else:
                        E1, J1 = subtree_EJ(make_pendant_star(p1, k1))
                    if p2 == 0:
                        E2, J2 = subtree_EJ(make_star(k2))
                    else:
                        E2, J2 = subtree_EJ(make_pendant_star(p2, k2))

                    for Ea, Ja, Eb, Jb, pa, ka, pb, kb in [
                        (E1, J1, E2, J2, p1, k1, p2, k2),
                        (E2, J2, E1, J1, p2, k2, p1, k1)
                    ]:
                        result = compute_rescue_ratio(Ea, Ja, Eb, Jb)
                        if result is not None:
                            ratio, T1, t23, w_k = result
                            if n not in best_overall or ratio < best_overall[n][0]:
                                best_overall[n] = (ratio, f"PS({pa},{ka})+PS({pb},{kb})")

    for n in sorted(best_overall.keys()):
        r, desc = best_overall[n]
        print(f"n={n:2d}: ratio={float(r):10.6f} = {r}  [{desc}]")

    # 2. Double star + path
    print()
    print("=== Double star + path (or pendant star) ===")
    best_ds = {}
    for k1 in range(1, 12):
        for k2 in range(1, 12):
            sz1 = 2 + k1 + k2
            E1, J1 = subtree_EJ(make_double_star(k1, k2))
            for p2 in range(0, 12):
                for k_arm in range(1, 12):
                    if p2 == 0:
                        sz2 = 1 + k_arm  # star
                        E2, J2 = subtree_EJ(make_star(k_arm))
                    else:
                        sz2 = p2 + 1 + k_arm  # pendant star
                        E2, J2 = subtree_EJ(make_pendant_star(p2, k_arm))
                    n = 2 + sz1 + sz2
                    if n > 35:
                        continue
                    for Ea, Ja, Eb, Jb in [(E1, J1, E2, J2), (E2, J2, E1, J1)]:
                        result = compute_rescue_ratio(Ea, Ja, Eb, Jb)
                        if result is not None:
                            ratio, T1, t23, w_k = result
                            desc = f"DS({k1},{k2})+PS({p2},{k_arm})"
                            if n not in best_ds or ratio < best_ds[n][0]:
                                best_ds[n] = (ratio, desc)

    for n in sorted(best_ds.keys()):
        r, desc = best_ds[n]
        print(f"n={n:2d}: ratio={float(r):10.6f} = {r}  [{desc}]")

    # 3. Merge all
    for n, (r, desc) in best_ds.items():
        if n not in best_overall or r < best_overall[n][0]:
            best_overall[n] = (r, desc)

    print()
    print("=== Overall minimum per n ===")
    for n in sorted(best_overall.keys()):
        r, desc = best_overall[n]
        print(f"n={n:2d}: ratio={float(r):10.6f} = {r}  [{desc}]")


if __name__ == '__main__':
    main()
