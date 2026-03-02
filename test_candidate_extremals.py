"""
Test specific tree families for the Karlin rescue ratio.

The extremal is always s=2, ell=1, k=2 at a support vertex.
We construct trees as: root r with 1 leaf + 2 non-leaf subtrees T1, T2.

The tree has n = 1 + 1 + |T1| + |T2| vertices (root + leaf + subtrees).

Candidate subtrees to test:
  - Path P_m (m vertices)
  - Star S(1^k) = K_{1,k}
  - Caterpillar
  - Double star
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


def path_EJ(m):
    """E and J polynomials for a path P_m rooted at one endpoint.

    Path: v0 - v1 - ... - v_{m-1}, rooted at v0.
    """
    if m == 1:
        return [1], [1]
    # Build bottom-up
    # dp[v][0] = E_v, dp[v][1] = J_v (= dp[v][1]/x for include-root)
    E = [1]  # leaf
    J = [1]  # leaf
    for _ in range(m - 1):
        # v has one child c with (E_c, J_c) = (E, J)
        I_c = _polyadd(E, [0] + J)
        new_E = I_c  # E_v = product of I over children = I_c
        new_J = E    # J_v = product of E over children = E_c
        E, J = new_E, new_J
    return E, J


def star_EJ(k):
    """E and J for star K_{1,k} rooted at center.

    Center has k leaf children. E = (1+x)^k, J = 1.
    """
    E = [1]
    for _ in range(k):
        E = _polymul(E, [1, 1])
    J = [1]
    return E, J


def broom_EJ(arm_len, num_leaves):
    """E and J for a broom: path of length arm_len with num_leaves leaves at the end.

    Structure: v0 - v1 - ... - v_{arm_len-1} - center
    where center has num_leaves leaf children.
    Rooted at v0 (or equivalently, the path endpoint).
    """
    # Start at center (has num_leaves leaves)
    E, J = star_EJ(num_leaves)
    # Build up the path
    for _ in range(arm_len):
        I_c = _polyadd(E, [0] + J)
        new_E = I_c
        new_J = E
        E, J = new_E, new_J
    return E, J


def compute_rescue_ratio(E1, J1, E2, J2, ell=1):
    """Compute the Karlin rescue ratio at k=2 for a support vertex
    with ell leaves, subtree 1 (E1,J1) and subtree 2 (E2,J2).

    Returns (ratio_as_Fraction, Term1, t23) or None if no rescue needed.
    """
    # E_acc after processing leaf: (1+x)^ell
    E_acc = [1]
    for _ in range(ell):
        E_acc = _polymul(E_acc, [1, 1])
    J_acc = [1]

    # Process first subtree
    g0, h0 = E1, J1
    I_c0 = _polyadd(g0, [0] + h0)
    E_acc = _polymul(E_acc, I_c0)
    J_acc = _polymul(J_acc, g0)

    # Check second subtree at k=2
    g1, h1 = E2, J2
    A = _polymul(E_acc, g1)
    B = _polymul(E_acc, h1)
    C = _polymul(J_acc, g1)

    k = 2
    Ck = _coeff(C, k)
    Ckm = _coeff(C, k-1)
    Ckp = _coeff(C, k+1)
    Bk = _coeff(B, k)
    Bkm = _coeff(B, k-1)

    star_margin = Bk * Ck - Bkm * Ckp
    if star_margin >= 0:
        return None  # ⋆ doesn't fail

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


def main():
    print("=== Testing path + path families ===")
    print("Root has 1 leaf + path P_a + path P_b (a ≤ b)")
    print(f"{'n':>4s} {'a':>3s} {'b':>3s} {'ratio':>12s} {'float':>10s} {'T1':>10s} {'t23':>10s} {'w_k':>6s}")
    print("-" * 65)

    best_per_n = {}
    for total_sub in range(4, 30):  # total subtree size = a + b
        for a in range(2, total_sub - 1):
            b = total_sub - a
            if b < a:
                continue
            n = 2 + a + b  # root + leaf + two subtrees
            E1, J1 = path_EJ(a)
            E2, J2 = path_EJ(b)
            # Try both orderings
            for Ea, Ja, Eb, Jb, aa, bb in [(E1, J1, E2, J2, a, b), (E2, J2, E1, J1, b, a)]:
                result = compute_rescue_ratio(Ea, Ja, Eb, Jb)
                if result is not None:
                    ratio, T1, t23, w_k = result
                    if n not in best_per_n or ratio < best_per_n[n][0]:
                        best_per_n[n] = (ratio, aa, bb, T1, t23, w_k, 'P+P')

    for n in sorted(best_per_n.keys()):
        r, a, b, T1, t23, w_k, fam = best_per_n[n]
        print(f"{n:4d} {a:3d} {b:3d} {str(r):>12s} {float(r):10.6f} {T1:10d} {t23:10d} {w_k:6d}")

    print()
    print("=== Testing path + star families ===")
    print("Root has 1 leaf + path P_a + star K_{1,b}")
    print(f"{'n':>4s} {'a':>3s} {'b':>3s} {'ratio':>12s} {'float':>10s}")
    print("-" * 50)

    best_ps = {}
    for a in range(2, 20):
        for b in range(2, 20):
            n = 2 + a + 1 + b  # root + leaf + path_a + star(center + b leaves)
            E1, J1 = path_EJ(a)
            E2, J2 = star_EJ(b)
            for Ea, Ja, Eb, Jb, aa, bb in [(E1, J1, E2, J2, a, b), (E2, J2, E1, J1, b, a)]:
                result = compute_rescue_ratio(Ea, Ja, Eb, Jb)
                if result is not None:
                    ratio, T1, t23, w_k = result
                    if n not in best_ps or ratio < best_ps[n][0]:
                        best_ps[n] = (ratio, aa, bb, T1, t23, w_k, 'P+S')

    for n in sorted(best_ps.keys()):
        r, a, b, T1, t23, w_k, fam = best_ps[n]
        print(f"{n:4d} {a:3d} {b:3d} {str(r):>12s} {float(r):10.6f}")

    print()
    print("=== Testing path + broom families ===")
    print("Root has 1 leaf + path P_a + broom(arm_len, num_leaves)")
    print(f"{'n':>4s} {'a':>3s} {'arm':>4s} {'lv':>3s} {'ratio':>12s} {'float':>10s}")
    print("-" * 55)

    best_pb = {}
    for a in range(2, 15):
        for arm in range(1, 10):
            for lv in range(2, 10):
                broom_size = arm + 1 + lv  # path vertices + center + leaves
                n = 2 + a + broom_size
                if n > 30:
                    continue
                E1, J1 = path_EJ(a)
                E2, J2 = broom_EJ(arm, lv)
                for Ea, Ja, Eb, Jb in [(E1, J1, E2, J2), (E2, J2, E1, J1)]:
                    result = compute_rescue_ratio(Ea, Ja, Eb, Jb)
                    if result is not None:
                        ratio, T1, t23, w_k = result
                        if n not in best_pb or ratio < best_pb[n][0]:
                            best_pb[n] = (ratio, a, arm, lv, T1, t23, w_k, 'P+B')

    for n in sorted(best_pb.keys()):
        r, a, arm, lv, T1, t23, w_k, fam = best_pb[n]
        print(f"{n:4d} {a:3d} {arm:4d} {lv:3d} {str(r):>12s} {float(r):10.6f}")

    print()
    print("=== Overall minimum per n (across all families) ===")
    all_best = {}
    for n in set(list(best_per_n.keys()) + list(best_ps.keys()) + list(best_pb.keys())):
        candidates = []
        if n in best_per_n:
            candidates.append(best_per_n[n])
        if n in best_ps:
            candidates.append(best_ps[n])
        if n in best_pb:
            candidates.append(best_pb[n])
        if candidates:
            all_best[n] = min(candidates, key=lambda x: x[0])

    for n in sorted(all_best.keys()):
        entry = all_best[n]
        r = entry[0]
        fam = entry[-1]
        if fam == 'P+P':
            print(f"n={n:2d}: ratio={float(r):10.6f} = {r}  [{fam}: P_{entry[1]}+P_{entry[2]}]")
        elif fam == 'P+S':
            print(f"n={n:2d}: ratio={float(r):10.6f} = {r}  [{fam}: P_{entry[1]}+S_{entry[2]}]")
        elif fam == 'P+B':
            print(f"n={n:2d}: ratio={float(r):10.6f} = {r}  [{fam}: P_{entry[1]}+Broom({entry[2]},{entry[3]})]")


if __name__ == '__main__':
    main()
