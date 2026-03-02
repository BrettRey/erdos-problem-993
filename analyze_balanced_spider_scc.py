"""Analyze the balanced spider family S(2^s) for SCC 4-term decomposition.

The balanced spider S(2^s) has s arms of length 2 from a center, plus one leaf
attached to the center (making it a support vertex). This family achieves the
infimum of the diag/|cross| ratio.

For the tight case: root at the center (support vertex), one P_3 non-leaf child,
one big subtree. We compute exact formulas.
"""

import sys
sys.path.insert(0, '/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993')
from indpoly import _polymul, _polyadd


def coeff(poly, k):
    if poly is None or k < 0 or k >= len(poly):
        return 0
    return poly[k]


def balanced_spider_polys(s):
    """Compute DP polys for balanced spider S(2^s) = star with s arms of length 2.

    n = 2s + 1 vertices: center + s*(path-endpoint + intermediate).
    Root at center. Each arm contributes:
      E_c = 1 + 2x (independent sets in P_2 not including arm root)
      J_c = 1 + x  (independent sets in P_2 including arm root)
      I_c = 1 + 3x + x^2
    """
    # E_c = (1+2x), J_c = (1+x) for each arm
    # E(center) = ∏ I_c = (1+3x+x²)^s
    # J(center) = ∏ E_c = (1+2x)^s
    E = [1]
    J = [1]
    I_c = [1, 3, 1]
    E_c = [1, 2]
    for _ in range(s):
        E = _polymul(E, I_c)
        J = _polymul(J, E_c)
    return E, J


def scc_terms(e, b, k):
    """Compute SCC = e_{k+1}*b_k - e_k*b_{k+1}."""
    return coeff(e, k+1)*coeff(b, k) - coeff(e, k)*coeff(b, k+1)


def four_term_decomposition(E_acc, J_acc, P, Q, R):
    """Compute the 4-term decomposition SCC = QQ + QR + RQ + RR.

    e_old = (1+x)I_old = (1+x)(E_acc + x*J_acc)
    b_old = E_acc

    e_new = Q * e_old + x(1+x)R * b_old
    b_new = P * b_old
    """
    I_old = _polyadd(E_acc, [0] + J_acc)
    e_old = _polyadd(I_old, [0] + I_old)

    eQ = _polymul(e_old, Q)
    EQ = _polymul(E_acc, Q)
    ER = _polymul(E_acc, R)
    x1xER = _polyadd([0] + ER, [0, 0] + ER)  # x(1+x)*ER
    xER = [0] + list(ER)

    E_new = _polymul(E_acc, P)
    J_new = _polymul(J_acc, Q)
    I_new = _polyadd(E_new, [0] + J_new)
    e_new = _polyadd(I_new, [0] + I_new)

    max_k = max(len(e_new), len(E_new)) - 1

    results = []
    for k in range(max_k):
        qq = coeff(eQ, k+1)*coeff(EQ, k) - coeff(eQ, k)*coeff(EQ, k+1)
        qr = coeff(eQ, k+1)*coeff(xER, k) - coeff(eQ, k)*coeff(xER, k+1)
        rq = coeff(x1xER, k+1)*coeff(EQ, k) - coeff(x1xER, k)*coeff(EQ, k+1)
        rr = coeff(x1xER, k+1)*coeff(xER, k) - coeff(x1xER, k)*coeff(xER, k+1)
        cross = qr + rq
        diag = qq + rr
        scc = diag + cross
        results.append((k, qq, qr, rq, rr, cross, diag, scc))
    return results


def main():
    print("Balanced spider S(2^s) SCC 4-term analysis")
    print("=" * 70)

    # The tightest configuration: root at center of S(2^s), with one P_3 non-leaf child
    # and (s-1) other arms folded into one big subtree.
    # But actually, the simplest model: root at center with ALL s arms as separate factors.

    # For the observed tight case, E_acc = [1,3,1] (= I(P_3)) and the big factor
    # is the product of (s-1) arms. Let's profile this.

    print("\n--- Direct: root at center of S(2^s), all arms are factors ---")
    print(f"{'s':>3} {'n':>4} {'min_ratio':>10} {'at_k':>5} {'QQ':>8} {'|cross|':>8} {'RR':>8}")

    for s in range(2, 25):
        n = 2*s + 1
        # E_acc starts at [1], J_acc at [1]
        # Each factor: P = [1,3,1], Q = [1,2], R = [1,1]
        E_acc = [1]
        J_acc = [1]
        P = [1, 3, 1]
        Q = [1, 2]
        R = [1, 1]

        min_ratio = None
        min_info = None

        for stage in range(1, s+1):
            results = four_term_decomposition(E_acc, J_acc, P, Q, R)
            for k, qq, qr, rq, rr, cross, diag, scc in results:
                if cross < 0:
                    ratio = diag / abs(cross)
                    if min_ratio is None or ratio < min_ratio:
                        min_ratio = ratio
                        min_info = (k, stage, qq, abs(cross), rr, scc)
            E_acc = _polymul(E_acc, P)
            J_acc = _polymul(J_acc, Q)

        if min_ratio is not None:
            k, stage, qq, ac, rr, scc = min_info
            print(f"{s:3d} {n:4d} {min_ratio:10.6f} {k:5d} {qq:8d} {ac:8d} {rr:8d}")

    # Now profile the observed tight case: one P_3 factor + one big factor
    print("\n--- Observed tight case: 1st factor = single arm (P_3), 2nd = (s-1) arms ---")
    print(f"{'s':>3} {'n':>4} {'min_ratio':>10} {'at_k':>5} {'QQ':>8} {'|cross|':>8} {'RR':>8}")

    for s in range(2, 25):
        n = 2*s + 1

        # First factor: single P_2 arm
        # E_acc after stage 1: [1, 3, 1], J_acc: [1, 2]
        E_acc = [1, 3, 1]
        J_acc = [1, 2]

        # Second factor: (s-1) arms folded into one subtree
        # IS poly = (1+3x+x²)^(s-1), exclude = (1+2x)^(s-1), include = (1+x)^(s-1)
        P2 = [1]
        Q2 = [1]
        R2 = [1]
        for _ in range(s-1):
            P2 = _polymul(P2, [1, 3, 1])
            Q2 = _polymul(Q2, [1, 2])
            R2 = _polymul(R2, [1, 1])

        results = four_term_decomposition(E_acc, J_acc, P2, Q2, R2)

        min_ratio = None
        min_info = None
        for k, qq, qr, rq, rr, cross, diag, scc in results:
            if cross < 0:
                ratio = diag / abs(cross)
                if min_ratio is None or ratio < min_ratio:
                    min_ratio = ratio
                    min_info = (k, qq, abs(cross), rr, scc)

        if min_ratio is not None:
            k, qq, ac, rr, scc = min_info
            print(f"{s:3d} {n:4d} {min_ratio:10.6f} {k:5d} {qq:8d} {ac:8d} {rr:8d}")

    # B1/B2 decomposition for the balanced spider
    print("\n--- B1/B2 2-term analysis for S(2^s) ---")
    print(f"{'s':>3} {'n':>4} {'B1_B2_ratio':>12} {'B1':>12} {'B2':>12} {'SCC':>10}")

    for s in range(2, 30):
        n = 2*s + 1
        # All arms as individual factors
        E_acc = [1]
        J_acc = [1]
        P = [1, 3, 1]
        Q = [1, 2]
        R = [1, 1]

        min_ratio = None
        min_info = None

        for stage in range(1, s+1):
            I_old = _polyadd(E_acc, [0] + J_acc)
            e_old = _polyadd(I_old, [0] + I_old)

            # B1 = Δ_k(e_new, E_old*Q), B2 = Δ_k(e_new, x*E_old*R)
            E_new = _polymul(E_acc, P)
            J_new = _polymul(J_acc, Q)
            I_new = _polyadd(E_new, [0] + J_new)
            e_new = _polyadd(I_new, [0] + I_new)

            EQ = _polymul(E_acc, Q)
            ER = _polymul(E_acc, R)
            xER = [0] + list(ER)

            max_k = max(len(e_new), len(E_new)) - 1
            for k in range(max_k):
                b1 = coeff(e_new, k+1)*coeff(EQ, k) - coeff(e_new, k)*coeff(EQ, k+1)
                b2 = coeff(e_new, k+1)*coeff(xER, k) - coeff(e_new, k)*coeff(xER, k+1)
                scc = b1 + b2
                if b2 < 0:
                    ratio = b1 / abs(b2)
                    if min_ratio is None or ratio < min_ratio:
                        min_ratio = ratio
                        min_info = (k, stage, b1, b2, scc)

            E_acc = E_new
            J_acc = J_new

        if min_ratio is not None:
            k, stage, b1, b2, scc = min_info
            print(f"{s:3d} {n:4d} {min_ratio:12.6f} {b1:12d} {b2:12d} {scc:10d}")

    # Exact formulas for large s
    print("\n--- Checking B1 = 4s²+5s-5, B2 = -2s²+s+5 for S(2^s) ---")
    for s in range(2, 30):
        n = 2*s + 1
        E_acc = [1]
        J_acc = [1]
        P = [1, 3, 1]
        Q = [1, 2]
        R = [1, 1]

        for stage in range(1, s+1):
            I_old = _polyadd(E_acc, [0] + J_acc)
            e_old = _polyadd(I_old, [0] + I_old)

            E_new = _polymul(E_acc, P)
            J_new = _polymul(J_acc, Q)
            I_new = _polyadd(E_new, [0] + J_new)
            e_new = _polyadd(I_new, [0] + I_new)

            EQ = _polymul(E_acc, Q)
            ER = _polymul(E_acc, R)
            xER = [0] + list(ER)

            max_k = max(len(e_new), len(E_new)) - 1

            if stage == s:
                # At the last stage (k=1 is where the tight case lives for uniform arms)
                # Actually let's find the tight k
                for k in range(max_k):
                    b1 = coeff(e_new, k+1)*coeff(EQ, k) - coeff(e_new, k)*coeff(EQ, k+1)
                    b2 = coeff(e_new, k+1)*coeff(xER, k) - coeff(e_new, k)*coeff(xER, k+1)
                    if b2 < 0:
                        b1_formula = 4*s*s + 5*s - 5
                        b2_formula = -(2*s*s - s - 5)
                        match_b1 = "OK" if b1 == b1_formula else f"MISMATCH ({b1} vs {b1_formula})"
                        match_b2 = "OK" if b2 == b2_formula else f"MISMATCH ({b2} vs {b2_formula})"
                        print(f"s={s:3d}: k={k}, B1={b1}, B2={b2}, formula B1={b1_formula} {match_b1}, formula B2={b2_formula} {match_b2}")
                        break

            E_acc = E_new
            J_acc = J_new


if __name__ == '__main__':
    main()
