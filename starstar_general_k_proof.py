"""
Attempt algebraic proof of w_k >= 0 for star+star at general k.

For the star+star family at k >= 3 (where B_k = binom(a+1,k)):

w_k = dk(A,C) + B_k*C_k - B_{k-1}*C_{k+1}

where:
  C_k = binom(s,k),  s = a+b
  A_{k+1} = binom(s+1,k+1) + binom(b+1,k)
  A_k = binom(s+1,k) + binom(b+1,k-1)
  B_k = binom(a+1,k)  [for k >= 3]
  B_{k-1} = binom(a+1,k-1)  [for k >= 4]

We proved: dk(A,C) = c_k(C) + delta_k
where c_k(C) = binom(s,k)^2 - binom(s,k-1)*binom(s,k+1) >= 0 (Newton)
and delta_k = binom(b+1,k)*binom(s,k) - binom(b+1,k-1)*binom(s,k+1).

So w_k = c_k + delta_k + binom(a+1,k)*binom(s,k) - binom(a+1,k-1)*binom(s,k+1).

For k >= 3: w_k = c_k + dk(F, C) where F = x(1+x)^{b+1} + (1+x)^{a+1}.
"""

from math import comb
from fractions import Fraction


def compute_wk_exact(a, b, k):
    """Compute w_k exactly for star+star."""
    s = a + b

    Ck = comb(s, k)
    Ckm = comb(s, k - 1) if k >= 1 else 0
    Ckp = comb(s, k + 1)

    # A = (1+x)^{s+1} + x(1+x)^{b+1}
    Ak = comb(s + 1, k) + (comb(b + 1, k - 1) if k >= 1 else 0)
    Akp = comb(s + 1, k + 1) + comb(b + 1, k)

    # B = (1+x)^{a+1} + x + x^2
    Bk = comb(a + 1, k) + (1 if k == 1 else 0) + (1 if k == 2 else 0)
    Bkm = (comb(a + 1, k - 1) if k >= 1 else 0) + (1 if k == 2 else 0) + (1 if k == 3 else 0)

    dk_AC = Akp * Ck - Ak * Ckp
    star_margin = Bk * Ck - Bkm * Ckp

    wk = dk_AC + star_margin
    ck = Ck * Ck - Ckm * Ckp  # LC gap of C

    # delta_k
    delta_k = comb(b + 1, k) * Ck - (comb(b + 1, k - 1) if k >= 1 else 0) * Ckp

    return {
        'wk': wk,
        'dk_AC': dk_AC,
        'star_margin': star_margin,
        'ck': ck,
        'delta_k': delta_k,
        'decomp_check': ck + delta_k - dk_AC,  # should be 0
    }


def analyze_decomposition():
    """Check whether c_k alone dominates w_k."""
    print("=== Decomposition analysis: w_k = c_k + delta_k + star_margin ===\n")
    print("For k >= 3, star_margin = dk((1+x)^{a+1}, (1+x)^s)")
    print("and delta_k = dk(x(1+x)^{b+1}, (1+x)^s)\n")

    # Track: what fraction of w_k comes from c_k vs correction?
    min_ratio_ck_wk = None
    worst = None

    for a in range(1, 50):
        for b in range(a, 50):
            s = a + b
            for k in range(0, s):
                r = compute_wk_exact(a, b, k)
                assert r['decomp_check'] == 0, f"Decomp failed at ({a},{b},{k})"
                assert r['wk'] >= 0, f"NEGATIVE w_k at ({a},{b},{k}): {r['wk']}"

                if r['ck'] > 0 and r['wk'] > 0:
                    ratio = Fraction(r['ck'], r['wk'])
                    if min_ratio_ck_wk is None or ratio < min_ratio_ck_wk:
                        min_ratio_ck_wk = ratio
                        worst = (a, b, k, r)

    print(f"Min c_k/w_k ratio: {float(min_ratio_ck_wk):.6f} = {min_ratio_ck_wk}")
    print(f"  at (a,b,k) = ({worst[0]},{worst[1]},{worst[2]})")
    r = worst[3]
    print(f"  w_k = {r['wk']}, c_k = {r['ck']}, delta_k = {r['delta_k']}, star_margin = {r['star_margin']}")
    print(f"  correction = delta_k + star_margin = {r['delta_k'] + r['star_margin']}")
    print()

    # Check: is c_k >= |correction| always?
    max_correction_ratio = None
    for a in range(1, 50):
        for b in range(a, 50):
            s = a + b
            for k in range(0, s):
                r = compute_wk_exact(a, b, k)
                correction = r['delta_k'] + r['star_margin']
                if correction < 0 and r['ck'] > 0:
                    ratio = Fraction(-correction, r['ck'])
                    if max_correction_ratio is None or ratio > max_correction_ratio:
                        max_correction_ratio = ratio
                        worst_cr = (a, b, k, r)

    if max_correction_ratio is not None:
        print(f"Max |correction|/c_k ratio: {float(max_correction_ratio):.6f} = {max_correction_ratio}")
        print(f"  at (a,b,k) = ({worst_cr[0]},{worst_cr[1]},{worst_cr[2]})")
        r = worst_cr[3]
        print(f"  correction = {r['delta_k'] + r['star_margin']}, c_k = {r['ck']}")
        if max_correction_ratio < 1:
            print("  ==> c_k ALWAYS dominates |correction| ==> w_k >= 0 PROVED!")
    else:
        print("Correction is NEVER negative ==> w_k >= c_k >= 0 trivially!")


def check_alternative_decomp():
    """Try: w_k = dk(A+xB, (1+x)C) or similar."""
    print("\n=== Alternative decomposition ===")
    print("Checking w_k = dk(E_full, J_full) where E_full = A + xB, J_full = C\n")

    # Actually, w_k in the incremental identity relates E_{new} and J_{new}
    # E_new = I_{root} = E_acc · I_t where I_t = g + xh = (1+x)^b + x = (1+x)^b + x
    # J_new = E_acc_prev · g (but E_acc_prev ≠ E_acc)
    # Actually, the E≽J condition we want is: dk(E_root, J_root) >= 0

    # For the star+star tree rooted at the support vertex:
    # Children: leaf, star(a), star(b)
    # E_root = (1+x) · [(1+x)^a + x] · [(1+x)^b + x]
    # J_root = 1 · (1+x)^a · (1+x)^b = (1+x)^s

    # So we want dk(E_root, J_root) >= 0 for all k = 0, ..., mode(I_root)-1
    # where I_root = E_root + x·J_root

    for a, b in [(1, 4), (3, 7), (5, 9), (10, 20)]:
        s = a + b

        # Compute E_root and J_root
        # E_root = (1+x) * ((1+x)^a + x) * ((1+x)^b + x)
        # J_root = (1+x)^s

        # (1+x)^a + x
        Ia = [comb(a, k) for k in range(a + 1)]
        Ia[1] += 1  # add x

        # (1+x)^b + x
        Ib = [comb(b, k) for k in range(b + 1)]
        Ib[1] += 1  # add x

        # Multiply: E_root_partial = Ia * Ib
        deg_partial = len(Ia) + len(Ib) - 2
        E_partial = [0] * (deg_partial + 1)
        for i, ai in enumerate(Ia):
            for j, bj in enumerate(Ib):
                E_partial[i + j] += ai * bj

        # E_root = (1+x) * E_partial
        E_root = [0] * (len(E_partial) + 1)
        for i, ei in enumerate(E_partial):
            E_root[i] += ei
            E_root[i + 1] += ei

        # J_root = (1+x)^s
        J_root = [comb(s, k) for k in range(s + 1)]

        # I_root = E_root + x*J_root
        deg_I = max(len(E_root), len(J_root) + 1)
        I_root = [0] * deg_I
        for i in range(len(E_root)):
            I_root[i] += E_root[i]
        for i in range(len(J_root)):
            I_root[i + 1] += J_root[i]

        # Find mode
        mode = max(range(len(I_root)), key=lambda k: I_root[k])

        print(f"(a,b) = ({a},{b}), s={s}, n={s+4}, mode={mode}")
        print(f"  {'k':>3s} {'dk(E,J)':>12s} {'E_{k+1}':>10s} {'J_k':>10s} {'E_k':>10s} {'J_{k+1}':>10s}")

        for k in range(min(mode + 2, len(J_root))):
            Ek = E_root[k] if k < len(E_root) else 0
            Ekp = E_root[k + 1] if k + 1 < len(E_root) else 0
            Jk = J_root[k] if k < len(J_root) else 0
            Jkp = J_root[k + 1] if k + 1 < len(J_root) else 0

            dk_EJ = Ekp * Jk - Ek * Jkp
            print(f"  {k:3d} {dk_EJ:12d} {Ekp:10d} {Jk:10d} {Ek:10d} {Jkp:10d}")

        print()


if __name__ == '__main__':
    analyze_decomposition()
    check_alternative_decomp()
