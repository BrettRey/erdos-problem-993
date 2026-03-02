"""
Prove that the Karlin rescue ratio > 1 for all star+star pairs (a,b).

The rescue ratio at k=2 for a support vertex with 1 leaf + star(a) + star(b) is:

    ratio(a,b) = N(a,b) / D(a,b)

where (from the closed-form derivation):
    C = (1+x)^s  (s = a+b)
    A = (1+x)^{s+1} + x(1+x)^{b+1}   [E_acc times g1]
    B = (1+x)^{a+1} + x(1+x)          [E_acc times h1 = [1]]

At k=2:
    C_2 = binom(s,2),  C_1 = s,  C_3 = binom(s,3)
    A_2 = binom(s+1,2) + binom(b+1,2)
    A_3 = binom(s+1,3) + binom(b+1,3)
    B_2 = binom(a+1,2) + 1
    B_1 = (a+1) + 1 = a+2

We need:
    w_2 = dk(A,C) + (B_2*C_2 - B_1*C_3)/C_2 >= 0

i.e., Term1/C_2 + (Term2+Term3)/C_2 >= 0, which is w_2 >= 0.

Strategy: compute w_2 in closed form and prove it's positive for a,b >= 1 with a+b >= 5.
(We need s = a+b >= 5 so that ⋆ fails; for s <= 4 either no failure or trivially positive.)
"""

from fractions import Fraction
from math import comb
import sympy as sp


def exact_ratio(a, b):
    """Compute exact rescue ratio for star(a)+star(b) at k=2."""
    s = a + b

    # C = (1+x)^s
    C2 = comb(s, 2)
    C1 = s
    C3 = comb(s, 3)

    # A = (1+x)^{s+1} + x(1+x)^{b+1}
    # A_k = binom(s+1,k) + binom(b+1,k-1)  for k >= 1
    A2 = comb(s+1, 2) + comb(b+1, 1)  # binom(s+1,2) + (b+1)
    A3 = comb(s+1, 3) + comb(b+1, 2)

    # B = (1+x)^{a+1} + x(1+x) = (1+x)^{a+1} + x + x^2
    B2 = comb(a+1, 2) + 1
    B1 = (a+1) + 1  # = a+2

    # Star margin (must be < 0 for ⋆ failure)
    star_margin = B2 * C2 - B1 * C3

    # dk(A,C) at k=2
    dk_AC = A3 * C2 - A2 * C3

    # Term1 = C_2 * dk_AC
    Term1 = C2 * dk_AC

    # t23 = C3*(B2*C1 - B1*C2) + B2*(C2*C2 - C1*C3)
    t23 = C3 * (B2 * C1 - B1 * C2) + B2 * (C2**2 - C1 * C3)

    # w_2 = (Term1 + t23) / C2
    if C2 == 0:
        return None
    w2 = Fraction(Term1 + t23, C2)

    if t23 >= 0 or Term1 <= 0:
        return None, star_margin, w2

    ratio = Fraction(Term1, -t23)
    return ratio, star_margin, w2


def verify_table():
    """Verify against the known exhaustive table."""
    known = {
        (1, 4): Fraction(10, 1),
        (1, 5): Fraction(7, 1),
        (2, 5): Fraction(43, 8),
        (2, 6): Fraction(19, 4),
        (3, 6): Fraction(59, 14),
        (3, 7): Fraction(75, 19),
        (4, 7): Fraction(26, 7),
        (4, 8): Fraction(32, 9),
        (5, 8): Fraction(100, 29),
        (5, 9): Fraction(10, 3),
        (6, 9): Fraction(125, 38),
        (6, 10): Fraction(147, 46),
    }

    print("=== Verification against exhaustive table ===")
    all_ok = True
    for (a, b), expected_ratio in sorted(known.items()):
        result = exact_ratio(a, b)
        ratio = result[0]
        w2 = result[2]
        ok = ratio == expected_ratio
        if not ok:
            all_ok = False
        print(f"  ({a},{b}): ratio={ratio} {'==' if ok else '!='} {expected_ratio}, w2={w2}")
    print(f"All match: {all_ok}")
    return all_ok


def symbolic_analysis():
    """Derive closed-form symbolic expression for w_2."""
    a, b = sp.symbols('a b', positive=True, integer=True)
    s = a + b

    C2 = s * (s - 1) // 2  # sympy integer division
    # Use sympy rational
    C2 = sp.Rational(1, 2) * s * (s - 1)
    C1 = s
    C3 = sp.Rational(1, 6) * s * (s - 1) * (s - 2)

    A2 = sp.Rational(1, 2) * (s + 1) * s + (b + 1)
    A3 = sp.Rational(1, 6) * (s + 1) * s * (s - 1) + sp.Rational(1, 2) * (b + 1) * b

    B2 = sp.Rational(1, 2) * (a + 1) * a + 1
    B1 = a + 2

    # dk(A,C) = A3*C2 - A2*C3
    dk_AC = sp.expand(A3 * C2 - A2 * C3)

    # Term1 = C2 * dk_AC
    Term1 = sp.expand(C2 * dk_AC)

    # t23 = C3*(B2*C1 - B1*C2) + B2*(C2^2 - C1*C3)
    t23 = sp.expand(C3 * (B2 * C1 - B1 * C2) + B2 * (C2**2 - C1 * C3))

    print("\n=== Symbolic expressions ===")

    # w2 = (Term1 + t23) / C2
    w2_num = sp.expand(Term1 + t23)
    print(f"w2 * C2 = Term1 + t23 = {sp.factor(w2_num)}")

    # The ratio
    ratio_expr = sp.Rational(-1, 1) * Term1 / t23  # ratio = Term1/|t23| = -Term1/t23 when t23 < 0
    ratio_simplified = sp.simplify(ratio_expr)
    print(f"ratio = {ratio_simplified}")

    # Factor w2_num
    w2_factored = sp.factor(w2_num)
    print(f"w2_num factored = {w2_factored}")

    # Check: substitute specific values
    print(f"\n  Check (5,9): w2_num = {w2_num.subs([(a,5),(b,9)])}")
    print(f"  Check (5,9): C2 = {C2.subs([(a,5),(b,9)])}")
    print(f"  Check (5,9): w2 = {w2_num.subs([(a,5),(b,9)])} / {C2.subs([(a,5),(b,9)])} = {sp.Rational(w2_num.subs([(a,5),(b,9)]), C2.subs([(a,5),(b,9)]))}")

    # Try to prove w2_num > 0 for a,b >= 1, a+b >= 5
    # Express in terms of a and s = a+b, then analyze
    print("\n=== Positivity analysis ===")

    # Substitute b = s - a and simplify
    s_sym = sp.Symbol('s', positive=True, integer=True)
    w2_in_s = w2_num.subs(b, s_sym - a)
    w2_in_s = sp.expand(w2_in_s)
    w2_in_s_collected = sp.collect(w2_in_s, a)
    print(f"w2_num(a, s-a) = {w2_in_s_collected}")

    # Factor
    w2_in_s_factored = sp.factor(w2_in_s)
    print(f"w2_num(a, s-a) factored = {w2_in_s_factored}")

    # Also compute Term1 and t23 separately
    Term1_simplified = sp.factor(Term1)
    t23_simplified = sp.factor(t23)
    print(f"\nTerm1 factored = {Term1_simplified}")
    print(f"t23 factored = {t23_simplified}")

    return w2_num, Term1, t23


def scan_w2():
    """Scan w2 for all a,b to verify positivity."""
    print("\n=== Scanning w2 for a=1..100, b=1..100 ===")
    min_w2 = None
    min_ab = None
    neg_count = 0
    total = 0

    for a in range(1, 101):
        for b in range(a, 101):
            result = exact_ratio(a, b)
            if result[0] is None:
                continue
            w2 = result[2]
            total += 1
            if w2 < 0:
                neg_count += 1
                print(f"  NEGATIVE w2 at ({a},{b}): w2={w2}")
            if min_w2 is None or w2 < min_w2:
                min_w2 = w2
                min_ab = (a, b)

    print(f"Total ⋆-failing pairs: {total}")
    print(f"Negative w2: {neg_count}")
    print(f"Minimum w2: {min_w2} at {min_ab}")

    # Find minimum w2 (absolute) for small a,b where ⋆ fails
    print("\n=== Smallest w2 values ===")
    results = []
    for a in range(1, 50):
        for b in range(a, 50):
            result = exact_ratio(a, b)
            if result[0] is None:
                continue
            w2 = result[2]
            results.append((w2, a, b, result[0]))

    results.sort()
    for w2, a, b, ratio in results[:20]:
        print(f"  w2={float(w2):10.4f} ({w2}) at ({a},{b}), ratio={float(ratio):.6f}")


def main():
    verify_table()
    symbolic_analysis()
    scan_w2()


if __name__ == '__main__':
    main()
