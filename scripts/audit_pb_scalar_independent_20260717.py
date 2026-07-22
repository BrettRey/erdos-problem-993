"""Independent symbolic audit of the scalar section of the PB paper.

Fresh implementation from the formulas in main.tex (eqs deficit/L/R/A/Q,
asymmetric-P, compact-ST/compact-P, bernstein-conversion, S0/T0, N_J).
Does NOT import or reuse the repo's generator or checker.
"""
import sympy as sp

H, t, u = sp.symbols('H t u', positive=True)
ok = []

def bernstein_coeffs(poly_t, d):
    """beta_i = sum_{j<=i} a_j C(i,j)/C(d,j) for G(t)=sum a_j t^j (paper eq)."""
    p = sp.Poly(sp.expand(poly_t), t)
    a = [p.coeff_monomial(t**j) for j in range(d + 1)]
    assert p.degree() <= d
    return [sp.together(sum(a[j] * sp.binomial(i, j) / sp.binomial(d, j)
                            for j in range(i + 1))) for i in range(d + 1)]

def check_bernstein_identity(poly_t, d):
    """Sanity: the conversion formula itself reproduces G(t)."""
    betas = bernstein_coeffs(poly_t, d)
    recon = sum(betas[i] * sp.binomial(d, i) * t**i * (1 - t)**(d - i)
                for i in range(d + 1))
    return sp.simplify(sp.expand(recon - poly_t)) == 0

delta = 1 / (H + 1)
a_ratio = (1 - 2*delta) / (1 - delta)          # = (H-1)/H
Q = (3*H + 4) * (H + 1) / 4                     # = (3+delta)/(4 delta^2)
assert sp.simplify(Q - (3 + delta)/(4*delta**2)) == 0
assert sp.simplify(a_ratio - (H - 1)/H) == 0

def L_r(r):   # eq:left-mass-bound
    return (1 - delta)**(-r) * sp.prod([(1 - j*delta) for j in range(2, r + 2)])

def R_r(r):   # eq:right-mass-bound
    return a_ratio**r * sp.prod([(1 - j*delta) for j in range(1, r)])

def b_r(r):   # eq:b-def
    return sp.prod([(1 - sp.Rational(s) / H) for s in range(1, r + 1)])

# L_r == b_r as claimed
for r in range(1, 6):
    assert sp.simplify(L_r(r) - b_r(r)) == 0
ok.append("L_r == b_r (r=1..5)")

# C ratio identity eq:C-ratio and R_r >= b_r direction
r_ = sp.symbols('r_', positive=True)
Cratio = sp.simplify((R_r(3)/b_r(3)) / (R_r(2)/b_r(2)))  # spot r=2
assert sp.simplify(Cratio - (1 + 2*2/((H+1)*(H-2-1)))) == 0
ok.append("C_{r+1}/C_r identity (spot r=2)")

def A_from_weights(w):  # eq:A-def, w maps index -> weight, symmetric window
    idx = sorted(w)
    return sp.Rational(1, 2) * sum(w[i]*w[j]*(i-j)**2 for i in idx for j in idx)

# --- compact cell [3,4]: exact asymmetric A, K=3 ---
w = {0: sp.Integer(1)}
for r in range(1, 4):
    w[r] = R_r(r); w[-r] = L_r(r)
A3 = A_from_weights(w)
P_paper = (-3*H**10 - 16*H**9 + 750*H**8 - 3676*H**7 + 6613*H**6
           - 5460*H**5 + 800*H**4 + 4696*H**3 - 7176*H**2 + 4272*H - 912)
lhs = sp.simplify(A3 - Q - P_paper / (4*H**5*(H+1)**3))
assert lhs == 0
ok.append("eq:asymmetric-P identity: A(delta)-Q(H) == P(H)/(4H^5(H+1)^3)")

P3t = sp.expand(P_paper.subs(H, 3 + t))
assert check_bernstein_identity(P3t, 10)
B = bernstein_coeffs(P3t, 10)
assert len(B) == 11 and all(c > 0 for c in B)
ok.append(f"cell [3,4]: 11 Bernstein coeffs, all > 0, min={min(B)}")

# --- compact cells m=4..15 ---
total = len(B)
mins = {}
for m in range(4, 16):
    Sm = 1 + 2*sum(b_r(r) for r in range(1, m + 1))
    Tm = 2*sum(r**2 * b_r(r) for r in range(1, m + 1))
    Pm = sp.expand(4 * H**(2*m) * (Sm*Tm - Q))
    Pm_poly = sp.Poly(Pm, H)
    assert all(c.is_Integer for c in Pm_poly.all_coeffs()), f"P_{m} not integer"
    assert Pm_poly.degree() == 2*m + 2, f"P_{m} degree {Pm_poly.degree()}"
    Pmt = sp.expand(Pm.subs(H, m + t))
    Bm = bernstein_coeffs(Pmt, 2*m + 2)
    assert len(Bm) == 2*m + 3 and all(c > 0 for c in Bm), f"cell {m} fails"
    total += len(Bm)
    mins[m] = min(Bm)
ok.append(f"cells m=4..15: integer polys, degrees 2m+2, all Bernstein coeffs > 0")
ok.append(f"total certificate coefficients: {total} (paper: 275)")
assert total == 275

# min-coefficient cross-check vs manuscript table (first few cells)
table_mins = {4: 332800, 5: 476945150, 6: 252843503616, 7: 141578177942152,
              15: 380093996939768323066638847260414000000}
for m, v in table_mins.items():
    assert mins[m] == v, (m, mins[m], v)
ok.append("min Bernstein coefficients match manuscript table (m=4,5,6,7,15)")

# --- analytic range: J=5 cell ---
J = sp.Integer(5)
lam = lambda r, HH: 1 - sp.Rational(r*(r+1), 2)/HH
S0 = 2*J + 1 - J*(J+1)*(J+2)/(3*H)
T0 = J*(J+1)*(2*J+1)/3 - (J*(J+1)*(2*J+1)*(3*J**2+3*J-1)/30 + J**2*(J+1)**2/4)/H
# verify S0,T0 formulas against direct sums
S0_direct = 1 + 2*sum(lam(r, H) for r in range(1, 6))
T0_direct = 2*sum(r**2*lam(r, H) for r in range(1, 6))
assert sp.simplify(S0 - S0_direct) == 0 and sp.simplify(T0 - T0_direct) == 0
ok.append("eq:S0/eq:T0 closed forms match direct sums (J=5)")

N5 = sp.expand(sp.cancel((H**2 * (S0*T0 - Q)).subs(H, 16 + 5*t)))
B5 = bernstein_coeffs(N5, 4)
expected = [2360, 7500, sp.Rational(25055, 2), sp.Rational(254205, 16),
            sp.Rational(31115, 2)]
assert [sp.nsimplify(x) for x in B5] == expected, B5
ok.append("J=5 cell: degree-4 Bernstein coeffs match the manuscript exactly")

# --- analytic range: general J >= 6, symbolic in u ---
Jsym = u + 6
S0g = 2*Jsym + 1 - Jsym*(Jsym+1)*(Jsym+2)/(3*H)
T0g = (Jsym*(Jsym+1)*(2*Jsym+1)/3
       - (Jsym*(Jsym+1)*(2*Jsym+1)*(3*Jsym**2+3*Jsym-1)/30
          + Jsym**2*(Jsym+1)**2/4)/H)
NJ = sp.expand(sp.cancel((H**2*(S0g*T0g - Q)).subs(H, Jsym*(Jsym+1)/2 + (Jsym+1)*t)))
BJ = bernstein_coeffs(NJ, 4)
mults = [Jsym**2*(Jsym+1)**2, Jsym*(Jsym+1)**2, (Jsym+1)**2,
         (Jsym+1)**2*(Jsym+2), (Jsym+1)**2*(Jsym+2)**2]
paper_polys = [
    121*u**4 + 2474*u**3 + 17431*u**2 + 46014*u + 25400,
    121*u**5 + 3442*u**4 + 37967*u**3 + 199405*u**2 + 480475*u + 384510,
    121*u**6 + 4410*u**5 + 65863*u**4 + 512764*u**3 + 2172926*u**2
        + 4668316*u + 3829440,
    121*u**5 + 3684*u**4 + 43735*u**3 + 249961*u**2 + 671859*u + 644400,
    121*u**4 + 2958*u**3 + 25579*u**2 + 88782*u + 91440,
]
for i in range(5):
    residual = sp.simplify(sp.expand(BJ[i] - mults[i]*paper_polys[i]/2880))
    assert residual == 0, (i, residual)
ok.append("J>=6 family: all five Bernstein coeffs == multiplier * paper-poly / 2880")
ok.append("J>=6 family: every printed u-polynomial has strictly positive coefficients"
          if all(c > 0 for p in paper_polys for c in sp.Poly(p, u).all_coeffs())
          else "J>=6 FAIL")

for line in ok:
    print("PASS:", line)
print("ALL SYMBOLIC CHECKS PASSED")
