"""Large-h trend of L_centre on S(h,m,2) by the closed form (flint fmpz_poly, exact):
  Ic=(1+x)^2+x, A=Ic^m+x(1+x)^{2m}, B=Ic^m, I=A^h+xB^h, Jc=Ic^{hm}=B^h, Jhub=(1+x)^{2m}A^{h-1}
  L_centre(k) = D_c/(h+1) + h D_hub/(m+2)   (exact sign test with integers)
Reports the violating k-range and F_centre = (k+1) L / (T2_c/(h+1) + h T2_hub/(m+2)) max (float diagnostic).
Only u = centre is examined (a violation at one u refutes R1).
Usage: python trend_closed_form.py > data/trend_closed_form.jsonl"""
import json, time
import flint
X = flint.fmpz_poly([0, 1]); O = flint.fmpz_poly([1, 1])
def co(p, k):
    return int(p[k]) if 0 <= k <= p.degree() else 0
for m in (10, 11, 12, 16, 24):
    for h in (30, 50, 100, 200, 400):
        t0 = time.time()
        Ic = O * O + X
        B = Ic ** m
        A = B + X * O ** (2 * m)
        I = A ** h + X * B ** h
        Jc = B ** h
        Jh = O ** (2 * m) * A ** (h - 1)
        n = 1 + h * (3 * m + 1)
        alpha = I.degree()
        lo, q = -((-n) // 4), -((-(2 * alpha - 1)) // 3)
        viol = []
        bestF = (-1e300, None)
        W1, W2 = m + 2, h + 1
        for k in range(lo, q + 1):
            a, b = k * co(I, k - 1), (k + 1) * co(I, k)
            T1c, T2c = a * co(Jc, k), b * co(Jc, k - 1)
            T1h, T2h = a * co(Jh, k), b * co(Jh, k - 1)
            LM = (T1c - T2c) * W1 + h * (T1h - T2h) * W2      # = L * (h+1)(m+2)
            G2 = T2c * W1 + h * T2h * W2
            if LM > 0:
                viol.append(k)
            F = (k + 1) * LM / G2
            if F > bestF[0]:
                bestF = (F, k)
        print(json.dumps({"m": m, "h": h, "n": n, "alpha": alpha, "window": [lo, q], "n_viol_levels": len(viol),
                          "viol_k_range": [min(viol), max(viol)] if viol else None,
                          "F_centre_max": bestF[0], "F_at_k": bestF[1], "secs": round(time.time() - t0, 1)}), flush=True)
