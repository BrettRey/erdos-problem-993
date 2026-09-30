"""Float diagnostic (not a proof): predicted h->infinity limit of F_centre on S(h,m,2).
As h grows the centre's free mass vanishes, so sign L_centre -> sign D_hub, and in the
grand-canonical limit (k = h*mu(lam)) F -> max over the window of  e(lam) - 1, where
e = dln(pi)/dln(mu), pi = P(hub free) = (1+lam)^{2m}/A, A = ((1+lam)^2+lam)^m + lam (1+lam)^{2m},
mu = lam A'/A (mean size of one hub-star), window top mu = (2/3)(2m+1), bottom mu ~ (3m+1)/4."""
import math
def comp(m, lam):
    A = ((1 + lam) ** 2 + lam) ** m + lam * (1 + lam) ** (2 * m)
    return A
def lnpi(m, lam):
    return 2 * m * math.log1p(lam) - math.log(comp(m, lam))
def mu(m, lam, eps=1e-6):
    return lam * (math.log(comp(m, lam * (1 + eps))) - math.log(comp(m, lam * (1 - eps)))) / (2 * eps * lam)
def elasticity(m, lam, eps=1e-5):
    dlnpi = (lnpi(m, lam * (1 + eps)) - lnpi(m, lam * (1 - eps)))
    dlnmu = (math.log(mu(m, lam * (1 + eps))) - math.log(mu(m, lam * (1 - eps))))
    return dlnpi / dlnmu
for m in (9, 10, 11, 12, 16, 24):
    top = (2 / 3) * (2 * m + 1)
    lo_, hi_ = 1e-3, 1e3
    for _ in range(200):
        mid = math.sqrt(lo_ * hi_)
        if mu(m, mid) < top: lo_ = mid
        else: hi_ = mid
    lam_top = lo_
    grid = [lam_top * (i / 400) for i in range(40, 401)]
    best = max((elasticity(m, l) - 1, l) for l in grid if mu(m, l) >= (3 * m + 1) / 4)
    print(f"m={m}: lambda at window top {lam_top:.4f}, e-1 at top {elasticity(m, lam_top) - 1:.4f}, max_(window) e-1 = {best[0]:.4f} at lambda {best[1]:.4f}")
