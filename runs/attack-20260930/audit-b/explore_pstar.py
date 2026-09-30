# DIAGNOSTIC (floats): true supremum of F_p(lam,y) = y q^p / log(1+r), r = lam e^{-y}, q=r/(1+r)
# over y>=0 at lam = lam_max (the worst case, see REPORT), versus p; and the crude interpolation
# used in the paper's Lemma 4.2 proof.
import math
def Fsup(p, lam=12.0, N=200000, Y=60.0):
    best=0; ybest=None
    for i in range(1,N+1):
        y = Y*i/N
        r = lam*math.exp(-y); q=r/(1+r)
        v = y*q**p/math.log1p(r)
        if v>best: best=v; ybest=y
    return best,ybest
for lam in [12.0, 10.0, 9.0, 8.0, 6.0]:
    print("lam",lam)
    for p in [1.34,1.40,1.45,1.5,1.52,1.55,1.6,1.7,1.8,1.9,2.0]:
        b,yb=Fsup(p,lam)
        print(f"  p={p:.2f} sup F={b:.5f} at y={yb:.3f}  a=2-2/p={2-2/p:.4f}")
