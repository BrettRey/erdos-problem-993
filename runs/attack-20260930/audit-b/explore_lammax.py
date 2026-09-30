# DIAGNOSTIC (floats): paper Section 7 bound E_lam X >= alpha*log(1+lam)/Q*(lam), eq (7.4);
# which lam_max makes the ratio exceed 2/3 (what Theorem 1.1 needs) vs 17/25 (paper) vs 64/95 (Lean);
# and the true Lemma 4.2 threshold p*(lam_max).
import math
def Qstar(l): return math.sqrt((1+l)/l)*(3*math.log(2)/math.sqrt(2)+2*math.log((math.sqrt(l)+math.sqrt(1+l))/(1+math.sqrt(2))))
def ratio(l): return math.log(1+l)/Qstar(l)
def Fsup(p, lam, N=60000, Y=40.0):
    best=0
    for i in range(1,N+1):
        y=Y*i/N; r=lam*math.exp(-y); q=r/(1+r)
        v=y*q**p/math.log1p(r)
        best=max(best,v)
    return best
def pstar(lam):
    lo,hi=1.0,2.0
    for _ in range(40):
        mid=(lo+hi)/2
        if Fsup(mid,lam)<1: hi=mid
        else: lo=mid
    return hi
for l in [6,7,8,9,9.5,10,10.5,11,11.5,12]:
    ps=pstar(l)
    print(f"lam={l:5.1f} ratio={ratio(l):.5f}  (>2/3? {ratio(l)>2/3}; >0.68? {ratio(l)>0.68})  p*={ps:.4f} a*={2-2/ps:.4f} 1-a*={2/ps-1:.4f}")
# crude: where ratio crosses 2/3
lo,hi=1.0,12.0
for _ in range(60):
    m=(lo+hi)/2
    if ratio(m)>2/3: hi=m
    else: lo=m
print("ratio=2/3 at lam=",hi, " p* there:",pstar(hi))
# paper-proof interpolation limit: b_min solving 12 e^{-1-3b/2}=b
lo,hi=0.9,1.0
for _ in range(80):
    m=(lo+hi)/2
    if 12*math.exp(-1-1.5*m)>m: lo=m
    else: hi=m
bmin=hi; A=2*math.sqrt(12)/math.e
smax=abs(math.log(bmin))/(2*(math.log(A)+abs(math.log(bmin))))
print("paper Lemma4.2: b_min=",bmin," A=2sqrt12/e=",A," s_max=2-p <",smax," 1-a <",smax/(2-smax))
