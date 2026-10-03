"""Exact checks of the external review's claims (2026-10-03)."""
from fractions import Fraction as F
import random, math
random.seed(993)

def pmf(ps):
    f=[F(1)]
    for p in ps:
        g=[F(0)]*(len(f)+1)
        for k,x in enumerate(f):
            g[k]+=x*(1-p); g[k+1]+=x*p
        f=g
    return f
def deltas(f):
    e=[F(0)]+f+[F(0)]
    return [1-e[k]*e[k+2]/e[k+1]**2 for k in range(len(f))]

viol={"recip":0,"support_D":0,"support_mu":0,"mode_lower":0,"mode_upper":0,"tail_new":0,"tail_pitman_claim":0}
maxVdc=F(0); minVdc=None; tested=0
for trial in range(3000):
    n=random.randint(2,36); mode=random.random(); ps=[]
    for _ in range(n):
        if mode<0.3: d=random.choice([7,13,29,61]); x=F(random.randint(1,d-1),d)
        elif mode<0.6:
            x=F(random.choice([1,2,3]),random.choice([50,100,400]))
            if random.random()<0.5: x=1-x
        else: x=F(random.randint(1,999),1000)
        ps.append(x)
    V=sum(p*(1-p) for p in ps)
    f=pmf(ps); dl=deltas(f); N=len(f)-1
    # reciprocal deficit, all adjacent pairs
    for k in range(N):
        if abs(1/dl[k+1]-1/dl[k])>1: viol["recip"]+=1
    if V<1: continue
    tested+=1
    mu=sum(ps)
    D=next(k for k in range(1,N+1) if f[k]<f[k-1]); c=D-1
    for k in range(N+1):
        if dl[k] < 1/(4*V+abs(k-D)): viol["support_D"]+=1
        if dl[k] < 1/(4*V+abs(k-mu)+2): viol["support_mu"]+=1
    if dl[c] < 1/(4*V+1): viol["mode_lower"]+=1
    if not V*dl[c] < 2: viol["mode_upper"]+=1
    maxVdc=max(maxVdc,V*dl[c]); minVdc=V*dl[c] if minVdc is None else min(minVdc,V*dl[c])
    # tails
    for r in range(1,N-D+1):
        ratio=f[D+r]/f[D]
        b_new=F(1)
        for j in range(1,r+1): b_new*= (4*V-1)/(4*V+j-1)
        if ratio>b_new: viol["tail_new"]+=1
        b_p=F(1)
        for j in range(1,r+1): b_p*= V/(V+j-1)
        if ratio>b_p: viol["tail_pitman_claim"]+=1
print("laws with V>=1:",tested)
print("violations:",viol)
print("V*delta_c range:", float(minVdc), float(maxVdc))

# ULC counterexample g_{m+j} ∝ C(2m,m+j) q^|j|
q=F(1,2)
for m in [5,20,80,200]:
    w=[math.comb(2*m,m+j)*q**abs(j) for j in range(-m,m+1)]
    Z=sum(w); g=[x/Z for x in w]
    # ULC(2m): g_k/C(2m,k) log-concave
    h=[g[k]/math.comb(2*m,k) for k in range(2*m+1)]
    ulc=all(h[k]**2>=h[k-1]*h[k+1] for k in range(1,2*m))
    D=m+1
    dD=1-g[D-1]*g[D+1]/g[D]**2
    mean=sum(k*g[k] for k in range(2*m+1)); var=sum((k-mean)**2*g[k] for k in range(2*m+1))
    print(f"ULC m={m}: ULC={ulc} D_check={g[m+1]<g[m] and g[m]>=g[m-1]} delta_D={float(dD):.5f} formula={(2*m+1)/(m*(m+2)):.5f} V={float(var):.4f} V*delta_D={float(var*dD):.4f}")
