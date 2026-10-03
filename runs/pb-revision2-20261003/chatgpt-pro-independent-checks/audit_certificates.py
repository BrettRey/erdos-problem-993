#!/usr/bin/env python3
"""Independent Section 4 audit, rebuilt from PDF formulas; Python stdlib only.
No supplied supplement was available or used. Exact fractions and sparse polynomials.
"""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import json

# Four indeterminates H,J,t,u. Sparse monomial maps to exact rational coefficient.
NV=4
ZERO=(0,)*NV
class P:
    def __init__(self, terms=0):
        if isinstance(terms,P): self.a=dict(terms.a)
        elif isinstance(terms,dict): self.a={e:F(c) for e,c in terms.items() if c}
        else: self.a={} if not terms else {ZERO:F(terms)}
    def __add__(self, other):
        other=P(other); out=dict(self.a)
        for e,c in other.a.items():
            out[e]=out.get(e,F(0))+c
            if not out[e]: del out[e]
        return P(out)
    __radd__=__add__
    def __neg__(self): return P({e:-c for e,c in self.a.items()})
    def __sub__(self, other): return self+-P(other)
    def __rsub__(self,other): return P(other)+-self
    def __mul__(self,other):
        other=P(other); out={}
        for e,c in self.a.items():
            for f,d in other.a.items():
                g=tuple(x+y for x,y in zip(e,f))
                out[g]=out.get(g,F(0))+c*d
        return P(out)
    __rmul__=__mul__
    def __truediv__(self, other): return self*F(1,other)
    def __pow__(self,n):
        assert isinstance(n,int) and n>=0
        out=P(1); a=self
        while n:
            if n&1: out=out*a
            a=a*a; n//=2
        return out
    def __eq__(self,other): return self.a==P(other).a
    def degree(self,index): return max((e[index] for e in self.a),default=-1)
    def coeff(self,index,n):
        out={}
        for e,c in self.a.items():
            if e[index]==n:
                f=list(e); f[index]=0; f=tuple(f)
                out[f]=out.get(f,F(0))+c
        return P(out)
    def scalar(self):
        assert not self.a or set(self.a)=={ZERO},self.a
        return self.a.get(ZERO,F(0))
    def subs(self,replacements):
        vals=[P(replacements.get(i,X[i])) for i in range(NV)]
        pows=[[P(1)] for _ in range(NV)]
        for i in range(NV):
            for n in range(1,self.degree(i)+1): pows[i].append(pows[i][-1]*vals[i])
        out=P(0)
        for e,c in self.a.items():
            z=P(c)
            for i,n in enumerate(e): z=z*pows[i][n]
            out=out+z
        return out

def var(i):
    e=[0]*NV; e[i]=1
    return P({tuple(e):F(1)})
X=[var(i) for i in range(NV)]
H,J,t,u=X

def poly_coeff(coefficients,var):
    return sum((P(c)*var**i for i,c in enumerate(coefficients)),P(0))

def bernstein(g,index,d):
    assert g.degree(index)<=d
    beta=[sum((g.coeff(index,j)*F(comb(i,j),comb(d,j)) for j in range(i+1)),P(0)) for i in range(d+1)]
    # Reconstruct in the Bernstein basis; this catches a conversion implementation error.
    reconstruction=sum((beta[i]*comb(d,i)*X[index]**i*(1-X[index])**(d-i) for i in range(d+1)),P(0))
    assert reconstruction==g
    return beta

def coeffs(g,index,d=None):
    if d is None: d=g.degree(index)
    return [str(g.coeff(index,k).scalar()) for k in range(d+1)]

checks=[]
def check(name,condition):
    if not condition: raise AssertionError(name)
    checks.append(name)

def product(seq):
    z=P(1)
    for y in seq: z=z*y
    return z

results={'audit_basis':'Rebuilt independently from formulas in the revised PDF; no supplement supplied or run.',
         'arithmetic':'Python standard library fractions.Fraction; symbolic polynomial identities; no floating point or numerical grids.',
         'cells':{}}
Q4=(3*H+4)*(H+1)

# K=3. Common denominator Z=H^3(H+1)^2 for all seven weights.
Z=H**3*(H+1)**2
weights={0:Z}
for r in range(1,4):
    weights[-r]=product(H-s for s in range(1,r+1))*H**(3-r)*(H+1)**2
    weights[r]=(H-1)**r*product(H+1-j for j in range(1,r))*H**(3-r)*(H+1)**(3-r)
Snum=sum(weights.values(),P(0))
Tnum=sum((r*r*w for r,w in weights.items()),P(0))
Mnum=sum((r*w for r,w in weights.items()),P(0))
Anum=Snum*Tnum-Mnum**2
# Independently recover the pairwise expression from all unordered pairs.
pairnum=sum((weights[i]*weights[j]*(i-j)**2 for i in range(-3,4) for j in range(i+1,4)),P(0))
check('K=3 moment numerator equals direct pairwise numerator',Anum==pairnum)
Pprinted=poly_coeff([-912,4272,-7176,4696,800,-5460,6613,-3676,750,-16,-3],H)
check('Eq.4.6 numerator/denominator identity',4*Anum-Q4*Z**2==H*(H+1)*Pprinted)
check('Eq.4.6 denominator cancellation',Z**2==H*(H+1)*H**5*(H+1)**3)
betas=bernstein(Pprinted.subs({0:3+t}),2,10)
check('K=3 all eleven Bernstein coefficients strictly positive',all(b.scalar()>0 for b in betas))
results['cells']['3']={'degree':10,'count':len(betas),'bernstein':[str(b.scalar()) for b in betas],
                       'minimum':str(min(b.scalar() for b in betas))}

# For cell m, H^m Sm and H^m Tm have integer coefficients.
allbeta=[b.scalar() for b in betas]
for m in range(4,16):
    B=P(1); sn=H**m; tn=P(0)
    for r in range(1,m+1):
        B=B*(H-r)
        term=B*H**(m-r)
        sn=sn+2*term; tn=tn+2*r*r*term
    pm=4*sn*tn-Q4*H**(2*m)
    check(f'cell {m} numerator has integer coefficients',all(c.denominator==1 for c in pm.a.values()))
    check(f'cell {m} numerator degree {2*m+2}',pm.degree(0)==2*m+2)
    bs=bernstein(pm.subs({0:m+t}),2,2*m+2)
    check(f'cell {m} all {len(bs)} Bernstein coefficients strictly positive',all(b.scalar()>0 for b in bs))
    allbeta += [b.scalar() for b in bs]
    results['cells'][str(m)]={'degree':2*m+2,'count':len(bs),'bernstein':[str(b.scalar()) for b in bs],
                             'minimum':str(min(b.scalar() for b in bs)),
                             'H_power_coefficients_ascending':coeffs(pm,0)}
check('finite-cell total count is 275',len(allbeta)==275)

# Summation formulas checked as polynomial antidifferences, with base values.
sigma3=J**2*(J+1)**2/4
sigma4=J*(J+1)*(2*J+1)*(3*J**2+3*J-1)/30
sumrr=J*(J+1)*(J+2)/3
sum2r2=J*(J+1)*(2*J+1)/3
for name,formula,increment in [('sigma3',sigma3,J**3),('sigma4',sigma4,J**4),
                              ('sum r(r+1)',sumrr,J*(J+1)),('sum 2r^2',sum2r2,2*J**2)]:
    check(name+' antidifference identity',formula-formula.subs({1:J-1})==increment)
    check(name+' zero base value',formula.subs({1:0})==0)
# H*S_tilde and H*T_tilde derived directly from the lambda sums.
HS=(2*J+1)*H-sumrr
HT=sum2r2*H-(sigma3+sigma4)
NJ=HS*HT-H**2*Q4/4
Nprinted=-F(3,4)*H**4-F(7,4)*H**3+(J*(J+1)*(2*J+1)**2/3-1)*H**2 \
        -((2*J+1)*(sigma3+sigma4)+J**2*(J+1)**2*(J+2)*(2*J+1)/9)*H \
        +J*(J+1)*(J+2)*(sigma3+sigma4)/3
check('general quartic NJ expansion in Section4.2',NJ==Nprinted)

b5=bernstein(NJ.subs({0:16+5*t,1:5}),2,4)
b5printed=[F(2360),F(7500),F(25055,2),F(254205,16),F(31115,2)]
check('J=5 all five printed quartic Bernstein coefficients', [b.scalar() for b in b5]==b5printed)
check('J=5 all five coefficients strictly positive',all(b.scalar()>0 for b in b5))
results['J5_bernstein']=[str(b.scalar()) for b in b5]

# Triangular cell in J, followed by the shift J=u+6.
bJ=bernstein(NJ.subs({0:J*(J+1)/2+(J+1)*t}),2,4)
mu=[J**2*(J+1)**2,J*(J+1)**2,(J+1)**2,(J+1)**2*(J+2),(J+1)**2*(J+2)**2]
pi_lists=[
 [25400,46014,17431,2474,121],
 [384510,480475,199405,37967,3442,121],
 [3829440,4668316,2172926,512764,65863,4410,121],
 [644400,671859,249961,43735,3684,121],
 [91440,88782,25579,2958,121]]
results['J_ge_6_factorizations']=[]
for i,cs in enumerate(pi_lists):
    lhs=bJ[i].subs({1:u+6})
    rhs=mu[i].subs({1:u+6})*poly_coeff(cs,u)/2880
    check(f'J>=6 beta_{i}=mu_{i} pi_{i}(u)/2880 identity',lhs==rhs)
    check(f'J>=6 pi_{i} all coefficients strictly positive',all(c>0 for c in cs))
    results['J_ge_6_factorizations'].append({'i':i,'pi_power_coefficients_ascending':cs,
                                          'beta_u_power_coefficients_ascending':coeffs(lhs,3)})

# Check shared finite-cell endpoints as algebraic S,T numerators, using zero bm.
for m in range(4,16):
    check(f'cell {m} added endpoint weight b_{m}({m}) vanishes',product(H-s for s in range(1,m+1)).subs({0:m})==0)

results['passed_checks']=checks
results['number_of_passed_checks']=len(checks)
results['total_finite_Bernstein_coefficients']=len(allbeta)
results['all_finite_Bernstein_coefficients_positive']=all(b>0 for b in allbeta)
results['smallest_finite_Bernstein_coefficient']=str(min(allbeta))
base=Path(__file__).with_suffix('')
base.with_suffix('.json').write_text(json.dumps(results,indent=2)+'\n')
log=[]
log.append('Independent exact-rational Section 4 certificate reconstruction: PASS')
log.append('No supplied supplement was available or used.')
log.append(f'Passed {len(checks)} exact algebraic assertions.')
log.append('Eq. (4.6): printed degree-10 numerator matches unsymmetrized K=3 A-Q exactly.')
log.append('Cells, degree, count, minimum positive Bernstein coefficient:')
for m,item in results['cells'].items():
    log.append(f'  {m}: degree={item["degree"]}, count={item["count"]}, min={item["minimum"]}')
log.append('Finite coefficients: 11+264=275; all strictly positive.')
log.append('J=5 coefficients: '+', '.join(results['J5_bernstein']))
log.append('J>=6: all five symbolic mu_i*pi_i(u)/2880 identities match; every pi_i coefficient positive.')
log.append('All Bernstein conversions also checked by exact symbolic reconstruction.')
log.append('All four sum formulas checked by exact polynomial finite differences and zero base case.')
log.append('Full exact coefficients and assertions: '+str(base.with_suffix('.json')))
text='\n'.join(log)+'\n'
base.with_suffix('.log').write_text(text)
print(text,end='')
