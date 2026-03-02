#!/usr/bin/env python3
"""
Test: LC(E) + LC(J) + J<=E + E≽J prefix  ==>  STP2(E,J)?
FINAL VERSION: strictly positive entries (all >= 1).
Reduced exhaustive range; bigger random phase.
"""
import random, time, sys

def is_lc(s):
    for k in range(1, len(s)-1):
        if s[k]*s[k] < s[k-1]*s[k+1]: return False
    return True

def smode(s):
    mx = max(s)
    for i,v in enumerate(s):
        if v == mx: return i
    return 0

def pfx_dom(E, J):
    m = smode(E); n = min(len(E), len(J))
    for k in range(min(m, n-1)):
        if k+1>=n: break
        if E[k+1]*J[k] < E[k]*J[k+1]: return False
    return True

def stp2(E, J):
    n = min(len(E), len(J)); f = []
    for k in range(1, n-1):
        if J[k+1]*E[k-1] > J[k]*E[k]: f.append(k)
    return f

def show(E, J, fails, idx):
    m = smode(E)
    print(f"  #{idx}: E={E}, J={J}, mode(E)={m}")
    for k in fails:
        l=J[k+1]*E[k-1]; r=J[k]*E[k]
        where = "TAIL" if k>m else ("PFX" if k<m else "MODE")
        print(f"       k={k}[{where}]: J[{k+1}]*E[{k-1}]={J[k+1]}*{E[k-1]}={l} > J[{k}]*E[{k}]={J[k]}*{E[k]}={r} (Δ={l-r})")
    frats = [f"{J[k]/E[k]:.4f}" for k in range(len(E))]
    print(f"       J/E: {frats}")
    sys.stdout.flush()

def exh4(ma=6):
    print(f"=== EXHAUSTIVE+ len=4, a≤{ma} ===", flush=True)
    tot=p=sf=0; ex=[]
    for a in range(1,ma+1):
        for b in range(1, a*a+1):
            bc = b*b//a
            for c in range(1, bc+1):
                E=[1,a,b,c]
                for j1 in range(1,a+1):
                    for j2 in range(1,b+1):
                        if j1*j1<j2: continue
                        for j3 in range(1,c+1):
                            if j2*j2<j1*j3: continue
                            J=[1,j1,j2,j3]; tot+=1
                            if not pfx_dom(E,J): continue
                            p+=1
                            f=stp2(E,J)
                            if f:
                                sf+=1
                                if len(ex)<20: ex.append((list(E),list(J),f))
    print(f"  pairs={tot:,} passed={p:,} STP2_fails={sf}", flush=True)
    return p,sf,ex

def exh5(ma=3):
    print(f"\n=== EXHAUSTIVE+ len=5, a≤{ma} ===", flush=True)
    tot=p=sf=0; ex=[]
    for a in range(1,ma+1):
        for b in range(1, a*a+1):
            bc=b*b//a
            for c in range(1, bc+1):
                bd=c*c//b if b>0 else 0
                for d in range(1, bd+1):
                    E=[1,a,b,c,d]
                    if not is_lc(E): continue
                    for j1 in range(1,a+1):
                        for j2 in range(1,b+1):
                            if j1*j1<j2: continue
                            for j3 in range(1,c+1):
                                if j2*j2<j1*j3: continue
                                for j4 in range(1,d+1):
                                    if j3*j3<j2*j4: continue
                                    J=[1,j1,j2,j3,j4]; tot+=1
                                    if not pfx_dom(E,J): continue
                                    p+=1
                                    f=stp2(E,J)
                                    if f:
                                        sf+=1
                                        if len(ex)<20: ex.append((list(E),list(J),f))
        print(f"  a={a}/{ma}, pairs={tot:,}, passed={p:,}, fails={sf}", flush=True)
    print(f"  TOTAL: pairs={tot:,} passed={p:,} STP2_fails={sf}", flush=True)
    return p,sf,ex

def gen_pos_lc(length, mx):
    E=[1, random.randint(2,mx)]
    for i in range(2,length):
        bnd=E[i-1]*E[i-1]//E[i-2]
        if bnd<1: return E
        E.append(random.randint(1,bnd))
    return E

def rnd_phase(N):
    print(f"\n=== RANDOM POSITIVE ({N:,} trials) ===", flush=True)
    p=sf=0; ex=[]
    for t in range(N):
        ln=random.randint(4,10)
        mx=random.choice([5,10,30,100,500,2000,10000])
        E=gen_pos_lc(ln,mx)
        if len(E)<4: continue
        mE=smode(E)
        st=random.choice(["b","np","st","md"])
        J=[1]
        for k in range(1,len(E)):
            ek=E[k]
            if st=="b":
                J.append(random.randint(1,ek))
            elif st=="np":
                if k<=mE and E[k-1]>0:
                    bnd=max(1,min(ek*J[k-1]//E[k-1],ek))
                    J.append(random.randint(max(1,bnd-max(1,bnd//5)),bnd))
                else:
                    J.append(random.randint(max(1,ek*2//3),ek))
            elif st=="st":
                if k<=mE and E[k-1]>0:
                    bnd=max(1,min(ek*J[k-1]//E[k-1],ek))
                    J.append(random.randint(max(1,bnd//2),bnd))
                else:
                    J.append(random.randint(max(1,ek*3//4),ek))
            elif st=="md":
                if k<=mE:
                    J.append(random.randint(max(1,ek-3),ek))
                else:
                    J.append(random.randint(max(1,ek//2),ek))
        if len(J)!=len(E): continue
        if not is_lc(J): continue
        if any(J[k]>E[k] for k in range(len(E))): continue
        if not pfx_dom(E,J): continue
        p+=1
        f=stp2(E,J)
        if f:
            sf+=1
            if len(ex)<30: ex.append((list(E),list(J),f))
        if t%2_000_000==0 and t>0:
            print(f"  ...{t:,} trials, {p:,} passed, {sf} fails", flush=True)
    print(f"  passed={p:,} STP2_fails={sf}", flush=True)
    return p,sf,ex

def sharp_drop_phase(N):
    print(f"\n=== SHARP-DROP E ({N:,} trials) ===", flush=True)
    p=sf=0; ex=[]
    for _ in range(N):
        ln=random.randint(5,10)
        peak=random.randint(50,5000)
        mp=random.randint(1,ln-2)
        E=[1]
        for k in range(1,mp+1):
            if k==1: E.append(random.randint(2,peak))
            else:
                bnd=E[k-1]*E[k-1]//E[k-2]
                lo=E[k-1]+1
                if lo>bnd: break
                E.append(random.randint(lo,bnd))
        if len(E)<mp+1: continue
        for k in range(mp+1,ln):
            bnd=E[k-1]*E[k-1]//E[k-2]
            if bnd<1: break
            hi=max(1,min(bnd,E[k-1]//2))
            E.append(random.randint(1,hi))
        if len(E)<4 or not is_lc(E): continue
        mE=smode(E)
        J=[1]
        for k in range(1,len(E)):
            if k<=mE and E[k-1]>0:
                bnd=max(1,min(E[k]*J[k-1]//E[k-1],E[k]))
                J.append(random.randint(max(1,bnd*3//4),bnd))
            else:
                J.append(random.randint(max(1,E[k]*4//5),E[k]))
        if len(J)!=len(E) or not is_lc(J): continue
        if any(J[k]>E[k] for k in range(len(E))): continue
        if not pfx_dom(E,J): continue
        p+=1
        f=stp2(E,J)
        if f:
            sf+=1
            if len(ex)<20: ex.append((list(E),list(J),f))
    print(f"  passed={p:,} STP2_fails={sf}", flush=True)
    return p,sf,ex

def main():
    random.seed(42); t0=time.time()
    print("="*70)
    print("TESTING (POSITIVE): {LC(E)+LC(J)+J<=E+E≽J pfx} => STP2(E,J)?")
    print("="*70, flush=True)

    AE=[]; GP=GF=0

    for func,args in [(exh4,(6,)),(exh5,(3,))]:
        p,f,ex = func(*args)
        GP+=p; GF+=f; AE.extend(ex)

    p,f,ex = rnd_phase(10_000_000)
    GP+=p; GF+=f; AE.extend(ex)

    p,f,ex = sharp_drop_phase(5_000_000)
    GP+=p; GF+=f; AE.extend(ex)

    el=time.time()-t0
    print("\n"+"="*70)
    print("GRAND SUMMARY (ALL POSITIVE)")
    print(f"  Passed all 4 conditions: {GP:,}")
    print(f"  STP2 failures:           {GF:,}")
    if GP>0: print(f"  Failure rate:            {GF/GP*100:.4f}%")
    print(f"  Elapsed:                 {el:.1f}s")

    if AE:
        AE.sort(key=lambda x:(len(x[0]),sum(x[0])))
        uniq=[]; seen=set()
        for E,J,f in AE:
            k=(tuple(E),tuple(J))
            if k not in seen: seen.add(k); uniq.append((E,J,f))
        print(f"\n  COUNTEREXAMPLES ({len(uniq)} unique, first 25):\n")
        for i,(E,J,f) in enumerate(uniq[:25],1):
            show(E,J,f,i)
        # location
        pc=tc=mc=0
        for E,J,fs in uniq:
            m=smode(E)
            for k in fs:
                if k<m: pc+=1
                elif k==m: mc+=1
                else: tc+=1
        print(f"\n  FAILURE LOCATIONS ({len(uniq)} unique):")
        print(f"    PREFIX (k<mode): {pc}")
        print(f"    AT MODE (k=mode): {mc}")
        print(f"    TAIL (k>mode): {tc}")
    else:
        print("\n  NO COUNTEREXAMPLES. Implication holds for positive sequences.")
    print("="*70, flush=True)

if __name__=="__main__": main()
