#!/usr/bin/env python3
"""Direct approach: bound excess(total) via coefficient ratios.

At a merge total = dp0 + dp1, with mode(total) at position M:
  total[M] ≥ total[M-1]  (M is the mode)
  i.e. dp0[M] + dp1[M] ≥ dp0[M-1] + dp1[M-1]

Since dp1[k] = dp1_shift[k-1] where dp1_shift = prod(dp0[c_i]):
  dp1[M] = dp1_shift[M-1]
  dp1[M-1] = dp1_shift[M-2]

So: dp0[M] + dp1_shift[M-1] ≥ dp0[M-1] + dp1_shift[M-2]

Rearranging: dp0[M] - dp0[M-1] ≥ dp1_shift[M-2] - dp1_shift[M-1]

If dp0 is "decreasing" at M (dp0[M] < dp0[M-1]), then dp1_shift must be 
"increasing enough" at M-1 to compensate.

But can we use this to bound the excess directly?

ALTERNATIVE APPROACH: The Cauchy-Davenport-like argument.

For two nonneg sequences a,b with supports [0,α₁] and [0,α₂]:
  mode(a+b) ∈ {mode(a), mode(a)+1, ..., mode(a)+mode(b)+1, ...}?
  Actually mode(a+b) is just wherever a+b peaks.

SIMPLEST APPROACH: Show that for the total IS poly, 
  mode ≤ floor(mean + 1-ε) for some ε depending on tree structure.

Let me try: for the IS poly of ANY graph (not just trees),
  we need i_M ≥ i_{M-1}  (M is the mode)
  and i_M ≥ i_{M+1}

The mean is: μ = Σk i_k / Σi_k

If mode M > ceil(μ), then M ≥ μ + 1, so there's a lot of "mass below M".
This means sum_{k≤M-1} i_k is large relative to sum_{k≥M} i_k.
But i_M is the largest coefficient...

KEY: For a discrete distribution on {0,...,d} with mode M:
  μ = Σk p_k where p_k = i_k / Σi_k
  M = argmax p_k

  μ = M + Σk (k-M) p_k = M + Σ_{k≠M} (k-M) p_k

  So μ - M = Σ_{k<M} (k-M) p_k + Σ_{k>M} (k-M) p_k
            = negative_part + positive_part

  μ - M = -Σ_{k<M} (M-k) p_k + Σ_{k>M} (k-M) p_k

  For μ - M > -1, we need:
    Σ_{k>M} (k-M) p_k > Σ_{k<M} (M-k) p_k - 1

  i.e., the right tail weighted sum exceeds the left tail weighted sum minus 1.

  Since p_M ≥ p_k for all k: Σ_{k<M} (M-k) p_k ≤ Σ_{k<M} (M-k) p_M 
                                                    = p_M * M(M-1)/2 ... no, 
                                                    = p_M * Σ_{j=1}^{M} j 
                                                    = p_M * M(M+1)/2
  
  This bound is too loose. Let me think differently.

Actually, the simplest observation: 
  Σ_{k=0}^{α} p_k = 1
  p_M = max(p_k)

  mean = Σk p_k
  Σ_{k=0}^{α} k p_k = mean
  
  Since p_M is the max: p_M ≥ 1/(α+1) (pigeonhole).
  
  Now: p_M + p_{M-1} ≥ p_M ≥ 1/(α+1)
  
  And the mass below M: Σ_{k<M} p_k = S_L
  The mass above M: Σ_{k>M} p_k = S_R
  Mass at M: p_M = 1 - S_L - S_R

  mean = Σ_{k<M} k p_k + M p_M + Σ_{k>M} k p_k
       ≥ 0 * S_L + M * p_M + (M+1) * S_R
       = M * p_M + (M+1) * S_R
       = M * (1 - S_L - S_R) + (M+1) * S_R
       = M - M*S_L - M*S_R + M*S_R + S_R
       = M - M*S_L + S_R

  So: mean ≥ M - M*S_L + S_R
  
  For excess = M - mean ≤ M*S_L - S_R
  We want excess < 1, so: M*S_L - S_R < 1
  i.e. M*S_L < 1 + S_R

  Hmm, this requires S_L < 1/M + S_R/M, which isn't guaranteed.

Let me just compute and track S_L, S_R for IS polys.
"""

import subprocess
import math
from indpoly import independence_poly
from graph6 import parse_graph6

def main():
    print("=" * 70)
    print("  MASS DECOMPOSITION: S_L, p_M, S_R")
    print("=" * 70)
    print()
    
    print(f"{'n':>3} {'excess':>7} {'M':>3} {'α':>3} {'S_L':>6} {'p_M':>6} {'S_R':>6} "
          f"{'M*S_L':>7} {'S_R':>6} {'bound':>7}")
    print("-" * 65)
    
    for n in range(5, 23):
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        lines = [l for l in proc.stdout.split(b"\n") if l]
        
        max_excess = -float('inf')
        best = None
        
        for line in lines:
            tn, adj = parse_graph6(line)
            poly = independence_poly(tn, adj)
            total = sum(poly)
            mean = sum(k * poly[k] for k in range(len(poly))) / total
            M = max(range(len(poly)), key=lambda k: poly[k])
            excess = M - mean
            
            if excess > max_excess:
                max_excess = excess
                alpha = len(poly) - 1
                pM = poly[M] / total
                SL = sum(poly[k] for k in range(M)) / total
                SR = sum(poly[k] for k in range(M+1, len(poly))) / total
                best = (M, alpha, SL, pM, SR)
        
        M, alpha, SL, pM, SR = best
        bound = M * SL - SR  # excess ≤ M*SL - SR
        print(f"{n:3d} {max_excess:7.4f} {M:3d} {alpha:3d} {SL:6.4f} {pM:6.4f} {SR:6.4f} "
              f"{M*SL:7.4f} {SR:6.4f} {bound:7.4f}")
    
    print()
    print("The bound M*S_L - S_R is always < 1? Let's check:")
    print("(This would give a proof via mass decomposition)")

if __name__ == "__main__":
    main()
