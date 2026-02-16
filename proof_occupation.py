#!/usr/bin/env python3
"""WHY does gap=+1 force excess(dp₀) ≤ 0?

SETUP (single child c):
  dp₀[v] = total[c]
  dp₁[v] = x · dp₀[c]
  total[v] = total[c] + x · dp₀[c]

Gap=+1 means: mode(total[v]) = mode(total[c]) + 1.
Let M = mode(total[c]).

total[v][M+1] ≥ total[v][M]
total[c][M+1] + dp₀[c][M] ≥ total[c][M] + dp₀[c][M-1]

Since M = mode(total[c]): total[c][M] ≥ total[c][M+1]
So: dp₀[c][M] - dp₀[c][M-1] ≥ total[c][M] - total[c][M+1] ≥ 0

I.e., dp₀[c] is INCREASING at position M: dp₀[c][M] > dp₀[c][M-1].
So mode(dp₀[c]) ≥ M+1 (since dp₀[c] hasn't peaked yet at M).

Actually, let me be precise: dp₀[c][M] > dp₀[c][M-1] means the 
sequence is still increasing going from M-1 to M. This means 
mode(dp₀[c]) ≥ M. Actually it means mode(dp₀[c]) > M-1, so ≥ M.

But wait: dp₀[c][M] > dp₀[c][M-1] just says position M is larger 
than M-1. The mode could still be at M (if dp₀[c] peaks there).

OK so we have: mode(dp₀[c]) ≥ M.

Now total[c] = dp₀[c] + dp₁[c], with mode M.
dp₀[c] has mode ≥ M.
dp₁[c] = x · stuff, so dp₁[c][k] is shifted up by 1.

Since total[c] = dp₀[c] + dp₁[c] has mode M, and dp₀[c] has mode ≥ M:
The mode didn't move or went DOWN when adding dp₁[c].

total[c][M] = dp₀[c][M] + dp₁[c][M] (maximum)
total[c][M+1] = dp₀[c][M+1] + dp₁[c][M+1] ≤ total[c][M]

If mode(dp₀[c]) = M: then dp₀[c][M] ≥ dp₀[c][M+1], and dp₁[c] 
just doesn't push the mode up.

If mode(dp₀[c]) > M: then dp₀[c][M] < dp₀[c][M+1] (going up still).
But total[c][M+1] ≤ total[c][M], so:
dp₀[c][M+1] + dp₁[c][M+1] ≤ dp₀[c][M] + dp₁[c][M]
dp₀[c][M+1] - dp₀[c][M] ≤ dp₁[c][M] - dp₁[c][M+1]
LHS > 0 (dp₀ increasing), so dp₁[c][M] > dp₁[c][M+1].
This means dp₁[c] is DECREASING at M+1, so mode(dp₁[c]) ≤ M.

INTERESTING: dp₁[c] has mode ≤ M while dp₀[c] has mode ≥ M.
The two "fight" at position M, with dp₁[c] pushing mass below M.

NOW: the key is to show excess(dp₀[v]) = excess(total[c]) ≤ 0.
excess(total[c]) = M - mean(total[c]).

mean(total[c]) = (|dp₀[c]| · mean(dp₀[c]) + |dp₁[c]| · mean(dp₁[c])) 
                 / (|dp₀[c]| + |dp₁[c]|)

Since mode(dp₀[c]) ≥ M and by induction excess(dp₀[c]) < 1:
  mean(dp₀[c]) > mode(dp₀[c]) - 1 ≥ M - 1.

Since mode(dp₁[c]) ≤ M and dp₁[c] = x · stuff:
  mean(dp₁[c]) ≥ 1 (at least 1, since dp₁[c][0] = 0).
  Actually, mean(dp₁[c]) = 1 + mean(P_c) where P_c = ∏ dp₀[gc].
  So mean(dp₁[c]) ≥ 1.

So: mean(total[c]) = weighted average of mean(dp₀[c]) and mean(dp₁[c])
                   ≥ weighted average of (M-1) and 1... this isn't tight enough.

BETTER: Let w = |dp₁[c]| / |total[c]| (weight of dp₁[c] in total).
mean(total[c]) = (1-w)·mean(dp₀[c]) + w·mean(dp₁[c])

For mean(total[c]) ≥ M:
(1-w)·mean(dp₀[c]) + w·mean(dp₁[c]) ≥ M

Since mean(dp₀[c]) ≥ M-1 (by induction, excess < 1):
(1-w)(M-1) + w·mean(dp₁[c]) ≥ M
M - 1 - wM + w + w·mean(dp₁[c]) ≥ M
w·(mean(dp₁[c]) + 1 - M) ≥ 1
w ≥ 1/(mean(dp₁[c]) + 1 - M)

Hmm, for M large and mean(dp₁[c]) ≈ M (roughly), this gives w ≥ 1/1 = 1,
which is too strong.

Actually wait, mean(dp₁[c]) = 1 + mean(P_c). And P_c = ∏ dp₀[grandchildren].
mean(P_c) can be comparable to M for large trees.

Let me try a DIFFERENT approach. Instead of proving excess(total[c]) ≤ 0 
outright (which is ALMOST true but has 0.03% exceptions with excess up 
to 0.07), let me show that the correction ALWAYS exceeds excess(dp₀),
using the structural constraints more precisely.

THE CORRECTION at the parent merge:
  correction = w_parent · (1 + mean(dp₀[c]) - mean(total[c]))

where w_parent = |dp₀[c]| / (|total[c]| + |dp₀[c]|).

Hmm wait, for single child:
  dp₀[v] = total[c]
  dp₁[v] = x · dp₀[c]
  
  w = |x·dp₀[c]| / (|total[c]| + |dp₀[c]|) = |dp₀[c]| / (|total[c]| + |dp₀[c]|)
  
  μ_dp₁ = 1 + mean(dp₀[c])
  μ_dp₀ = mean(total[c])
  
  correction = w · (1 + mean(dp₀[c]) - mean(total[c]))

The condition excess < 1 requires:
  correction > excess(dp₀[v]) = excess(total[c]) = M - mean(total[c])

So: w · (1 + mean(dp₀[c]) - mean(total[c])) > M - mean(total[c])

Let μ = mean(total[c]). Then:
  w · (1 + mean(dp₀[c]) - μ) > M - μ

Now: mean(dp₀[c]) = (μ·|total[c]| - mean(dp₁[c])·|dp₁[c]|) / |dp₀[c]|
Wait, total[c] = dp₀[c] + dp₁[c], so:
  μ = (|dp₀[c]|·mean(dp₀[c]) + |dp₁[c]|·mean(dp₁[c])) / |total[c]|

And mean(dp₁[c]) = 1 + mean(P_c).

Let me define:
  a = |dp₀[c]|, b = |dp₁[c]|, s = a + b = |total[c]|
  α = mean(dp₀[c]), β = 1 + mean(P_c) = mean(dp₁[c])
  μ = (aα + bβ)/s

The condition becomes:
  w · (1 + α - μ) > M - μ
  
  w = a/(s + a) = a/(2a + b)  [since s = a+b, dp₀[c] has total a]
  
  1 + α - μ = 1 + α - (aα + bβ)/s = 1 + α(1 - a/s) - bβ/s
            = 1 + αb/s - bβ/s = 1 + b(α - β)/s = 1 - b(β - α)/s
  
  Wait, if β > α (dp₁ has higher mean): 1 + α - μ = 1 - b(β-α)/s < 1.
  If β < α: 1 + α - μ > 1.
  
  But β = 1 + mean(P_c) and α = mean(dp₀[c]).
  β - α = 1 + mean(P_c) - mean(dp₀[c]).
  
  Now at the CHILD: mode(dp₀[c]) ≥ M and mode(dp₁[c]) ≤ M.
  And total[c] has mode M.
  
  This isn't leading to a clean bound. Let me try the PROBABILISTIC 
  interpretation instead.

PROBABILISTIC VIEW:
  I(T; λ) = Σ_S λ^|S| where S runs over all IS.
  At λ=1: each IS gets weight 1.
  The IS-size distribution has mean μ = I'(1)/I(1) and mode M.
  
  CLAIM: μ ≥ M - 1, i.e., the expected IS size is within 1 of the mode.
  
  Intuitively: the IS distribution is "concentrated" enough that 
  the mode can't be far from the mean.
  
  For trees, the IS distribution is a CONVOLUTION of local contributions 
  (from the belief propagation). Each vertex contributes roughly 
  P(v) ≈ λ/(1+λ) to the mean, and the mode of each local contribution 
  is 0 or 1. The sum of these has mode close to the sum of their means.

Let me just investigate what P(v) values look like at the tightest cases
and see if there's a pattern we can exploit.
"""

import subprocess
from graph6 import parse_graph6

def poly_mode(p):
    return max(range(len(p)), key=lambda k: p[k])

def poly_mean(p):
    s = sum(p)
    if s == 0: return 0
    return sum(k*p[k] for k in range(len(p))) / s

def polymul(a, b):
    n = len(a) + len(b) - 1
    r = [0]*n
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            r[i+j] += ai*bj
    return r

def polyadd(a, b):
    n = max(len(a), len(b))
    r = [0]*n
    for i in range(len(a)): r[i] += a[i]
    for i in range(len(b)): r[i] += b[i]
    return r

def compute_is_poly(adj, n):
    """Compute IS polynomial using subgraph removal."""
    # Simple recursive with memoization on vertex subsets
    from indpoly import independence_poly
    return independence_poly(n, adj)

def main():
    print("=" * 70)
    print("  OCCUPATION PROBABILITIES AT TIGHTEST CASES")
    print("=" * 70)
    print()
    
    # For each tree, compute the IS poly, find tie points, 
    # and compute occupation probabilities.
    
    from indpoly import independence_poly
    
    tightest_cases = []
    
    for n in range(5, 18):
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        lines = [l for l in proc.stdout.split(b"\n") if l]
        
        for line in lines:
            tn, adj = parse_graph6(line)
            poly = independence_poly(tn, adj)
            
            mode = poly_mode(poly)
            mean = poly_mean(poly)
            excess = mode - mean
            
            if excess > 0.5:
                tightest_cases.append((excess, n, poly, adj, line))
    
    tightest_cases.sort(reverse=True)
    
    print(f"  Top 10 tightest cases (highest excess at λ=1):")
    print(f"  {'excess':>7} {'n':>3} {'mode':>4} {'mean':>7} {'poly'}")
    
    for excess, n, poly, adj, line in tightest_cases[:10]:
        mode = poly_mode(poly)
        mean = poly_mean(poly)
        # Compute occupation prob of each vertex
        # P(v) = derivative contribution
        # mean = Σ_v P(v) where P(v) = i_v(T;1)/I(T;1) and 
        # i_v = (independence poly of T with v forced in) = I(T-N[v];1)
        
        # Compute P(v) for each v
        pvs = []
        for v in range(n):
            # I(T-N[v]) at λ=1
            # N[v] = v and its neighbors
            nv = {v} | set([u for u in range(n) if v in adj[u] or u in adj[v]])
            # Actually adj is adjacency list from parse_graph6
            # adj[v] is the list of neighbors of v
            nv = {v} | set(adj[v])
            remaining = [u for u in range(n) if u not in nv]
            
            if not remaining:
                # T - N[v] is empty, I = 1
                pv = 1.0 / sum(poly)
            else:
                # Build subgraph
                remap = {u: i for i, u in enumerate(remaining)}
                new_n = len(remaining)
                new_adj = [[] for _ in range(new_n)]
                for u in remaining:
                    for w in adj[u]:
                        if w in remap:
                            new_adj[remap[u]].append(remap[w])
                
                sub_poly = independence_poly(new_n, new_adj)
                pv = sum(sub_poly) / sum(poly)
            
            pvs.append(pv)
        
        print(f"  {excess:7.4f} {n:3d} {mode:4d} {mean:7.4f} P(v)_sorted={sorted(pvs, reverse=True)[:5]}")
    
    print()
    
    # KEY OBSERVATION: If ALL P(v) are in (0, 1/2], then μ = Σ P(v) ≤ n/2.
    # And mode ≤ α(T) where α(T) is the IS number.
    # For trees: α(T) ≈ n/2 (since trees are bipartite).
    
    # But more importantly: for the HARD-CORE MODEL at general λ,
    # the BP recursion gives exact P(v) values. These satisfy:
    # P(v) = λ / (1 + λ + λ·Σ_{u∈N(v)} P(u)/(1-P(u)))... 
    # Actually the exact recursion for trees is:
    # 
    # For a tree rooted at root, DFS from leaves:
    # P(leaf) = λ/(1+λ)
    # P(v) = λ · ∏_{children c of v} (1 - P(c)) / 
    #         (∏_{children c} (1-P(c)) + λ · ∏_{children c} (1-P(c)))
    # Actually: P(v) = λ·B(v) / (A(v) + λ·B(v))
    # where A(v) = ∏ (A(c) + λB(c)) = ∏ Z(T_c) and B(v) = ∏ A(c) = ∏ Z_0(T_c)
    # Hmm, think of it differently.
    # 
    # Z_0(v) = ∏_c Z(c)  (v not in IS, children unconstrained)
    # Z_1(v) = λ · ∏_c Z_0(c)  (v in IS, children forced out)
    # Z(v) = Z_0(v) + Z_1(v)
    # P(v) = Z_1(v) / Z(v) = λ · ∏_c Z_0(c) / (∏_c Z(c) + λ · ∏_c Z_0(c))
    #       = λ / (∏_c Z(c)/Z_0(c) + λ) = λ / (∏_c (1 + Z_1(c)/Z_0(c)) + λ)
    
    # Since P(c) = Z_1(c)/Z(c), and Z_1(c)/Z_0(c) = P(c)·Z(c)/Z_0(c):
    # Hmm, let's define R(c) = Z_1(c)/Z_0(c) = P(c)/(1-P(c)) (odd ratio).
    # Then P(v) = λ / (∏_c (1 + R(c)) + λ)
    
    # For a leaf: R = λ, P = λ/(1+λ), 1+R = 1+λ.
    # For v with children c₁,...,cₖ:
    # P(v) = λ / (∏(1+R(cᵢ)) + λ)
    
    # The mean IS size = Σ_v P(v).
    # The mode... is harder to relate to the P(v)'s.
    
    # BUT: for a tree at fugacity λ=1:
    # P(v) = ∏_c Z_0(c) / (∏_c Z(c) + ∏_c Z_0(c))
    #       = ∏_c (Z_0(c)/Z(c)) / (1 + ∏_c (Z_0(c)/Z(c)))
    #       = ∏_c (1-P(c)) / (1 + ∏_c (1-P(c)))
    
    # Wait, that's at λ=1:
    # P(v) = 1 · ∏_c (1-P(c)) / (1 + ∏_c (1-P(c))) = ∏(1-P(c)) / (1+∏(1-P(c)))
    # So P(v) = Q(v) / (1+Q(v)) where Q(v) = ∏_c (1-P(c)).
    # P(v) < 1/2 iff Q(v) < 1 iff ∏(1-P(c)) < 1.
    # Since 0 < 1-P(c) < 1, Q(v) < 1 iff at least one child exists.
    # So P(v) < 1/2 for any non-leaf vertex! ✓
    # And P(leaf) = 1/(1+1) = 1/2 at λ=1.
    
    # So for trees at λ=1: ALL occupation probabilities are ≤ 1/2.
    # μ = Σ P(v) ≤ n/2.
    # And α(T) ∈ {⌈n/2⌉} for paths, up to n-1 for stars.
    
    # Hmm, this doesn't directly constrain mode vs mean.
    
    print("  INSIGHT: For trees at λ=1, P(v) ≤ 1/2 for all v.")
    print("  Mean = Σ P(v) ≤ n/2.")
    print()
    print("=" * 70)

if __name__ == "__main__":
    main()
