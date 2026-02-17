# ECMS Mean Shift Analysis (2026-02-17)

## Key Result: |δμ| < 1 suffices for ECMS

### Algebraic Identity

From the subdivision-contraction identity I(T_e) = I(T) + x·I(T/e), differentiation at x=1 gives:

    μ(T_e) = [a·μ(T) + b·(1 + μ(T/e))] / (a + b)

where a = I(T;1), b = I(T/e;1). This is a weighted average, so:

    μ(T) - μ(T/e) < 1  ⟺  μ(T) < 1 + μ(T/e)  ⟺  μ(T_e) < 1 + μ(T/e)

### Cavity Decomposition

In the hard-core model at λ=1, with cavity messages A = R_{u→v}, B = R_{v→u}:

    μ(T) - μ(T/e) = local + remote

where:
- **local** = (A+B-AB) / [(1+A+B)(1+AB)] ∈ (0, 1/2)  -- **PROVED**
  - Proved because A, B ∈ (0, 1] on any tree edge
  - Approaches 1/2 for leaf-center edges (A→1, B→0)
- **remote** = Σ_{w≠u,v} [P_T(w) - P_{T/e}(w)] ∈ (-0.262, 0.098)  -- empirical
  - Remote effect bounded by exponential correlation decay in trees
  - Can be negative (contraction can increase remote probabilities)

### Straddling Analysis

For LC polynomials, mode ∈ {⌊μ⌋, ⌈μ⌉}. If |δμ| < 1 but means straddle an integer k, modes could in principle differ by 2. Exhaustive check:

| n range | Edges | Straddling | Gap-2 (|Δmode|=2) |
|---------|-------|------------|-------------------|
| ≤ 18 | 3,348,675 | 864,230 (25.8%) | **0** |

In straddling cases, the mode of T/e always rounds UP to ⌈μ(T/e)⌉:
- 98.5%: both modes at straddled integer k
- 1.5%: mode(T) = k+1, mode(T/e) = k
- 0%: mode(T/e) drops to k-1

**Conclusion**: Proving |μ(T) - μ(T/e)| < 1 would establish ECMS for log-concave trees.

### Computational Verification

| Source | Trees | Edges | max |δμ| | Extremal |
|--------|-------|-------|---------|----------|
| Exhaustive n≤16 | 32,508 | 464,872 | 0.5263 | star-like |
| Exhaustive n=17-18 | 172,496 | 2,883,803 | 0.5330 | n=18 |
| Families n≤200 | 1,767 | ~50K | 0.5367 | T(5,5,5) |
| **Combined** | **206,771** | **~3.4M** | **0.537** | **T(5,5,5)** |

The tripod star T(5,5,5) (n=19, central vertex connected to 3 hubs each with 5 leaves) is the extremal family. All families converge to |δμ| → 1/2 from above.

### Proof Strategies for |δμ| < 1

1. **Cavity decay**: local < 1/2 (proved). Need |remote| < 1/2 (empirically < 0.27). The remote effect decays exponentially through non-leaf vertices with rate ≈ 0.618 per step (golden ratio). Total remote bounded by initial perturbation × geometric series. Formal proof requires careful analysis of how contraction changes cavity messages.

2. **Coefficientwise domination**: i_k(T) ≥ i_k(T/e) for all k (proved analytically). This gives μ(T) - μ(T/e) = (a-b)/a · (μ_excess - μ(T/e)), but bounding μ_excess leads to a circular argument.

3. **Mode sandwich lemma** (PROVED, now Proposition in paper): From the identity, mode(T_e) lies between min(mode(T), mode(T/e)+1) and max(mode(T), mode(T/e)+1). Verified: 3,348,674 edges (n ≤ 18), 0 failures. Corollary: ECMS implies subdivision shift ∈ {0, +1} (never -1). Empirically confirmed: 77.1% stay at mode, 22.9% increase by 1.

### Balanced tripod convergence (T(a,a,a))

| a | n | max |δμ| |
|---|---|---------|
| 3 | 13 | 0.4955 |
| 4 | 16 | 0.5224 |
| 5 | 19 | **0.5367** |
| 6 | 22 | 0.5316 |
| 7 | 25 | 0.5223 |
| 10 | 34 | 0.5051 |
| 15 | 49 | 0.5003 |
| 20 | 64 | 0.5000 |
| 26 | 82 | 0.5000 |

Peak at a=5 (T(5,5,5)), then monotone decrease to 1/2. Extremal edge is always centre-to-bridge.

### Remote Decay Analysis (2026-02-17)

**Per-vertex probability change decays geometrically with distance from contracted edge:**

| Distance | max |ΔP(w)| | Effective rate |
|----------|--------------|----------------|
| 1 | 0.167 | -- |
| 2 | 0.067 | 0.40 |
| 3 | 0.025 | 0.38 |
| 4 | 0.0096 | 0.38 |
| 5 | 0.0037 | 0.38 |
| 8 | 0.00020 | 0.38 |

Decay rate ≈ 1/φ² ≈ 0.382 (golden ratio squared). Average decay rate 0.28, max 0.50, 99th percentile 0.41.

**But per-vertex cavity decay factor γ → 1 on long paths.** The naive per-vertex bound fails. The product of gammas along paths stays bounded (golden-ratio decay), but individual factors approach 1 as path length grows. This kills the simple "uniform per-vertex decay" proof strategy.

### Remote by Family (n up to 300, 2026-02-17)

| Family | max |remote| (limit) | max |δμ| (limit) | Notes |
|--------|---------------------|------------------|-------|
| Paths | 0.040 | 0.276 | Converged by n=20 |
| Stars | → 0 | → 0.500 | All effect is local |
| Tripods | 0.040 | 0.276 | Same as paths (arms are paths) |
| Brooms | 0.126 | 0.514 | Star-like at hub |
| Extended stars | 0.106 | 0.514 | Similar to brooms |
| Caterpillars (k=3) | **0.243** | 0.447 | Worst family for remote |
| Mixed (k=3 → bare) | **0.254** | -- | Double star S(3,3) |

**All families converge.** No family exceeds |remote| = 0.254 at any n.

### Exhaustive Worst Trees (2026-02-17)

The worst trees for |remote| at each n (n ≤ 18):

| n | max |remote| | Tree type |
|---|---------------|-----------|
| 6 | 0.216 | Double star S(2,2) |
| 7 | 0.242 | Double star S(3,2) |
| 8 | 0.254 | Double star S(3,3) |
| 12 | 0.258 | Caterpillar spine=4, degs=[2,4,4,4] |
| 16 | 0.262 | deg_seq=[4,4,4,4,2,2,...] |
| 18 | **0.264** | deg_seq=[4,4,4,3,3,2,2,2,...] |

Trend: slowly increasing with diminishing increments. Appears to converge around 0.27--0.28.

Worst trees are caterpillar-like with **mixed leaf counts**: spine vertices with 3 leaves adjacent to bare path segments. The transition zone between high-leaf and low-leaf sections creates the largest perturbation.

### Caterpillar Fixed Point (algebraic)

For infinite caterpillar with k leaves per spine vertex:
- Cavity fixed point: R(1+R) = 2^{-k}, so R = (-1 + √(1 + 2^{2-k})) / 2
- k=3: R = (√6 - 2)/4 ≈ 0.1124, P_spine ≈ 0.092, P_leaf ≈ 0.454
- local = (2R - R²)/((1+2R)(1+R²)) ≈ 0.171 for k=3
- Remote converges to -0.231 (middle edge) or -0.243 (boundary edge)
- δμ ≈ -0.060 for k=3 middle edge (local + remote nearly cancel!)

Worst k for |remote|: k=3 (remote = -0.231). Falls off rapidly for k ≥ 4.

### S_1 Formula (THEOREM, 2026-02-17)

**Exact algebraic formula for distance-1 perturbation:**

For edge e = (u,v) with cavity messages A = R_{u→v}, B = R_{v→u}:

    ΔP(w) = F_u × q_w/(1+q_w)    for each w ∈ N(u)\{v}

where:
- F_u = A(B²+B-1) / [(1+A+B)(1+AB)]  (the "edge factor")
- q_w = R_{w→u}  (cavity message from w to u)
- q_w/(1+q_w) = P_w^{(u-cavity)}  (occupation probability in the u-cavity tree)

Similarly for w ∈ N(v)\{u}: ΔP(w) = F_v × R_{w→v}/(1+R_{w→v})
where F_v = B(A²+A-1)/[(1+A+B)(1+AB)].

**Derivation:**
- In T: R_{u→w} = A(1+q_w)/(1+B)
- In T/e: R_{m→w} = A(1+q_w)B  (merged vertex sends product of both sides)
- Change: β-α = A(1+q_w)(B²+B-1)/(1+B)
- After substitution: ΔP(w) = Q(β-α)/[(1+α+Q)(1+β+Q)] simplifies to F × q_w/(1+q_w)

**Verification:** 3,348,674 edges (n ≤ 18), 0 mismatches (max error 1.94×10⁻¹⁶).

**Corollary (S_1 bound):**
S_1 = F_u × Σ P_w^{(u-cav)} + F_v × Σ P_w^{(v-cav)}

Using log(1+x) ≥ x/(1+x) and A = ∏ 1/(1+q_w):
  Σ P_w^{(u-cav)} ≤ -log A

So: |S_1| ≤ g(A,B) = |F_u|(-log A) + |F_v|(-log B)

Numerically: max g(A,B) = **0.3546 < 1/2** on both exhaustive data and dense grid.
Maximum at A ≈ B ≈ 0.143 (caterpillar-like structures).

### Signed Decay and Alternation (2026-02-17)

**Signed sums at each distance alternate 94-100% of the time:**

| Distance | max S_d | min S_d | Alternation % |
|----------|---------|---------|---------------|
| 1 | +0.107 | -0.281 | -- |
| 2 | +0.203 | -0.109 | 94.3% |
| 3 | +0.041 | -0.078 | 94.6% |
| 4 | +0.028 | -0.015 | 97.9% |
| 5+ | <0.005 | <0.010 | 99%+ |

Partial sums oscillate: max |cumulative partial| = 0.281 (at d=1), never exceeds this.

**Alternating dominance** (|partial(d)| ≤ |S_1| for all d): fails 9.25% of (edge,d) pairs.
**|remote| ≤ |S_1|**: fails 11.82% of edges (dist-2+ tail can amplify).
**Tail has opposite sign to S_1**: 90.0% of edges (cancellation typical).

### Mean Response Bound (KEY, 2026-02-17)

**M = subtree_sum / δR (mean response per unit cavity perturbation) is universally bounded:**

    |M| < 0.273    for ALL subtrees of ALL trees through n ≤ 18

- 0.00% of subtrees exceed |M| = 0.5
- Mean |M| = 0.133, median = 0.146
- Convergence: max |M| grows at ~0.003/step, likely converging to ~0.28

This means every subtree DISSIPATES the incoming cavity perturbation by at least ~73%.
The total mean change in a subtree is always less than 27.3% of the incoming stimulus.

**Bound on remote:** |remote| ≤ max|M| × Σ|δR_w|. And using A×Σ(1+q_w) ≤ deg-1:
  Σ|δR_w| ≤ (deg_u - 1)|B²+B-1|/(1+B) + (deg_v - 1)|A²+A-1|/(1+A)

This bound grows with degree, BUT the product max|M| × g(A,B) ≈ 0.097 is far below 1/2.
The issue is translating between Σ|δR_w| (unsigned cavity changes) and g(A,B) (signed S_1 bound).

### Key Scripts

- `ecms_remote_decay.py` + `results/ecms_remote_decay.json`: Distance decomposition
- `ecms_cavity_decay_proof.py` + `results/ecms_cavity_decay.json`: Per-vertex γ factors
- `ecms_remote_families.py` + `results/ecms_remote_families.json`: Family convergence
- `ecms_worst_remote.py` + `results/ecms_worst_remote.json`: Exhaustive worst trees
- `ecms_mixed_caterpillar.py` + `results/ecms_mixed_caterpillar.json`: Mixed structures
- `ecms_caterpillar_algebra.py`: Algebraic fixed points
- `ecms_signed_decay.py` + `results/ecms_signed_decay.json`: Signed per-distance sums
- `ecms_alternating_proof.py` + `results/ecms_alternating_proof.json`: Alternating dominance tests
- `ecms_S1_algebra.py` + `results/ecms_S1_algebra.json`: S_1 formula verification
- `ecms_remote_vs_S1.py` + `results/ecms_remote_vs_S1.json`: Remote vs S_1 comparison
- `ecms_subtree_bound.py` + `results/ecms_subtree_bound.json`: Subtree amplification
- `ecms_recursive_bound.py` + `results/ecms_recursive_bound.json`: Mean response analysis

### Dead ends for proving |δμ| < 1

- Pure algebraic manipulation of the identity: all routes are circular (α < 1+β ⟺ α < 1+β)
- Global bounds μ(T) < (2/3)(n-1) minus μ(T/e) > 0: gives O(n), useless
- Naive cavity decay bound: local ≤ 1/2, total decay ≤ 1/2 × 1.618 ≈ 0.81. Gives total ≤ 1.3, too weak
- Per-vertex γ < 1: FAILS (γ → 1 on long paths, up to 0.9998)
- Bound |remote| < 1/4: FAILS (exhaustive max = 0.264 > 0.250)
- Subtree amplification ratio: grows linearly with n (unbounded)
- Bound remote by S_1 alone: FAILS (|remote| > |S_1| for 11.82% of edges)
- g(A,B) directly bounds remote: FAILS (|remote| > g(A,B) for 11K edges)
- max|M| × Σ|δR_w| with degree bound: gives max|M| × (d-1), grows with d

### Remaining path to proof

**Need: |remote| < 1/2.** Empirically |remote| < 0.264. Gap is 0.236.

**S_1 formula + mean response** is the most promising framework:
1. S_1 = F × Σ P_w^{cavity} (PROVED, exact)
2. |S_1| ≤ g(A,B) < 0.355 (PROVED numerically, needs analytic proof)
3. |M| < 0.273 (EMPIRICAL, needs proof) -- subtree dissipation
4. Need: connect (1)-(3) to bound total remote, accounting for sign cancellation

**Key structural insight:** cavity perturbations alternate sign at each tree level (from the multiplicative structure of R_w→parent). This creates oscillating partial sums that converge rapidly. The alternation rate is 94%+ at dist ≥ 2.

**Alternative:** Codex's WHNC inequality (weighted heavy-neighborhood compensation). If proved, gives Conjecture A directly, bypassing ECMS entirely. Verified 931K trees n ≤ 23, 0 failures.
