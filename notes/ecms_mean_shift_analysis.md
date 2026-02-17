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

### Key Scripts

- `ecms_remote_decay.py` + `results/ecms_remote_decay.json`: Distance decomposition
- `ecms_cavity_decay_proof.py` + `results/ecms_cavity_decay.json`: Per-vertex γ factors
- `ecms_remote_families.py` + `results/ecms_remote_families.json`: Family convergence
- `ecms_worst_remote.py` + `results/ecms_worst_remote.json`: Exhaustive worst trees
- `ecms_mixed_caterpillar.py` + `results/ecms_mixed_caterpillar.json`: Mixed structures
- `ecms_caterpillar_algebra.py`: Algebraic fixed points

### Dead ends for proving |δμ| < 1

- Pure algebraic manipulation of the identity: all routes are circular (α < 1+β ⟺ α < 1+β)
- Global bounds μ(T) < (2/3)(n-1) minus μ(T/e) > 0: gives O(n), useless
- Naive cavity decay bound: local ≤ 1/2, total decay ≤ 1/2 × 1.618 ≈ 0.81. Gives total ≤ 1.3, too weak (but not by much!)
- Per-vertex γ < 1: FAILS (γ → 1 on long paths)
- Bound |remote| < 1/4: FAILS (exhaustive max = 0.264 > 0.250)

### Remaining path to proof

**Need: |remote| < 1/2.** Empirically |remote| < 0.264. Gap is large.

Best prospect: **path-level decay product.** While individual γ → 1, the product along any path from the contracted edge to a distant vertex decays at rate ~(1/φ²)^d ≈ 0.382^d per hop. The total remote contribution at distance d is bounded by initial_perturbation × 0.382^d × (branching at d). Need to show branching doesn't overwhelm decay.

Alternative: **Codex's WHNC inequality** (weighted heavy-neighborhood compensation). If proved, gives Conjecture A directly, bypassing ECMS entirely.
