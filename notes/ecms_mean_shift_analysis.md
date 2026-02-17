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

### Dead ends for proving |δμ| < 1

- Pure algebraic manipulation of the identity: all routes are circular (α < 1+β ⟺ α < 1+β)
- Global bounds μ(T) < (2/3)(n-1) minus μ(T/e) > 0: gives O(n), useless
- Naive cavity decay bound: local ≤ 1/2, total decay ≤ 1/2 × 1.618 ≈ 0.81. Gives total ≤ 1.3, too weak (but not by much!)
