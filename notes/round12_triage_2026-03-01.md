# Round 12 Triage (2026-03-01)

**Models:** Codex 5.3 (1 instance), GPT 5.2 Pro (1 instance)
**Missing:** GPT P2 (STP2 algebraic proof) — output not received

## Codex P1: STP2 + Prefix E≽J + Karlin Split

### STP2 (CONFIRMED, 0 failures)
- Exhaustive n≤18: 921,265 factor instances, 4,868,953 inequality comparisons, **0 failures**
- Partial n=19-20: 714,243 more instances, 6,013,829 comparisons, **0 failures**
- Tightest slack: 5 at (q,r)=([1,3,1],[1,2,1]), (m,n)=(1,0)
- 50,926 distinct factor pairs tested

### Prefix E≽J at n=23-25 (sampled, not exhaustive)
- 2000 trees per n, 26,350 total support vertices
- **0 prefix failures, 0 full-range failures**
- Consistent with prefix E≽J being universal (first full-range failure at n=32)

### Karlin split prefix vs tail (exhaustive n≤18)
- 9,458,502 total (step,k) positions evaluated
- **w_k < 0**: 0 occurrences (prefix AND tail) — w_k ≥ 0 everywhere through n≤18
- Negative correction events: 67,409 (prefix), 167,814 (tail)
- **Worst prefix ratio**: main/|corr| = 1.4 (main=91, corr=-65, w=26)
- **Worst tail ratio**: main/|corr| = 4.4375
- **Prefix is HARDER** than tail for the Karlin-vs-correction comparison

### n=32 Karlin split at failing step
- STP2 holds on both factors (P_2 factor: vacuous; broom factor: 66 checks, 0 fails)
- At k=15 (the failure): **main_15 = 4694, corr_15 = -8192**
- Correction EXCEEDS Karlin: ratio = 4694/8192 = 0.573
- This is why w_15 = -3498 < 0: the Karlin term alone cannot rescue the tail

## GPT P3: s≥2 Completion Strategy

### THE KEY NEW FINDING: Shifted Ratio Condition ⋆

Using the curvature-augmented identity (same as Round 10 but now REDUCED):

```
C_k · w_k = C_k · d_k(A,C) + C_{k+1} · d_{k-1}(B,C) + B_k · c_k(C)
```

where A = E_acc·g, B = E_acc·h, C = J_acc·g.

Term1 ≥ 0 (Karlin), Term3 ≥ 0 (curvature bonus). So w_k ≥ 0 ⟺ Term2 + Term3 ≥ 0.

GPT's algebra: Term2 + Term3 ≥ 0 simplifies (mixed terms cancel!) to:

**⋆: B_k / B_{k-1} ≥ C_{k+1} / C_k**

This is a "shifted ratio dominance" condition. It says: B's growth rate at position k-1 beats C's growth rate at position k.

### Why ⋆ is better than "Karlin dominates correction"
- ⋆ IGNORES Term1 entirely — it only needs Term2 + Term3 ≥ 0
- The Karlin term is a FREE bonus on top
- ⋆ is a single one-line ratio inequality (vs the 2D CB expansion)
- ⋆ is local (depends on 4 consecutive coefficients of B and C)

### s=1 mechanism clarified
- The correction C_k is NOT automatically ≥ 0 from "h has nonneg coefficients"
- Tree-realizability is essential: (g,h) share the same child factor structure
- For s=1, any negativity occurs only beyond the prefix (consistent with SCC failing at n=28)

### Child ordering matters for ⋆
- B = E_acc · h where h = J_last (the last child processed)
- Ordering affects which child's J appears in B
- Choosing the "easiest" child last (smallest J influence) maximizes B_k/B_{k-1}
- Scalar proxy: order by decreasing μ(I_i) - μ(E_i) or J_i(1)/E_i(1)

### s=2 proof skeleton
1. Curvature identity (exact algebra) ✓
2. Term1 ≥ 0 (Karlin) ✓
3. Term3 ≥ 0 (LC of C) ✓
4. Need: ⋆ on prefix for at least one ordering
5. ⋆ with s=2 definitions: B = (1+x)^ℓ · I_1 · J_2, C = E_1 · E_2

**Shifted-growth lemma** (the NEW proof target):
For k < mode(I(T)), show that
```
((1+x)^ℓ · I_1 · J_2)_k / ((1+x)^ℓ · I_1 · J_2)_{k-1} ≥ (E_1·E_2)_{k+1} / (E_1·E_2)_k
```
for at least one ordering of children 1, 2.

### Product closure conjecture
- Use Karlin's CB expansion: d_k(P*g, Q*g) = Σ_{i<j} (P_j Q_i - P_i Q_j) · M_{i,j,k}(g)
- Each PF2 convolution amplifies by τ_k(g) depending on g's Toeplitz minors
- Ratio gap should INCREASE with each child processed

### Absolute margin = 1 strategy
- w_k is integer, so w_k > 0 ⟹ w_k ≥ 1
- Need: classify when both Karlin and curvature can simultaneously vanish
- Then strict positivity + integrality + finite base cases closes it

## Strategic Assessment

### The proof now has a CLEAN target: the shifted ratio condition ⋆

**For each s≥2 support vertex, for k < mode, for at least one child ordering:**
```
B_k · C_k ≥ B_{k-1} · C_{k+1}
```

This is a single determinantal inequality. It's testable, local, and has a probabilistic interpretation (comparison of consecutive hazard rates).

### ⋆ Computational Test (COMPLETED)

**Result: ⋆ alone is NOT sufficient, but the full 3-term identity IS.**

⋆ condition (B_k·C_k ≥ B_{k-1}·C_{k+1}):
- 2,440,182 incremental steps checked (all s≥2 SVs, n≤18)
- **84,025 step-level prefix ⋆ failures** (3.4%)
- **1,801 SVs** (0.44%) have no ordering satisfying ⋆ on full prefix

Full curvature-augmented identity (Term1+Term2+Term3 ≥ 0):
- **0 failures** across 2,440,182 steps, ALL orderings
- Karlin term (Term1) rescues every ⋆ failure
- Min ratio Term1/|Term2+Term3| = 3.33 (Karlin always 3.3x bigger than deficit)
- Min w_k = 0 (never negative)
- **Ordering freedom NOT needed**: w_k ≥ 0 at EVERY ordering

Per-s breakdown (⋆ alone):
- s=2: 325K SVs, 1389 fail (0.43%)
- s=3: 73K SVs, 373 fail (0.51%)
- s=4: 10K SVs, 36 fail (0.33%)
- s=5: 1293 SVs, 3 fail (0.23%)
- s≥6: 0 fails

### Key Insight: 3-Term = 2-Term (Algebraic Equivalence)

The curvature-augmented 3-term identity is just the 2-term identity multiplied by C_k:

2-term: w_k = d_k(E_acc·g, J_acc·g) + (B_k·C_k - B_{k-1}·C_{k+1})
3-term: C_k·w_k = Term1 + Term2 + Term3

Multiply 2-term by C_k and expand Term2+Term3:
C_k·(B_k·C_k - B_{k-1}·C_{k+1}) = B_k·C_k² - B_{k-1}·C_k·C_{k+1}
                                   = B_k·(C_k² - C_{k-1}·C_{k+1}) + C_{k+1}·(B_k·C_{k-1} - B_{k-1}·C_k)
                                   = Term3 + Term2  ✓

So the "curvature rescue" is an artifact of the split. The REAL proof target is simply:
**Karlin(k) + B_k·C_k - B_{k-1}·C_{k+1} ≥ 0** (the ⋆ quantity IS the correction)

### Updated Priority Ranking
1. **Prove w_k ≥ 0 for s≥2** (full 3-term identity): the proof target is Term1+Term2+Term3 ≥ 0, not ⋆ alone
2. **STP2 algebraic proof** (still waiting on GPT P2 output)
3. **Quantify: what fraction of w_k ≥ 0 comes from Karlin vs curvature?**
4. **Understand failure trees**: at ⋆-failing SVs, what structural property requires the Karlin rescue?
