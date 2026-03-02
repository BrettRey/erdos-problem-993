# GPT 5.2 Round 5 Triage (2026-02-28)

## Instance assignments
- **Instance 1**: Incremental factor induction (fresh redirect from dead SV1)
- **Instance 2**: Cauchy-Binet curvature domination (has binomial-minor context)
- **Instance 3**: TN₃/block-Toeplitz construction (continuation from Round 4)

## Convergent findings (all three instances agree)

### 1. Leaf-factor case is PROVED

Both Instances 1 and 2 independently proved the leaf case. Instance 1's proof is cleaner:

**Proof (Instance 1).** Attaching a leaf means I_new = (1+x)I_old, E_new = E_old. Then:
- a'_k = e_k (coefficients of (1+x)I_old)
- b'_k = b_k (unchanged)
- d'_k = e_{k+1}b_k - e_kb_{k+1} = Δ_k (old)

So the new unsmoothed d' IS the old SCC quantity Δ_k, which is ≥ 0 by hypothesis. Then
applying the binomial upgrade lemma (LR-dominance + LC of E ⟹ (1+x) LR-dominance)
gives Strong C for the new pair. QED.

### 2. Off-diagonal cross-minors are the sole obstruction

All three instances converge: the product-level SCC decomposes into:
- **Diagonal terms**: weighted sums of factor-level adjacent minors d^{(m)}_t
  (controllable via factor-level SCC)
- **Off-diagonal terms**: non-adjacent minors K_{u,v} = p_{u+1}q_v - p_u q_{v+1}
  (the obstruction, can be negative)

Instance 2's formula (3.3): Δ⁺_k = Σ_{u,v} e_{k-u} b_{k-v} · K_{u,v}

Instance 2's bivariate identity (3.5): F_{eP,bQ} = F_{e,b} · PQ + eb · F_{P,Q}

### 3. Tree-realizability is essential (no purely axiomatic closure)

Confirmed by all instances and by our counterexample verification. Generic pairs
satisfying J≤E + LC + Cond C can still fail product closure (2 failures / 50K).

## Instance-specific contributions

### Instance 1: LC brick wall (CRITICAL INSIGHT)

**Claim:** Any tree F can appear as dp[v][0] = I(F; x) by attaching a new vertex v to F.
So once LC fails for some tree (n=26), you CANNOT assume LC of E^{(k)} = ∏I_{c_i}
uniformly. The LC hypothesis has to be either local-in-k or eliminated entirely.

**Impact:** Kills the naive incremental approach (Prompt 1) as stated. The proof of closure
MUST use the full Strong Condition C inequality, not just "curvature bonuses from LC."

**Status:** This is logically correct but less devastating than it sounds, because:
1. At n≤22 (our verification range), LC holds at all intermediate stages (89.9M checks)
2. LC failures at n=26 are exactly 2 trees out of 279M — extremely rare
3. E^{(k)} at a support vertex is ∏I_{c_i} where each I_{c_i} is a subtree IS poly.
   The n=26 LC failures might not be realizable as such products at the relevant stages.

To determine: Can the n=26 non-LC trees actually appear as factors in an intermediate
product at a support vertex? If not, LC of E^{(k)} might still hold universally at
support vertices even though global LC fails.

### Instance 2: Kernel form + bivariate antisymmetrization

**Clean deliverables:**
1. W_{I,E}(x,y) = E(x)I(y) - I(x)E(y) bivariate formulation
2. Two-factor decomposition: W_{I,E} = W_{I₁,E₁}·(E₂·I₂) + (I₁·E₁)·W_{I₂,E₂}
3. Cauchy-Binet for product curvature c_k^{prod} (fully rigorous, c_k^{prod} ≥ 0)
4. Off-diagonal reduction to "tree-only bound" on long minors

**The missing bound (Instance 2's "tree-only bound"):**
For tree-derived factors, the long minors L_{p,q} = b_p a_q - a_p b_q should satisfy:
|L_{p,q}| ≤ Σ_{t=p}^{q-1} (a_t/b_t) · c_{t+1}

### Instance 3: Semigroup encoding + TN₃ obstruction diagnosis

**Clean deliverables:**
1. Matrix encoding: U(I,E,J) = [[I, xJ], [0, E]] (upper triangular 2×2 polynomial matrix)
2. Exact multiplicativity: U(triple₁) · U(triple₂) = U(product triple)
3. Toeplitz lift: T(U₁U₂) = T(U₁) · T(U₂)

**TN₃ obstruction (honest diagnosis):**
- Block-triangular structure forces all 3×3 minors to factor as b_{k-1} × (2×2 minor)
- No 3×3 minor of T(U) or T(Ũ) equals b_{k-1}·Δ_k as a polynomial identity
- Need: 3-dimensional representation ρ of the 2×2 upper-triangular semigroup
  such that Toeplitz minors of T(ρ(U)) encode Strong C

**Suggested direction:** Sym² representation (symmetric square of 2×2 matrices → 3×3)
- Multiplicative by construction
- Right dimension for 3×3 minors
- Unresolved: basis/normalization to make Toeplitz minor = Δ_k

## MY COMPUTATIONAL FINDING (during triage)

**Long minor curvature bound verified (n≤14, 4.66M long minors):**

For ALL tree-derived factors (P=I_c, Q=E_c):
- 23.1% of long minors L_{p,q} are negative (1.076M / 4.66M)
- max |L_{p,q}| / sum(c_t for t in [p+1,q]) = **0.500** (exactly 1/2!)
- ZERO cases where curv_sum = 0 but L < 0

**This means:** |L_{p,q}| ≤ (1/2) · Σ_{t=p+1}^{q} c_t for ALL tree-derived factors.

This is exactly Instance 2's "tree-only bound" with an explicit constant of 1/2!
The curvature budget is ALWAYS at least 2× the negative long minor.

**Significance:** If this 1/2 bound can be proved (or even just verified at larger n),
it would close the product closure proof via Instance 2's framework:
1. Diagonal: regroup via factor-level SCC (✓)
2. Off-diagonal: bound by curvature with factor 1/2 (✓ computationally)
3. Product curvature via Cauchy-Binet absorbs the debt (✓ proven)

## UPDATED: Factor-level curvature budget FAILS (verified n≤19)

The "1/2 bound" |L_{p,q}| ≤ (1/2)·Σ c_t was verified at n≤14, but **breaks** at n=15:

| n | max ratio |L|/Σc | inf cases |
|---|-----------|------------|
| ≤14 | 0.500 | 0 |
| 15 | 0.656 | 0 |
| 17 | 0.848 | 0 |
| 18 | 0.859 | 0 |
| 19 | **1.067** | 0 |

At n=19: L_{8,10} = -192, Σ c_t = 180, ratio = 1.067. The factor-level curvature
is INSUFFICIENT to control the long minors.

**All violations share the same structure:**
- Tree deg seq: [4, 3, 3, 3, 2, 2, 2, 2, ...]
- Span 2 (q-p=2), at the tail (near independence number)
- Subtree is nearly the full tree (n-1 vertices)
- n=19 violation: P=[1,18,136,565,...,34,1], Q=[1,17,121,470,...,18,1]

**Key conclusion:** Instance 2's "tree-only bound" (|L_{p,q}| ≤ Σ c_t) is FALSE
in its literal form. The control must come from the PRODUCT-LEVEL weighting
(e_{k-u}·b_{k-v} in the kernel sum), not from factor-level curvature alone.

But product-level SCC holds at ALL 89.9M intermediate stages (n≤22, 0 failures).
So the proof must work at the product level, using the accumulated structure.

## Action items (revised)

1. **Round 6 prompts**: 3 new prompts incorporating the curvature budget failure
2. **Product-level kernel analysis**: profile Σ e_{k-u}b_{k-v}·K_{u,v} directly
3. **Instance 3 Sym² follow-up**: explicit computation of 3×3 Toeplitz minors
4. **Check whether product-level K_{u,v} weighting always saves**: the e_{k-u}b_{k-v}
   weights must concentrate on indices where K_{u,v} ≥ 0
