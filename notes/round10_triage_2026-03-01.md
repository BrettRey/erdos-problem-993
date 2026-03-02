# Round 10 Triage (2026-03-01)

**Models:** GPT 5.2 Pro (3 instances), Codex 5.3 (2 instances)

## Prompt 1: W-Form Scan + s>=2 Profiling

**GPT P1 results (exhaustive n<=18, partial n=21-22):**
- W-form: 0 failures, 14.3M checks (n<=18 exhaustive)
- Partial n=21-22: 0 failures (~400K trees scanned)

**s>=2 extremal profile (n<=18):**

| s | instances | neg Term2 rate | min ratio | min margin |
|---|-----------|---------------|-----------|------------|
| 2 | 325,087 | 55.6% | 1.264 | 1 |
| 3 | 73,510 | 56.1% | 1.320 | 1 |
| 4 | 10,786 | 54.2% | 1.406 | 1 |
| 5 | 1,293 | 51.1% | 1.552 | 1 |
| 6 | 141 | 46.8% | 1.552 | 1 |
| 7 | 14 | 35.7% | 3.000 | 1 |
| 8 | 1 | 0% | - | 1 |

Key: min ratio INCREASES with s (opposite of s=1 behavior). Absolute margin = 1 universally.

**Factor-level invariants at T_{3,4} (n=28):**
- Factor E_t LC: YES
- Factor SCC: YES
- Factor leaf-aug: YES
- Factor E_t >= J_t: **NO** (d_12 = -2036)

So the factor-level invariant is {LC, SCC, J<=E}, NOT {E>=J}. The overall E>=J at the support vertex EMERGES from the product structure despite factor-level E>=J failure.

**E>=J at all rootings of T_{3,4}:**
- v (top root, deg 5): HOLDS (min minor = 0 at k=13)
- x-vertices (support, deg 2): HOLDS (min minor = 0 at k=14)
- w-vertices (non-support, deg 5): **FAILS** (d_13 = -4064)
- y-vertices (leaves, deg 1): **FAILS** (d_13 = -1972)

E>=J is strictly a support-vertex property.

## Prompt 2: W-Form Algebraic Structure

**Key identity (GPT P2):**
```
W_k = C_k * Δ_k(E_new, C) - B_k * c_k(C)
```
So W_k >= 0 is STRONGER than E_new >= C (ratio dominance). Curvature B_k*c_k(J) is a built-in rescue.

**s=2 linearity decomposition:**
W_k[f] is linear in f, so for s=2 with f = (1+x)^ℓ*(E_1 + xJ_1):
```
W_k = W_k[(1+x)^ℓ E_1] + W_k[x(1+x)^ℓ J_1]
```
Reduces to two sub-problems with known (f,g) structure.

**Cauchy-Binet double-sum (W-double):**
```
W_k = Σ_{i,j} f_i g_j [J_k(q_{k+1-i}q_{k-j} - q_{k-i}q_{k+1-j})
                       + J_{k+1}(r_{k-i}q_{k-1-j} - r_{k-1-i}q_{k-j})]
```

**Precise obstruction:** The mixed bracket L_{i,j} = r_{k-i}q_{k-1-j} - r_{k-1-i}q_{k-j} is NOT sign-definite. No amount of factor-level constraints controls it pointwise.

**Aggregated condition:** Summing over j with weights g_j gives:
```
Σ_j g_j L_{i,j} = r_{k-i}J_{k-1} - r_{k-1-i}J_k
```
So the aggregated obstruction reduces to: does r_{k-i}/r_{k-1-i} <= J_k/J_{k-1} hold? This is where LC of J (monotone J_k/J_{k+1} ratios) enters.

**Approach ranking:** D (transport) > C (s=2 recursion) > B (CB lens) > A (factor-pair).

## Prompt 3: Correction Term Bounds

**GPT P3 explicit s=2 extremals:**
- Smallest margin: 7-vertex tree (leaf + P2 arm + K_{1,2} arm), W_2 = 33 (Term1=36, Term2=-3)
- Smallest ratio: balanced double-star at n=20 (L=8 leaves each arm), ratio = 1.6635 at k=2
- Closed form: ratio = ((2L-1)(3L²+5L+4)) / (4(L-1)(L²+2L-4)) → 3/2 as L→∞

**Dead approaches (confirmed by GPT P3 + Codex P3):**
- B = A - f*(q-r): d_{k-1}(f*(q-r), C) < 0 on tree data
- Consecutive Karlin: J_k*K_k >= J_{k+1}*K_{k-1} FALSE on tree data

**Key structural insight (GPT P3):**
q - r has zero constant term (divisible by x) because q-r comes from expanding at least one xJ_ui choice. This is a real tree-specific constraint.

## Codex Results

**Codex P1:** Created `round10_wform_verification.py` (25KB), running n=21-22 with 8 workers. Results pending.

**Codex P3:** Confirmed GPT's findings independently:
- B=A-f*(q-r) fails
- Consecutive Karlin fails
- s=2 star-arm extremal: ratio 1.6635 at n=20

## Strategic Assessment

The proof is now cleanly structured:

1. **Target:** E >= J at support vertices (PROVED for s=1, OPEN for s>=2)
2. **Framework:** W-form induction (Karlin + correction)
3. **Karlin term:** PROVED non-negative
4. **Correction term:** bounded by Karlin empirically, no algebraic proof
5. **s=1 handles 62.8%** of support vertices. s=2 handles ~30%. s>=3 is rare.
6. **The obstruction** is precisely the mixed bracket L_{i,j} in the CB expansion
7. **The resolution** likely involves monotone transport using LC of J

The s≥2 ratio being bounded away from 1 (→ 3/2 for s=2, higher for s≥3) is actually GOOD news: unlike s=1 (ratio → 1), the s≥2 case has a gap that might be provable.

## Round 10 Addendum: Detailed Algebra (GPT P2/P3 full outputs)

### Full CB expansion (GPT P2, equation W-double)

```
W_k = Σ_{i,j} f_i g_j [J_k · K_{i,j} + J_{k+1} · L_{i,j}]
```

where:
- K_{i,j} = q_{k+1-i}·q_{k-j} - q_{k-i}·q_{k+1-j} (Toeplitz 2×2 minor of q = E_t, antisymmetric)
- L_{i,j} = r_{k-i}·q_{k-1-j} - r_{k-1-i}·q_{k-j} (mixed bracket, NOT antisymmetric)

Sign control for the Karlin term: for i<j, (f_i g_j - f_j g_i) ≤ 0 (from f≽g) and K_{i,j} ≤ 0 (from q PF2), giving each CB summand ≥ 0.

### Curvature identity (GPT P2 + P3, independently derived)

```
C_k · d_k(E_new, C) = W̃_k + B_k · c_k(C)
```

where W̃_k = C_k·d_k(A,C) + C_{k+1}·d_{k-1}(B,C) and c_k(C) = C_k² - C_{k-1}·C_{k+1} ≥ 0.
So W̃_k ≥ 0 is STRONGER than the target d_k(E_new,C) ≥ 0. The curvature B_k·c_k(C) is a built-in rescue.

### STP2 condition (GPT P2)

```
r_{m+1} · q_n ≤ r_m · q_{n+1}   for all m > n ≥ 0
```

Would resolve the mixed bracket obstruction entirely. Probably FALSE for trees, needs testing.

### Explicit extremals (GPT P3)

Smallest s=2 margin: 7-vertex tree (1 leaf + P2 arm + K_{1,2} arm), at k=2:
Term1 = 36, Term2 = -3, W = 33. Full polynomials:
- A = (1,5,9,7,2), B = (1,3,2), C = (1,3,3,1)

Balanced double-star closed form at k=2:
ratio = ((2L-1)(3L²+5L+4)) / (4(L-1)(L²+2L-4)) → 3/2 as L→∞

### Tree-specific constraint (GPT P3)

q - r has zero constant term (divisible by x): because q-r comes from expanding at least one xJ_ui choice in the product formula for E_t vs J_t.

### T_{3,3} broom (n=22) factor data (GPT P1)

Factor at x-support vertex: E_t ≽ J_t FAILS at k=9 with d_9 = -252.
Factor polynomials:
- E_t = (1, 19, 155, 716, 2073, 3916, 4843, 3794, 1717, 348, 4)
- J_t = (1, 16, 108, 404, 924, 1338, 1220, 668, 196, 24, 1)

## Diagonal/Cross Decomposition (Post-Round-10 Discovery)

**NEW** decomposition of the P2 minor w_k:

```
w_k = Σ_{i,j} E_acc_i · J_acc_j · φ(i,j)
    = D_k (diagonal, i=j) + X_k (cross, i≠j)
```

**D_k** = Σ_i E_acc_i · J_acc_i · d_{k-i}(f, g)   (weighted sum of factor LR minors)
**X_k** = w_k - D_k                                  (off-diagonal interaction)

### Verification (n ≤ 18, 410K support vertices, s ≥ 2):
- D_k ≥ 0 in ALL cases (0 failures)
- X_k ≥ 0 in ALL cases (0 failures)
- Both parts independently non-negative!

### Verified at T_{3,4} broom (n=28, s=1):
- D_k ≥ 0 and X_k ≥ 0 at all k
- Cross-term dominates: X_7 = 2598976972 >> D_7 = 58794144

### Interpretation (HPC "projection"):
Individual factor-level LR minors d_p(I_t, E_t) can be negative (~14% of support vertices).
But the aggregate product structure (weights E_acc_i · J_acc_j) ensures both the weighted
diagonal sum AND the off-diagonal cross-term are non-negative. This is the "emergent
property" that projects from the product mechanism.

## Next Steps (Round 11)

1. **Extend** diagonal/cross verification to n=20 (confirm both D≥0 and X≥0)
2. **Test** aggregated condition (Agg) and STP2
3. **Prove** D_k ≥ 0 or X_k ≥ 0 algebraically (either suffices for W≥0)
4. **s=2 binomial case**: f^(0) = (1+x)^ℓ · g gives A = (1+x)^ℓ · C — exploit binomial smoothing
5. SCC search background job: 0 fails through n=24 (39.3M trees), running n=25
