# Subdivision Lemma: Corrected Findings (2026-02-15)

## Bug Fix

The previous analysis used the wrong formula: A_old = x²R_uR_v + xP_uP_v (an extra factor of x). The correct formula is:

**A(x) = P_uP_v + xR_uR_v**

where P_u = IS poly of u-side with u included, R_u = IS poly with u excluded.

Derivation: I(T') = I(T) + A where T' subdivides edge uv with new vertex w:
- w not in S: IS of T minus edge uv = I(T) + P_uP_v (adding {u,v both in} case)
- w in S: forces u,v out, contributing x·R_u·R_v

Verified: A_theory = I(T') - I(T) matches direct computation for 10,018 edges (n≤12), 0 mismatches.

## Key Results (all with CORRECT formula)

### A(x) is always unimodal and log-concave

Verified for 9,071,864 edge subdivisions (all trees through n=19). Zero failures.

### I(T') is always unimodal and log-concave (through n=19)

Zero failures in 9M+ checks. (LC fails at n=26 but unimodality never fails.)

### A ≤ (1+x)I coefficientwise

Verified through n=16 (464,871 edges), 0 violations.

### mode(A) - d(I) distribution

| diff | count | percent |
|------|-------|---------|
| +0 | 2,412,894 | 26.6% |
| +1 | 6,650,565 | 73.3% |
| +2 | 8,405 | 0.09% |

A's mode is always at or past I's mode. The gap (d_A = d+2) is very rare.

### Combined tail: I(T') nonincreasing from d+1

Verified: 0 failures in 9M+ checks. I(T') is nonincreasing at every position k ≥ d+1.

### d(I') never moves earlier

d(I(T')) ≥ d(I(T)) always. d(I') = d(I) in 77.6%, d(I') = d(I)+1 in 22.4%.

### B = (1+x)I - A is always unimodal

Zero failures in 9M+ checks. B is not always LC (64 rare failures).

## Conditional Proof of Subdivision Lemma

**Theorem.** If the following conditions hold for tree T with edge uv:
- (C1) A = P_uP_v + xR_uR_v is unimodal
- (C2) first_descent(A) ≥ first_descent(I(T)) = d

Then I(T') is unimodal.

**Proof.** Let d_A = first_descent(A).

*Case d_A = d (26.6% of edges):*
For k < d: both I and A nondecreasing → I+A nondecreasing.
For k ≥ d: both I and A nonincreasing → I+A nonincreasing.
I+A is unimodal with peak at d. ✓

*Case d_A = d+1 (73.3% of edges):*
For k ≤ d-1: both ascending → I+A ascending.
For k ≥ d+1: both descending → I+A descending.
At k = d: one uncertain position. Sign pattern: +...+, ?, -...-.
Regardless of ?, the sequence has at most one transition from nondecreasing to nonincreasing. Unimodal. ✓

*Case d_A = d+2 (0.09% of edges):*
For k ≤ d-1: both ascending → I+A ascending.
For k ≥ d+2: both descending → I+A descending.
Two uncertain positions {d, d+1}. Sign pattern: +...+, ?₁, ?₂, -...-.
Unimodal iff no valley (no sign pattern (-, +) at (?₁, ?₂)).

Empirical verification: all 8,405 gap cases have sign pattern (+, -). The combined tail condition (I+A nonincreasing from d+1) always holds, giving ?₂ ≤ 0 and thus ruling out valleys.

**Verified:** All conditions (C1, C2) plus combined tail hold for 9,071,864 edge subdivisions (all trees n ≤ 19). Zero failures.

### Gap case analysis (d_A = d+2)

Sign patterns at (?₁, ?₂):
- (+, -): 8,405/8,405 = 100%

This remarkable rigidity means: in EVERY gap case, I+A peaks at exactly d+1 (one position past I's peak). The peak never stays at d or shifts further.

Combined tail margins at d+1:
- min margin: 8
- p10: 590
- median: 1,427
- max: 5,863

A_rise / I_drop ratio at d+1:
- max: 0.164 (A's rise is at most 16.4% of I's drop)
- median: 0.035

### What's needed for a full proof

The conditional proof (C1 + C2 → unimodal for 99.91% of cases) is a genuine theorem. To complete the proof:

1. **Prove C1**: A = P_uP_v + xR_uR_v is unimodal. Both P_uP_v and xR_uR_v are products of LC polynomials (hence LC), but their sum is not automatically unimodal.

2. **Prove C2**: first_descent(A) ≥ d. This says A is nondecreasing throughout I's ascending phase.

3. **Close the 0.09% gap**: Either prove d_A ≤ d+1 (false: 8,405 counterexamples) or prove the combined tail at d+1 when d_A = d+2. The algebraic bound A ≤ (1+x)I combined with LC of I is NOT sufficient (A[d+1] < 2I[d+2] in all gap cases). A tighter structural argument is needed.

### Failed algebraic approach for the gap

To prove (I+A)[d+2] ≤ (I+A)[d+1] when d_A = d+2:
- Need: A[d+2] - A[d+1] ≤ I[d+1] - I[d+2]
- Using A ≤ (1+x)I: get A[d+2] ≤ I[d+2] + I[d+1]
- This gives: need A[d+1] ≥ 2I[d+2]
- But A[d+1] < 2I[d+2] in ALL 8,405 gap cases (tightest: A[d+1]/(2I[d+2]) = 0.49)
- So the A ≤ (1+x)I bound alone cannot prove the combined tail.

## Computational Certificates (CORRECTED)

| Check | Edges | n range | Failures |
|-------|-------|---------|----------|
| A(x) unimodal | 9,071,864 | 3–19 | 0 |
| A(x) LC | 9,071,864 | 3–19 | 0 |
| I(T') unimodal | 9,071,864 | 3–19 | 0 |
| I(T') LC | 9,071,864 | 3–19 | 0 |
| A ≤ (1+x)I | 464,871 | 3–16 | 0 |
| d(A) ≥ d(I) | 9,071,864 | 3–19 | 0 |
| d(I') ≥ d(I) | 9,071,864 | 3–19 | 0 |
| Combined tail | 9,071,864 | 3–19 | 0 |
| B=(1+x)I-A unimodal | 9,071,864 | 3–19 | 0 |
| B=(1+x)I-A LC | 9,071,864 | 3–19 | 64 |
| Valley risk (d_A=d+2) | 8,405 | 7–19 | 0 |

### Extended certificates (n ≤ 20)

| Check | Count | n range | Failures |
|-------|-------|---------|----------|
| Delta bound (A ascending through d) | 24,710,099 edges | 3–20 | 0 |
| xRR ascending before d(I) | 24,710,099 edges | 3–20 | 0 |
| P_u[i] ≤ R_u[i-1] | 24,710,099 edges | 3–20 | holds |
| PP[k] ≤ RR[k-2] | 24,710,099 edges | 3–20 | holds |
| ECMS: \|mode(I(T))-mode(I(T/e))\| ≤ 1 | 24,710,099 edges | 3–20 | 0 |
| δ(T,v) ≤ 1 (vertex removal) | 26,056,121 pairs | 3–20 | 0 |
| \|δ(T,v)\| ≤ 1 | 26,056,121 pairs | 3–20 | 0 |
| I(T')=I(T)+xI(T/e) identity | 66,697 edges | 3–14 | 0 |
| Tight ratio A[k+1]/A[k] ≥ 1.125 | 3,348,674 edges | 3–18 | 0 |
| Tightest position at k=d-1 | 3,348,674 edges | 3–18 | 100% |

## Subdivision-Contraction Identity (KEY DISCOVERY 2026-02-15)

**Theorem.** For any tree T and edge e = (u,v):

    I(T_e) = I(T) + x · I(T/e)

where T_e is the subdivision of e (insert new vertex) and T/e is the contraction
of e (merge u,v into one vertex).

**Proof:** Decompose IS of T_e based on whether the new vertex w is in S:
- w ∉ S: edges u-w, w-v don't block u or v. S restricted to V(T) is any IS of T
  plus possibly {u,v both in} (allowed since u,v not adjacent in T_e). Contribution: I_u·I_v.
- w ∈ S: forces u,v out. Contribution: x·R_u·R_v.
So I(T_e) = I_u·I_v + x·R_u·R_v.

Meanwhile, A = I(T_e) - I(T) = P_uP_v + xR_uR_v (derived earlier).
And P_uP_v + xR_uR_v = x(x·∏R_i + ∏I_i) = x·I(T/e) where T/e is the contraction.

**Verified computationally:** 66,697 edges (n ≤ 14), 0 mismatches.

**Consequences:**
1. The delta bound (A ascending through d) is equivalent to mode(I(T/e)) ≥ mode(I(T)) - 1.
   Edge contraction drops IS mode by at most 1.
2. mode(A) - d = mode(I(T/e)) - mode(I(T)) + 1. The distribution {0, 1, 2} maps to
   contraction mode change of {-1, 0, +1}.
3. Unimodality of I(T_e) reduces to: sum of I(T) and x·I(T/e) is unimodal.
   Since both are tree IS polynomials (unimodal by induction hypothesis), and their
   modes are within 2 of each other, the sum is automatically unimodal when modes
   differ by ≤ 1 (99.91% of cases).

**Edge Contraction Mode Stability (ECMS):** |mode(I(T)) - mode(I(T/e))| ≤ 1.
Verified: 24,710,099 edges (n ≤ 20). Zero violations.

Distribution of mode(I(T/e)) - mode(I(T)):
- -1: 26.6% (contraction drops mode)
- 0: 73.3% (contraction preserves mode)
- +1: 0.09% (contraction raises mode)

## New Properties Discovered (2026-02-15)

### xRR Ascending Property (STRONG)

**Statement:** For every edge (u,v) of tree T, x·R_u·R_v never descends before d(I(T)).
Equivalently: mode(R_u·R_v) ≥ d(I(T)) - 1.

Verified for 24,710,099 edges (n ≤ 20). Zero failures.

This is arguably the cleanest formulation of C2. It says the "root-excluded product" xRR
is always ascending throughout the ascending phase of I(T).

mode(RR) vs d(I) distribution (n ≤ 20):
- mode(RR) = d+1: 39,555 (0.16%)
- mode(RR) = d: 16,229,025 (65.68%)
- mode(RR) = d-1: 8,441,519 (34.16%)

So mode(RR) ∈ {d-1, d, d+1}. The property mode(RR) ≥ d-1 is tight
(34% of edges hit the lower bound).

### Vertex Removal Mode Stability: |δ(T,v)| ≤ 1

**Statement:** For any tree T and vertex v, |mode(I(T)) - mode(IS(T-v))| ≤ 1.

δ(T,v) distribution (n ≤ 20, 26,056,121 (tree,vertex) pairs):
- δ = -1: 194,170 (0.75%)
- δ = 0: 20,260,669 (77.76%)
- δ = +1: 5,601,282 (21.50%)

Max |δ| = 1 through n=20. Zero violations of δ ≥ 2.

Note: δ < 0 means removing a vertex INCREASES the mode. Example: K_{1,3} center
removal gives δ = -1 (mode goes from 1 to 2).

Closed neighbourhood removal δ_N(T,v) = mode(I(T)) - mode(IS(T-N[v])) can be
much larger (up to +10 at n=20), so the single-vertex bound is sharp.

### Decomposition I(T') = I_u·I_v + x·R_u·R_v

This is an equivalent (and cleaner) form of the subdivision formula. Both summands
are log-concave (products of LC polynomials; prepending 0 preserves LC).

mode(I_u·I_v) vs mode(x·R_u·R_v) distribution (n ≤ 20):
- gap = -2: 49,925 (0.20%)
- gap = -1: 13,336,766 (53.97%)
- gap = 0: 11,320,164 (45.81%)
- gap = +1: 3,244 (0.01%)

For 99.8% of edges, modes differ by ≤ 1. For nonneg unimodal sequences with
modes within 1, their sum is automatically unimodal. The 0.2% with gap = ±2
need an additional argument (all verified to be unimodal).

### Failed Proof Approaches (do NOT revisit)

1. **Product mode domination**: mode(∏R_c) ≤ mode(∏I_c) is FALSE (40,355 failures,
   n ≤ 18). Despite R_c ≤ I_c coefficientwise, the mode of the product can exceed.

2. **δ(T,v) ≥ 0 assumption**: FALSE. δ can be -1 (0.53% of cases). Removing a vertex
   CAN increase the IS polynomial mode.

3. **Inductive proof of δ ≤ 1**: Tried via factorwise product domination of subtree
   polynomials. Breaks because mode(∏R_c) > mode(∏I_c) is possible.

4. **Fibonacci sufficient condition** (RR[k] ≥ RR[k-1] + PP[k]): FALSE (2,368,507
   failures at n ≤ 20). Too crude a bound on the PP drop vs RR rise.

5. **"Modes within 1" approach**: |mode(II) - mode(xRR)| ≤ 1 would give automatic
   unimodality of I(T'), but this FAILS for 0.2% of edges (gap = 2).

## Implications

### For the paper
The subdivision lemma has a clean conditional proof covering 99.91% of cases, with three verifiable conditions. The remaining 0.09% is handled by an empirically verified combined tail condition. This is a strong computational result suitable for Experimental Mathematics.

### For a full proof
The most promising routes:
1. **Prove A unimodal algebraically** using the product structure A = P_uP_v + xR_uR_v
2. **Prove d_A ≥ d** — perhaps using Newton's inequality on the component polynomials
3. **Prove combined tail** — needs a bound tighter than A ≤ (1+x)I; perhaps the factored structure of A provides additional constraints

### For minimal counterexample reduction
If proved, any minimal counterexample to unimodality has no degree-2 vertices (is homeomorphically irreducible). Combined with finite kernel verification, this would reduce the conjecture to a finite problem.
