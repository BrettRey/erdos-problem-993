# Round 14 Triage (2026-03-01)

**Models:** Codex 5.3 (1 instance), GPT 5.2 Pro (2 instances)

## Instance 1 (Codex): W-Form Ratio Profiling

### w_k ≥ 0 verification
- **Exhaustive n=9-17**: 0 failures. n=17: 48,629 trees, 27,811 SVs with s≥2, 707K checks.
- **Partial n=18**: 50K trees (~40%), 884K checks, 0 failures.
- **Sampled n=21-24**: 5K trees each, ~170-200K checks each, 0 failures.
- Total: 0 failures across all checks.

### Min W ratio by (n, s) — exhaustive n=9-17

| n  | s=2    | s=3    | s=4   | s=5   |
|----|--------|--------|-------|-------|
| 9  | 10.000 | —      | —     | —     |
| 10 | 5.625  | —      | —     | —     |
| 11 | 4.667  | 16.250 | —     | —     |
| 12 | 4.684  | 5.273  | —     | —     |
| 13 | 4.214  | 4.747  | —     | —     |
| 14 | 3.947  | 4.662  | 7.400 | —     |
| 15 | 3.714  | 4.518  | 6.385 | —     |
| 16 | 3.556  | 4.429  | 6.772 | 8.824 |
| 17 | 3.448  | 4.530  | 6.628 | 8.030 |

**Key**: s=2 always has the lowest ratio. Global min always at s=2.

### Star+star extremal extended to n=24

| n  | a | b  | ratio  |
|----|---|----|--------|
| 21 | 6 | 11 | 3.167  |
| 22 | 7 | 11 | 3.105  |
| 23 | 7 | 12 | 3.076  |
| 24 | 8 | 12 | 3.043  |

Limit: 1+√3 ≈ 2.732. **s≥2 ratio stays bounded away from 1.**

### Absolute margin: min w_k often = 0

| n  | s=2 | s=3 | s=4 | s=5 |
|----|-----|-----|-----|-----|
| 9  | 35  | —   | —   | —   |
| 10 | 0   | —   | —   | —   |
| 11 | 0   | 4   | —   | —   |
| 12 | 0   | 6   | —   | —   |
| 13 | 0   | 0   | —   | —   |
| 14 | 0   | 0   | 8   | —   |
| 15 | 0   | 0   | 0   | —   |
| 16 | 0   | 0   | 0   | 6   |
| 17 | 0   | 0   | 0   | 0   |

**min w_k = 0 frequently.** Absolute margin is NOT a useful robustness signal.
Ratio-based argument is the viable approach.

### Chain STP2: 0 failures
- Both STP2(E_t, J_t) and STP2(I_t, E_t) checked on all directed-edge factors.
- Exhaustive n=9-17, partial n=18 (50K trees): **0 failures**.
- Tightest slack: 5 (at n=9).

## Instance 2 (GPT): Cauchy-Binet Pairwise Domination — **DEAD**

### CB expansion derived correctly
- correction(k) = Σ_{i,j} E_i J_j · M_{i,j}(k) where M_{i,j}(k) = h_{k-i}g_{k-j} - h_{k-1-i}g_{k+1-j}
- Antisymmetric bracket L_{ij} = (A) - (B) where both (A) ≤ 0 and (B) ≤ 0 under STP2
- L_{ij} has **indeterminate sign** (difference of two ≤ 0 quantities)

### **Pairwise domination FAILS**
- |L_{ij}| ≤ K_{ij} is FALSE in general
- Combined K_{ij} + L_{ij} ≥ 0 is also NOT sign-definite
- **CB pairwise approach is conclusively dead**

### Star+star case works because h=[1] makes L extremely sparse
- Only nonzero at boundary indices (a=0,1)
- Explains why ratio ≈ 3 in extremal family

### Recommendation
- CB analysis is diagnostic, not sufficient for proof
- Suggests: derive 1-dimensional inequality for s=2 eliminating CB sums

## Instance 3 (GPT): STP2 Tree Induction — MAJOR CLARIFICATIONS

### Important correction: STP2 equivalence was WRONG
**STP2(f,g):** g(m+1)/g(m) ≤ f(n+1)/f(n) for m > n.
NOT equivalent to "g(k)/f(k) nonincreasing" (as previously stated in memory).

### Diagonal reduction lemmas (VALUABLE)
- **Lemma A**: Under LC(f), STP2(f,g) reduces to diagonal: g(k+1)f(k-1) ≤ g(k)f(k)
- **Lemma B**: Under LC(g), STP2(f,g) reduces to diagonal: g(n+2)f(n) ≤ g(n+1)f(n+1)
- These collapse the two-parameter inequality to one-parameter.

### Single-child step: CANNOT bootstrap from STP2(E,J) alone
- Deriving STP2(I_c, E_c) from STP2(E_c, J_c) + LC(E_c) is FALSE.
- Counterexample: E=(1,2,4,7), J=(1,1,1). STP2(E,J) holds but STP2(I,E) fails at (2,1).
- **Fix**: Use chain STP2 as the IH. Then single-child step is trivial.

### Multi-child step: THE ENTIRE REMAINING DIFFICULTY
- Approach A (incremental bracket): structurally aligned with STP2, but needs "off-diagonal domination"
- Approach B (Toeplitz/Schur): σ(F) = (F-1)/x shift operator, M(F) = [F, 0; σ(F), 1] matrix.
  M(FG) = M(F)·M(G) — makes shift multiplicative! Promising but incomplete.
- Approach C (ratio characterization): correct target form but not a closure engine

### STP2 does NOT follow from abstract axioms
Explicit counterexample: E=(1,19,77,72,50,22), J=(1,17,21,24,25).
Satisfies: LC, J≤E, deg(J)=deg(E)-1, prefix E≽J. But STP2 FAILS at (m,n)=(3,2).
**Tree-specific structure essential.**

### Chain STP2 is the right invariant
- STP2(E→J) and STP2(I→E) hold (0 failures)
- STP2(J→E) and STP2(E→I) FAIL
- Directional asymmetry matches DP arrows
- Chain STP2 = minimal closed-looking package

### σ/M(F) matrix trick — most promising new idea
- σ(FG) = σ(F)·G + σ(G) (exact identity, requires F(0)=G(0)=1)
- 2×2 matrix M(F) = [F, 0; σ(F), 1] satisfies M(FG) = M(F)·M(G)
- Could encode STP2 as TN2/TP2 condition on block matrix
- Product closure would follow from "multiply TN2 matrices → TN2"
- **Status**: promising direction, not yet completed

### Main obstacle identified
Tree-specific "off-diagonal domination" principle for convolution minors.
STP2 alone insufficient for product closure. Need to exploit shared-factor structure.

## Strategic Assessment

### What's now PROVED or CONFIRMED
1. w_k ≥ 0 at all s≥2 steps (exhaustive n≤17, sampled to n=24)
2. Star+star extremal (confirmed n≤24, ratio → 1+√3 > 1)
3. Chain STP2 (both directions, 0 failures n≤18)
4. Diagonal reduction lemmas for STP2 under LC
5. s=1 case proved

### What's DEAD
1. **CB pairwise domination** (Instance 2 conclusive)
2. **STP2 via abstract axioms** (counterexample, Instance 3)
3. **Absolute margin arguments** (min w_k = 0, Instance 1)
4. **Single-child bootstrap from STP2(E,J) alone** (counterexample, Instance 3)

### The SINGLE remaining gap
**Chain STP2 multi-child closure**: Given chain STP2 for each child, show
STP2(∏I_i, ∏E_i) holds.

Generic product closure FAILS (Round 13 counterexample). Tree-realizable closure
needs the shared-factor structure (I_i = E_i + xJ_i).

### Most promising next directions (ranked)

1. **σ/M(F) matrix encoding** (GPT Instance 3): Make shift multiplicative via 2×2 matrices.
   If STP2 can be encoded as a TN2 condition on block matrices built from M(E), M(I),
   closure follows from Karlin's TN2 product theorem.

2. **1-dimensional s=2 inequality** (GPT Instance 2 recommendation): Derive a closed-form
   inequality for the s=2 case that eliminates CB sums, working directly with ratio profiles.

3. **Ratio bound for s≥2**: Since ratio → 1+√3 > 1 (not → 1), a clean ratio-based proof
   may be possible. The extremal star+star is fully characterized; general trees are "further
   from the boundary."
