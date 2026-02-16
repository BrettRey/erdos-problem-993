# Theoretical Attack Strategies for Erdős #993 (UPDATED)

## Problem Summary

The conjecture: The independent set sequence of every tree is unimodal.
- Verified through n=26 (447M trees)
- 2 log-concavity failures at n=26 (both subdivided stars)
- Brooms get closest to violation (nm → 1 as s → ∞)

## CRITICAL CORRECTION

**Brooms are already solved!** 
- Broom B(p,s) with s ≥ 2 is a spider
- Li-Li-Yang-Zhang (2025) proved all spiders are log-concave
- Therefore all brooms are unimodal

**The actual concerning configuration: SUBDIVIDED STARS (SST)**
- The two LC failures at n=26 are subdivided stars
- These are NOT spiders (they have multiple degree-3+ vertices)
- Galvin constructed families of SSTs with LC failures

## What's Been Tried

### 1. Broom Unimodality - NOW SOLVED
- **Status**: Covered by spider theorem!
- Li-Li-Yang-Zhang proved spiders are log-concave
- Brooms are spiders for s ≥ 2, paths for s ≤ 1

### 2. Subdivision Lemma - MOST PROMISING
- **Status**: Empirical hold through n=19 (5.7M subdivisions)
- **Goal**: Prove edge subdivision preserves unimodality
- If true: minimal counterexamples have no degree-2 vertices
- **Key data**: No counterexample found in exhaustive checks
- See `notes/subdivision_lemma.md` for full details

### 3. Two-Branch Vertex Class (C2)
- **Status**: Verified through n=24, 0 failures
- 196,635 trees tested
- Worst LC ratio: 0.846

### 4. Tight Mode Analysis
- Characterized when mode(I(T-w)) = d(I(T))
- Key invariants identified

## High-Value Theoretical Attacks

### Attack A: Complete Subdivision Proof (HIGHEST PRIORITY)

The subdivision lemma states:
> Let T' be obtained by subdividing any edge of T. 
> If I(T) is unimodal, then I(T') is unimodal.

**Why it matters**:
- If true, any minimal counterexample has no degree-2 vertices
- This dramatically reduces the search space
- Would explain why subdivided stars are the obstruction

**What's needed**:
- Prove tail ratio-monotonicity for forest polynomials
- Or: prove a first-difference dominance inequality
- Or: finite kernel check for leaf-light cases

**Empirical status**:
- No counterexample: n≤19 exhaustive (5.7M checks)
- Ratio-tail conditions hold in all checked cases
- Multiple proof paths attempted (see notes)

### Attack B: Two-Branch Proof (C2 Class)

- Trees with ≤2 branch vertices
- Could use transfer-matrix / inductive proof
- 196,635 tested, 0 failures

### Attack C: Mode Descent Analysis
- Characterize when violations can occur
- Sharp peak condition: d(I) = mode + 1

## Updated Priority Order

| Priority | Attack | Status | Impact |
|----------|--------|--------|--------|
| 1 | **Subdivision lemma** | 80% empirical | HIGHEST |
| 2 | C2 class proof | Exploratory | HIGH |
| 3 | Mode descent | Conceptual | MEDIUM |

## Key Files

- `notes/subdivision_lemma.md` - The main attack
- `notes/broom_unimodality_proof.md` - Now superseded by spider theorem
- `notes/tight_mode_analysis.md` - Invariant characterization
- `notes/steedman_invariants.md` - Blocked approach
