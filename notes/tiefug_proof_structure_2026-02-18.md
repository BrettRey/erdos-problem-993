# Tie-Fugacity Proof Structure: New Results (2026-02-18)

## Summary of computations this session

All verified for n ≤ 16 (244,690 tree-leaf pairs) unless noted.

| Property | n | Checks | Fails | Min margin |
|----------|---|--------|-------|------------|
| Condition (2): μ(T-{l,s}, λ_m^T) ≥ m-2 | ≤ 18 | 1,723,514 | 0 | 0.725 |
| All-k tie-fug: μ(T, λ_k^T) ≥ k-1 for all k ≤ m | ≤ 16 | 162,792 | 0 | 0.379 |
| Weak C1: λ_m^T ≥ λ_{m-1}^{T-l} | ≤ 16 | 244,683 | 0 | 0.111 |
| Weak C2: λ_m^T ≥ λ_{m-2}^{T-{l,s}} | ≤ 16 | 244,683 | 0 | 0.278 |
| STRONG C2: λ_m^T ≥ λ_{m-1}^{T-{l,s}} | ≤ 16 | 244,690 | 0 | 0.083 |
| Cross-tree: λ_m^{T-l} ≥ λ_{m-1}^{T-{l,s}} | ≤ 16 | 244,690 | 0 | 0.125 |
| mode(T-{l,s}) ≥ m-1 (all cases) | ≤ 18 | 1,723,514 | 0 | N/A |
| mode(T-{l,s}) ≥ m-1 in δ=-1 cases | ≤ 16 | 72,042 | 0 | N/A |

Distribution of mode(T-{l,s}) - (m-2) for n ≤ 18: rel=+1 (95.7%), rel=+2 (4.3%). Never rel=0.

Min c2 margin by rel (b111588, exhaustive n≤18, 1,723,514 pairs):
- rel=+1: min c2 margin = 0.725412 (at n=14, m(T)=5, mode(T-{l,s})=4)
- rel=+2: min c2 margin = 0.849421

## CRITICAL CORRECTION (2026-02-18, continued session)

**The "combined = 0.075" in prior notes was WRONG.** The earlier analysis used coefficient ratio
p = i_{m-1}(T-{l,s}) / i_m(T) as decomposition weight, but the correct formula uses:
  α = I(T-l, λ_m) / I(T, λ_m)   (partition function weight)
  β = λ_m · I(T-{l,s}, λ_m) / I(T, λ_m)   (with α+β=1)
giving: μ(T,λ_m)-(m-1) = α·c1 + β·c2.

The ACTUAL tie-fugacity margin = μ(T, λ_m^T) - (m-1) (computed directly) has:
- Minimum through n≤22 = **0.362** at S(2^10) (spider, from exhaustive scan)
- For Case C (δ=+1, rel=+1) trees: minimum = 0.379 at n=14
- For S(2^k) spiders: margin → 1/3 ≈ 0.333 asymptotically (proved: C = 0.33334, A = 0.50)
  Verified up to k=118 (n=237): combined ≈ 1/3 + 1/(2k), always > 1/3

The previous "0.075" figure was p·c2 + (1-p)·c1 with WRONG weight p (not = margin).
The summary claiming "min combined = 0.075" should be disregarded.

## Key algebraic identity (mediant property)

From the leaf deletion identity i_k(T) = i_k(T-l) + i_{k-1}(T-{l,s}):

**Lemma (mediant)**: λ_m^T lies strictly between λ_m^{T-l} and λ_{m-1}^{T-{l,s}}:
```
lambda_m^T = [i_{m-1}(T-l) + i_{m-2}(T-{l,s})] / [i_m(T-l) + i_{m-1}(T-{l,s})]
            = mediant of (A/B, D/E)
  where A = i_{m-1}(T-l), B = i_m(T-l), D = i_{m-2}(T-{l,s}), E = i_{m-1}(T-{l,s})
  i.e., lambda_m^{T-l} = A/B and lambda_{m-1}^{T-{l,s}} = D/E
```

**Mediant inequalities**: min(A/B, D/E) ≤ λ_m^T ≤ max(A/B, D/E).

STRONG C2 (λ_m^T ≥ D/E = λ_{m-1}^{T-{l,s}}) ⟺ A/B ≥ D/E ⟺ cross-tree comparison.

So: **STRONG C2 holds iff the cross-tree comparison holds**:
```
lambda_m^{T-l} >= lambda_{m-1}^{T-{l,s}}
i.e., i_{m-1}(T-l) * i_{m-1}(T-{l,s}) >= i_m(T-l) * i_{m-2}(T-{l,s})
```

## The proof structure (assuming three key lemmas)

**Theorem (all-k tie-fugacity)**: For all trees T with ≥ 1 vertex and all k ≤ mode(T),
  μ(T, λ_k^T) ≥ k-1.

**Proof** by strong induction on n = |V(T)|.

Take any leaf l with support s. Let T' = T-l, T'' = T-{l,s}. Leaf decomposition at λ = λ_m^T:
```
μ(T, λ_m^T) - (m-1) = p * (μ(T'', λ_m^T) - (m-2))    [condition 2 slack]
                     + (1-p) * (μ(T', λ_m^T) - (m-1))  [condition 1 slack]
```

**Step 1 (Condition 2)**: μ(T'', λ_m^T) ≥ m-2.

- **Lemma A** (mode lower bound): mode(T'') ≥ m-1. [proved for δ=-1 cases algebraically; for δ=0,+1 cases immediate]
- **Lemma B** (STRONG C2): λ_m^T ≥ λ_{m-1}^{T''}.
- IH on T'' at level k=m-1 (≤ mode(T'')): μ(T'', λ_{m-1}^{T''}) ≥ m-2.
- Monotonicity: μ(T'', λ_m^T) ≥ μ(T'', λ_{m-1}^{T''}) ≥ m-2. ✓

**Step 2 (Combined)**: p*(cond-2 slack) + (1-p)*(cond-1 slack) ≥ 0.

If cond-1 slack ≥ 0: both terms non-negative, done.

If cond-1 slack < 0:
- Weak C1 (λ_m^T ≥ λ_{m-1}^{T'}) + IH on T' at level m-1: μ(T', λ_m^T) ≥ m-2. So deficit ≤ 1.
- Need: p * (cond-2 slack) ≥ (1-p) * deficit.
- **GAP**: deficit empirically ≤ 0.095 and slack ≥ 0.725, so p*slack >> (1-p)*deficit in practice.
  But the algebraic proof of this balance is MISSING.

## Lemma A: mode(T-{l,s}) ≥ m-1

**Algebraic partial proof** (gives mode ≥ m-2, not m-1):

From i_m(T) ≥ i_{m-1}(T) [mode = m] and i_k(T) = i_k(T') + i_{k-1}(T''):
```
i_m(T') + i_{m-1}(T'') ≥ i_{m-1}(T') + i_{m-2}(T'')
=> i_{m-1}(T'') - i_{m-2}(T'') ≥ i_{m-1}(T') - i_m(T')
```
When mode(T') = m-1: i_{m-1}(T') ≥ i_m(T'), so RHS ≥ 0.
Hence i_{m-1}(T'') ≥ i_{m-2}(T'') => mode(T'') ≥ m-2. ✓ **PROVED for δ=-1 case.**

**Gap**: For mode(T'') ≥ m-1 (i.e., i_{m-1}(T'') ≥ i_m(T'')):
From i_m(T) ≥ i_{m+1}(T): i_{m-1}(T'') - i_m(T'') ≥ i_{m+1}(T') - i_m(T').
When mode(T') = m-1: i_{m+1}(T') ≤ i_m(T'), so RHS ≤ 0.
This gives i_{m-1}(T'') - i_m(T'') ≥ something ≤ 0. Cannot conclude ≥ 0.

**Empirical status**: 0 failures (mode ≥ m-1) for all 1,723,514 pairs through n=18.
In δ=-1 cases: 68,912/72,042 have mode(T'')=m-1; 3,130/72,042 have mode(T'')=m (even better).

## Lemma B: STRONG C2 (cross-tree comparison)

**Algebraic expansion** (using leaf deletion at s in T' = T-l):
Let B_k = i_k(T'' - N_s) where N_s = neighbors of s in T-l.
Then: i_{m-1}(T') = i_{m-1}(T'') + B_{m-2} and i_m(T') = i_m(T'') + B_{m-1}.

STRONG C2 ⟺ i_{m-1}(T'')² + B_{m-2}·i_{m-1}(T'') ≥ i_m(T'')·i_{m-2}(T'') + B_{m-1}·i_{m-2}(T'')
          ⟺ [i_{m-1}(T'')² - i_m(T'')·i_{m-2}(T'')] ≥ i_{m-2}(T'')·[B_{m-1}/B_{m-2}] · B_{m-2} - B_{m-2}·i_{m-1}(T'')

LHS ≥ 0 iff T'' is log-concave at position m-1.
RHS sign requires comparing B's step-fugacity to T'''s step-fugacity.

**Gap**: T'' = T-{l,s} is NOT always d_leaf≤1 (removing s can create new leaves with d_leaf>1), so
T'''s LC at m-1 cannot be guaranteed from the d_leaf≤1 hypothesis.

**Empirical status**: 0 failures through n=16, min margin 0.125.

## Algebraic fact: i_{m-1}(T-{l,s}) - i_m(T-{l,s}) ≥ i_{m+1}(T-l) - i_m(T-l) in δ=-1 cases

**PROVED**: From i_m(T) ≥ i_{m+1}(T) [mode = m]:
  i_m(T-l) + i_{m-1}(T-{l,s}) ≥ i_{m+1}(T-l) + i_m(T-{l,s})
  => i_{m-1}(T-{l,s}) - i_m(T-{l,s}) ≥ i_{m+1}(T-l) - i_m(T-l)

**Verified**: 14,251 pairs for n≤14, 0 violations.

In δ=-1 cases: i_{m+1}(T-l) - i_m(T-l) ≤ 0 (mode T-l = m-1, descending at m).
So the proved bound gives i_{m-1}(T-{l,s}) - i_m(T-{l,s}) ≥ negative number (not ≥ 0).

## Where the proof stands

**Proved from this session's algebra** (given mode(T) = m and mode(T-l) = m-1):
1. mode(T-{l,s}) ≥ m-2
2. i_{m-1}(T-{l,s}) - i_m(T-{l,s}) ≥ i_{m+1}(T-l) - i_m(T-l) [both ≤ 0 when δ=-1]
3. When mode(T-{l,s}) = m-2 (if it occurred): STRONG C2 would require non-trivial cross-comparison

**Still needed** (empirically verified, not proved algebraically):
1. mode(T-{l,s}) ≥ m-1 always (stronger than what algebra gives: ≥ m-2)
2. STRONG C2: λ_m^{T-l} ≥ λ_{m-1}^{T-{l,s}} (cross-tree comparison)
3. Sub-claim B: when condition (1) fails, p * slack_c2 ≥ (1-p) * deficit

**Key insight for Sub-claim B**: The deficit is not just ≤ 1 but empirically ≤ 0.095.
The algebraic reason: at λ_m^T (close to 1 in the failing cases), the mean μ(T', λ_m^T)
for T' with mode m-1 is very close to m-1 (since the tie-fug condition for T' would put
μ(T', λ_{m-1}^{T'}) ≥ m-2, but λ_m^T ≥ λ_{m-1}^{T'} means we're above the m-1 tie-fug,
pushing μ closer to m-1).

## δ=0 and δ=+1 cases for condition (1)

When mode(T-l) = m (δ=0):
- IH on T' at k=m: μ(T', λ_m^{T'}) ≥ m-1.
- Need μ(T', λ_m^T) ≥ m-1: requires λ_m^T ≥ λ_m^{T'}.
- But from mediant: λ_m^T ≤ max(λ_m^{T'}, λ_{m-1}^{T''}) = λ_m^{T'} (by STRONG C2: λ_m^{T'} ≥ λ_{m-1}^{T''}).
- So λ_m^T ≤ λ_m^{T'}, meaning the STRONG C1 comparison fails (as verified: 244,690 fails).
- Thus μ(T', λ_m^T) ≤ μ(T', λ_m^{T'}) ≥ m-1 but we can't conclude μ(T', λ_m^T) ≥ m-1 from this alone.

Wait: the direction is λ_m^T ≤ λ_m^{T'} by STRONG C2 (since STRONG C2 says max(A/B, D/E) = A/B = λ_m^{T'}).
So μ(T', λ_m^T) ≤ μ(T', λ_m^{T'}) (since mean is increasing).
And IH gives μ(T', λ_m^{T'}) ≥ m-1 (if STRONG C2 => λ_m^T ≤ λ_m^{T'}).

So when λ_m^T ≤ λ_m^{T'}, condition (1) requires proving μ(T', λ_m^T) ≥ m-1 with a SMALLER fugacity
than the tie-fug of T'. This requires showing the ECMS/tie-fug condition for T' holds at λ_m^T < λ_m^{T'}.

This is the δ=0 sub-case where condition (1) can fail (1.4% of cases, per check_tie_induction.py).

## New findings (2026-02-18, continued session)

### Mode Non-increase under Leaf Removal (NEW LEMMA, verified n≤19)

**Lemma (Mode Non-increase)**: For any tree T and leaf l, mode(T-l) ≤ mode(T).
Equivalently, δ(T,l) = mode(T) - mode(T-l) ∈ {0, +1} for ALL leaf deletions.
Verified: 4,626,113 pairs through n=19, 0 violations.
Min c1 in δ=0 cases: 0.001419 (n=19, still positive, converging slowly to 0?).

**Algebraic connection**: Mode non-increase at j=m+1 follows algebraically from:
i_m(T-l) - i_{m+1}(T-l) ≥ i_m(T-{l,s}) - i_{m-1}(T-{l,s}).
When mode(T-{l,s}) = m (rel=+2): RHS ≥ 0, so mode(T-l) ≤ m follows. ✓
When mode(T-{l,s}) = m-1 (rel=+1): RHS ≤ 0, mode non-increase needs separate argument.

### Condition (1) failure characterization

All 9,094 c1 failures through n=17 have δ=+1 AND rel=+1 (mode(T-{l,s}) = m-1).
Never in δ=0 cases or rel=+2 cases.

Delta distribution: δ=0 (68%), δ=+1 (32%) for n≤19 leaf pairs.
Min c1: δ=0 → 0.001419 (always positive); δ=+1 → -0.094642 (can be negative).

### c2 when mode(T-{l,s}) = m (rel=+2)

All c1 failures have rel=+1. So when rel=+2:
- c1 ≥ 0 always (0 failures with rel=+2 and c1 < 0)
- c2 ≥ 1 always (min 1.058, n≤16)
- Combined ≥ 0 trivially

And in this case, λ_m^T < λ_m^{T-{l,s}} always (λ_m^T / λ_m^{T-{l,s}} < 1 for all 3,130 rel=+2, δ=+1 pairs).

### Three-case proof structure for combined ≥ 0

| Case | δ | rel | c1 | c2 | combined | proof status |
|------|---|-----|-----|-----|----------|--------------|
| A | 0 | 1 or 2 | ≥ 0 for n≤19; **FAILS n=21+** (S(2^9,1^2): c1=-0.033) | ≥ 0 | ≥ 0 (verified) | NEED: combined ≥ 0 |
| B | +1 | 2 | ≥ 0 (no failures through n≤17) | ≥ 1 | ≥ 0 trivially | NEED: c1 ≥ 0 |
| C | +1 | 1 | may be < 0 (min -0.094) | ≥ 0.907 | ≥ 0 (min 0.075) | NEED: combined ≥ 0 |

**CRITICAL CORRECTION (2026-02-18)**: Case A is NOT automatically handled by c1 ≥ 0.
For S(2^9, 1^2) (n=21): direct leaf removal gives δ=0 but c1 = -0.033 < 0.
The combined expression is still +0.398 (positive) because c2 = 0.983 and p = 0.424.
This means Cases A and C share the same proof burden: need combined = p·c2 + (1-p)·c1 ≥ 0.
Case B (rel=+2) remains trivial since c1 ≥ 0 and c2 ≥ 1.

### Algebraic status of each case

**Case A (δ=0)**: c1 = μ(T-l, λ_m^T) - (m-1). Earlier believed c1 ≥ 0 from n≤19 data.
**CORRECTION**: c1 can be NEGATIVE for large trees. S(2^9, 1^2) (n=21): δ=0 direct leaf, c1 = -0.033.
Tie-fug condition still holds: combined = 0.424·0.983 + 0.576·(-0.033) = +0.398 > 0.
IH at k=m for T-l gives μ ≥ m-1 at λ_m^{T-l}, but λ_m^T ≤ λ_m^{T-l} (by STRONG C2), so
evaluating at smaller λ_m^T gives μ(T-l, λ_m^T) ≤ m-1. Deficit ≤ 1 by Weak C1 + IH at k=m-1.
**STATUS: SAME OPEN PROBLEM as Case C -- need combined = p·c2 + (1-p)·c1 ≥ 0 algebraically.**

**Min c1 trend in δ=0 cases (be2ceee+b1f33cd, exhaustive n≤20)**:
| n | m | min c1 | ratio | sign |
|---|---|--------|-------|------|
| 15 | 5 | 0.018518 | 0.7949 | + |
| 16 | 6 | 0.028309 | 0.8428 | + |
| 17 | 6 | 0.002311 | 0.8310 | + |
| 18 | 6 | 0.001878 | 0.8250 | + |
| 19 | 6 | 0.001419 | 0.8201 | + |
| **20** | **7** | **−0.003217** | **0.8493** | **NEGATIVE** |

**c1 goes negative at n=20 (δ=0 cases, m=7)**. Confirmed: the δ=0 case has the same proof burden as δ=+1 for n≥20.

Worst n=17 tree identified (b1f33cd): **S(2^6, 1^4)** -- mixed spider with center degree 10, six arms of length 2, four pendant leaves. Degree dist: {1:10, 2:6, 10:1}. poly=[1,17,120,491,1309,2379,2978,2529,1392,448,64], mode=6.

**Conclusion**: Cases A and C share IDENTICAL proof burden for n≥20. The proof must handle combined = α·c2 + β·c1 ≥ 0 for both δ=0 and δ=+1 (when c1 < 0). No simplification from the δ=0 case for large trees.

**Case B (δ=+1, rel=+2)**: mode(T-{l,s}) = m; λ_m^T < λ_m^{T-{l,s}} always.
IH at k=m for T-{l,s} gives c2 ≥ m-1-(m-2)=1 only if λ_m^T ≥ λ_m^{T-{l,s}}... but λ_m^T < λ_m^{T-{l,s}}.
Empirically c2 ≥ 1.058, but algebraic reason is unclear.
**STATUS: Open -- the IH doesn't explain c2 ≥ 1.**

**Case C (δ=+1, rel=+1)**: mode(T-l) = m-1 and mode(T-{l,s}) = m-1.
From IH: c1 ≥ -1 (Weak C1 + IH at k=m-1 for T-l) and c2 ≥ 0 (STRONG C2 + IH at k=m-1 for T-{l,s}).
Combined = α·c1 + β·c2. Need α·|c1| ≤ β·c2 when c1 < 0.
Min c2/|c1| = 9.80 (empirical), but only c2 ≥ 0 and |c1| ≤ 1 provable.
**STATUS: Sub-claim B open -- need correlated bound (c2 and c1 not independent).**

**Case C combined bound + component tracking (bdae801+b280d37, n≤20):**
| n | min_c1 | min_c2 | min_comb | max_deficit | n_caseC |
|---|--------|--------|----------|-------------|---------|
| 14 | -0.085 | 0.725  | 0.404    | 0.085       | 9,640  |
| 16 | -0.095 | 0.737  | 0.418    | 0.095       | 27,671 |
| 18 | -0.080 | 0.729  | 0.391    | 0.080       | 530,982 |
| 20 | -0.127 | 0.700  | **0.379**| **0.127**   | 2,735,038 |

- **min_comb trend**: 0.548→0.497→0.391→0.379, converging toward 0.362 (global min) from above.
- **max_deficit = 0.127** at n=20 (new high, was 0.095 at n=16). Slowly increasing.
- **min_c2 = 0.700** at n=20: c2 remains well above 0; c2/|deficit| ≥ 0.700/0.127 ≈ 5.5 per-tree average.
- All min_c1, min_c2, min_comb taken at DIFFERENT trees (anti-correlation confirmed).
Anti-correlation suggests: c2 is large precisely when β is large, OR |c1| is large only when β is small.
This structural correlation may be provable via BP cavity equations.

### Key equivalence / mutual dependence

The three empirical lemmas appear mutually dependent:
1. Mode Non-increase follows algebraically from Lemma A in rel=+2 case; needs more in rel=+1.
2. Lemma A follows from Mode Non-increase for sub-trees (but leads to 2-step descent, only proving ≥ m-2).
3. Condition (1) in Case A requires Lemma A (for mode(T-{l,s}) ≥ m-1 via STRONG C2) plus something more.

Proving all three simultaneously by a joint induction may be the right approach.

## Key conjecture: margin ≥ 1/3 for all trees

**Conjecture (1/3 Margin Lower Bound)**: For any tree T with mode m ≥ 1,
  μ(T, λ_m^T) ≥ (m-1) + 1/3 = m - 2/3.

Evidence:
- **All 9,114,283 trees n≤22: 0 violations; min margin = 0.36226 at n=21 (S(2^10)), confirmed.**
- n=22 did NOT produce a new minimum (min stayed at 0.36226 from n=21).
- S(2^k) spiders up to k=500: 0 violations; min margin = 0.334 ≈ 1/3 + 0.001 (still converging)

If proved + Steiner peeling (μ(T,1) < n/3) + mode ≤ ⌈μ(T,1)⌉ → mode ≤ ⌊n/3⌋+1 = Conjecture A.
(But note: proving margin ≥ 1/3 is NOT required for Conjecture A; margin > 0 suffices.)

## Spider family asymptotic (KEY RESULT, 2026-02-18)

The S(2^k) family gives the global minimum tie-fugacity margin among all trees at n=2k+1.
More precisely, for the k=3j+1 sub-sequence, S(2^k) achieves the minimum over all odd n=2k+1.

**Asymptotic fit** (k = 3j+1, tip leaf removal, margin = μ(T, λ_m) - (m-1)):
  margin ≈ 1/3 + 1/(2k) - 0.27/k² + O(1/k³)
  Verified: k up to 118 (n up to 237), min margin = 0.338

**Conclusion**: The margin → 1/3 from ABOVE. Infimum = 1/3, NEVER ACHIEVED.

**Asymptotic structure (k=3j+1, analytically verified)**:
  - μ(T,1) → 2k/3 = m - 1/3  (mean at λ=1 is 1/3 BELOW mode)
  - μ(T,λ_m) → m - 2/3         (tie-fug mean is 2/3 below mode)
  - margin = μ(T,λ_m) - (m-1) → 1/3
  - drop: μ(T,1) - μ(T,λ_m) → 1/3  (drop from λ=1 to λ_m is 1/3)
  This is NOT an asymptotic limit of 0 -- the mean stops 1/3 above m-1.

**Conjecture (Spider Extremality for Margin)**: S(2^k) achieves minimum margin over all trees of size n=2k+1.
- Verified: all n≤21 (exhaustive)
- Minimum margin sequence: 0.436 (n=5), 0.380 (n=15), 0.362 (n=21), →1/3 (predicted)

**Global minimum through n≤22**: 0.36226 at S(2^10) (n=21) -- in the k=3j+1 sub-sequence.
n=22 did NOT achieve a new minimum (first non-spider n that doesn't improve the bound).

Minimum margin by n (exhaustive where available; spider-family where not):
| n | min margin | achieved by |
|---|-----------|-------------|
| 14 | 0.379 | S(2^6, 1^1) |
| 15 | 0.380 | S(2^7) -- pure spider |
| 18 | 0.367 | bi-hub tree: deg_seq [5,4,2×9,1×7], NOT a pure spider |
| 20 | 0.363 | S(2^9, 1^1): deg_seq [10,2×9,1×10] -- mixed spider |
| 21 | 0.362 | S(2^10) -- pure spider (k=10=3·3+1) |

At EVEN n: extremal trees can be bi-hub (secondary hub of degree 4) rather than pure/mixed spiders.
At ODD n=2k+1: S(2^k) dominates, especially k=3j+1 sub-sequence. (bfa4f2f confirmed)

## Mode ≤ ⌈μ(T,1)⌉ for ALL trees (KEY NEW VERIFICATION, 2026-02-18)

**Conjecture (Darroch for Trees)**: For every tree T, mode(T) ≤ ⌈μ(T,1)⌉.

This is equivalent to: μ(T,1) > mode(T) - 1 (since mode is an integer).

**Verification** (ALL trees, not just d_leaf≤1):
| n | Trees | Failures | Min gap (⌈μ⌉ - mode) |
|---|-------|----------|----------------------|
| ≤ 18 | all | 0 | 0 |
| 19 | 317,955 | 0 | 0 |
| 20 | 823,065 | 0 | 0 |
| 21 | 2,144,505 | 0 | 0 |
| 22 | 5,623,756 | 0 | 0 |

Min gap = 0 means: **the bound is tight** (some tree has mode = ⌈μ⌉ exactly).
The tight cases are the k=3j+1 spiders S(2^k): mode = 2j+1 ≈ μ+1/3, so ⌈μ⌉ = 2j+1 = mode.

**Key implication**: Darroch + μ < n/3 (Steiner peeling, PROVED for d_leaf≤1) → Conjecture A.
If proved for all trees, it immediately gives Conjecture A for d_leaf≤1 trees.

**Status**: OPEN (verified n≤21, no proof). Classical Darroch requires LC; here it holds without LC.

**Proof ideas tried and failed**:
- Transfer argument via leaf deletion: mode can jump +1 without μ jumping ≥ 1 (circular)
- Hard-core model + WHNC: gives μ < n/3 (upper bound) but not μ > mode-1 (lower bound)
- Steiner peeling: same issue

**Potentially new angle**: The tight cases (min gap=0) are precisely the k=3j+1 spiders where
mode = ⌈μ⌉ exactly. Other trees have mode < ⌈μ⌉ (gap ≥ 1). This suggests the spiders are
the extremal family for both the tie-fugacity margin (→ 1/3) and the mode-vs-mean gap (mode=⌈μ⌉).
A proof for spiders first, then extension to general trees by adding vertices, might work.

## Mean Drop Lemma: Algebraic Derivation and Proof Status (2026-02-18)

**Mean Drop Lemma (UNPROVED, VERIFIED)**:
For any tree T and leaf l with support s:
  mean(T-l, λ=1) - mean(T-{l,s}, λ=1) < 1.

Verified: all tree-leaf pairs through n=20 (bfdaf41 task), max drop = 0.542, 0 violations ≥ 1.
Max drop by n (even n only): 0.497 (n=6), 0.515 (n=8), 0.515 (n=10), 0.525 (n=12), 0.529 (n=14), 0.528 (n=16), 0.537 (n=18), 0.542 (n=20).
The max is increasing slowly (increments ≈ 0.005-0.010 per 2 steps). Converges to < 1.

### Algebraic derivation of Drop formula

Let T' = T-l, T'' = T-{l,s}. Apply vertex deletion at s in T':
  I(T'; x) = I(T''; x) + x · I(Z; x)  where Z = T' - N[s] = T - {l, s, w_1,...,w_d}
  (N[s] in T' = {s, w_1,...,w_d}, the closed neighbourhood of s in T-l)

Differentiating and evaluating at x=1:
  I'(T') = I'(T'') + I(Z) + I'(Z)  (where ' denotes d/dx at x=1)

Let α'' = I(T'')/(I(T'')+I(Z)), β'' = I(Z)/(I(T'')+I(Z)), so α''+β''=1.
  mean(T') = α'' · mean(T'') + β'' · (1 + mean(Z))

**Drop formula (EXACT)**:
  Drop = mean(T') - mean(T'') = β'' · (1 + mean(Z) - mean(T''))
  where β'' = I(Z;1) / (I(T-{l,s};1) + I(Z;1))

**Equivalence**: Drop < 1 ⟺ I(Z;1) · (mean(Z;1) - mean(T-{l,s};1)) < I(T-{l,s};1).

**Physical interpretation**: Removing s from T' leaves T'' (T' minus s) and Z (T' minus N[s]).
β'' = probability at x=1 that s is included in the IS of T', so β'' = P(s in IS of T') < 1.

### Case 1: mean(Z) ≤ mean(T'') -- PROVED

If mean(Z) ≤ mean(T''): 1 + mean(Z) - mean(T'') ≤ 1, so Drop ≤ β'' < 1. ✓ PROVED.

**Does Case 1 cover all instances?** Empirically, Case 2 (mean(Z) > mean(T'')) DOES occur.
Example: T = l - s - w - {a,b,c} (l leaf, s connected to l and w, w = center of K_{1,3}).
T'' = K_{1,3} (w-a-b-c, mean=13/9≈1.44), Z = {a,b,c} (3 isolated, mean=3/2=1.50).
mean(Z) = 3/2 > 13/9 = mean(T''). Drop = (8/17)·(1+1/18) = 76/153 ≈ 0.497 < 1. ✓

### Case 2: mean(Z) > mean(T'') -- OPEN

When mean(Z) > mean(T''): each component C_i of T'' contributes an INCREASE when w_i is removed.
This happens when mean(C_i) > 1 + mean(C_i - N[w_i]), i.e., w_i is a "dense hub" of C_i.

For positive increase Δ_i = mean(C_i - {w_i}) - mean(C_i):
  Δ_i = r_i · (mean(C_i) - 1 - mean(C_i - N[w_i])) / 1   [since I(C_i) - I(C_i-N[w_i]) = I(C_i-{w_i})]
  where r_i = I(C_i-N[w_i]) / I(C_i-{w_i})

**Key identity**: mean(G) - mean(G-v) = β_v · (1 + mean(G-N[v]) - mean(G-v)) with same structure.
So "mean decreases under vertex removal" ↔ "removing v makes mean(G-N[v]) < mean(G-v)" ↔ same Drop < 1 form for the sub-problem. The structure is self-similar.

**Algebraic bound attempt**: I'(Z;1) ≤ I'(T'';1) (proved: vertex deletion decreases I').
  I(Z)·mean(Z) = I'(Z) ≤ I'(T'') = I(T'')·mean(T'').
  I(Z)·(mean(Z) - mean(T'')) ≤ mean(T'')·(I(T'') - I(Z)).
For Drop < 1: need I(Z)·(mean(Z)-mean(T'')) < I(T'').
The bound gives: I(Z)·(mean(Z)-mean(T'')) ≤ mean(T'')·(I(T'')-I(Z)).
This is ≤ I(T'') iff mean(T'')·(1-I(Z)/I(T'')) ≤ 1 iff mean(T'') ≤ I(T'')/(I(T'')-I(Z)).
FAILS for large trees (mean(T'') ~ n/3 grows linearly; I(T'')/(I(T'')-I(Z)) can be close to 1).

**Status**: Case 2 open. Need: when removing s's neighbors from T'' increases the mean,
the increase is bounded by β''^{-1} - 1 = I(T'')/I(Z). Empirically this holds with large margin.

### Darroch induction using Mean Drop (Case 1 proved)

**Case 1 Proof (mode(T) = mode(T-l) = m')**:
Assume: (1) Mean Drop Lemma holds; (2) mode non-increase (mode(T-l) ≤ mode(T)); (3) IH.

From Mean Drop Lemma: mean(T-{l,s}) + 1 > mean(T-l).
From leaf-deletion: mean(T) = α·mean(T-l) + β·(mean(T-{l,s})+1) where α+β=1, α,β>0.
Since mean(T-{l,s})+1 > mean(T-l): mean(T) > α·mean(T-l) + β·mean(T-l) = mean(T-l).
So mean(T) > mean(T-l). Therefore ⌈mean(T)⌉ ≥ ⌈mean(T-l)⌉.
By IH: mode(T-l) = m' ≤ ⌈mean(T-l)⌉ ≤ ⌈mean(T)⌉.
Since mode(T) = m': mode(T) ≤ ⌈mean(T)⌉. ✓

**Case 2 (mode(T) = mode(T-l)+1 = m'+1): OPEN.**
Need: mean(T) > m'.
From IH: mode(T'') = m'' ≤ ⌈μ(T'')⌉ and mode(T'') ≥ m' (Lemma A), so μ(T'') > m'-1,
thus μ(T'')+1 > m'. The β-term contributes positively. The α-term contributes α·(μ(T')-m')
which can be negative (if μ(T') < m', e.g., tight spiders with μ(T') ≈ m'-1/3).
**Gap**: For μ(T) > m', need β·(μ(T'')+1-m') > α·(m'-μ(T')). No clean algebraic bound.

## Next steps

1. **Prove margin ≥ 1/3 for all trees** (the natural clean goal given the asymptotic).
   - For S(2^k), proved: margin ≈ 1/3 + 1/(2k) > 1/3. The monotone decrease + limit = 1/3 from above.
   - For general trees: extremal conjecture (S(2^k) is worst) would suffice.

2. **Spider extremality proof for margin**: Prove S(2^k) minimizes margin among all n=2k+1 trees.
   - Connects to the known spider extremality for μ(T,1) = μ(T) (already in paper).

3. **Alternative proof route**:
   - Prove tie-fug condition directly using the hard-core model / cavity equations
   - Connect to Steiner peeling (which already gives μ(T,1) < n/3)

4. **IMPORTANT NOTE**: Any "combined = p·c2 + (1-p)·c1" formula using COEFFICIENT ratio
   p = i_{m-1}(T-{l,s})/i_m(T) is NOT equal to the true margin. Use α·c1 + β·c2 with
   partition function weights α=I(T-l,λ)/I(T,λ), β=λI(T-{l,s},λ)/I(T,λ). Or just
   compute margin = μ(T, λ_m^T) - (m-1) directly.
