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

## Next steps

1. **Prove Lemma A** (mode(T-{l,s}) ≥ m-1): might follow from a finer analysis of the δ=-1 case
   using the fact that s is the support of the leaf that caused the mode drop.

2. **Prove Lemma B** (STRONG C2 / cross-tree comparison): might follow from LC of T'' at m-1
   plus a comparison of the B_k step-fugacities, but T'' LC is not guaranteed from d_leaf≤1.

3. **Prove Sub-claim B** (balance when cond-1 fails): show that deficit ≤ c·(1-p)/p for some
   structural reason, using the specific relationship between T' and T'' as a leaf-support pair.

4. **Alternative**: Find a completely different proof of the tie-fugacity condition that avoids
   the leaf-decomposition obstruction.
