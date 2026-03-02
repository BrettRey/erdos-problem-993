> **PARTIALLY OBSOLETE (2026-03-01):** The "Condition C" (SCC) product closure path discussed here fails at n=28. See `notes/scc_false_n28_2026-03-01.md`. Scan data remains valid; proof strategy is dead.

# Round 4 Scan Results (2026-02-28)

## Setup

Three scans run in a single pass over all 9,114,283 trees with n <= 22 (same dataset as P* scan).
Script: `scan_p2_round4.py`. Runtime: 1456s (24 min) with 8 workers.

## Scan 1: SV1 (Single Sign Variation of LR minors)

**Conjecture:** d_k = a_{k+1}*b_k - a_k*b_{k+1} has at most 1 sign change on prefix {0,...,m-1}.

**Result: FALSE.** 9,790 trees have > 1 sign change (all have exactly 2; at the support-vertex level, max 3 changes seen). First failure at n=10.

SV1 failure counts by n:

| n | Failures | Max changes |
|---|----------|-------------|
| 10 | 1 | 2 |
| 11 | 1 | 2 |
| 12 | 2 | 2 |
| 13 | 8 | 2 |
| 14 | 13 | 3 |
| 15 | 18 | 3 |
| 16 | 52 | 3 |
| 17 | 115 | 3 |
| 18 | 201 | 3 |
| 19 | 476 | 3 |
| 20 | 1,209 | 3 |
| 21 | 2,274 | 3 |
| 22 | 5,420 | 3 |

Note: d_k < 0 at some support vertex is common (8.4M of 59.9M checks, ~14%). The sign pattern can be ++-+ or ++-++-, but P2 still holds because the (1+x) smoothing compensates.

**Verdict: Kill the SV1 approach. Do NOT send to GPT 5.2.**

## Scan 2: Canonical Max-ell Rooting

**Conjecture B2:** P* holds at r* = argmax_v ell(v) (vertex with most leaf neighbors).

**Result: 0 FAILURES.** Clean across all 9.1M trees.

ell distribution:

| max ell | Trees | Fraction |
|---------|-------|----------|
| 1 | 403,400 | 4.4% |
| 2 | 3,163,272 | 34.7% |
| 3 | 3,060,965 | 33.6% |
| 4 | 1,505,580 | 16.5% |
| 5 | 608,243 | 6.7% |
| 6+ | 372,823 | 4.1% |

**Key observation:** 95.6% of trees have max ell >= 2, meaning E = (1+x)^2 * A with double smoothing. Only 4.4% of trees are "skinny" (max ell = 1).

**Implications for GPT 5.2 Prompt 2:**
- For ell >= 2, the double (1+x)^2 smoothing may be strong enough to close P2 from A >= B + B LC alone.
- For ell = 1, the single (1+x) smoothing suffices computationally but needs the full Delta identity argument.

## Scan 3: Condition C (Delta Identity)

**Condition C:** b_{k-1} * Delta_k = b_{k-1}*d_k + b_k*d_{k-1} + a_{k-1}*c_k >= 0.

**Result: 0 FAILURES.** Clean across all 59,916,124 support-vertex checks, all trees n <= 22.

This is the integer version of the decomposition:
    Delta_k = d_k + (b_k/b_{k-1})*d_{k-1} + (a_{k-1}/b_{k-1})*c_k

where c_k = b_k^2 - b_{k-1}*b_{k+1} (LC gap of B).

**This is the strongest result.** It says: at every support vertex of every tree, P2 follows from this 3-term decomposition. Proving Condition C algebraically would close the entire P2 problem.

### What Condition C requires for an algebraic proof:

The three terms in b_{k-1}*Delta_k are:
1. **b_{k-1}*d_k**: current LR minor (can be negative)
2. **b_k*d_{k-1}**: previous LR minor weighted by b_k (memory term, can be negative)
3. **a_{k-1}*c_k**: LC gap of B weighted by a_{k-1} (curvature bonus, always nonneg when B is LC)

Since a_{k-1} >= b_{k-1} (A >= B coefficientwise), the curvature bonus is at least b_{k-1}*c_k = b_{k-1}(b_k^2 - b_{k-1}*b_{k+1}).

**Product closure is the key challenge:** A and B are products of factor pairs (I_c, E_c), and the LR minors of products involve cross terms. The question for GPT 5.2 Prompt 3 is whether Condition C propagates through multiplication.

## Factor-Level Analysis (follow-up scans)

### Factor-level ratio dominance: FAILS

496,210 of 1,614,001 factors (~31%) have some d_k^(c) < 0. First example: P_3 factor I_c = 1+3x+x^2, E_c = 1+2x+x^2, with d_1 = -1.

### Factor-level Condition C: 0 FAILURES

11,941,358 factor pairs (I_c, E_c) from all trees n <= 20, every factor satisfies Condition C. Also: every factor E_c is LC.

### Product closure of Condition C: 0 FAILURES

- 1,184 unique factor pairs from n <= 12
- 701,520 exhaustive pairwise products: ALL satisfy Condition C
- 19,535 unique factors from n <= 15, 100K random products: ALL pass

**Conclusion:** Condition C holds at every level of the induction:
1. At the leaf level (trivially: I = 1+x, E = 1)
2. At each subtree factor level (verified)
3. At the product level (verified)
4. At the root level (verified)

### The proof structure

To prove unimodality:
1. Prove Condition C holds for each factor pair (I_c, E_c). This is an inductive claim: it holds for factors one level down.
2. Prove product closure: if (I_1, E_1) and (I_2, E_2) satisfy Condition C with E_1, E_2 log-concave and I_i >= E_i coefficientwise, then (I_1*I_2, E_1*E_2) also satisfies Condition C.
3. Use Condition C + the identity to get P2 at support vertices.
4. P3 is already proved. P⋆ = P2 + P3 implies unimodality.

**The key algebraic challenge is step 2.** The product closure question requires understanding how the 3-term sum b_{k-1}*d_k + b_k*d_{k-1} + a_{k-1}*c_k behaves under Cauchy products.

## Summary

| Conjecture | Status | Next step |
|------------|--------|-----------|
| SV1 | **FALSE** | Kill (don't send prompt) |
| Canonical max-ell (B2) | **VERIFIED** n<=22 | Send Prompt 2 |
| Condition C | **VERIFIED** n<=22 | Send Prompt 3 (highest priority) |

**Priority ranking for GPT 5.2 round 4:**
1. Prompt 3 (Condition C + product closure) -- highest value, directly closes P2 if proved
2. Prompt 2 (canonical max-ell) -- useful structural insight, may simplify the ell>=2 case
3. Prompt 1 (SV1) -- CANCEL, conjecture is false
