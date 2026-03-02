# E ≽ J Universal Ratio Dominance (2026-03-01)

## Discovery

**E ratio-dominates J at ALL k, not just the prefix, at every support vertex of every tree.**

Previously, P2 was verified only in the prefix (k ≤ mode-1). The new result:

| Metric | Value |
|--------|-------|
| P2 (E ≽ J) fails at ANY k | **0** |
| Trees checked | 9.1M+ (n ≤ 22) |
| Support vertices | 59.9M+ |
| Total checks | 907M+ |

E ≽ (1+x)J fails ~12% of checks, so ratio dominance does NOT extend to (1+x)J.

## Definition

E ≽ J means: E_{k+1}·J_k ≥ E_k·J_{k+1} for ALL k ≥ 0.
Equivalently: the ratio J_k/E_k is nonincreasing (for E_k > 0).

At a support vertex r:
- E = dp[r][0] = (1+x)^ℓ · ∏ I(T_{c_j})
- J = dp[r][1]/x = ∏ dp[c_j][0]

In IS counting terms: E_k = #{IS of size k excluding r}, J_k = #{IS of size k+1 including r}.

## 2-Term SCC Decomposition

**Identity (PROVED):** Δ_k = c_k + LR_k

where:
- Δ_k = ((1+x)I)_{k+1}·E_k - ((1+x)I)_k·E_{k+1} (SCC)
- c_k = E_k² - E_{k-1}·E_{k+1} (LC gap of E, always ≥ 0)
- LR_k = E_k·(J_k+J_{k-1}) - E_{k+1}·(J_{k-1}+J_{k-2}) (LR minor of E vs (1+x)J)

**Proof:** (1+x)I = (1+x)E + x(1+x)J. By linearity of minors:
Δ_k = Δ_k((1+x)E, E) + Δ_k(x(1+x)J, E) = c_k + LR_k. □

### Statistics (n ≤ 20)

| Metric | Value |
|--------|-------|
| Identity failures | **0** |
| c_k < 0 | **0** |
| LR_k < 0 | 12.3% of checks |
| SCC < 0 | **0** |
| Min utilization (c_k+LR_k)/c_k | 0.200 |
| Avg c_k/|LR_k| when LR<0 | 19.8 |

Tightest: star K_{1,19} at k=1, c_k=190, LR_k=-152, Δ_k=38.

### Pendant-Star Extremal Analysis

The tightest tree for both SCC and E ≽ J correction ratio is the **pendant-star**:
leaf u₁ -- support vertex r -- star center c -- {m leaves}, where m = n-3.

At the single incremental stage (root at r), at k=1:
- main = m(m+1)/2 (Karlin main part, from LC of (1+x)^{m+1})
- corr = m(3-m)/2 (correction, negative for m ≥ 4)
- total = 2m (absolute margin)
- **ratio = (m+1)/(m-3) = (n-2)/(n-6) → 1** as n → ∞

Verified: EXACT match at every n from 7 to 20 (minimum over all 823K trees).

**Implications:**
- **Ratio-based bounds are DEAD** (ratio → 1, not bounded away from 1)
- **Absolute margin grows linearly** (2m = 2(n-3) → ∞)
- Same extremal family as SCC 2-term decomposition (same c_k, LR_k formulas)
- Ratio decays as 1 + 4/(n-6), polynomially (not exponentially)

## Incremental Preservation

E^{(t)} ≽ J^{(t)} holds at EVERY incremental product stage.

| Metric | Value |
|--------|-------|
| E^{(t)} ≽ J^{(t)} fails | **0** (11.9M stages, n ≤ 20) |
| Karlin main part A ≽ J^{(t)} | **0** failures |
| Correction xB negative | 961K stages (~8%) |
| Factor I_t ≽ E_t | FAILS (~30%) |
| Factor E_t ≽ J_t | FAILS (~14%) |
| Min main/|corr| (when corr < 0) | 1.286 = 9/7 (pendant-star n=20) |
| Histogram: ratio > 10 | 95.2% of neg-corr checks |
| Histogram: ratio 1-1.5 | 7 checks (0.0%) |

### Inductive Structure

Stage 0: E^{(0)} = (1+x)^ℓ ≽ J^{(0)} = [1]. Trivially holds.

Stage t: E^{(t)} = E^{(t-1)}·I_t, J^{(t)} = J^{(t-1)}·E_t where I_t = E_t + x·J_t.

**Karlin step:** A = E^{(t-1)}·E_t. By induction E^{(t-1)} ≽ J^{(t-1)} and E_t is PF2 (LC).
By Karlin's TP2 theorem: A = E^{(t-1)}·E_t ≽ J^{(t-1)}·E_t = J^{(t)}.

**Decomposition:**
Δ_k(E^{(t)}, J^{(t)}) = Δ_k(A, J^{(t)}) + Δ_k(x·E^{(t-1)}·J_t, J^{(t)})
- First term ≥ 0 (Karlin). PROVED.
- Second term can be negative. Needs compensation.
- 0 total failures: first term ALWAYS compensates.

### What's needed for complete proof

Prove: Δ_k(A, J^{(t)}) ≥ |Δ_k(x·E^{(t-1)}·J_t, J^{(t)})| when second term < 0.

This has the SAME structure as the SCC proof (nonneg main + negative correction).

## Combined Invariant

The following properties hold at every support vertex and every incremental stage:
1. E is LC (products preserve LC)
2. J ≤ E coefficientwise (proved: J^{(t)} = J^{(t-1)}·E_t ≤ E^{(t-1)}·E_t ≤ E^{(t)})
3. E ≽ J (verified, 0 failures)
4. SCC ≥ 0 (verified, 0 failures)

Properties 1-2 are PROVED algebraically. Properties 3-4 are verified but need proofs.

## s=1 Reduction (NEW)

When support vertex r has exactly one non-leaf child c (s=1, ℓ ≥ 1 leaves):

**E≽J follows from SCC of subtree at c** via Karlin + transitivity:
1. SCC at c: (1+x)·I_c ≽ E_c
2. Karlin: (1+x)^ℓ·I_c ≽ (1+x)^{ℓ-1}·E_c
3. (1+x)^{ℓ-1}·E_c ≽ E_c (trivial)
4. Transitivity: E = (1+x)^ℓ·I_c ≽ E_c = J. QED.

**Statistics (n ≤ 18):**
- s=1: 692,736 vertices (62.8%)
- s≥2: 410,848 vertices (37.2%)
- 0 E≽J failures at s=1; 0 SCC failures at subtree

The pendant-star is s=1, so its E≽J (the tightest case) is provable from star SCC.
Proving SCC universally gives E≽J for 63% of cases for free.

## Key Open Question

Both E ≽ J and SCC reduce to the same algebraic gap:
"a nonneg main term (from Karlin/LC) compensates a negative correction (from x·J factor)"

Can these be proved simultaneously? The combined invariant {LC, J≤E, E≽J, SCC} might be
the right inductive hypothesis if the corrections in (3) and (4) can be bounded together.

## Connection to Previous Decompositions

The 2-term SCC decomposition Δ_k = c_k + LR_k is related to the 3-term identity:
b_{k-1}·Δ_k = b_{k-1}·d_k + b_k·d_{k-1} + a_{k-1}·c_k

The 2-term form is CLEANER (1 negative term vs 2) but loses the b_{k-1} weight.
The 3-term form has the a_{k-1}/b_{k-1} amplification factor (always ≥ 1).

## Correction Ratio Profile (n ≤ 20)

| n | min ratio | formula (n-2)/(n-6) |
|---|-----------|---------------------|
| 7 | 5.000 | 5.000 |
| 10 | 2.000 | 2.000 |
| 15 | 1.444 | 1.444 |
| 20 | 1.286 | 1.286 |

Histogram of main/|corr| when corr < 0 (4.2M checks):
- ratio > 10: 95.2%, ratio 1-1.5: 0.0% (7 checks), ratio < 1: 0 (impossible)
- Pendant-star at k=1 is the ONLY minimizer family

## Diagonal Convolution Structure

The incremental step on (E, J) is DIAGONAL:
```
E_new = I_t * E_old,  J_new = E_t * J_old
```
M = diag(I_t, E_t) in convolution sense.

**Transitivity approach**: if I_t ≽ E_t and E_old PF2, then E_new ≽ E_t*E_old ≽ J_new.
BUT: I_t ≽ E_t FAILS at ~30% of factors (≠ E_t ≽ J_t, which fails ~14%).
Product structure must rescue: aggregate products always satisfy E ≽ J.

## Reformulation

SCC = c_k + LR_k can be written as:
Δ_k = E_k·S_k - E_{k+1}·S_{k-1}

where S_k = I_k + J_k = E_k + J_{k-1} + J_k = E_k + [(1+x)J]_k.

Equivalently: S_k/S_{k-1} ≥ E_{k+1}/E_k (S decreases no faster than E).

For the star: S = (1+x)^m + (1+x)·1 at the relevant indices. The condition reduces
to LC of E + boundary corrections at k=1,2.
