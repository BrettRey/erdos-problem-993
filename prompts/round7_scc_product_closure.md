# Round 7: Proving SCC Product Closure

## The Problem

We're trying to prove that the independence polynomial of every tree is unimodal (Erdős Problem #993). We've reduced the problem to showing that **Strong Condition C (SCC)** is preserved under incremental products at support vertices.

## Setup

Root a tree T at a support vertex r (adjacent to at least one leaf). The DP gives:
- E = dp[r][0] (exclude-root polynomial), J = dp[r][1]/x (include-root polynomial)
- Define e = (1+x)I = (1+x)(E + xJ), b = E

**SCC states:** Δ_k(e, b) ≥ 0 for all k, where Δ_k(f,g) = f_{k+1}·g_k - f_k·g_{k+1}.

Equivalently: the ratio e_k/b_k is nondecreasing (when b_k > 0).

**SCC implies unimodality** via the identity b_{k-1}·Δ_k = b_{k-1}·d_k + b_k·d_{k-1} + a_{k-1}·c_k where d_k = e_k - b_k - b_{k-1}, c_k = b_k² - b_{k-1}·b_{k+1} (LC gap), a_k = I_k = e_k - b_{k-1}.

**SCC is verified computationally:** 930M+ checks through n=22, ZERO failures.

## Incremental Product Structure

At a support vertex r with ℓ leaf neighbors and non-leaf children c₁,...,c_s, the polynomials are built by multiplying in one factor at a time:

```
E_acc = ∏ I(T_{c_i}),  J_acc = ∏ E(T_{c_i})
```

Each factor: P = I_c = E_c + xJ_c, Q = E_c, R = J_c (so P = Q + xR).

The transformation at each step:
```
E_new = E_acc · P,     J_new = J_acc · Q
e_new = Q · e_old + x(1+x)R · b_old    [PROVED identity]
b_new = P · b_old
```

This is a convolution matrix M = [[Q, x(1+x)R], [0, P]] acting on the pair (e, b).

**All three properties SCC, LC(E), J≤E are verified at EVERY intermediate stage** (930M checks, 0 failures).

## What's Proved

### Karlin's TP2 Closure Theorem
If A is PF2 (nonneg + LC) and f ≥_{lr} g, then A*f ≥_{lr} A*g.

### SCC of Leaf Factor
For a single leaf factor: P = 1+x, Q = 1, R = 1. Then SCC = Δ_k((1+x)I, E) = c_k(E) (LC gap). Since E is LC at all stages (verified, and PF2 products remain PF2), leaf factors trivially preserve SCC.

## Decompositions (all verified through n=22)

### B decomposition (best 2-term)
```
SCC = B1 + B2
B1 = Δ_k(e_new, E_old·Q)     ← ALWAYS ≥ 0 (930M checks, 0 failures)
B2 = Δ_k(e_new, x·E_old·R)   ← can be negative
Min B1/|B2| when B2 < 0: 2.41 (infimum = 2, balanced spider)
```

### B1 sub-decomposition (Karlin proof)
```
B1 = T1 + T2
T1 = Δ_k(E_old·(1+x)I_c, E_old·E_c)  ← PROVED ≥ 0 by Karlin
     (E_old is PF2, (1+x)I_c ≥_{lr} E_c by SCC of factor subtree)
T2 = Δ_k(x(1+x)·J_old·E_c, E_old·E_c) ← can be negative
Min T1/|T2| when T2 < 0: 3.66 (n=22)
```

### 4-term bilinear decomposition
```
SCC = QQ + QR + RQ + RR
QQ = Δ_k(e_old·Q, E_old·Q)           ← PROVED ≥ 0 (Karlin + SCC_old)
QR = Δ_k(e_old·Q, x·E_old·R)         ← usually negative
RQ = Δ_k(x(1+x)·E_old·R, E_old·Q)   ← usually positive
RR = c_{k-1}(E_old·R)                ← PROVED ≥ 0 (LC closure)
```

Cross = QR + RQ: positive 98.7% of the time. When cross < 0:
- Min (QQ+RR)/|cross| = 2.09 (infimum = 2)
- RR does 96% of the work (QQ only 4% in tight cases)

## The Gap

We have **two provably nonneg terms** (QQ from Karlin, RR from LC) and need them to dominate the cross terms when the cross is negative. The ratio (QQ+RR)/|cross| ≥ 2.09 with infimum 2.

### Why the infimum is 2

The factor of 2 comes from e = (1+x)I. The (1+x) multiplier makes e "twice as heavy" as b at the mode. Formally, Δ_k((1+x)f, f) = c_k(f) (the LC gap), so the (1+x) contributes one full LC-gap worth of positivity.

For the balanced spider S(2^s) at k=1:
- B1 = 4s² + 5s - 5
- B2 = -2s² + s + 5
- SCC = 2s² + 6s
- B1/|B2| → 2 as s → ∞

### What doesn't work

1. **Cauchy-Schwarz** (QQ·RR ≥ QR·RQ): FALSE (min ratio 0.6)
2. **Additive cone closure**: SCC cone is NOT closed under addition of pairs
3. **Factor-level curvature budget**: |L_{p,q}| ≤ Σc_t fails at n≥15
4. **Generic product closure** (without tree-realizability): FALSE
5. **HWZZ partial synchronicity**: FALSE for tree (I,E) pairs

### What might work

1. **The (1+x) amplification**: e = (1+x)I gives one extra LC-gap of margin. This is what makes the infimum 2 rather than 1. Can we formalize "the (1+x) factor always provides enough margin"?

2. **Weighted Turán inequality**: Is there a weighted version of SCC that is EXACTLY preserved (zero gap)? Then actual SCC would have margin ≥ the weight excess.

3. **Semigroup cone**: M = [[Q, x(1+x)R], [0, P]] acts by convolution. Can we show the SCC cone is M-invariant using properties of Q, R, P beyond just nonnegativity/LC?

4. **Two-step approach**: Step 1: convolve with Q (preserves SCC by Karlin). Step 2: add the R-correction. Can the R-correction be bounded relative to the Q-step output?

5. **Coefficient-level analysis at the mode**: The tight cases all have k ≈ mode(E_old·P). Near the mode, the ratio e_k/b_k = (1+x)·(1 + xJ_k/(E_k))... what structural constraints apply here?

## Specific Questions

1. For the 4-term decomposition, is there a TP2/total positivity result that bounds |cross| ≤ c·diag with c < 1 for structured (e, b) pairs?

2. The upper-triangular structure M = [[Q, x(1+x)R], [0, P]] with positive (off-)diagonal entries: are there cone-preservation results for such semigroups in the TP2 literature?

3. Since Q = E_c and R = J_c satisfy Q ≥ R coefficientwise (J ≤ E for tree DPs), does this constraint on the matrix entries help? The (1,2) entry x(1+x)R is "smaller" than x(1+x)Q.

4. The fact that SCC, LC(E), and J≤E ALL hold at every intermediate stage (not just the final product) is extremely strong. Can we use all three as a JOINT inductive hypothesis?

## Data Available

- 930M+ checks through n=22 (all 9.1M trees), zero SCC failures
- Min ratios: B1/|B2| = 2.41, diag/|cross| = 2.09, T1/|T2| = 3.66
- Balanced spider exact formulas: B1 = 4s²+5s-5, B2 = -2s²+s+5
- All tight cases: 2 factors, k ≈ 10, stage 2
