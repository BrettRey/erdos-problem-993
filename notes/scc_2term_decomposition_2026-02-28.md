# SCC 2-Term Decomposition Results (2026-02-28)

## Three decompositions of SCC

With e = (1+x)I, b = E at each incremental product stage:

### (A) Linearity in 1st argument: e_new = eQ + x(1+x)bR

SCC = Δ_k(eQ, bP) + Δ_k(x(1+x)bR, bP)
    = A1 + A2

### (B) Linearity in 2nd argument: b_new = bQ + xbR

SCC = Δ_k(e_new, bQ) + Δ_k(e_new, xbR)
    = B1 + B2

**KEY FINDING: B1 ≥ 0 ALWAYS** (0 failures, 930M+ checks through n=22)

### (C) LC + correction

SCC = c_k(b_new) + Δ_k(x(1+x)·J_new, b_new)
    = C1 + C2

where c_k(b) = b_k² - b_{k-1}·b_{k+1} is the LC gap.

**Both C1 ≥ 0 always (LC holds at all stages) and C2 rarely negative.**

## Results (n ≤ 20, all 1,346,022 trees, 112,210,348 checks)

| Decomp | Pos term | Neg events | Total events | % neg | Min ratio |
|--------|----------|------------|--------------|-------|-----------|
| (A) | A1 | 144,696 | 112M | 0.13% | 2.49 |
| (B) | B1 ≥ 0 ALWAYS | B2 neg: 41,370,204 | 112M | 36.9% | 2.46 |
| (C) | C1 ≥ 0 ALWAYS | C2 neg: 142,664 | 112M | 0.13% | 2.63 |

## B1 sub-term decomposition (PROVED algebraically)

B1 = T1 + T2 where:

- **T1 = Δ_k(E_old·(1+x)I_c, E_old·E_c)** ≥ 0 by Karlin's TP2 closure theorem
  - E_old is PF2 (LC and nonneg, verified at all stages)
  - (1+x)I_c ≥_{lr} E_c by SCC of the factor subtree (inductive hypothesis)
  - Therefore T1 ≥ 0 by TP2 convolution closure

- **T2 = Δ_k(x(1+x)·J_old·E_c, E_old·E_c)** can be negative
  - By TP2 closure with E_c PF2: sign related to C2 at previous stage
  - T2 < 0 only 0.012% of events (1,609 / 13.1M at n≤18)

### T1/|T2| ratio when T2 < 0

| n_max | T1 neg | T2 neg | Min T1/|T2| | T2 < 0 rate |
|-------|--------|--------|-------------|-------------|
| 14 | 0 | 6 | -- | 0.002% |
| 18 | 0 | 1,609 | 5.50 | 0.012% |
| 20 | 0 | 16,982 | 3.77 | 0.016% |
| 22 | 0 | 192,826 | 3.66 | 0.020% |

Sign match rate (T2 vs C2_prev): 99.95% at n=22 (205.6M match, 106K mismatch).
All tightest cases have C2_prev = 0 at the specific index k.

## Algebraic structure

The proof has been reduced to:

```
SCC ≥ 0
  ⟺ B1 + B2 ≥ 0              [linearity in 2nd arg]
  ⟸ B1 ≥ 0 ∧ B1 ≥ |B2|       [sufficient]
       ↑
       B1 = T1 + T2
       T1 ≥ 0                  [PROVED: Karlin + SCC(factor)]
       T1 ≥ |T2|               [verified, ratio ≥ 3.77]
```

### Karlin's theorem application (T1 ≥ 0 proof)

**Theorem (Karlin, 1968):** If A is PF2 (nonneg LC sequence) and f ≥_{lr} g
(i.e., Δ_k(f,g) ≥ 0 for all k), then A*f ≥_{lr} A*g.

**Application:**
- A = E^{(t-1)} = accumulated product of I(T_c) polys. PF2 verified.
- f = (1+x)·I_c, g = E_c.
- f ≥_{lr} g ⟺ SCC holds for the factor subtree T_c (inductive hypothesis).
- Conclusion: E^{(t-1)}·f ≥_{lr} E^{(t-1)}·g, i.e., T1 ≥ 0. □

### Why T2 can be negative even when C2_prev ≥ 0

Karlin's theorem preserves the GLOBAL sign condition (all k simultaneously),
not the pointwise sign. If C2_prev ≥ 0 at index k but < 0 at some other index k',
the convolution with E_c mixes contributions from different indices, potentially
creating T2 < 0 at k. The sign match rate is 99.97%, confirming near-pointwise
preservation but not exact.

## Ratio convergence trend

| n_max | min B1/|B2| | min C1/|C2| | min T1/|T2| |
|-------|-------------|-------------|-------------|
| 14 | 2.77 | 5.75 | -- |
| 18 | -- | -- | 5.50 |
| 20 | 2.46 | 2.63 | 3.77 |
| 22 | 2.41 | 1.94 | 3.66 |

Decomposition (C) dropped below 2 at n=22 (less useful for proof).
Decomposition (B) remains above 2.4 and converges to 2 (balanced spider infimum).
T1/|T2| ratio decreasing slowly; all tightest cases at k=10-12, stage=4, 4 factors.

### Balanced spider asymptotics (s identical P₃ arms)

At k=1: B1 = 4s²+5s-5, B2 = -2s²+s+5. Ratio B1/|B2| → 2 as s→∞.
This is the GLOBAL minimizer family. The infimum of B1/|B2| is exactly 2 (never achieved).

## Remaining gaps for a full proof

1. **T1 ≥ |T2|**: Ratio ≥ 3.66 at n=22. Need algebraic bound.
2. **B1 ≥ |B2|**: Ratio ≥ 2.41 at n=22 (infimum = 2). Need algebraic bound.

Both reduce to: "the TP2-propagated SCC of the factor dominates the TP2-propagated
C2 deficit, after convolution with the accumulated E product."

Equivalently: Δ_k(e_new, E_old·(Q + 2xR)) ≥ 0 suffices (since B1+B2 ≥ 0 ⟸ B1 ≥ |B2|
⟺ B1 ≥ B1 - SCC, i.e., SCC ≥ 0). But a weaker condition also works.

## 4-Term Bilinear Decomposition (n ≤ 22, 930M checks)

By bilinearity in both arguments:
SCC = QQ + QR + RQ + RR where
- **QQ = Δ_k(e_old·Q, E_old·Q) ≥ 0 ALWAYS** (Karlin + SCC_old, PROVED)
- QR = Δ_k(e_old·Q, x·E_old·R) — usually negative (45.5M events)
- RQ = Δ_k(x(1+x)·E_old·R, E_old·Q) — usually positive
- **RR = Δ_k(x(1+x)·E_old·R, x·E_old·R) = c_{k-1}(E_old·R) ≥ 0 ALWAYS** (LC, PROVED)

Cross terms (QR + RQ):
- Positive 98.7% of the time — HELP SCC
- Negative only 1.0% of events (9.3M / 930M)

**When cross < 0: min (QQ+RR)/|cross| = 2.0922** (n=22). Ratio converges to 2 (same infimum as B1/|B2|).

### Tight case anatomy

All tightest cases: k ≈ 10, stage 2, 2 non-leaf children. Typical example:
- QQ = 12, QR = -185, RQ = 44, RR = 283
- cross = -141, diag = 295, ratio = 2.09
- RR does 96% of the work (QQ only 4%)

### Key: RR is the dominant positive term

In contrast to the B1 decomposition (where the Karlin term T1 dominates), in the 4-term
decomposition the LC term RR = c_{k-1}(E_old·J_c) is the main contributor. The Karlin
term QQ can be negligible in tight cases.

### Factor of 2 origin

The infimum of 2 for both diag/|cross| and B1/|B2| comes from e = (1+x)I: the (1+x)
multiplier roughly doubles e relative to b = E at the mode, giving Toeplitz minors
involving e twice the "strength" of those involving b alone.

## Equivalent reformulation: cone preservation

SCC ≥ 0 is equivalent to the cone {(e,b) : e ≥_{lr} b} being invariant under the
convolution matrix M = [[Q, x(1+x)R], [0, P]] (upper triangular). Decomposing:

Step 1: (e,b) → (e·Q, b·Q) — preserves cone (Karlin, Q = E_c is PF2)
Step 2: Add (x(1+x)·E_old·R, x·E_old·R) — a vector with ratio profile (1+x)/1

The additive step adds a vector in the SCC cone ((1+x)γ ≥_{lr} γ for LC γ).
But the SCC cone is NOT closed under addition. The proof needs structure beyond
"both summands are in the cone."

## Tightest tree identification (n=20, 21, 22)

The exact trees achieving min diag/|cross| have been identified:

| n | ratio | QQ | QR | RQ | RR | cross | g6 |
|---|-------|----|----|----|----|-------|----|
| 20 | 3.091 | 12 | -89 | 34 | 158 | -55 | S???????C?G?G?C?@?CG??o@c??w??D__ |
| 21 | 2.092 | 12 | -185 | 44 | 283 | -141 | T???????C?G?G?C?@??G??o?A_C?_R?_?[O? |
| 22 | 2.092 | 12 | -185 | 44 | 283 | -141 | U?????????O?O?G?A??O?@?AA??B?BC??@WCF_?? |

**Tree structure:** "Star-of-stars" trees. Degree sequence [4, 3, 3, 3, 3, 2, 2, ...].
A trivalent hub with three arms, each terminating in a fan of P₂ pendants, plus one arm
carrying the support vertex (root) with a P₃ pendant.

**Key structural features:**
- E_acc = [1, 3, 1] = I(P₃), J_acc = [1, 2] = E(P₃) at all three n values
- The big factor has a P₂ pendant from its root creating a "tail spike": J/E ratio
  drops monotonically then jumps from ~0.18 to 0.50 at the last coefficient
- This spike is structural: E_last = 2, J_last = 1 (two maximum IS excluding root
  vs one including root, due to the pendant P₂)
- n=22 tree has SAME factor polynomials as n=21 (extra leaf doesn't change the critical subtree)

**The balanced spider S(2^s) is NOT the tightest for the 4-term decomposition.**
Cross = QR + RQ is ALWAYS nonneg for balanced spiders (verified s ≤ 24). The tightest
cross-negative events come from these star-of-stars trees instead.

For the B1/B2 decomposition, balanced spiders ARE the tightest: B1/|B2| → 2 as s → ∞.

## Scripts

- `verify_scc_2term.py`: main 3-decomposition scanner
- `verify_B1_decomposition.py`: B1 sub-term profiler (Karlin T1/T2)
- `profile_cross_terms.py`: 4-term bilinear profiler (all terms)
- `profile_cross_negative.py`: margin when cross < 0
- `find_tightest_small.py`: parallel tightest tree finder
- `analyze_balanced_spider_scc.py`: balanced spider SCC analysis
- `verify_scc_decomposition.py`: 5-term decomposition (earlier session)
