# GPT 5.2 Pro × Codex 5.3 × Claude: Id1 Synthesis (2026-03-01)

## Overview

GPT 5.2 Pro contributed three algebraic identities that recast the E ≽ J correction term.
Claude independently verified all three computationally. Combined with Codex 5.3's Round 8
findings, this gives the cleanest formulation of the remaining gap.

## Identity 1 (Id1): Correction = Memory + Curvature

For any sequences U, V with V_k > 0:

```
U_k V_k - U_{k-1} V_{k+1} = (V_{k+1}/V_k) · d_{k-1}(U,V) + (U_k/V_k) · c_k(V)
```

where d_{k-1}(U,V) = U_k V_{k-1} - U_{k-1} V_k (LR minor) and c_k(V) = V_k² - V_{k-1}V_{k+1} (LC gap).

**Proof:** One-line cancellation (verified as algebraic identity).

**Applied to** U = B = E_old · J_t, V = J^(t) = J_old · E_t:
- c_k(J^(t)) ≥ 0 always (J^(t) is PF2)
- So d_{k-1}(B, J^(t)) is the **sole source of negativity** in the correction

**Verified:** 17,239,045 checks (n ≤ 18), 0 identity failures.

## Integer Form (IF): Three-Term Stage Inequality

Multiply through by J_k:

```
J_k · Δ_k(A,J)  +  J_{k+1} · d_{k-1}(B,J)  +  B_k · c_k(J)  ≥  0
  Term 1 (≥0)        Term 2 (can be neg)        Term 3 (≥0)
```

- Term 1: Karlin main part (Δ_k(A,J) ≥ 0 by TP2)
- Term 2: Mixed LR minor (negative in 22.75% of checks)
- Term 3: Curvature bonus from LC of J (always ≥ 0, never zero when Term 2 < 0)

**Verified:** 0 failures through n = 20 (all trees, all support-vertex rootings).

## Weaker Form (W): Two-Term Sufficiency

Even without the curvature bonus, the Karlin term alone suffices:

```
J_k · Δ_k(A,J) + J_{k+1} · d_{k-1}(B,J) ≥ 0     (W)
```

**Verified:** 0 failures through n = 20.

### W Margin Profile

| n | min W ratio | s=1 extremal |
|---|-------------|--------------|
| 10 | 1.556 | ℓ=7 |
| 14 | 1.320 | ℓ=11 |
| 18 | 1.225 | ℓ=15 |
| 20 | 1.195 | ℓ=17 |

**s=1 formula:** min W ratio = ℓ(ℓ+1)/(ℓ-1)² → 1 as ℓ→∞

**Absolute margins:**
- W margin (s=1 extremal) = ℓ(3ℓ-1)/2 → ∞
- IF margin (s=1 extremal) = 2ℓ² → ∞

**By s-value (n ≤ 20):**
- s=1: min ratio 1.195 (pendant-star extremal)
- s=2: min ratio 1.664
- s=3: min ratio 1.783
- s=4: min ratio 1.742

## Leaf-Augmentation LR Dominance

GPT 5.2 Pro's conceptual insight: g_k = c_k(E) + d_k(I,E) is the LR minor of (I+xE, E).

**Tree interpretation:** I+xE = IS poly of tree with one extra leaf at root.
So g_k ≥ 0 means: "leaf-augmented tree ratio-dominates root-deleted forest."

**Verified:** 0 failures through n = 22 (ALL rootings, not just support).

## Binomial Smoothing Lemma (PROVED)

**Lemma:** (1+x)B ≽ B iff B is LC.

**Proof:** d_k((1+x)B, B) = (B_{k+1}+B_k)B_k - (B_k+B_{k-1})B_{k+1} = B_k² - B_{k-1}B_{k+1}.

**Consequence:** At a support vertex with ℓ leaves and non-leaf product (A, B):
E = (1+x)^ℓ A, J = B. If SCC holds ((1+x)A ≽ B) and B is LC, then:
E = (1+x)^{ℓ-1}·(1+x)A ≽ (1+x)^{ℓ-1}·B ≽ B = J.

So **E ≽ J reduces to SCC product closure for the non-leaf factors.**

**J is LC at all rootings:** Verified 0 failures through n ≥ 22.
Also PROVED: J is a product of PF2 factors (each E_c is PF2), and PF2 is closed under convolution.

## Stronger Form: I+xE ≽ (1+x)E

**FAILS** (many failures starting at n=3). Too strong.

## Comparison of Decompositions

| Decomposition | Positive terms | Negative term | Min ratio (n=20) |
|---------------|---------------|---------------|------------------|
| One-step (Instance 3) | D_k (Karlin) | T_k | 1.286 (s=1) |
| SCC 2-term | c_k(E) | LR_k | 1.94 |
| SCC B1+B2 | B1 (Karlin) | B2 | 2.41 |
| IF (GPT 5.2) | J_k·Δ_k + B_k·c_k(J) | J_{k+1}·d_{k-1}(B,J) | — |
| W (GPT 5.2) | J_k·Δ_k(A,J) | J_{k+1}·d_{k-1}(B,J) | 1.195 |

All decompositions have ratio → 1 (pendant-star extremal) and absolute margin → ∞.

## Proof Architecture (Updated)

```
PROVED: P3 (tail domination) — leaf-swap injection
PROVED: J is PF2 (product of PF2 factors)
PROVED: Binomial smoothing ((1+x)B ≽ B for LC B)
PROVED: Karlin main term (Δ_k(A,J) ≥ 0 by TP2)
PROVED: s=1 case (SCC + transitivity, 63% of support vertices)

VERIFIED (0 fails):
  - Leaf-augmentation: I+xE ≽ E at ALL rootings (n ≤ 22)
  - IF/W form: 3-term and 2-term integer forms (n ≤ 20)
  - E ≽ J at ALL support vertices (907M+ checks, n ≤ 22)
  - SCC at ALL support vertices (930M+ checks, n ≤ 22)

OPEN: SCC product closure for s ≥ 2 non-leaf factors
  Equivalently: (1+x)·∏I_t ≽ ∏E_t (the single remaining gap)
  All decompositions agree: tree-realizability is essential, ratio → 1 but margins → ∞
```

## Next Targets

1. **Cauchy-Binet for W:** Write J_k·Δ_k(A,J) and J_{k+1}·d_{k-1}(B,J) as sums over
   pairs (i<j) with weights (f_j g_i - f_i g_j) ≥ 0 times kernel minors.
   Show kernel minors satisfy pointwise inequality from factor SCC.

2. **Leaf-augmentation algebraically:** Prove I+xE ≽ E for all tree rootings.
   This is g_k = c_k(E) + d_k(I,E) ≥ 0 — the "leaf-augmented tree dominates
   root-deleted forest" statement. Equivalent to Instance 2's target.

3. **Tree-realizable cone lemma (Instance 1):** Characterize the cone of
   tree-realizable (E_c, J_c) pairs and show the correction bound holds on it.

## Files

- `verify_id1_decomposition.py` — Id1 identity and IF/W verification (COMPLETE)
- `profile_id1_margin.py` — W margin profiling by s and n (COMPLETE)
- `verify_leaf_augmentation.py` — Leaf-augmentation at all rootings (RUNNING)
- `verify_onestep_decomposition.py` — One-step identity (Instance 3, COMPLETE)
- `find_s1_minimizer.py` — Pendant-star extremal analysis (COMPLETE)
