# Round 8 Codex 5.3 Triage (2026-03-01)

## Instance 1: Inductive Bound (correction vs main) — COMPLETE (121K tokens)

**Key finding: 6 abstract constraints are INSUFFICIENT, but tree-realizability saves it.**

- Synthetic violations exist where all 6 constraints hold but `main + corr < 0`
  - Example: `E_t = [1,1,1], J_t = [1]` (never occurs as tree DP factor)
- Exhaustive cross-product: 2,973 tree-realized (U,V) x 1,184 tree-realized (Q,R) = **3.52M checks, 0 violations**
- Conclusion: tree-realizability provides the missing structure beyond the 6 listed constraints

**Correction ratio by s-value (n ≤ 20, UPDATED):**

| s | vertices | checks | neg_corr | min D/|T| | Status |
|---|----------|--------|----------|-----------|--------|
| 1 | 4,955,759 | 68.8M | 2.3M | 1.2857 (= 9/7) | PROVED (SCC reduction) |
| 2 | 2,411,269 | 53.0M | 1.5M | 1.8000 (= 9/5) | Key battleground |
| 3 | 578,941 | 17.1M | 379K | 3.0000 | Comfortable |
| 4 | 90,160 | 3.3M | 52K | 4.1075 | Comfortable |

Note: min ratios DECREASED from n≤18 to n≤20. The s=1 pendant-star extremal gives ratio = (n-2)/(n-6) → 1. **Ratio-based bounds are dead.** Must use absolute margin (total = 2m → ∞).

**Factor ordering for s=2 (n ≤ 18):**

- ∃ factor order with correction ≥ 0 at ALL k: **99.1%** (322,165/325,087)
- BOTH orders nonneg: **85.8%** (278,772/325,087)
- 0.9% of s=2 vertices have negative correction in BOTH orderings — but still main+corr ≥ 0
**Best-order analysis for s=2 (n ≤ 18):**

- Worst "best-order" min ratio = **3.905** (n=18). Both orderings give same ratio.
- Choosing the optimal factor ordering improves the margin dramatically vs arbitrary ordering (3.9 vs 2.0).
- But this only helps within a fixed tree; the absolute margin for E ≽ J is what matters.

**Synthetic counterexample (all 6 constraints satisfied but E ≽ J fails):**

- E^(t-1) = J^(t-1) = [1,1], E_t = [1,1,1], J_t = [1]
- At k=2: main = 0, corr = -1, total = -1. FAILS.
- But [1,1,1] is NOT a tree-realizable E_c (no tree has dp[c][0] = [1,1,1])
- Conclusion: **tree-realizable cone lemma** is the missing proof ingredient

**Instance 1 final line count:** 1701 lines, 121K tokens

## Instance 2: Curvature Rescue (T3 ≥ |T1|)

**Key finding: circularity is REAL, modified target works.**

- `c_k ≥ |d_k| ⟺ d_{k-1} ≥ 0` is **FALSE** (15,028 counterexamples at n≤16)
- `d_k < 0` at k=1 CAN happen when (1+x)^ℓ included (n=5 broom)
- Both d_k and d_{k-1} negative simultaneously: **26.7%** of all 907M checks
- Forward induction on d_k ≥ 0 is impossible

**New algebraic target:** `g_k = c_k + d_k = b_k·a_k - b_{k+1}·a_{k-1} ≥ 0`
- Exact identity: `T3 - |T1| = b_{k-1}·g_k + j_{k-2}·c_k`
- g_k ≥ 0 + amplification gives T3 ≥ |T1| directly
- g_k ≥ 0 means: E_k·I_k ≥ E_{k+1}·I_{k-1} (ratio dominance of I vs E at adjacent indices)

## Instance 3: Combinatorial/Structural Approaches (COMPLETE, 129K tokens)

**Key finding: one-step decomposition is the right lens, not factorwise dominance.**

- Y_k = (E·R)_k / (J·Q)_{k+1} nondecreasing: FAILS ~8% of stages (not usable)
- s=2, stage-2 min ratio = 3.333 at n=18 (171K negative events)
- s=1 is the global minimizer at ratio 1.333 (already PROVED)
- Literal combinatorial injection φ: A_k × B_{k+1} → A_{k+1} × B_k too rigid; suggests fractional matching/Strassen coupling

**One-step algebraic decomposition:**

Write A = E_old·E_t, B = E_old·J_t, C = J_old·E_t, so E_new = A + xB, J_new = C. Then:

```
Δ_k(E_new, J_new) = Δ_k(A,C) + (B_k·C_k - B_{k-1}·C_{k+1})
                   = D_k (≥0, Karlin) + T_k (can be neg)
```

**Ratio criterion:** Define α_k = A_{k+1}/A_k, β_k = B_k/B_{k-1}, γ_k = C_{k+1}/C_k. Then:

```
γ_k - β_k ≤ (A_k / B_{k-1}) · (α_k - γ_k)
```

- Karlin gives α_k ≥ γ_k (RHS positive)
- SCC at c_t should bound β_k-vs-γ_k deficit
- This is the exact necessary and sufficient condition for each inductive step

**s=2 lifting route:** Use SCC at c_2 to control the ratio criterion. Empirically very favorable (min main/|corr| = 2.0 overall, 3.33 at genuine stage-2).

## Cross-Instance Synthesis

1. **s=1 PROVED** (SCC reduction, 63% of support vertices)
2. **s≥2 has 2x+ better margins** (min ratio ≥ 2.0 vs 1.33)
3. **Generic product closure FALSE** (synthetic counterexamples)
4. **Tree-realizability essential** (3.52M cross-product, 0 violations)
5. **Best algebraic targets**:
   - g_k = c_k + d_k ≥ 0 (Instance 2: ratio dominance of I vs E)
   - Ratio criterion γ_k - β_k ≤ (A_k/B_{k-1})(α_k - γ_k) (Instance 3: one-step induction)
6. **Three convergent views**: all instances agree the proof must use tree-specific structure at s≥2, with Karlin providing the main term and SCC bounding the correction

## Independent Verification (Claude, 2026-03-01)

**One-step decomposition identity: CONFIRMED.** 2,021,980 checks (n=3..16), 0 identity failures.

**Karlin main term D_k ≥ 0: CONFIRMED.** 0 failures across all checks.

**Correction term profile (n=3..18):**

| n | T<0 count | T total | T<0 % | s=1 ratio | s=2 ratio |
|---|-----------|---------|-------|-----------|-----------|
| 10 | 83 | 2,877 | 2.9% | 2.000 | — |
| 14 | 5,076 | 165,103 | 3.1% | 1.500 | — |
| 16 | 41,013 | 1,302,831 | 3.1% | 1.400 | 3.556 |
| 18 | 335,299 | 10,453,779 | 3.2% | 1.333 | 3.333 |

**s=1 ratio formula CONFIRMED:** min D/|T| = **(n−2)/(n−6) → 1** as n→∞.
This matches Instance 1's pendant-star formula exactly. At s=1 the one-step decomposition
IS the full E≽J gap (only one non-leaf child), so the formulas coincide.

**s≥2 ratios much more comfortable:** s=2 ≈ 3.3, s=3 ≈ 4.2, s=4 ≈ 4.9 at n=18.
These are BETTER than the product-level ratios from Instance 1 because issues haven't compounded.

## Remaining Gap

Find tree-specific structural property implying `main + corr ≥ 0` at s≥2.
Margins large (factor 2+), but the right algebraic lemma hasn't been identified.
Most promising: show SCC implies the ratio criterion γ_k - β_k ≤ (A_k/B_{k-1})(α_k - γ_k).

## Files

- `results/round8_codex_instance1.md` (1701 lines, 121K tokens, complete)
- `results/round8_codex_instance2.md` (1426 lines, complete)
- `results/round8_codex_instance3.md` (1376 lines, 129K tokens, complete)
- `verify_onestep_decomposition.py` (full verification n≤20, 4 instances still running)
- `verify_onestep_quick.py` (quick verification n=17-18, COMPLETE)
