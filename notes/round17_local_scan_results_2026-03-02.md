# Round 17 Local Scan Results (2026-03-02)

## Overview

Three local scans ran while R17 prompts were dispatched to GPT/Codex.
All completed successfully with exact integer arithmetic.

---

## Scan 1: CB Pairwise Decomposition (test_cb_pairwise.py)

**Scope**: n=3..15, 13,186 trees, 59,061 support-vertex rootings, 24,652 steps (t >= 2),
with boundary-inclusive CB indices (includes index -1).

| Metric | Checks | Failures | Rate | First failure |
|--------|--------|----------|------|---------------|
| X_k >= 0 | 253,819 | 2,015 | 0.79% | n=13, k=1, val=-3 |
| S(i,j,k) >= 0 | 6,996,985 | 295,355 | 4.22% | n=7, i=-1, j=1, k=1 |
| F(i,j) >= 0 | 652,441 | 370,015 | 56.7% | n=6, i=-1, j=1 |
| STP2(P,Q) >= 0 | 4,102,212 | 0 | 0% | (none) |

### Key findings

1. **X_k >= 0 is FALSE**. The cross sum in D_k + X_k = Λ_k^{new} can be negative.
   This kills the "prove X_k >= 0" approach that was the R16/R17 target.

2. **F(i,j) >= 0 is FALSE** at nearly half the cases. The symmetrized bracket
   A(i)B(j) + A(j)B(i) - A(i-1)B(j+1) - A(j-1)B(i+1) does NOT have a definite sign.
   Schur-concavity of the 2x2 permanent for LC sequences is dead.

3. **S(i,j,k) >= 0 is FALSE** at 3% of cases. Pairwise (above + below diagonal)
   symmetric sums don't individually dominate.

4. **STP2(P,Q) = STP2(I_c, E_c) is TRUE universally**. Zero failures across 4.10M checks.
   Whenever Δ_{i,j}(A,B) < 0 (i > j), the complementary factor minor Δ_{k-i,k-j}(P,Q) >= 0.

### Implications for R17 GPT prompts

- **Prompt 2** (prove X_k >= 0): Target is DEAD. GPT will discover this is false.
  The relevant question is now: prove D_k + X_k >= 0 holistically.

- **Prompt 1** (Codex scan): Our scan already covers Tasks 1-3 of this prompt.
  Task 4 (STP2(I,E)) is covered by Scan 2 below.

---

## Scan 2: STP2(I,E) Derivation Check (test_stp2_ie_derivation.py)

**Scope**: n=3..18, 205,002 trees, 3,553,676 rootings (ALL rootings, not just support),
36,709,721 total (m, root) checks.

### Main result: STP2(I,E) holds at ALL rootings of ALL trees n <= 18.

Zero failures. This confirms STP2(I,E) is universal, not just at support vertices.

### Correction term analysis

The STP2(I,E) diagonal condition: c_m(E) + [E(m)J(m-1) - E(m+1)J(m-2)] >= 0

The correction E(m)J(m-1) - E(m+1)J(m-2) is negative in 1,104,396 cases:
- At support vertices: 1,088,362 (vast majority)
- At non-support vertices: 16,034

### Rescue ratios (c_m / |correction| when correction < 0)

| Scope | Min ratio | Tree (n=18) | Root | m |
|-------|-----------|-------------|------|---|
| Overall | 2.625 | Q????A?O@?A?A?@??OEA?KO?N?? | root=7 (non-support) | 9 |
| Support only | 3.842 | same tree | root=0 (support) | 9 |
| Non-support only | 2.625 | same tree | root=7 | 9 |

### Key findings

1. **LC gap always rescues**: min ratio 2.625 means c_m(E) > 2.6 |correction|.
   Substantial headroom at all rootings.

2. **Correction is negative at support vertices too**: the chain argument from R17 Prompt 3
   (E ≽ J + LC(E) + LC(J) => correction >= 0) does NOT work. The ratio chain breaks.

3. **Same extremal tree** for both support and non-support worst cases (n=18 tree).

---

## Scan 3: CB Rescue Profile (test_cb_rescue_profile.py, corrected)

**Scope**: n=3..17, 81,135 trees, 412,437 support-vertex rootings, 185,447 steps (t >= 2),
2,109,316 total X_k checks.

### Critical correction (index boundary)

Initial profiling omitted the CB boundary index `j = -1` in cross-term gap grouping.
After fixing the index window, reconstruction sanity holds exactly:

- `X_k = sum_g gap_sum(g)` in **all 2,109,316 checks** (0 failures).

This invalidates the earlier provisional claim that cumulative gap sums stay nonnegative
for all gaps in `X_k < 0` cases.

### Main quantitative results

- `X_k < 0` in 25,001 / 2,109,316 checks (1.1853%); first appears at n=13.
- `D_k / |X_k|` when `X_k < 0`:
  - min = 2412/527 = 4.576850
  - median = 45/4 = 11.25
  - p95 = 139
- `D_k >= |X_k|` in **100%** of `X_k < 0` cases (25,001 / 25,001).

### Gap-sign profile

- `gap_sum(1) >= 0` in all checks (0 negatives).
- Negative rates by gap:
  - g=2: 27.76%
  - g=3: 32.30%
  - g=4: 20.86%
  - g=5: 9.33%
  - g=6: 2.07%
  - g=7: 0.12%
  - g=8: 0.0009%
  - g >= 9: 0 negatives observed.

### Prefix-cumulative test (X_k < 0 cases only)

- `sum_{g<=1} gap_sum(g) >= 0`: 100%
- `sum_{g<=2} gap_sum(g) >= 0`: 11.48%
- `sum_{g<=3} gap_sum(g) >= 0`: 0.272%
- `sum_{g<=4} gap_sum(g) >= 0`: 0%

So the naive prefix-cumulative rescue does **not** work with correct indexing.

---

## Strategic Assessment

### What's confirmed (all 0 failures)
- Λ_k^{new} >= 0 (ladder-minor preservation): inductively at all child steps
- STP2(I,E) at ALL rootings (not just support vertices)
- STP2(P,Q) = STP2(I_c, E_c) universally
- D_k >= 0 (diagonal from IH)
- D_k + X_k >= 0 (the full CB sum)
- In all sampled `X_k < 0` events (n<=17), `D_k` dominates by factor >= 4.576.

### What's dead
- X_k >= 0 alone
- F(i,j) >= 0 (symmetrized bracket)
- S(i,j,k) >= 0 (pairwise symmetric sum)
- Correction >= 0 in STP2(I,E) derivation
- Naive prefix-cumulative gap rescue (`sum_{g<=G} gap_sum(g) >= 0` for all G)

### The remaining gap

**Prove D_k + X_k >= 0 holistically.** The diagonal D_k (nonneg from IH) must dominate
the occasionally-negative cross sum X_k.

Approaches to try:
1. **Averaged CB per-term bound** for gap-1 terms (j = i+1): proved T(i,j) >= 0.
   For larger gaps: need alternative argument.
2. **Weighted/telescoping gap aggregation** (not plain prefix sums).
3. **Quadratic form**: express D_k + X_k as a quadratic form in coefficient vectors
   and show PSD using both STP2 conditions.
4. **Direct proof of Λ_k^{new} >= 0** using Form 2 (factor-minor grouping):
   Σ_{i,j} A(i)B(j)·Δ_{k-i,k-j}(P,Q). Nonneg for i >= j (STP2(P,Q)).
   Need i < j terms to be dominated by i > j terms.
