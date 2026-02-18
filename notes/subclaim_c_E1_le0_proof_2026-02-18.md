# Sub-claim C: Proof of E₁ ≤ 0 for k ≡ 0,2 (mod 3), k ≥ 5 (2026-02-18)

This note closes the algebraic gap "E₁ ≤ 0 for k ≥ 5, k ≡ 0,2 (mod 3)" in the
Sub-claim C proof for mixed spiders S(2^k, 1^j).

Context: `notes/subclaim_c_algebra_attempt_2026-02-18.md`,
`notes/unit_leaf_c2_algebra_2026-02-18.md`.

---

## Setup and goal

For the j=1 branch of S(2^k, 1^j) at its tie-fugacity λ₁ = i_{m₁-1}/i_{m₁} and
mode m₁, define:

    A₁ = 2kλ₁/(1+2λ₁) + λ₁/(1+λ₁) - (m₁-1)
    B₁ = 1 + kλ₁/(1+λ₁) - (m₁-1)
    r₁ = λ₁(1+λ₁)^(k-1) / (1+2λ₁)^k
    E₁ = r₁(B₁-A₁)/(1+r₁)

Since r₁ > 0 and 1+r₁ > 0, we have E₁ ≤ 0 iff B₁ - A₁ ≤ 0.

**Goal:** Prove B₁ ≤ A₁ for all k ≡ 0,2 (mod 3), k ≥ 5.

---

## Step 1: Simplify B₁ - A₁

    B₁ - A₁ = [1 + kλ/(1+λ) - (m₁-1)] - [2kλ/(1+2λ) + λ/(1+λ) - (m₁-1)]
             = 1 + kλ/(1+λ) - 2kλ/(1+2λ) - λ/(1+λ)
             = 1 + (k-1)λ/(1+λ) - 2kλ/(1+2λ)

(dropping subscripts on λ₁ for clarity).

Combine over the common denominator (1+λ)(1+2λ):

    = [(1+λ)(1+2λ) + (k-1)λ(1+2λ) - 2kλ(1+λ)] / [(1+λ)(1+2λ)]
    = [1 + 3λ + 2λ² + (k-1)λ + 2(k-1)λ² - 2kλ - 2kλ²] / [(1+λ)(1+2λ)]
    = [1 + 3λ + 2λ² + kλ - λ + 2kλ² - 2λ² - 2kλ - 2kλ²] / [(1+λ)(1+2λ)]
    = [1 + 2λ - kλ] / [(1+λ)(1+2λ)]
    = [1 + λ(2-k)] / [(1+λ)(1+2λ)]

So:

    **B₁ - A₁ = [1 - (k-2)λ] / [(1+λ)(1+2λ)]**

Since (1+λ)(1+2λ) > 0, we have B₁ ≤ A₁ iff (k-2)λ ≥ 1, i.e.:

    **E₁ ≤ 0  ⟺  λ ≥ 1/(k-2).**

(For k ≥ 3, so k-2 ≥ 1 and the bound is positive.)

---

## Step 2: The b-bound on λ₁ from the unit-leaf c₂ proof

From `notes/unit_leaf_c2_algebra_2026-02-18.md`, for the j=1 branch:

    I_B(x) = (1+2x)^k,  b_t = C(k,t) · 2^t.

The b-bound states:

    λ₁ ≥ τ := b_{m₁-2} / b_{m₁-1} = [C(k, m₁-2)·2^{m₁-2}] / [C(k, m₁-1)·2^{m₁-1}]
             = C(k, m₁-2) / [2 · C(k, m₁-1)]
             = (m₁-1) / [2(k - m₁ + 2)].

This bound is proved in `unit_leaf_c2_algebra_2026-02-18.md` (log-concavity of b_t
plus e_r/e_{r-1} bound).

Mode formula (from MEMORY.md and mode_j1 in scripts): m₁ = ⌊(2k+3)/3⌋.

Computing τ by residue class:

**k ≡ 0 (mod 3), k = 3t:**
    m₁ = (6t+3)//3 = 2t+1, so m₁-1 = 2t, k-m₁+2 = 3t-(2t+1)+2 = t+1.
    τ = 2t / [2(t+1)] = t/(t+1).

**k ≡ 2 (mod 3), k = 3t+2:**
    m₁ = (6t+7)//3 = 2t+2, so m₁-1 = 2t+1, k-m₁+2 = (3t+2)-(2t+2)+2 = t+2.
    τ = (2t+1) / [2(t+2)].

(Note: the task description had τ = t/(t+2) for k ≡ 0 mod 3, which is incorrect.
The correct value is t/(t+1), as confirmed by the script via exact coefficient
computation.)

(Note: k ≡ 1 (mod 3) is excluded from Sub-claim C.)

---

## Step 3: Verify τ ≥ 1/(k-2) for k ≥ 5

We need τ(k-2) ≥ 1.

**Case k ≡ 0 (mod 3), k = 3t, t ≥ 1:**

    τ(k-2) = [t/(t+1)] · (3t-2) = t(3t-2)/(t+1).

Need t(3t-2)/(t+1) ≥ 1, i.e., t(3t-2) ≥ t+1, i.e., 3t²-3t-1 ≥ 0.

The discriminant is 9+12=21, positive root at t = (3+√21)/6 ≈ 1.26.
For integer t ≥ 2 (k ≥ 6): 3(4)-6-1=5 > 0. ✓
For t=1 (k=3): 3-3-1=-1 < 0. This is outside the theorem scope (k ≥ 5).

So for k ≡ 0 (mod 3), k ≥ 6 (t ≥ 2): E₁ ≤ 0 follows from the b-bound.

**Case k ≡ 2 (mod 3), k = 3t+2, t ≥ 1:**

    τ(k-2) = [(2t+1)/(2(t+2))] · (3t) = 3t(2t+1) / [2(t+2)].

Need 3t(2t+1)/[2(t+2)] ≥ 1, i.e., 3t(2t+1) ≥ 2(t+2), i.e., 6t²+3t ≥ 2t+4,
i.e., 6t²+t-4 ≥ 0.

The discriminant is 1+96=97, positive root at t = (-1+√97)/12 ≈ 0.737.
For integer t ≥ 1 (k ≥ 5): 6+1-4=3 > 0. ✓

So for k ≡ 2 (mod 3), k ≥ 5 (t ≥ 1): E₁ ≤ 0 follows from the b-bound. ✓

---

## Step 4: Base cases not covered by Step 3

The only remaining cases are:

- k=3 (k ≡ 0 mod 3, t=1): outside scope of the theorem (k ≥ 5 stated), but
  checked for completeness.
- k=5 (k ≡ 2 mod 3, t=1): covered by Step 3 above (6t²+t-4=3>0 at t=1). ✓

So the b-bound route actually covers ALL k ≥ 5 for both residue classes.

For k=3, τ(k-2) = [1/2]·1 = 1/2 < 1, so the b-bound does not give E₁ ≤ 0.
But the direct check (see verify_subclaim_c_E1_le0.py) confirms B₁-A₁ at k=3
is positive (E₁ > 0 there), which is consistent with Sub-claim C being proved
directly at k=3 by other means (the raw margin difference is still negative,
since A_gap dominates).

---

## Step 5: Summary of the proof of E₁ ≤ 0

**Theorem.** For all k ≡ 0,2 (mod 3) with k ≥ 5, and for the j=1 branch of
S(2^k, 1^j) at its tie-fugacity λ₁:

    B₁ - A₁ ≤ 0,  equivalently  E₁ ≤ 0.

**Proof.**

1. By direct algebra (Step 1): B₁ - A₁ = [1-(k-2)λ₁]/[(1+λ₁)(1+2λ₁)]. Since
   the denominator is positive, B₁ ≤ A₁ iff λ₁ ≥ 1/(k-2).

2. By the b-bound from unit_leaf_c2_algebra_2026-02-18.md: λ₁ ≥ τ =
   (m₁-1)/[2(k-m₁+2)].

3. It remains to show τ ≥ 1/(k-2), i.e., (k-2)τ ≥ 1.

   - For k ≡ 2 (mod 3), k=3t+2, t ≥ 1: (k-2)τ = 3t(2t+1)/[2(t+2)] ≥ 1 iff
     6t²+t-4 ≥ 0, true for all t ≥ 1 (since 6(1)²+1-4=3>0 and the quadratic
     is increasing for t ≥ 1).

   - For k ≡ 0 (mod 3), k=3t, t ≥ 2: (k-2)τ = t(3t-2)/(t+1) ≥ 1 iff
     3t²-3t-1 ≥ 0, true for all t ≥ 2 (since 3(4)-6-1=5>0 and the expression
     is increasing for t ≥ 1).

4. Base case k=3t, t=1 (k=3): not in scope (k ≥ 5 required). Verified directly
   (E₁ > 0 at k=3, but A_gap compensates and Sub-claim C still holds).
   Base case k=3t+2, t=1 (k=5): covered by case analysis (t ≥ 1 above).

Therefore E₁ ≤ 0 for all k ≥ 5 with k ≡ 0,2 (mod 3). □

---

## Step 6: Combined sufficient condition for Sub-claim C

With E₁ ≤ 0 established, Sub-claim C reduces to verifying:

    A_gap = A₀ - A₁ > -E₀.

This is verified exactly for k ≤ 500 in prove_subclaim_c_algebra.py (A_gap ≥ 1/(4k)
and A_gap > -E₀, both 0 failures). The gap A_gap ≥ 1/(4k) is a separate algebraic
bound not closed here; it is verified computationally.

---

## Verification

Script: `verify_subclaim_c_E1_le0.py`

Checks (all pass, verified 2026-02-18):
1. Polynomial inequalities (F3a) 6t²+t-4 ≥ 0 (t ≥ 1, k ≡ 2 mod 3) and
   (F3b) 3t²-3t-1 ≥ 0 (t ≥ 2, k ≡ 0 mod 3): min values 3 and 5. PASS.
2. b-bound λ₁ ≥ τ using exact Fraction coefficients for k ≤ 200: PASS.
   Minimum gap λ₁-τ ≈ 0.0110 at k=200.
3. (k-2)τ ≥ 1 for k ≡ 0,2 (mod 3), k ≥ 5, k ≤ 500: PASS.
   Minimum gap (k-2)τ-1 = 1/2 at k=5.
4. E₁ ≤ 0 directly via IS polynomial coefficients (exact Fraction), k ≤ 200: PASS.
   Algebraic identity B₁-A₁ = [1-(k-2)λ]/[(1+λ)(1+2λ)]: PASS (0 failures).
   Maximum E₁ (closest to 0): ≈ -1.7×10⁻⁸ at k=200.
