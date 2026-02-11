# Steedman-style attempt: rooted-pair invariants

Goal: find a rooted-tree invariant on the pair (P, Q) that is preserved
under the root-composition rules and forces unimodality of I = P + Q.

Notation for a rooted tree T (root r):
  P_T(x) = polynomial for independent sets in T with r excluded.
  Q_T(x) = polynomial for independent sets in T with r included.

These satisfy the standard DP:
  P_T = ∏_i (P_i + Q_i)
  Q_T = x ∏_i P_i
where i ranges over children of r, and (P_i, Q_i) are the child rooted pairs.

Define the two derived sequences:
  B_T(x) = Q_T(x) / x = ∏_i P_i   (drop the leading 0 coefficient)
  F_T(x) = P_T(x) - B_T(x)
  R_T(x) = Q_T(x) + B_T(x)
so that
  I_T(x) = P_T(x) + Q_T(x) = F_T(x) + R_T(x).

Combinatorial meanings:
  - B_T counts independent sets with r included (with the forced x removed),
    i.e., independent sets of T - N[r].
  - F_T counts independent sets with r excluded and at least one child root
    included (P_T minus the “all children excluded” case).
  - R_T is the adjacent-sum transform of B_T (see Lemma 1 below).

## Empirical invariants (checked by ad hoc scripts)

Tested on:
  - all non-isomorphic trees for n <= 10 (all roots),
  - 200 random trees for each n in {20,30,40,50} (all roots),
  - the two n = 26 log-concavity failures from results/analysis_n26.json
    (all roots).

Observed without a counterexample in these tests:
  E1. Coefficientwise dominance: P_T >= B_T, so F_T has nonnegative coeffs.
      (This is actually exact by algebra: ∏(P_i + Q_i) >= ∏ P_i.)
  E2. B_T is log-concave and unimodal for every tested root.
  E3. F_T is log-concave and unimodal for every tested root.
  E4. R_T is log-concave and unimodal for every tested root.

Note: E2–E4 are empirical only. E1 is a direct algebraic inequality.

## Lemma 1 (adjacent-sum preserves log-concavity)

Let b_0, b_1, ..., b_m be a nonnegative log-concave sequence with no
internal zeros. Define r_k = b_k + b_{k+1} for k = 0..m-1. Then r_k is
log-concave.

Proof sketch:
Let r_k = b_k + b_{k+1}. We need r_k^2 >= r_{k-1} r_{k+1}. Define
  a = b_{k-1}, b = b_k, c = b_{k+1}, d = b_{k+2}.
Log-concavity gives b^2 >= a c and c^2 >= b d, so ratios
  b/a >= c/b >= d/c.

Write
  r_{k-1} = a + b, r_k = b + c, r_{k+1} = c + d.
After algebra, the inequality r_k^2 >= r_{k-1} r_{k+1} reduces to
  (b+c)^2 >= (a+b)(c+d).
Using the ratio monotonicity above (b/a >= c/b >= d/c), one can show
  (b+c)^2 - (a+b)(c+d) >= 0 (e.g., reduce to the case b/a = c/b, d/c = c/b).

Consequently, if B_T is log-concave, then R_T = B_T + shift(B_T) is
log-concave.

## Why this might help

If E2 and E3 can be proved (B_T and F_T log-concave for all roots),
then R_T is log-concave by Lemma 1. The remaining step is to show that
the sum F_T + R_T is unimodal. This is not automatic: sums of log-concave
or even unimodal sequences can fail to be unimodal.

However, the special structure of (F_T, R_T) may allow a tailored argument:
  - F_T has zero constant term; R_T has positive constant term (R_0 = 1).
  - F_T and R_T arise from the same rooted tree decomposition and are not
    arbitrary log-concave sequences.
  - Empirically, I_T = F_T + R_T is unimodal for all tested cases.

Two possible proof targets:
  T1) Prove a stronger inequality linking F_T and R_T, e.g., a one-sided
      control on the first differences that forbids a “valley” in F_T + R_T.
  T2) Find a root-selection rule (e.g., choose r on a longest path or a
      specific centroid) so that the resulting pair satisfies an additional
      ordering (or difference) inequality, and show unimodality for that
      root. Since I_T is independent of the root, existence of *some* root
      with the property would be sufficient.

## Next concrete steps

1) Try to prove E2: B_T is log-concave for every rooted tree.
   - B_T = I(T - N[r]; x), so this would amount to log-concavity for every
     forest of the form T - N[r]. A weaker inductive claim may suffice.

2) Try to prove E3 by induction:
   - F_T = ∏(P_i + Q_i) - ∏ P_i
         = sum_{∅≠S} (∏_{i in S} Q_i)(∏_{j not in S} P_j).
   - If each term were log-concave and the sum were “compatible” (e.g.,
     supports nested or ratio-ordered), log-concavity would follow.

3) Empirically, for many trees there exists a root with mode(F_T) <= mode(R_T).
   If a general root-selection rule can guarantee this together with a
   monotonicity condition, unimodality of F_T + R_T may follow.

This note is exploratory; no claims here are used in the paper.
