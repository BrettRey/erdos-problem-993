# Primus rerun audit, 26 September 2026
<!-- SUMMARY: Audit of the second Primus run on the depth-three window brief, after the platform's announced changes: no advance on the open case; one correct but minor exclusion (spiders, via a published theorem); all checkable numbers reproduce; one statement contradicts the supplied data; no code returned · status: audited · updated: 2026-09-26 -->

**Input.** Same brief as the first run (`primus_handoff_2026-09-23.md`,
pasted as text), rerun by Brett after Primus said it had changed the
platform. **Output.** A 9-page report, "Near-top coefficients of tree
independence sequences in a bounded order window", in Brett's `pdf-inbox`
(`paper-2.pdf`, SHA-256 `e4120f2e…`; `paper-2.md`, SHA-256 `7d0df7df…`).
No code bundle: the report says "Artifacts are available from the authors on
request". Brett reports the run took hours.

## Scored against the criteria fixed before the rerun

1. **Does any new claim survive the three supplied test cases?** Mostly yes.
   There's no new conjecture this time. But Section 4 says the paper "makes
   no use of" the supplied `b_1 = 0` theorem because "the hypothesis b1 = 0
   does not occur among the trees we examine in the window". Both halves are
   false: the proof of its own main theorem (Theorem 21) uses that theorem for
   the `b_1 = 0` case, and the third supplied window tree, which the report
   quotes in Section 2 (`(b_2, b_3, b_4) = (1, 16, 4216)`), has `b_1 = 0`
   (packet `replay.py`: `b_0..b_4 = [0, 0, 1, 16, 4216]`). The error does no
   mathematical damage, but it is the same kind as the first run's: a claim
   contradicted by input the report itself quotes.
2. **Does its code recompute the counts or read stored numbers?** Not
   assessable. No code was returned.
3. **If there's no advance, does it say so?** Much better than the first run.
   It states that the residual case is open, that nothing bears on #993
   beyond the bounded inequality, and exactly what bound is missing (an upper
   bound on `b_4` in terms of the extendable counts). The abstract's first
   sentence ("We prove a log-concavity inequality … for trees of order
   between 33 and 38 …") overstates before the qualification arrives, and the
   unproved target is labelled "Theorem 1".

## Mathematical content

- **New and correct, but minor.** Spiders are removed from the open case by
  citing Li, Li, Yang and Zhang (arXiv:2501.04245), whose abstract says "all
  spiders have log-concave independence polynomials". Log-concavity at every
  index includes the target, so this is a valid exclusion. It removes a small
  class already settled in print.
- **Correct and elementary.** Lemma 6: an independent set of size `alpha - d`
  whose closed neighbourhood leaves at least `2d` vertices uncovered extends
  to a maximum independent set. It's the standard halving argument; brute
  force over all trees with `n <= 12` found 0 violations in 9,159 cases.
- **Unchanged.** The open case (`b_1 > 0`, high density, `D < 0`) is exactly
  where it was, minus spiders. The stated obstruction is the one the brief
  already gave.

## Checks run (`scripts/primus_rerun_checks_20260926.py`)

- Path margins for `P_33`–`P_38` (Example 15): all six match.
- The three order-8 trees said to satisfy the open-case conditions
  (Proposition 24): exactly three exist among the 23 trees of order 8, with
  the stated profiles and `D` values.
- Spider deficiency formula (Proposition 10): no violations over all spiders
  with 3–5 legs of lengths 1–7.
- Citations checked: Li–Li–Yang–Zhang's spider theorem (arXiv abstract);
  Kadrawi–Levit's Conjecture 5.1 on log-concavity breaking at `alpha - k`
  for arbitrary `k` (arXiv:2305.01784); Fang et al.'s `17 alpha/25` interval
  and p. 17 remark.

## Minor defects

- Proposition 11's proof miscounts: legs `(1, 1, 2)` have two odd legs, not
  one, and its "orders 5, 7, 9, … in particular 34, 36, 38" mixes parities.
  The conclusion still holds (for example, legs `(1, 2, 2m)` give `delta = 0`
  at every even order).
- Definition 20 restricts the open case to trees in the window, yet
  Proposition 24 calls it "non-empty" on the strength of trees with
  `alpha = 4`. Remark 25 concedes the point.
- Example 15 calls the path a two-legged spider, against the report's own
  definition (exactly one vertex of degree at least 3).
- Remark 22 contains a sentence that doesn't parse ("Since a pair of theorems
  could not both hold if …").
- Remark 23 reports a search of about 1,000 trees per order at `n = 34–36`;
  nothing to check it against without artifacts.

## Verdict

**No advance on the open case.** It's a clear improvement over the first run
in honesty and computational accuracy. One statement still contradicts the
supplied data, and without returned code the "recompute or read" question
can't be answered.
