# Single-index laws and an assembly lemma for the C2 base facts

Date: 2026-08-11. Companion script: `verify_c2_single_index_laws_20260811.py`;
results: `results/c2_single_index_laws_20260811.json`; head certificates:
`scripts/bottleneck_head_certificates_20260811.py`. Produced with substantive
generative-AI assistance (Fable session, replay-verified computations). All
computations exact integer or exact rational arithmetic; the lemma below is
proved, the laws are finite-corpus conjectures.

## Setting

`notes/c2_bounded_pendant_core_2026-08-08.md` reduces the C2 conjecture
(every tree with at most two branch vertices is log-concave) to four base
facts per arm pair `(a, b)`: with `Q_a = prod P_{a_i}`, `R_a = prod
P_{a_i-1}`,

```
B = Q_aQ_b + x(R_aQ_b + Q_aR_b)      (adjacent two-hub tree)
C = Q_aQ_b + xR_aR_b                 (spider on the combined arms)
```

the facts are: B log-concave, C log-concave, `C ~_p B`, `xC ~_p B`.

**Fact 2 is now a citation, not a target.** C is the independence polynomial
of a spider, and Li, Li, Yang, and Zhang prove all spiders strongly
log-concave (arXiv:2501.04245, verified against
`notes/literature/li_2026_unimodality_two_families.txt` line 115 and the
strategy map). The open facts are B's log-concavity and the two
partial-synchronization relations.

## Single-index forms and the three laws

For coefficient sequences `c = coeffs(C)`, `b = coeffs(B)`, define

```
U_k = c_k b_k - c_{k-1} b_{k+1}
V_k = c_k b_k - c_{k+1} b_{k-1}
W_k = 2 c_{k-1} b_k - c_k b_{k-1} - c_{k-2} b_{k+1}
```

(out-of-range coefficients read 0; `W` is the diagonal of the
partial-synchronization form for the pair `(xC, B)`).

Three laws, exceptionless over every instance computed:

- **LAW-V.** `V_k >= 0` for every arm pair and every index.
- **LAW-W.** `W_k >= 0` for every arm pair and every index.
- **LAW-U.** `U_k < 0` happens only at reflected depth `N - k <= 3`
  (`N = deg B`), and only for some arm pairs (72 of 6,708 structured
  cases; the deepest failure seen is depth 3, at `() x (7,15,5)`).

Evidence: 6,708 structured pairs (all small multisets to weight 24 plus
structured large shapes), 1,458 random pairs to total weight 200, and
hill-climbing from the tightest pairs (mutating arm lengths, counts, and
hub assignment, 150 accepted-descent steps per restart). The climbs
plateau at normalized margins near 0.035 without approaching 0. On the
bottleneck family `(1,1) x (2h+1,2h+1)` all three forms are clean for
every `h <= 80`, including `U`.

## Assembly lemma (proved)

**Lemma.** Let `B`, `C` have nonnegative coefficients with `deg C = deg B
= N` and no internal zeros on their supports, and let `C` be log-concave.
Write `alpha_k = c_k/c_{k-1}` on the positive support (log-concavity
makes `alpha` nonincreasing there).

1. If `V_n >= 0` and `U_m >= 0`, the partial-synchronization inequality
   for `(C, B)` holds at the pair `(m, n)` for every `m >= n`.
2. If `V_n >= 0`, `U_m >= 0`, and `m >= n + 1`, the
   partial-synchronization inequality for `(xC, B)` holds at `(m, n)`;
   at `m = n` it is exactly `W_n >= 0`.

*Proof.* For `(C, B)` at `(m, n)`, split the inequality
`c_m b_n + c_n b_m >= c_{m+1} b_{n-1} + c_{n-1} b_{m+1}` into
`c_m b_n >= c_{m+1} b_{n-1}` and `c_n b_m >= c_{n-1} b_{m+1}`. When the
right side of a split piece is zero it is trivial; otherwise all four
entries are positive and the pieces read `alpha_{m+1} <= b_n/b_{n-1}` and
`b_{m+1}/b_m <= alpha_n`. Since `m + 1 >= n + 1` and `alpha` is
nonincreasing, `alpha_{m+1} <= alpha_{n+1} <= b_n/b_{n-1}`, the last step
being `V_n >= 0`. Likewise `b_{m+1}/b_m <= alpha_m <= alpha_n`, the first
step being `U_m >= 0`. For `(xC, B)` at `(m, n)` with `m >= n + 1`, the
inequality is `c_{m-1} b_n + c_{n-1} b_m >= c_m b_{n-1} + c_{n-2}
b_{m+1}`; the same two-piece split needs `alpha_m <= b_n/b_{n-1}` (from
`alpha_m <= alpha_{n+1} <= b_n/b_{n-1}`, using `m >= n + 1` and `V_n`)
and `b_{m+1}/b_m <= alpha_{n-1}` (from `b_{m+1}/b_m <= alpha_m <=
alpha_{n-1}`, using `U_m`). At `m = n` the `(xC, B)` inequality is
`W_n >= 0` by definition. Support-boundary cases degenerate to products
with a zero factor and hold trivially. QED.

**Consequence.** C's log-concavity (cited) plus LAW-V plus LAW-W plus
`U_k >= 0` on the body reduce both open synchronization facts to the
`U`-failure strip: the only pairs `(m, n)` not covered are those whose
larger index `m` lies at reflected depth `<= 3`. The base relations
themselves never fail in the corpus (checked directly on 2,000 structured
pairs in the verifier), so the strip pairs hold for head-window reasons
that a proof must supply per family. This is the C2 instantiation of the
banded-synchronization-plus-head-reserve architecture from
`notes/partial_sync_obstruction_2026-08-11.md`.

## Head certificates for the bottleneck family

For arms `(1,1) x (2h+1,2h+1)` (the asymptotically tight family of the
c2 note), every reflected coefficient of `B` and `C` at fixed depth `d`
is an explicit polynomial in `h`, derived from `[x^k]P_n =
binom(n+1-k, k)` and validated against exact instances for `h <= 30`:
for example `B_{N-1} = 2h^2+5h+6`, `C_{N-1} = h^2+3h+4` (matching the c2
note), with closed forms computed to depth 8. From these,
`scripts/bottleneck_head_certificates_20260811.py` computes exact
polynomial-in-h margins with Sturm-certified real-root isolation
(sympy `real_roots`, exact arithmetic, per the root-finding rule):

- Turan gaps of `B` at reflected depths 1 through 7: positive leading
  coefficients and no real root at or above `h = 2`; depths 1 and 2 have
  no real roots at all. Positivity for every in-range integer `h`
  follows, with small-`h` boundary cases covered by the exact instance
  checks.
- Partial-synchronization margins for both `(C, B)` and `(xC, B)` at all
  head pairs with depths at most 4: every real root of every margin
  polynomial is negative, so every margin is positive for all `h >= 0`.
  The `(dm=1, dn=0)` margin reproduces the c2 note's `h + 2` exactly.

Combined with the lemma, an all-`h` proof of facts 3 and 4 for the
bottleneck family now needs only: `U_k, V_k, W_k >= 0` for all `(h, k)`
(three two-parameter single-index families; clean empirically to
`h = 80`). B's log-concavity additionally needs the body depths beyond 7.

## Bilinear localization

`U`, `V`, `W` are bilinear in `(C, B)`. Writing `C = g1 + g2` and
`B = g1 + h2 + h3` with `g1 = Q_aQ_b`, `g2 = xR_aR_b`, `h2 = xR_aQ_b`,
`h3 = xQ_aR_b`, each form is a sum of six generator-pair terms. Across
126 arm pairs: for `V` the terms `(g1,g1)`, `(g2,g1)`, `(g2,h2)`,
`(g2,h3)` are nonnegative at every index while `(g1,h2)`, `(g1,h3)` go
negative somewhere; for `U` the clean set is `(g1,*)` and the `(g2,*)`
terms are dirty; for `W` only `(g1,g1)` is clean. So the laws hold by
aggregate cancellation, not termwise positivity, matching the two-hub
kill-test's finding that the missing positivity is aggregate. The
decomposition tells a proof exactly which cross terms must be dominated
and suggests reduction of the clean terms to single-arm interlacing
facts via the Hu-Wang-Zhao-Zhao common-convolution closure.

## What this changes

1. The two synchronization base facts are no longer two-index problems:
   off a three-deep head strip they follow from three single-index laws
   and a cited theorem.
2. The remaining universal targets, in order of apparent difficulty:
   LAW-V and LAW-W (single-index, aggregate, adversarially rigid),
   the head strip per family (polynomial certificates as for the
   bottleneck family), and B's log-concavity (independent lane, still
   needing an aggregate Turan-gap argument).
3. The bottleneck family is nearly closed: three clean single-index
   families in `(h, k)` stand between the current state and an all-`h`
   theorem for the synchronization facts.

## Limitations

The three laws are conjectures backed by 8,166 exact instances and
adversarial search; none is proved. The head certificates prove
polynomial positivity for the stated depths and family only. The lemma's
support hypotheses (no internal zeros) hold for independence polynomials
here but are used without a general statement about degenerate supports.
The bilinear cleanliness table is empirical (126 pairs) and the
single-arm reduction of the clean terms is a suggestion, not a proof.
