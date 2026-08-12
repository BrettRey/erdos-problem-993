# TP2/likelihood-ratio toolkit for the LAW-V and LAW-W targets

Date: 2026-08-12. Produced with generative-AI assistance; every lemma
below is either quoted from a source on disk or proved inline and
verified computationally in exact arithmetic (verification commands at
the end). One item is from memory and marked as such.

## Sources now on disk

Hu, Wang, Zhao, Zhao, *Convolution preserves partial synchronicity of
log-concave sequences*, Math. Inequal. Appl. 20 (2017) 91-103, fetched
as `notes/literature/hu_wang_zhao_zhao_2017_mia_partial_sync.{pdf,txt}`
(ele-math open access). The campaign previously used its statements "as
published" without the text; exact statements now quotable:

- **Definition 2.1 (synchronized).** `a_{k-1}b_{k+1} <= a_k b_k` and
  `a_{k+1}b_{k-1} <= a_k b_k` for all k.
- **Definition 2.2 (weakly synchronized).**
  `a_{k-1}b_{k+1} + a_{k+1}b_{k-1} <= 2 a_k b_k` for all k.
- **Definition 2.4 (partially synchronized).**
  `a_m b_n + a_n b_m >= a_{m+1}b_{n-1} + a_{n-1}b_{m+1}` for all
  `m >= n`. (Matches the campaign's usage exactly.)
- **Lemma 3.4.** A ~p B (both log-concave) implies uA + vB log-concave
  for nonnegative u, v.
- **Theorem 3.5.** Pairwise partially synchronized families are closed
  under nonnegative linear combinations on both sides of ~p.
- **Theorem 3.6.** For A, B, C log-concave with no internal zeros:
  A ~p B implies A*C ~p B*C. (The common-convolution closure. Note the
  hypothesis on C: log-concave, no internal zeros.)
- Their Example 2.3 shows weak synchronicity is NOT preserved by
  convolution; partial synchronicity is the repair.

## The one-sided forms as consecutive TP2 conditions

For coefficient sequences B, C define, as in
`notes/c2_single_index_laws_2026-08-11.md`,

```
U_k = c_k b_k - c_{k-1} b_{k+1},   V_k = c_k b_k - c_{k+1} b_{k-1}.
```

**Reformulation (verified on 6,000 random exact instances).**
`V_k >= 0` for all k is equivalent to: every consecutive 2x2 minor of
the two-row array with first row C and second row xB (B shifted one
step right) is nonnegative. Likewise `U_k >= 0` for all k is the same
condition on the array with rows B and xC. In stochastic-order
language these are discrete likelihood-ratio (TP2) orderings between a
sequence and a shift of the other.

## Two inline lemmas

**Lemma K1 (PF2 = log-concave).** A nonnegative sequence with no
internal zeros is log-concave iff every 2x2 minor of its Toeplitz array
`(a_{i-j})` is nonnegative. Proof: the minors are
`a_r a_s - a_{r+1} a_{s-1}` for `r <= s`; log-concavity is the case
`r = s`, and the general case follows by chaining the nonincreasing
ratio sequence across the gap. (Polya-frequency context: Branden's
survey, Section on PF sequences, fetched full-text in the research
vault.)

**Lemma K2 (common-convolution closure of the one-sided forms).** If
the pair (C, B) satisfies `V_k >= 0` for all k, and c is log-concave
with no internal zeros, then (C*c, B*c) satisfies `V_k >= 0` for all
k; same for U. Proof sketch: by the reformulation, V is TP2 of the
2-row array [C; xB]. Convolving both rows by c multiplies the array on
the right by the Toeplitz matrix of c, which is TP2 by Lemma K1;
Cauchy-Binet writes every 2x2 minor of the product as a nonnegative
sum of products of 2x2 minors. Since `(xB)*c = x(B*c)`, the resulting
array is [C*c; x(B*c)]. Verified on 4,000 random exact instances
(products of random real-rooted polynomials). This is classical total-
positivity technology (Karlin's composition-formula genre); no novelty
is claimed, but the statement in LAW-V form was not previously in the
campaign's toolkit. It is the one-sided analogue of HWZZ Theorem 3.6.

## What K2 buys, and what it cannot buy

**Free reductions.** Any common polynomial factor of a pair can be
stripped or introduced without affecting the V/U property, provided the
factor is log-concave with no internal zeros. Applied to the bilinear
generator pairs of the C2 program: the dirty term V(Q_aQ_b, xR_aQ_b)
carries common factor Q_b, so its failure geometry is exactly the
failure geometry of the single-hub pair (Q_a, xR_a) transferred by K2.
Dirty-term analysis therefore reduces to single-hub pairs, one hub at a
time.

**The obstruction stays where it was.** Neither U, V, W nor TP2
orderings are closed under sums (mixtures). The generator-pair table in
`results/c2_single_index_laws_20260811.json` is the campaign's own
witness: termwise positivity fails while the aggregate holds. The
stochastic-orders literature knows the same phenomenon for
likelihood-ratio order under mixtures; the standard reference is
Shaked and Shanthikumar, *Stochastic Orders*, ch. 1.C (STATEMENT FROM
MEMORY, UNVERIFIED: no copy on disk; arXiv and Semantic Scholar were
rate-limited during this session; searches run:
`all:"likelihood ratio order" AND all:convolution`, blocked 429). Any
packet claim should cite HWZZ Theorem 3.5 (which IS on disk) for the
sum step instead, since partial synchronization is precisely the
condition that repairs mixture-closure for the symmetric form.

## Consequences for the proof queue

1. LAW-V for arm pairs sharing structure reduces via K2 to smaller
   pairs; the irreducible core is the single-hub pair (Q_a, xR_a)
   together with the aggregate (mixture) step.
2. The aggregate step is exactly where HWZZ Theorem 3.5 operates for
   the symmetric relation; a one-sided analogue of 3.5 with a
   compensation term is the missing lemma. The bilinear table says
   which cross terms the compensation must dominate.
3. LAW-W is the diagonal of the (xC, B) partial-sync form, so HWZZ
   machinery applies to it only through the full two-index relation;
   the packet should pose W's aggregate step as its own graded target.

## Verification

```
python3 - (inline scripts, session 2026-08-12): K2 on 4,000 random
exact instances; V<->TP2 reformulation on 6,000; V-preservation
corollary on 3,000. All passes; scripts embedded in the LAW-V packet
self-test (gpt_attack/law_v_packet_2026-08-12/).
```
