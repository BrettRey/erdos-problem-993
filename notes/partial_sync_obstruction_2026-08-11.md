# Where partial synchronization fails on the pattern-graph families

Date: 2026-08-11. Companion script: `verify_partial_sync_obstruction_20260811.py`;
results: `results/partial_sync_obstruction_20260811.json`. This note arose through
substantive generative-AI assistance. All computations are exact integer
arithmetic; nothing here is a proof beyond the finite instances checked.

## Status and scope

The double-broom log-concavity note carries a partial-synchronization
induction (Hu, Wang, Zhao, and Zhao, *MIA* 20 (2017), 91–103) through the
connector recurrence. This experiment asks where that machinery breaks on
the known non-log-concave families, using the positive two-term
decompositions supplied by Bautista-Ramos, Guillén-Galván, and
Gómez-Salgado (arXiv:2603.14204):

- k-direction, their Eq. (3), with pendant graph `P_2` rooted at a leaf:
  `F_k = P_1 · F_{k-1} + x · P_2^{k-1} · G`.
- m-direction, their Eq. (10):
  `U_m = (G\v) · Z_k(H)^m + x · (G\N[v]) · H^{km}`,
  with `Z_k(H) = H^k + x (H\w)^k`.

Both summands are products of log-concave polynomials, hence log-concave.
By the HWZZ sum lemma, if the two summands were partially synchronized the
sum would be log-concave; these families break log-concavity, so summand
synchronization must fail somewhere. The experiment locates the failure
set exactly.

Families checked (parameters in the JSON): double brooms `T_{5,5,t}`,
`t <= 25` (control); the Kadrawi–Levit structure `3,k,k` as
`(T_{1,3}:P_2)_k^{(2)}`, `k <= 24`; and the three m-direction families
`(T_{1,2}:S_{2,4})_2^{(m)}` (one break), `(T_{1,3}:S_{2,4})_2^{(m)}` (two
consecutive breaks), `(T_{1,7}:S_{2,5})_2^{(m)}` (three consecutive
breaks), `m <= 60` on the first and `m <= 30` / `m <= 24` on the others.

## Verification

Three independent checks per instance, all passing (134/134 closed-form
checks over 158 instances):

1. The decomposition polynomial equals the direct rooted-tree DP
   polynomial on the explicitly constructed tree (exact list equality).
2. The reflected head coefficients match the closed forms of
   arXiv:2603.14204 Corollaries 4.8, 4.12, 4.13, and 4.14 symbolically
   evaluated at each parameter, as do the degrees.
3. At `(T_{1,2}:S_{2,4})_2^{(9)}`, the computed Turán deficit at the break
   index is `-238214`, matching Corollary 4.12's constant term exactly.

The double-broom control shows zero violations of any relation at every
`t`, as the all-parameter theorem requires.

Independent local replay and audit, 2026-08-11: the externally supplied
artifacts (which covered `m <= 14` on the first m-direction family and
sampled parameters on the others) replayed byte-identically; the full
ranges above were then run locally. A separate audit script
(`scripts/audit_partial_sync_obstruction_20260811.py`) with freshly
written violation, Turán, and break computations reproduced every finding,
and brute-force subset enumeration (no DP) confirmed the independence
polynomials of an independently reconstructed `T_{5,5,2}`, the KL `k = 2`
tree, `S_{2,4}`, and `T_{1,2}`. The closed forms in check 2 were re-read
against the paper (`notes/literature/arxiv_2603_14204.txt`). The extension
corrected one threshold in F4 and exposed a small-`m` exception in F5,
both now stated as computed.

## Findings

**F1 (single-pair obstruction, KL family).** For `(T_{1,3}:P_2)_k^{(2)}`,
the entire violation set is one pair, `(m,n) = (alpha-1, alpha-1)`, at
every `k` from 2 to 24. The whole failure of the induction is a single
diagonal inequality at reflected index 1.

**F2 (synchronization fails before log-concavity does).** The KL summand
pair already fails at `k = 2` and `k = 3`, where the polynomial is still
log-concave; LC first breaks at `k = 4`. Likewise `(T_{1,2}:S_{2,4})`
fails synchronization from `m = 3` while LC survives to `m = 9`, and the
two- and three-break families fail synchronization already at `m = 2`,
the smallest instance run, against LC thresholds of 16 and 13 for their
eponymous break counts. Synchronization is strictly stronger than the
conclusion it feeds, and the parameter gap measures the slack available
to a weakened relation.

**F3 (thin diagonal band anchored at the top corner).** In every
m-direction instance the violation set is a band hugging the diagonal:
max `|m-n|` grows roughly like `m/2` while `alpha` grows like `10m`, so
the band's angular width relative to `alpha` shrinks (width 26 at
`alpha = 603`). The pair `(alpha-1, alpha-1)` is violated in every
non-clean instance; the pair `(alpha, alpha)` never is.

**F4 (band depth grows; final-third containment is a small-m artifact).**
The deepest violated reflected index grows linearly in `m`: for
`(T_{1,2}:S_{2,4})`, depth/alpha runs 0.25 (m=14), 0.42 (m=30), 0.47
(m=40), 0.51 (m=50), 0.54 (m=60), with marginal ratio ~0.7. Measured
against the Levit–Mandrescu decreasing zone `k >= ceil((2 alpha - 1)/3)`
(the convention in `paper/main.tex`), the band touches the zone boundary
exactly at `m = 20` (depth 68 = boundary 68) and escapes it at every
`m >= 21` (depth 74 against boundary 71 at `m = 21`).

**F5 (break structure vs violation structure).** In the mature regime the
LC breaks sit in the contiguous reflected top blocks {1}, {1,2}, {1,2,3}
in the one-, two-, and three-break families, while the sync-violation
band is far deeper. The anchor at reflected index 1 is not universal,
though: in the three-break family the block is {2} alone for
`3 <= m <= 7`, becoming {1,2} at `m = 8` and {1,2,3} at `m = 13`,
because Corollary 4.14's reflected-index-1 deficit stays negative below
`m = 8`. Deep-band diagonal deficits are absorbed by the summands' own
Turán surpluses; only at the head does the deficit win. Nothing in the
corpus produces violations that could seed a rebound with rise distance
greater than 1, consistent with the campaign's no-cross-gap-rebound law.

## Mechanism

For any decomposition `F = A + B`, exactly

```
Turan(F)_j = Turan(A)_j + Turan(B)_j + S_j,
S_j = 2 a_j b_j - a_{j-1} b_{j+1} - a_{j+1} b_{j-1},
```

where `S_j` is the diagonal (`m = n = j`) partial-synchronization form.
Verified exactly at all indices of `(T_{1,2}:S_{2,4})_2^{(9)}`. The
families break LC precisely where the head diagonal deficit `-S_j`
outgrows the two summand surpluses; at the break index of that instance
the three terms are `3945235 + 1089228 - 5272677 = -238214`. The
off-diagonal synchronization inequalities never fail outside the thin
band, and the diagonal ones fail long before they matter.

## Implications

1. Full partial synchronization is the wrong invariant to carry through
   pattern recurrences: it dies at `k = 2` on the KL family. What
   survives everywhere in the corpus is synchronization off a thin
   diagonal band anchored at the top corner.
2. The candidate inductively stable strengthening is therefore two-part:
   banded partial synchronization (all pairs with `m - n` above a small
   threshold, or both indices below the head window) plus a quantitative
   head condition bounding the diagonal deficit `-S_j` by the summand
   Turán surpluses. The second part is a reserve inequality at a
   descent (the same shape as the Poisson-binomial
   `V · Delta_eff >= 1/4` lane), which suggests the existing
   reserve toolkit is the right instrument for the head window.
3. For the double-broom extension (endpoint bristles to general two-hub
   trees): the c2-bounded-pendant-core certificate plus a banded-sync
   induction with a head reserve is the concrete next target. The
   failure geometry says exactly which inequalities the head condition
   must cover: the reflected diagonal indices 1 through O(m), nothing
   off-diagonal beyond width O(m)/alpha.

## Limitations

Finite corpus: five families, bounded parameters. The band-geometry
claims (F3–F5) are exceptionless observations over that corpus, not
theorems; the depth growth in F4 has no proved asymptote. The
partial-synchronization relation is taken as stated in the double-broom
note; HWZZ symmetry and closure claims are used as published, not
re-proved here. The `3*` (base `T_K`) families and Galvin's `T_{m,t}` in
the m-direction were not run; the drivers accept them with a one-line
builder addition.
