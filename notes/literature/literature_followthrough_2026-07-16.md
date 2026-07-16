# Literature follow-through and certified family stress test

Date: 2026-07-16

## Outcome

The high-value actions from the literature intake are implemented. The manuscript now cites the three directly relevant 2026 papers, the probability note attributes both of its known analytic inputs correctly, and the spectral notes no longer treat a fixed positive-axis sector as a viable universal conjecture.

The new certified stress corpus contains 41 source-defined tree families. Exact integer formulas and recurrences were replayed against the generic tree DP, and every root-ranking decision was made with Arb intervals rather than floating-point roots. The computations support a clean separation:

- Jerrum--Patel's positive-axis mechanism defeats a global fixed sector, but its accessible positive-axis roots are remote from the minimum-modulus root.
- Near-dominant non-real roots in the same trees remain cusp-like: at height 8 the closest non-real pair has modulus ratio 1.41113 and negative-axis deviation 0.04758 radians.
- The Bautista recurrence loci are also populated by roots well beyond the near-dominant window in the tested range.
- None of this proves the surviving cusp envelope or a saddle-sector bridge. It does remove two misleading shortcuts: a universal sector and an inference from an unranked limiting locus to dangerous near-dominant roots.

## Concurrent-work boundary

Claude Code owns `DECISIONS.md` and `gpt_attack/bridge_window_unimodality/` during this run. Those paths were not edited, read for reuse, reformatted, or included in any patch. The scope reservation and both Jerrum--Patel source corrections were sent through shared memory so that the live bridge-window run could account for them without changing its frozen packet.

Codex's write scope was limited to the manuscript bibliography/prose, literature and spectral notes, and the standalone family harness, tests, and result certificate. A final process and modification-time check showed no overlapping concurrent writes.

## Manuscript and provenance changes

`paper/main_v2.tex` and `paper/references.bib` now include:

- Bautista-Ramos, Guillén-Galván, and Gómez-Salgado (2026), including the journal DOI, for the recurrence unification and one-to-five consecutive log-concavity breaks;
- Hibi, Kara, and Vien (2026) for the symmetric-unimodal tree constructions;
- Levit and Kadrawi (2026), Lemma 2.16, as a contemporaneous version of the adjacent-mode bridge;
- source-specific wording for the Ramos--Sun order-27 discrepancy. The text records the abstract's endpoint, the contrary exhaustive certificate, and the absence of a public order-27 witness without calling the endpoint a typo.

Jerrum--Patel was not forced into the manuscript because the current paper does not state the spectral-sector program. Its effect belongs in the spectral working notes unless that program is later added to the paper.

The finite Poisson--binomial note now identifies:

- Hillion--Johnson (2016), Theorem A.2 and Corollary A.3, equations (78)--(79), as the primary source of the cubic coefficient inequalities;
- Bobkov--Marsiglietti--Melbourne (2022), Theorem 1.1 and Corollary 3.2, as the source of the lattice max-atom/variance inequality.

The remaining novelty claim is deliberately narrower: no screened source gave the same finite first-strict-descent bound, while the endpoint-aware propagation, two-sided modal windows, and propagation-to-variance synthesis appear new. Confidence is moderate--high, not a claim of historical priority. A MathSciNet/zbMATH pass and an expert check remain appropriate before publication.

## Jerrum--Patel source audit

Lemma 17 of arXiv:2510.01466v2 states positive-axis accumulation for complete binary trees with every edge subdivided a fixed even number of times. Its proof deserves a caveat:

1. the normal-family step directly yields contiguous subsequences of the periodic branching word, hence phase-truncated balanced trees rather than plainly the exact family in the lemma statement;
2. the invocation says `Delta=2` and `H_2`, while Buys's notation appears to require `Delta=3` and `H_3=<f_1,f_2>`.

The broader maximum-degree-3 phase-truncated family still appears to support the positive-axis accumulation conclusion, so a universal fixed sector is not tenable. The exact uniformly subdivided formulation should not be cited as unproblematic without clarification from the authors or a repaired proof.

The harness keeps the statement and proof-supported lanes separate instead of silently identifying them.

## Certified stress corpus

| Lane | Cases | Vertex orders | Maximum tree degree | Maximum polynomial degree |
| --- | ---: | ---: | ---: | ---: |
| Jerrum--Patel stated exact family | 11 | 7--511 | 3 | 341 |
| Jerrum--Patel proof phase truncations | 9 | 60--187 | 3 | 104 |
| Bautista P2-pendant circle families | 18 | 26--142 | 35 | 72 |
| Bautista consecutive-break families | 3 | 177--315 | 17 | 164 |

Across all 41 cases:

- the closed form or compressed recurrence equals the generic tree-DP polynomial exactly;
- every applicable Bautista recurrence residual is exactly zero;
- every polynomial is unimodal and squarefree;
- the minimum-modulus root is certified unique, simple, negative, and real;
- root multiplicities sum to the polynomial degree;
- every modulus threshold is unambiguous at 192-bit Arb precision.

### Positive-axis versus near-dominant roots

For the unsplit binary family, the minimum positive-axis angle decreases from 2.18381 radians at height 2 to 1.10264 at height 8. The height-8 root has modulus ratio 20.40894 relative to the minimum root. This is consistent with the asymptotic scale separation between the positive neutral activity 4 and the negative threshold -4/27, whose ratio is 27.

The nine tested endpoint phases of five `k=1` periods are shallower. Their smallest positive-axis angle is 1.62114 radians, at modulus ratio 118.16224. They do not numerically reach the source's asymptotic accumulation regime.

The same height-8 unsplit polynomial also contains a much closer non-real pair, at ratio 1.41113, but its negative-axis deviation is only 0.04758 radians. Thus the finite spectrum exhibits both phenomena at once: a remote positive-axis root and a near-dominant cusp root.

### Bautista recurrence families

For the three `k=32` P2-pendant cases, the closest non-real roots have modulus ratios 2.17870--2.19833 and negative-axis deviations 1.04423--1.04621 radians. The roots closest to the limiting circle `|z+1/3|=1/3` have residuals between 1.95e-5 and 2.67e-4, but their modulus ratios are 7.58530--9.98530. The limiting circle is therefore real in the expected numerical sense but remote from the minimum-root scale in these cases.

The three consecutive-break threshold cases reproduce the source exactly:

| Parameters `(ell,n,k,m)` | Order | Polynomial degree | Exact LC-failure indices | Closest non-real ratio | Negative-axis deviation |
| --- | ---: | ---: | --- | ---: | ---: |
| `(2,4,2,9)` | 177 | 93 | 92 | 1.60400 | 0.17488 |
| `(3,4,2,16)` | 312 | 164 | 162, 163 | 1.94143 | 0.31030 |
| `(7,5,2,13)` | 315 | 164 | 161, 162, 163 | 1.70631 | 0.24334 |

All three remain unimodal. Their roots closest to the source-specific equimodular locus occur at ratios 5.55467--7.01376, again distinct from the closest non-real pair.

## Artifacts and replay

- Harness: `scripts/stress_literature_root_families.py`
- Exact non-root tests: `test_literature_root_families.py`
- Arb certificate: `results/literature_root_stress_20260716.json`
- Certificate SHA-256: `5e206c9b1572277ba4e380cd0d5e17454374e61d681452de82c2acc82f48e5fe`
- Environment recorded in the JSON: Python 3.14.6, python-flint 0.9.0, initial Arb precision 192 bits.

Exact generation command:

```bash
/tmp/erdos993-flint-venv/bin/python scripts/stress_literature_root_families.py \
  --output results/literature_root_stress_20260716.json \
  --precision 192
```

The certificate was generated twice. After removing `generated_at`, elapsed-time fields, per-case root timings, and the newly added environment-version field, the two JSON trees were byte-for-byte identical.

## Validation

- `python3 test_all.py`: 55 tests passed.
- `python3 -m unittest -v test_literature_root_families.py`: 3 tests passed.
- `python3 -m py_compile scripts/stress_literature_root_families.py test_literature_root_families.py`: passed.
- Full XeLaTeX/biber/XeLaTeX/XeLaTeX build: passed; 27-page PDF, no undefined citations, and no biber warnings.
- Citation-key audit: no missing, unused, or duplicate keys.
- House-style audit: the two findings introduced in the first manuscript patch were corrected. The remaining 51 findings predate this work or are mathematical/code-subscript false positives.
- The only remaining overfull box is the pre-existing long artifact-path paragraph in the exhaustive-search section.

## Program-level verdict

The global positive-axis sector is closed as a proof route. The live spectral target should be stated only in saddle-relevant, modulus-ranked terms and should allow the sector width to shrink with size or critical regime. The most defensible surviving target is the near-dominant cusp envelope, coupled to a quantified saddle-sector profile; it remains a conjectural bridge input.

The immediate high-leverage follow-ups are external rather than another broad search: clarify the Jerrum--Patel proof scope, run a publication-grade database check for the finite Poisson--binomial novelty claim, and compare Will Blair's rooted-bush witnesses against the existing hard-family corpus before launching duplicate enumeration.
