# Supplementary audit and replay manifest

This release candidate supports the computations in Proposition 3.1 of *Variance and local log-concavity of Poisson–binomial laws*. For 3 < H ≤ 16, the manuscript's Section 4.1 reduces Proposition 3.1 by monotonicity to 32 exact rational inequalities, which `scripts/verify_pb_compact_monotone.py` checks. As an independent second check of that range, the archive certifies the strict positivity of 275 exact Bernstein coefficients in thirteen cells. For H ≥ 16 a separate program checks, in exact arithmetic, every symbolic expansion printed in the manuscript's subsection "The range H ≥ 16". Proposition 3.1 is also formalized in full in Lean, and so is Theorem 1.1, including the reduction of Sections 2–3 and the cited inequalities it uses (both below).

The archive also contains the complete conditional Lean project described in the manuscript's supplement description. In the manuscript's numbering, that project formalizes Lemma 2.1 (the bound (2.8) and the endpoint exclusion), the bound q_D ≥ a, and the deduction of (1.6) from (1.3); its Lean names are `curvature_propagation`, `endpoint_exclusion`, `crossing_ratio_lower_bound`, `raw_drop_ge_effective`, and `raw_quarter_of_effective`. It assumes the recurrence (2.4). The full formalization of Theorem 1.1 below supersedes it, except for the two-line deduction of (1.6) from (1.3). The Lean comments predate the manuscript's terminology. "Curvature" there means the normalized deficit δ (not the log-curvature 𝒞_k). The recurrence they call "verified" enters as the hypothesis `hstep`. "Reserve", "raw" and "effective" are working names for the quantities in (2.8), (1.6) and (1.3).

The Lean project for Proposition 3.1 (`formalization/pb_scalar_inequality_aristotle_result/`) was produced against an earlier draft, in which Section 4.1 proved the range 3 < H ≤ 16 by Bernstein expansions. The manuscript then adopted the monotonicity route that this project's `PBScalar/Compact.lean` uses. Comments in that project's `README.md`, `ARISTOTLE_SUMMARY.md`, `Compact.lean`, `Statement.lean`, `PaperIdentities.lean` and `PaperCells/Defs.lean` that describe "the paper's" route, the identity A − Q = P(H)/(4H⁵(H+1)³), or the polynomials P_m therefore describe that earlier draft. In the current manuscript the Bernstein certificate is the supplement's second check of 3 < H ≤ 16, and the theorem statements and definitions are unchanged. The files are kept as Aristotle returned them; `LOCAL_REPLAY.md` records the replay and this change.

## Archive contents

| Archive path | Role |
|---|---|
| `scripts/verify_universal_pb_finite_bernstein.py` | SymPy generator; reconstructs the source identities and writes the compact and full certificates |
| `scripts/check_universal_pb_finite_bernstein_certificate.py` | Independently implemented standard-library checker |
| `scripts/verify_pb_compact_monotone.py` | Standard-library check of the 32 rational inequalities of Section 4.1, of the printed rounded-down ratios, and of the remark that the symmetrized bound fails at H = 7/2 |
| `scripts/verify_pb_cue_threshold.py` | Standard-library check of Example 1.6: V_1996 < 1 < V_1998 with rational bounds on π (Machin's formula), and the N = 10⁴ values |
| `scripts/verify_pb_large_h_range.py` | SymPy check of every expansion in the range H ≥ 16 (closed forms, the quartic N_J, the J = 5 coefficients, and β_i = μ_i π_i(u)/2880 identically in u) |
| `scripts/build_poisson_binomial_supplement.py` | Deterministic archive builder and post-build verifier |
| `results/universal_pb_finite_bernstein_certificate_2026-07-10.json` | Compact summary certificate |
| `results/universal_pb_finite_bernstein_full_certificate_2026-07-16.json` | Full exact certificate |
| `formalization/pb_effective_drop_aristotle/` | Complete conditional Lean project, including `PBReserve/Core.lean`, build metadata, proof context, and provenance |
| `formalization/pb_deduction_aristotle_result/` | Complete Lean project proving Theorem 1.1 (`PBDeduction.theorem_1_1`) from the probability-generating polynomial, including Sections 2–3, the Hillion–Johnson cubic inequalities, strict log-concavity and the maximal-mass bound; it contains a copy of `PBScalar/`. Its `PROOF_CONTEXT.md` is omitted from the archive because it quotes Hillion and Johnson's Appendix A verbatim |
| `formalization/pb_scalar_inequality_aristotle_result/` | Complete Lean project proving Proposition 3.1 (`PBScalar.scalar_inequality`), with Aristotle's summary, the local replay record `LOCAL_REPLAY.md`, build metadata, prompt and proof context |
| `CERTIFICATE.md` | This audit and replay manifest |
| `LICENSE` | MIT license |
| `MANIFEST.sha256` | SHA-256 digest of every other archive member; generated when the archive is built |

The checker does not import the generator, SymPy, or formula strings from the certificate. It implements the two scalar-window formulas directly with `fractions.Fraction`, reconstructs the certified power-basis numerators by cross multiplication, converts them to the Bernstein basis, verifies strict positivity, and recomputes all per-cell and aggregate digests.

## Reference environment

- Python 3.14.6
- SymPy 1.14.0 for the generator and the H ≥ 16 check
- Python standard library only for the checker and archive builder
- Lean 4.28.0 and Mathlib 4.28.0 for the three Lean projects

All six Python programs require Python 3.10 or later because their source uses modern type-annotation syntax. The exact replays below were completed with the reference versions listed above.

## Exact replay

Run from the repository root:

```bash
python3 scripts/verify_universal_pb_finite_bernstein.py \
  --out /tmp/universal_pb_finite_bernstein_summary.json \
  --full-out /tmp/universal_pb_finite_bernstein_full.json

cmp /tmp/universal_pb_finite_bernstein_summary.json \
  results/universal_pb_finite_bernstein_certificate_2026-07-10.json
cmp /tmp/universal_pb_finite_bernstein_full.json \
  results/universal_pb_finite_bernstein_full_certificate_2026-07-16.json

python3 -I -S scripts/check_universal_pb_finite_bernstein_certificate.py \
  results/universal_pb_finite_bernstein_full_certificate_2026-07-16.json
```

The two `cmp` commands must produce no output and return status 0. The checker must report `"status": "passed"`, 13 cells, 275 coefficients, and payload digest

```text
64cacd6c220fc3b67250de8c761adbd4cf451fdf025e1c6579db5c052999629a
```

Check the 32 inequalities of Section 4.1 (range 3 < H ≤ 16) with:

```bash
python3 scripts/verify_pb_compact_monotone.py
```

It must print `ALL CHECKS PASSED` and exit with status 0.

Check the range H ≥ 16 with:

```bash
python3 scripts/verify_pb_large_h_range.py
```

It must print `ALL CHECKS PASSED` and exit with status 0.

Check the threshold in Example 1.6 with:

```bash
python3 scripts/verify_pb_cue_threshold.py
```

It must print `certified V_1996 < 1 < V_1998: True` and `ALL CHECKS PASSED`. Its decimal enclosures are rounded outward in exact integer arithmetic.

Build and verify the deterministic supplementary archive with:

```bash
python3 scripts/build_poisson_binomial_supplement.py
cd paper/poisson_binomial
shasum -a 256 -c poisson_binomial_certificate_supplement.zip.sha256
```

Build the conditional Lean project from the repository root with:

```bash
cd formalization/pb_effective_drop_aristotle
lake exe cache get
lake build
```

The expected result is a successful build of `PBReserve.Core` and `PBReserve`. The project contains no `sorry`, `axiom`, `admit`, or `implemented_by` declaration. The wrapper `PBReserve.lean` contains only `import PBReserve.Core`; the proved theorem bodies are in `PBReserve/Core.lean`, which therefore has to accompany the wrapper.

Build the Lean proof of Proposition 3.1 and print its axioms with:

```bash
cd formalization/pb_scalar_inequality_aristotle_result
lake exe cache get
lake build
lake env lean PBScalar/Axioms.lean
```

The build must succeed, and `PBScalar.scalar_inequality`, `PBScalar.scalar_inequality_compact` and `PBScalar.scalar_inequality_large` must each depend only on `propext`, `Classical.choice` and `Quot.sound`. The theorem `scalar_inequality` states Proposition 3.1 for `0 < δ < 1/4` with `K` characterized by `(K+1)δ < 1 ≤ (K+2)δ`; the definitions are in `PBScalar/Defs.lean`. Its proof of the range 3 < H ≤ 16 (`PBScalar/Compact.lean`) is the monotonicity argument of Section 4.1, and of H ≥ 16 (`PBScalar/Large.lean`) follows Section 4.2. `LOCAL_REPLAY.md` records the replay, the escape-hatch search, and the statement comparison.

Build the Lean proof of Theorem 1.1 and print its axioms with:

```bash
cd formalization/pb_deduction_aristotle_result
lake exe cache get
lake build
lake env lean PBDeduction/Axioms.lean
```

The build must succeed, and all eleven theorems listed there, including `PBDeduction.theorem_1_1`, must depend only on `propext`, `Classical.choice` and `Quot.sound`. In `theorem_1_1`, the mass function `pbPmf n p` is the coefficient sequence of `∏ i (1 − p i + p i X)`, `pbVar n p = ∑ i p i (1 − p i)`, `IsFirstDescent` characterizes D, and `deficit` is δ_k. The definitions are in `PBDeduction/Defs.lean`. `LOCAL_REPLAY.md` in that folder tabulates the correspondence with the paper's statement.

## Certificate summary (second check of 3 < H ≤ 16)

The minimum in the table is the minimum reduced rational Bernstein coefficient of the stated cleared numerator. Its magnitude depends on the chosen normalization; strict positivity is the invariant claim.

| Cell | Interval | Degree | Count | Minimum | Bernstein-vector SHA-256 | Cell-payload SHA-256 |
|---|---:|---:|---:|---:|---|---|
| `asymmetric_3_4` | [3,4] | 10 | 11 | 22272 | `e1e452c9c93c30cd33ff0470c6d81a4c0bfb79fe07c36459cd4daaa2ffb789d5` | `b3e30a488c2726bc94f2169875e9e8af1a22cc22290e8a4e531bb0247dd60b2e` |
| `symmetric_4_5` | [4,5] | 10 | 11 | 332800 | `441f603ac3f60919665f154801ae9fec150ba2dac02382ac870105f92067c356` | `a709542891e6e44853140520fe94208131dacdfe150f565bf05d35676494625f` |
| `symmetric_5_6` | [5,6] | 12 | 13 | 476945150 | `80a652bad6b13e69c4f3200f6140f0f7ceb8d681b433a817927fcfde6f4c80be` | `e2add7f41523bf783a0441b08a745431daf49c36684bc89284a6b939309fc649` |
| `symmetric_6_7` | [6,7] | 14 | 15 | 252843503616 | `e3ad5553e47317b267ea56a45a2d853029a9a602545ba7f2357462a617cd0d85` | `92fa3fcd6f370f2d21aa96b78498680b1197ef84925004ec98e0540688c31031` |
| `symmetric_7_8` | [7,8] | 16 | 17 | 141578177942152 | `6c774572e048423c105ff978b7042abb398369666ac5f6a6f903df5e8fe7f980` | `8c2282305e86e273a5070ab0169f6b5533c134f2eb041bfc6d5abad2a2da6083` |
| `symmetric_8_9` | [8,9] | 18 | 19 | 92316620768411648 | `e9b192d522d49b1e66d1948d4bb774cce9eab468e3bca0f4ed1b06a69d467169` | `c3c331d7e51679168450e1e0c1a1229e5212295cd9951813da6bffd2914b500c` |
| `symmetric_9_10` | [9,10] | 20 | 21 | 71284663153149460650 | `876f43f48d0c3c1ea3a7fa237ac1ea527a7855286ebc58ac891969ca880718d1` | `0143a7121b548a0b77a675c0894bc4bed57c1c9e55c3279f1ed6a8110ba5b792` |
| `symmetric_10_11` | [10,11] | 22 | 23 | 65053547560474378240000 | `b741906966aed963db8b6dd41732fa6fd1ffab2ef99e7491db3ae904fd033091` | `9a4a971579e0fcbf88b26ed35793a9c30fca4aebb3c94c55c1f0c6bb9792caac` |
| `symmetric_11_12` | [11,12] | 24 | 25 | 69645362929101036447979796 | `6c318dad1104904865f9aa3c865bad147fed1a827cd868b49c2bb5d66476f70a` | `dedb6dcb3cb94435e6158a37b424cd50bb376952951c80bcd239a4b40374240a` |
| `symmetric_12_13` | [12,13] | 26 | 27 | 86706679760909016273411637248 | `643c6c7cf85a5535537a9eb95135fb2a8954b4866b83fbd40fa41a1a24c41728` | `c0f55dbd54755c2105f7f0ae8b52edeb1d3e5190351490066682a2bb44ad006b` |
| `symmetric_13_14` | [13,14] | 28 | 29 | 124437022735783915561664544889462 | `64f915da71a6e95d98d0b4451704204ee7a447fd71e6090dfd134ab3231bae18` | `b86f0b0b68b64edd3b068f3c839abf4ffaec1751c07dbb3eb9999c7aadcb2369` |
| `symmetric_14_15` | [14,15] | 30 | 31 | 204172536352322626071425358548172800 | `e31285e840aa69d9882286637143a1c78da0a2ff88ca3527d82e8c1e04b7f234` | `c9417f52e15bb935c4315e0294b3cf853c029d4545e83fec84f6273e47089221` |
| `symmetric_15_16` | [15,16] | 32 | 33 | 380093996939768323066638847260414000000 | `0850f4a393088df01ab34888ea804ee55937c3a9de9948c7a8166b365f5beb4c` | `2e97cdbe380bf1a8661d87aaa22db2d1040a4f3aa25d1dccaee1dc04ac591d07` |

Aggregate certificate digests:

- Compact certificate combined-cell digest: `c920cc3bc11eb1564047645ef6b8dd4221efa834b439115bcbad6fb8fdfe4330`
- Full certificate canonical-payload digest: `64cacd6c220fc3b67250de8c761adbd4cf451fdf025e1c6579db5c052999629a`

Whole-file digests before packaging:

| Repository path | SHA-256 |
|---|---|
| `scripts/verify_universal_pb_finite_bernstein.py` | `7913038a93a18cc9df0fbe770238543d80655de82c35fac9a48247ed1f8e1b61` |
| `scripts/check_universal_pb_finite_bernstein_certificate.py` | `f5762de3d7990f82d62e63d5f7007b6f9ec62b60eea325f6b68354b34ee146a7` |
| `scripts/build_poisson_binomial_supplement.py` | `e8ed0fdbd9129ab01ccf5d9c4c58e4c0f79dda5abb1dcf631df6ae32ea7ae429` |
| `scripts/verify_pb_compact_monotone.py` | `60cd99b4cb159672a30fbf74d78e029b54ba575fc5b788d4ac6492ce5957cfb4` |
| `scripts/verify_pb_cue_threshold.py` | `6c3268b8e90143883b0f2decdc33e603f417c312c4d686d89b63b7c2815c4e99` |
| `scripts/verify_pb_large_h_range.py` | `0a5e18d2605f79b348891b249cb1b86899b13ce247c48cb061e9cd5e48c0a81d` |
| `results/universal_pb_finite_bernstein_certificate_2026-07-10.json` | `6b91554d9ab1f43151e36c94c5c8c427c7bb057130b7f39d233b14c7ab3860c6` |
| `results/universal_pb_finite_bernstein_full_certificate_2026-07-16.json` | `5fbe0570403d3e49161e60371a8208e916895f002b15f68602c12bce9ed3aa69` |
| `LICENSE` | `8c6cac9c3f9dc235a38e5700048e097286a3f1e2cf5797aeee4577e0ca6970f0` |

## Scope and archive status

The source repository is <https://github.com/BrettRey/erdos-problem-993>. The deterministic ZIP built by the command above is the submission-ready release candidate; its adjacent `.sha256` file records the archive digest. The authoritative copy should be deposited with the journal's supplementary material and, if desired, in a separate immutable repository when the manuscript is submitted. No immutable public DOI has yet been assigned to this Poisson–binomial supplement, and the Zenodo DOI associated with the separate Erdős-problem paper must not be used for it.
