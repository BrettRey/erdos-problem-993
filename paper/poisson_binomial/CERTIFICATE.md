# Supplementary audit and replay manifest

This release candidate supports the compact-range computation in Proposition 3.1 of *A variance-scaled Turán inequality at the first descent of a Poisson–binomial mass function*. It certifies the strict positivity of 275 exact Bernstein coefficients in thirteen scalar cells. It does not machine-verify the probabilistic reduction or Theorem 1.1.

The archive also contains the complete conditional Lean project discussed in the manuscript's acknowledgement. That project formalizes recurrence propagation, endpoint exclusion, the first-crossing ratio bound, and the raw-from-effective corollary. It assumes the normalized Hillion–Johnson recurrence and does not formalize Proposition 3.1 or Theorem 1.1.

## Archive contents

| Archive path | Role |
|---|---|
| `scripts/verify_universal_pb_finite_bernstein.py` | SymPy generator; reconstructs the source identities and writes the compact and full certificates |
| `scripts/check_universal_pb_finite_bernstein_certificate.py` | Independently implemented standard-library checker |
| `scripts/build_poisson_binomial_supplement.py` | Deterministic archive builder and post-build verifier |
| `results/universal_pb_finite_bernstein_certificate_2026-07-10.json` | Compact summary certificate |
| `results/universal_pb_finite_bernstein_full_certificate_2026-07-16.json` | Full exact certificate |
| `formalization/pb_effective_drop_aristotle/` | Complete conditional Lean project, including `PBReserve/Core.lean`, build metadata, proof context, and provenance |
| `CERTIFICATE.md` | This audit and replay manifest |
| `LICENSE` | MIT license |
| `MANIFEST.sha256` | SHA-256 digest of every other archive member; generated when the archive is built |

The checker does not import the generator, SymPy, or formula strings from the certificate. It implements the two scalar-window formulas directly with `fractions.Fraction`, reconstructs the certified power-basis numerators by cross multiplication, converts them to the Bernstein basis, verifies strict positivity, and recomputes all per-cell and aggregate digests.

## Reference environment

- Python 3.14.6
- SymPy 1.14.0 for the generator
- Python standard library only for the checker and archive builder
- Lean 4.28.0 and Mathlib 4.28.0 for the conditional formalization

All three Python programs require Python 3.10 or later because their source uses modern type-annotation syntax. The exact replays below were completed with the reference versions listed above.

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

## Certificate summary

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
| `scripts/build_poisson_binomial_supplement.py` | `d42d0bcb37bf45f9e34371d3bfddb38137efea380bec41187c70d99abdb060c2` |
| `results/universal_pb_finite_bernstein_certificate_2026-07-10.json` | `6b91554d9ab1f43151e36c94c5c8c427c7bb057130b7f39d233b14c7ab3860c6` |
| `results/universal_pb_finite_bernstein_full_certificate_2026-07-16.json` | `5fbe0570403d3e49161e60371a8208e916895f002b15f68602c12bce9ed3aa69` |
| `LICENSE` | `8c6cac9c3f9dc235a38e5700048e097286a3f1e2cf5797aeee4577e0ca6970f0` |

## Scope and archive status

The source repository is <https://github.com/BrettRey/erdos-problem-993>. The deterministic ZIP built by the command above is the submission-ready release candidate; its adjacent `.sha256` file records the archive digest. The authoritative copy should be deposited with the journal's supplementary material and, if desired, in a separate immutable repository when the manuscript is submitted. No immutable public DOI has yet been assigned to this Poisson–binomial supplement, and the Zenodo DOI associated with the separate Erdős-problem paper must not be used for it.
