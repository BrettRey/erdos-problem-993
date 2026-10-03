# Verifier return: new §4.1 (monotonicity route) and formalization claims

Saved by the parent session (Claude Code, Opus 5.5) on 2026-10-03, from the final message of a
fresh-context `responsibility-verifier` subagent (same model family as the parent, read-only tools).
The report below is verbatim.

Inspected manuscript: `paper/poisson_binomial/variance-local-log-concavity-poisson-binomial.tex`,
uncommitted working tree after the §4.1 switch (on top of commit 35a6fd0).

---

No RESPONSIBILITY NOTICE. The 32 inequalities are kernel-checked in Lean (`decide +kernel`, axioms limited to the standard three per LOCAL_REPLAY), which is machine evidence independent of this review. The new 4.1 route has had no cross-family review, but `submission/portal-fields-ECP-2026-10-03.md` line 137 already tells Brett that. I run on the parent's model family, so my agreement with LOCAL_REPLAY is not a second vote.

## Verdict
I found no mathematical error in the new §4.1, the restructured §4 preamble or §4.2. There are no dangling references. `Compact.lean` proves exactly what §4.1 describes. The claims about the formalizations do not overclaim. One wording problem in the supplement description is worth fixing before submission, and several files that ship in the zip now contradict the restructured paper.

## Directly verified (by reading and by hand arithmetic)

**§4.1 (lines 673-708):**
- `a = 1 − δ/(1−δ)`.
- For j ≤ K−1 and δ ≤ 1/(K+1), the factor 1−jδ ≥ 2/(K+1) > 0 and decreases in δ. So each R_r is nonnegative and nonincreasing.
- The identity holds: ∏_{s=1}^r (1 − sδ/(1−δ)) = ∏_{j=2}^{r+1}(1−jδ)/(1−δ)^r = L_r.
- Each factor lies in [0,1] exactly when (s+1)δ ≤ 1, and decreases in δ.
- Monotonicity of A_K in δ follows because A is a quadratic form with nonnegative coefficients, evaluated on nonnegative weights that are ordered componentwise. (3+δ)/(4δ²) decreases in δ.
- The endpoint chain at eq:endpoint-check is valid.
- m < H ≤ m+1 ⟺ 1/(m+2) ≤ δ < 1/(m+1), and on that range K = max{r : r < H} = m, which matches eq:K-a, since (r+1)δ < 1 ⟺ r < H.
- 3 < H ≤ 4 ⟺ 1/5 ≤ δ < 1/4, with K = 3. 20 × 1/400 = 1/20 = 1/4 − 1/5.
- The cells [1/5,1/4) and [1/(m+2),1/(m+1)) for m = 4..15 cover [1/17,1/4) exactly. 12 + 20 = 32.
- The remark checks out: at H = 7/2, ST = (1073/343)(1600/343) ≈ 14.59, below Q = 16.3125.

**Printed ratios, recomputed by hand, with their margin above the printed floor:**

| Case | Recomputed ratio | Margin |
|---|---|---|
| m = 4 | 30.5523159375 / 28.5 = 1.0720111 | 0.000011 |
| m = 5 | 1.43933 | ≈0.0003 |
| m = 6 | 1.76390 | ≈0.0009 |
| m = 15 | 777.18 / 221 ≈ 3.51664 | ≈0.0006 |
| K = 3, interval ending at 1/4 | (1082/81)/(129900/9801) = 1.0078676 | 0.000068 |

The K = 3 interval ending at 99/400 gives ≈1.0255, which supports "smallest on the last interval".

**§4.2 and the preamble:**
- Q(H) = (3H+4)(H+1)/4 under δ = 1/(H+1).
- eq:C-ratio is correct, including the difference 2r.
- A_sym = ST.
- The Bernstein conversion formula, S̃_J, T̃_J and the N_J expansion are correct.
- All five J = 5 coefficients match.
- β_0 and β_4 at J = 6 match (15557.5 and 99568). The π_i and μ_i agree character for character with `Large.lean`.
- b_r, Q, S, T and the Bernstein conversion are all defined before §4.2 uses them.
- Every `\ref`/`\eqref` resolves to an existing label.

**Lean against the manuscript:**
- `Compact.lean` uses unsymmetrized `A δ K`.
- `cell_K3` covers [1/5,1/4] in 20 steps of 1/400 with K = 3.
- `cell_K` covers [1/(m+2),1/(m+1)] with K = m, for m = 4..15.
- `compact_main` derives K ∈ {3..15} from (K+1)δ < 1 ≤ (K+2)δ.
- `Defs.lean` matches eqs. right-mass-bound, left-mass-bound and A-def. `Large.lean` follows §4.2: triangular cells, J = 5 direct, J ≥ 6 uniform.
- `PBReserve/Core.lean` (`curvature_propagation`, `endpoint_exclusion`, `crossing_ratio_lower_bound`, `raw_quarter_of_effective`) matches the supplement's sentence about the second formalization, including that it is conditional on the one-sided `hstep`.

**Check program:** `verify_pb_compact_monotone.py` uses the same lo, hi, K and target as the text, and computes A from the definitions as the half double sum.

**Acknowledgements:**
- "Supplied the monotonicity argument" is supported. The packet's PROMPT asks only for monotonicity in the weights; the δ-monotonicity plus endpoint check is new in Aristotle's return.
- The referee-reconstruction sentence is supported by `runs/pb-revision2-20261003/INPUT_chatgpt_pro_referee_report.md` lines 5 and 487 (all 275 coefficients, plus the π_i factorizations).

## Findings

1. **Supplement description, lines 845-848 (imprecision, the one a referee would trip on).** "The proof of Proposition 3.1 consists of the arguments of Section 4 together with the exact computations *these programs* perform" comes straight after the sentence about the 275-coefficient second check. The nearest antecedent is therefore that certificate's generator and checker, which the proof no longer needs.
   - Fix: "…together with the computations performed by the first two programs (the 32 inequalities of Section~\ref{sec:compact} and the expansions for $H\geq16$); the second-check certificate is not needed."

2. **Shipped supplement files contradict the current paper (imprecision, supplement only).** Lines 111-117 of `build_poisson_binomial_supplement.py` package every `.md` and `.lean` file in `formalization/pb_scalar_inequality_aristotle_result/`. These still describe the old §4.1:
   - `LOCAL_REPLAY.md` line 64 ("different route from the paper's")
   - `README.md` lines 38-39, 53-57 and 120-122 (eq. `asymmetric-P`, eq. `compact-P`, "printed identity A − Q = P(H)/…")
   - `ARISTOTLE_SUMMARY.md` lines 14 and 27
   - `Compact.lean` line 6, `Statement.lean` lines 17-24, `PaperIdentities.lean` lines 4-20, `PaperCells/Defs.lean` line 6

   Those labels and that identity are no longer in the manuscript.
   - Fix: don't edit Aristotle's files, because LOCAL_REPLAY's provenance claims rely on them being unchanged. Instead, add a dated paragraph to `CERTIFICATE.md` following its line-5 precedent, and append a dated note to `LOCAL_REPLAY.md`. Both should say the comments describe the earlier draft, that the 275-coefficient certificate is now the second check, and that §4.1 adopted the Lean route.
   - Consequence: the zip must be rebuilt, which changes the SHA in `portal-fields-ECP-2026-10-03.md` line 81.

3. **`LOCAL_REPLAY.md` lines 43 and 45 (imprecision).** They give 1.0079 and 3.517, while the paper prints the rounded-down values 1.0078 and 3.516. My values (1.0078676, 3.51664) show the two agree if the replay rounded to nearest. Add one line saying so.

4. **Line 659, "nonnegative" (fine-tier).** Line 714 divides by L_r in C_r = R_r/L_r, so it needs b_r > 0, which r ≤ K < H gives. Write "positive".

5. **Line 677, `ST` used before eq:S-T (fine-tier).** Add "(see \eqref{eq:S-T})". That label is otherwise unreferenced.

6. **Line 687 (fine-tier).** Write $0<\delta_-<\delta_+\leq1/(K+1)$.

7. **Lines 660-665 (fine-tier).** The partial-derivative argument only shows monotonicity on the nonnegative orthant, and §4.1 compares two weight vectors. A one-clause fix: "A is a polynomial in the weights with nonnegative coefficients, so A(v) ≤ A(u) whenever 0 ≤ v ≤ u componentwise". That is what Lean's `Aw_mono` proves.

8. **Lines 405-407 (fine-tier).** "Whose margins are printed" holds for the 12 cells but only the minimum is printed for the 20 K = 3 pieces. Optionally write "(for the 20 pieces at K = 3, the smallest)".

No error and no overclaim in: the last sentence of the abstract, "Section 4 is also formalized in Lean", the description of the two formalizations, or the Aristotle acknowledgement. Theorem 1.1 and Sections 2-3 are correctly described as not formalized.

## Inferred, not reproduced
- `lake build`, the `#print axioms` output and the escape-hatch grep: these come from LOCAL_REPLAY and README. I had no way to execute anything.
- Printed ratios for m = 7..14: only the program checks these (lines 65 and 76). The Lean kernel check proves the inequalities hold, not their third decimals.
- **Needed from the parent:** confirm that `python3 scripts/verify_pb_compact_monotone.py` printed `ALL CHECKS PASSED`. The script is untracked (`??`) but already hashed in CERTIFICATE.md line 145 and included in the zip, so it needs committing before `/ship`.

## What would reverse this
- The program failing on a printed ratio would make the corresponding number a text error; the proof would still stand.
- A Lean rebuild showing `sorry`, `ofReduceBool` or a changed statement.
- Evidence that `Defs.lean` differs from the packet's definitions.
