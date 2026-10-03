<!-- Saved verbatim by the parent session from the subagent's final message (Claude Opus via responsibility-reader, read-only; level-category-audit pass adapted to philosophical category errors at Brett's direction, 2026-10-03). Subagents cannot write report files. The parent then applied findings 1-11, the "Three steps" and "uniformly" accuracy items, and "at the rightmost mode". The Lean docstring points were handled by a terminology note in CERTIFICATE.md, leaving the frozen Aristotle output unedited. A new supplement script (scripts/verify_pb_cue_threshold.py) now certifies the N = 1998 threshold. -->

LEVEL AND CATEGORY AUDIT (read-only), variance-local-log-concavity-poisson-binomial.tex

**Version note.** The manuscript changed while I was auditing it. My first read had 836 lines. My second read had 863 lines, adding Table 1 at l.304-323 and an itemized AI acknowledgement at l.838-857. DECISIONS.md l.997-1008 records that change. All line numbers below refer to the 863-line version, and I audited the new material too. Concurrent sessions are normal under the house rule, so this is a note, not a RESPONSIBILITY NOTICE. If the file has moved again, the quoted text will locate each item.

**Method.** I adapted the registry's four checks (/Users/brettreynolds/projects/LLM-CLI-projects/passes/registry/level-category-audit.yaml) to a mathematics paper:
- The levels are the random variable W, its law, its pmf f, its generating polynomial, and derived quantities (V, δ_k, C_k, D, q_k).
- The kinds of epistemic standing are: hand proof, computer-generated exact data, printed hand-checkable data, Lean-formalized step, heuristic, and numerical observation.

Check 1 (target defined before mechanism) passes. δ_k (l.75-79), C_k (l.84-88) and D (l.117-119) are all defined before the theorem, and the proof outline (l.378-399) comes after it. Everything I found falls under checks 2 to 4 and your item 5.

## Ranked findings (most misleading first)

**1. The standing of the computer-assisted step is blurred in every summary (items 2 and 4)**
- l.36-37 (abstract): "reduce the main inequality to a one-variable inequality, which is proved in exact arithmetic by Bernstein expansions."
- l.393-394: "Section~\ref{sec:scalar} proves the one-variable inequality in exact arithmetic by Bernstein expansions"
- l.825-826 (supplement): "Together with the arguments of Section~\ref{sec:scalar}, the programs prove Proposition~\ref{prop:scalar}"
- l.826-828: "the reduction of Theorem~\ref{thm:main} to it ... is proved in the text."

Why this is a category error and not idiom:
- "In exact arithmetic" says what kind of arithmetic was used, not who or what did it. The abstract never says that a computer produced 275 coefficients the paper does not print. Only the keyword "computer-assisted proof" (l.16) says so.
- The sentence also merges two standings. For H≥16 the data are printed and can be checked by hand (five rationals for J=5, and π_i polynomials with positive coefficients). For 3<H≤16 the proof depends on unprinted, machine-computed coefficients.
- Programs compute and check. The proof is the argument: Bernstein basis nonnegativity plus positivity of the computed coefficients. The supplement sentence makes "the programs" the subject of "prove".
- "The reduction ... is proved" treats an argument as a proposition.
- The Section 4 body itself is careful: "A SymPy computation (in the supplement) gives" (l.688), l.729-733, l.812-814. The slip is only in the summaries.

Repairs:
- Abstract: "...which is proved by Bernstein expansions, with a computer-assisted step in exact rational arithmetic."
- l.393-396: say which ranges use the 275 machine-computed coefficients (3<H≤16, i.e. 1/17≤δ<1/4), which use the five printed ones (1/22≤δ≤1/17), and which use the polynomial-coefficient argument (δ≤1/22).
- Supplement: "The proof of Proposition 3.1 consists of the arguments of Section 4 together with the exact computations these programs perform; the deduction of Theorem 1.1 from Proposition 3.1 (Sections 2 and 3) is given in full in the text and uses no computation."

**2. Acknowledgement: "each was checked before it was adopted" (l.843-844) has no agent and no method (item 2). Brett has to decide the wording.**
- In an AI-use disclosure, who did the checking is exactly what the reader wants to know. The passive invites the reading "the author verified these proofs".
- The project records describe model-based and script-based checks:
  - DECISIONS.md l.950: "Its mathematics was treated as candidate material, and each claim was checked before use."
  - l.962-963: an exact-arithmetic script over 1,839 laws and a same-family verifier.
  - l.967: "every hand-check of the written proofs so far is Claude-family."
  - l.983: the ChatGPT Pro check was "the first non-Claude check of the written proofs."
- I found no record either way of Brett checking them himself.
- Repair: name the checker and method truthfully. For example, "each was checked, by independent model-based proof review and exact computation on test laws [and by the author, if true], before it was adopted."
- Same paragraph, l.846-847: "independently reconstructed the computations of Section 4 and Example 1.6". It was the model that did the reconstructing, not the report, and the sentence doesn't say the results agreed. DECISIONS.md l.983 says "All 275 coefficients equal ours exactly." Repair: "...and its independent reconstruction of the computations of Section 4 and Example 1.6 agreed exactly with ours."

**3. Example 1.6, l.346-348: "Meckes and Meckes ... show that $X_N$ is a sum of $N$ independent Bernoulli variables" (item 3: random variable vs law)**
The cited results are equalities in distribution, and HKPV show that the distinction matters:
- Meckes–Meckes Prop. 3(1), /Users/brettreynolds/projects/LLM-CLI-projects/literature/meckes_2017_empirical_spectral_measure_unitary.md l.37: "there are independent Bernoulli random variables ξ1,..., ξN such that N_A =d Σ ξj" (equality in distribution).
- HKPV Theorem 7, /Users/brettreynolds/projects/LLM-CLI-projects/literature/hough_2006_determinantal_processes_independence.md l.53: the count "has the distribution of a sum of independent Bernoulli(λk) random variables."
- HKPV Example 23, same file l.145: "the Ik's are not measurable w.r.t. the process X in general."

So the cited result is misreported as an identity of random variables. The mathematics is unaffected, because Theorem 1.1 depends only on the law, and l.367-368 does use the right type ("have the same law"). Repair: "show that $X_N$ has the law of a sum of $N$ independent Bernoulli variables."

**4. l.92-94: "we ask whether this scale holds without approximation for every Bernoulli sum: does the variance alone force $\delta_k\geq\kappa/V$ near the mode for a universal $\kappa>0$?" (items 1 and 4)**
- V is a number and forces nothing. More seriously, "alone" credits the variance with something the paper itself shows needs the Bernoulli-sum structure.
- Proposition 1.5 (l.257-265) shows that the variance together with ULC does not give the bound. The proof uses Hillion–Johnson cubic inequalities, which are specific to Bernoulli sums.
- "This scale holds" is a smaller slip of the same kind: a bound on that scale holds, not the scale.
- Repair: "we ask whether every Bernoulli sum satisfies $\delta_k\geq\kappa/V$ near the mode for a universal $\kappa>0$, with no dependence on $n$ or on the $p_i$." This also matches the next paragraph's contrast between n and V (l.96-98).

**5. Abstract l.31-32: "Ultra-log-concavity alone gives no such bound." (item 5: a negative result about one quantity reads as being about another)**
- The nearest antecedent is the two-sided bound on δ_c. Proposition 1.5 proves only the δ_D case. The body is correctly scoped at l.257 ("No variance-scaled lower bound on $\delta_D$ holds for ULC laws"), and DECISIONS.md l.987 records "The ULC claim is scoped to δ_D". That fix reached the body but not the abstract.
- The paper's own example cannot support the δ_c reading, because g^(m) has δ_m ≥ 3/4 at its mode (l.280).
- I believe the δ_c version is true. Hand sketch, not computed: take g_{m+j} ∝ C(2m,m+j)·r(j), with r=1 for |j|≤1 and 2^{-(|j|-1)} beyond. This is ULC(2m) with unique mode m, δ_m = 1/(m+1), and variance tending to 24/5. Even so, the paper doesn't prove it.
- Repair: "Ultra-log-concavity alone gives no variance-scaled lower bound on $\delta_D$", or move the sentence before the δ_c clause.

**6. Table 1 (l.304-323): an inequality sits in a column of values (items 2 and 3)**
- The caption says "Lower bounds on $\delta_D$" and the column is "Value at $D$". The Theorem 1.1 entry is "$\geq1/4$", which is a statement about δ_D, not the value of the bound. Every other row gives the bound's value, and here that value is exactly 1/(4V) = 1/4.
- The final row, "$\delta_D$ itself", is the target quantity rather than a lower bound. It is set off by a rule, but the caption covers it.
- Repairs: make the Theorem entry "$1/4$" (or "$1/(4V)=1/4$"). Change the caption to "Lower bounds on $\delta_D$, and $\delta_D$ itself, in the balanced family...".
- The Johnson row's "$O(m^{-2})$" is an upper estimate of a lower bound's value. That is enough to show it tends to zero, so it's acceptable.

**7. l.280-282: "a ULC law can change curvature abruptly at the mode, which the reciprocal bound \eqref{eq:reciprocal} forbids for Bernoulli sums." (items 1 and 3)**
- Curvature belongs to log g (C_k, l.86), not to the law.
- What the sentence exhibits jumping is 1/δ. The paper set up δ_k as a bounded transform of C_k at l.84-88, and the qualitative point carries over, but the levels are swapped without saying so.
- "The bound forbids" makes an inequality an agent.
- The jump is between m and m+1, not "at the mode".
- Repair: "For $g^{(m)}$, $1/\delta_m\leq4/3$ while $1/\delta_{m+1}\to\infty$, so $1/\delta$ changes by an unbounded amount in one step; by \eqref{eq:reciprocal}, this cannot happen for a Bernoulli sum."

**8. Abstract l.29-30: "Cubic inequalities of Hillion and Johnson then give $\delta_k\geq1/(4V+|k-D|)$ at every $k$ in the support" (item 4)**
- The support bound needs Theorem 1.1 plus the reciprocal bound (2.5). That bound in turn needs the cubic inequalities and the positivity of δ_k from (1.7) (l.482).
- "Then" partly signals the dependence, but the part is still the grammatical subject.
- Repair: "Combined with cubic inequalities of Hillion and Johnson, this gives...". The matching sentence in the introduction (l.183-186) is fine, because it says "so the bound at $D$ extends".

**9. Computed facts stated flatly, and their scripts are not in the archive (item 2)**
- l.369-370: "among even $N$ it first reaches $1$ at $N=1998$." This is a finite computation (V_1996 < 1 ≤ V_1998). By my hand estimate from the displayed asymptotic, the margins near the crossing are about 4e-5 and 6e-5. So the result can't be read off "$+O(N^{-2})$" without an explicit error constant, and no method is stated.
- DECISIONS.md l.962 names `scripts/verify_pb_corollaries_20261003.py`, and l.998 names `scripts/verify_pb_table1_20261003.py`. The archive table in /Users/brettreynolds/projects/LLM-CLI-projects/papers/queue/erdos-problem-993/paper/poisson_binomial/CERTIFICATE.md l.9-20 lists neither.
- The acknowledgement (l.846-847) mentions "computations of ... Example 1.6" that the reader cannot find in the archive.
- Repair: "a direct evaluation of the middle form, with rigorous bounds on π, gives $V_{1996}<1\leq V_{1998}$". Then add that script to the archive and to the supplement description.
- The Table 1 rates are derived in the text (l.289-302), so they need no change. The figures at l.372-373 are marked "about", which is fine.

**10. l.234: "applied to these tilted laws (which are Bernoulli sums)" (item 3, low)**
A law is not a sum. Elsewhere the paper keeps the types apart: "Bernoulli sum" for random variables (l.409, l.437-438) and "Poisson–binomial" for laws (l.57, l.138). Repair: "(which are again Poisson–binomial)".

**11. l.70-71: "so Newton's inequalities make the pmf log-concave" (item 1, borderline, low)**
Inequalities don't make anything; they establish a property. Repair: "so, by Newton's inequalities, the pmf is log-concave".

## The supplement artifact (low rank, but it ships with the paper)

/Users/brettreynolds/projects/LLM-CLI-projects/papers/queue/erdos-problem-993/formalization/pb_effective_drop_aristotle/PBReserve/Core.lean
- **l.11:** "a truthful, reusable Lean verification of the recurrence". The recurrence is not verified. It is the hypothesis `hstep` (l.75, l.117), and l.8-9 of the same header says so. This contradicts the paper's accurate "conditional on" (l.829). Item 2.
- **δ is called "curvature" throughout:** l.4, l.19 ("propagated upper bound for curvature"), l.66-67 ("If neighboring curvatures obey..."), and l.108 ("An endpoint has curvature one"). In the paper's terms the endpoint value of C_k is infinite and δ is 1 (l.84-88). Item 3.
- **Internal working names the body has dropped:** "Poisson-binomial reserve" (l.4), "universal finite Poisson-binomial effective-drop theorem" (l.6-7), and "raw"/"effective" in theorem names, which CERTIFICATE.md l.5 lists.
- Repairing only the docstrings changes no proofs, but the archive would have to be rebuilt and rehashed.

## Other observations (accuracy, not category errors)

- **l.397-399, "Three steps can lose a constant".** The list leaves out the replacement of b_r by λ_r and the truncation at J in §4.2 (l.745-749, l.765-767). By my hand asymptotics, that step costs about a factor of 4 for large H. Also, symmetrization isn't applied on 3<H≤4 (l.687-688). Repair: "Several steps can lose a constant, among them..."
- **l.639-640, "Section~\ref{sec:large} treats $H\geq16$ uniformly".** This is inconsistent with l.743-744: "The first interval ($J=5$) is checked directly, and the rest ($J\geq6$) uniformly in $J$."

## Judged to be acceptable mathematical idiom

- **"an explicit family of binomial laws shows that no constant above $1/3$ is possible" (l.27-28).** A witness family plus the computation in Proposition 1.2's proof is a proof that κ⋆ ≤ 1/3.
- **"Darroch's theorem ... locate[s] the mode" (l.110-114).** A theorem stating where every mode lies does locate it. The sentence then correctly says these results "by themselves give no lower bound".
- **"a small value of $\delta$ ... forces small values nearby" (l.378-381) and "a small deficit at $D$ forces small deficits nearby" (l.404).** "Forces" means logical necessity under a displayed constraint, and the forced quantities are values of the same sequence.
- **"Heuristically, ... the two bounds force $V$..." (l.390-392).** The heuristic standing is marked, and the rigorous version follows.
- **l.183-186, "imply that $1/\delta_k$ changes by at most one per lattice step".** "Changes" means variation in k. "Imply" leaves out positivity and the endpoint values, which l.482-487 supply.
- **Ordinary metonymy:** "Section~\ref{sec:reduction} checks" (l.61), "extends" (l.211), "Log-concavity says" (l.450), "Rotation by $\pi$ shows" (l.367), "The middle form shows" (l.368-369), "Determinantal point processes supply" (l.341).
- **Abstract l.25-26, "for a normal density of variance $V$ that curvature is $1/V$".** log φ is quadratic, so its unit-spacing second difference is exactly 1/V. The body (l.89-91) marks the comparison as a suggestion of scale.
- **Careful level marking worth keeping:** l.137-139, l.145-146, l.288-289, l.336-341, l.367-368, l.300-302 ("converges in law"), and l.437-442.
- **Title "local log-concavity".** δ_k quantifies the slack in a property every Poisson–binomial law has. That's a defensible compression for a title.
- **The Lean description at l.828-834 is accurately conditional.** "Of the deficit bound in Lemma 2.1" means an abstract analogue whose instantiation to the pmf is not formalized, and "for an abstract nonnegative sequence" signals that.
- **"Aristotle ... produced the Lean proofs" (l.853-855) is accurate.**

## Comparisons with prior work (item 5): evidenced negative

- **Search:** I grepped for "stronger|weaker|by contrast|whereas|against order|improv|sharper|only \(". The hits were l.34, l.333, l.371, and l.428 (internal, not a comparison).
- **Read in full:** l.96-114, Table 1 (l.304-323), l.370-373, and l.325-334.
- **Result:** every comparison is like with like.
  - The intro paragraph, Table 1, Example 1.6 and the abstract compare lower bounds on δ_k or δ_D.
  - l.325-334 compares Pitman's upper bound on f_{D+1}/f_D with (1.6) rewritten as an upper bound on the same ratio.
  - Dümbgen–Wellner's δ_k > 1/(k+1) is correctly said to follow from a ratio result, and Pitman's (20) is correctly described as ratio bounds that give a bound on δ_k.
- **The old ratio-vs-deficit comparison is gone** (DECISIONS.md l.986). What remains are findings 5 and 6.

**Filler verbs:** a grep for "underl|involv|regulat|encod|represent|drives|mechanism" found nothing.

## Remaining drift I judge defensible

- l.725-729 says all 275 coefficients are positive before saying who computed them, but the same paragraph discloses that (l.729-733).
- "At the mode" (l.214), where there can be two modes. "At the rightmost mode $c$" would be exact.
- The "forces" idiom above.

Brett should accept these explicitly, or have them repaired.

## Files

- Audited manuscript: /Users/brettreynolds/projects/LLM-CLI-projects/papers/queue/erdos-problem-993/paper/poisson_binomial/variance-local-log-concavity-poisson-binomial.tex
- Pass procedure: /Users/brettreynolds/projects/LLM-CLI-projects/passes/registry/level-category-audit.yaml
- Supplement manifest: /Users/brettreynolds/projects/LLM-CLI-projects/papers/queue/erdos-problem-993/paper/poisson_binomial/CERTIFICATE.md
- Lean file: /Users/brettreynolds/projects/LLM-CLI-projects/papers/queue/erdos-problem-993/formalization/pb_effective_drop_aristotle/PBReserve/Core.lean
- Project records: /Users/brettreynolds/projects/LLM-CLI-projects/papers/queue/erdos-problem-993/DECISIONS.md (l.950, l.962-967, l.983-987, l.997-1008)
- Sources checked for finding 3: /Users/brettreynolds/projects/LLM-CLI-projects/literature/meckes_2017_empirical_spectral_measure_unitary.md and /Users/brettreynolds/projects/LLM-CLI-projects/literature/hough_2006_determinantal_processes_independence.md
