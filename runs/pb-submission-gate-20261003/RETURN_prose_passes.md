<!-- Saved verbatim by the parent session from the subagent's final message (Claude Opus via responsibility-reader, read-only; five editorial passes on paper/poisson_binomial/main.tex, 730-line version, 2026-10-03). Subagents cannot write report files. After this report the parent applied 1.1-1.5, 2.2-2.14, 3.2-3.3, 4.1-4.7, 5.1-5.4 and 5.6, and added one-sentence openers to Sections 2 and 3 (5.7, first half); the Section 2 retitle was not made. See DECISIONS 2026-10-03. -->

EDITORIAL PASSES: paper/poisson_binomial/main.tex (ECP), read-only, 2026-10-03

Version: all line numbers refer to the 730-line /Users/brettreynolds/projects/LLM-CLI-projects/papers/queue/erdos-problem-993/paper/poisson_binomial/main.tex. main.log (L687, "11 pages") and main.pdf were built before the latest edits (the revision verifier saw 726 lines). Until you rebuild, the page count and the absence of overfull boxes are unverified.
No responsibility notice is needed.

What I read: main.tex in full; .claude/rules/writing-style.md; runs/pb-submission-gate-20261003/RETURN_referee_readability.md and RETURN_revision_verifier.md; ejpecp.cls (by grep); formalization/pb_effective_drop_aristotle/PBReserve/Core.lean L60–190; and the arXiv/preprint texts in notes/literature/ for the citation locators.

====================================================================
PASS 1. CONTRIBUTION ALIGNMENT
====================================================================

What the front matter claims, in my words:
- Main result: for any finite Bernoulli sum with variance V ≥ 1, the normalized log-concavity slack at the first-descent index D is at least 1/(4V).
- Consequences:
  - the pmf drops by a factor of at most 1−1/(4V) per step after D;
  - the same holds for differences of Bernoulli sums;
  - at the rightmost mode the slack is at least 1/(5V).
- Sharpness: no positive bound holds without a variance threshold, and a binomial family at V=1 caps the constant at 1/3.
- Proof route: Hillion–Johnson cubic inequalities, then deficit propagation, then mass lower bounds near the mode, then a variance lower bound, then the BMM maximal-mass bound. This leaves a one-variable inequality, proved in exact arithmetic.
- Machine checks: the programs cover only Proposition 3.1. A conditional Lean formalization covers parts of Lemma 2.1, the bound q_D ≥ a, and (1.5) from (1.3).

Does the body deliver it?
- Every front claim has a matching statement and proof:
  - Thm 1.1: Sections 2–4.
  - (1.4): Thm 1.1 plus Prop 1.2.
  - Cor 1.3: proof at L210–232.
  - Transfer to the mode: L234–238.
  - Threshold example: L118–120.
  - Degeneracy of the ULC and Johnson bounds: L251–277.
- Nothing proved in Sections 2–4 is missing from the front.
- Checked by hand and correct: the Prop 1.2 algebra, the Lemma 2.1 induction, q_D ≥ a, the R_r and L_r telescoping products, (4.3) in both forms, Q(H), the mode-transfer inequality, and P(3) = 22272.

1.1 [Medium] L36–38, abstract: "The resulting one-variable inequality is proved symbolically for large arguments and by exact Bernstein expansions on a compact range."
- The abstract never introduces H. The only variable the reader has is δ (or V), and the symbolic part covers δ ≤ 1/17, i.e. small δ. So "large arguments" names the wrong end of the range.
- Both ranges also use Bernstein expansions. Section 4.2 (L661–691) computes degree-four Bernstein coefficients: rational ones for J=5, polynomials in u for J≥6. So "symbolically" versus "Bernstein" is the wrong contrast.
- L296–298 repeats it: "symbolically for δ≤1/17, and by Bernstein expansions with positive rational coefficients for 1/17≤δ<1/4".
- Abstract replacement: "The resulting one-variable inequality is proved in exact arithmetic, by Bernstein expansions with rational coefficients for $1/17\leq\delta_D<1/4$ and uniformly in an integer parameter for $\delta_D\leq1/17$."
- L296–298 replacement: "Section 4 proves the one-variable inequality in exact arithmetic: by Bernstein expansions with positive rational coefficients for $1/17\leq\delta<1/4$, and for $\delta\leq1/17$ by Bernstein expansions whose coefficients are polynomials with positive coefficients in an integer parameter."

1.2 [Medium-low] L62–64: "Section~\ref{sec:reduction} checks that neither operation changes the quantities defined below."
- Section 2.1 (L311–313) says translation "shifts the first-descent index by $n_1$". D is one of the quantities defined below.
- Replacement: "Section~\ref{sec:reduction} checks that neither operation changes $V$ or the deficits defined below, and that translation shifts the first-descent index by the number of summands with $p_i=1$."

1.3 [Medium-low] L99–101: "Such a bound says that the pmf turns over at its peak at a rate set by the variance alone: Corollary~\ref{cor:tail} converts it into a relative drop of at least $1/(4V)$ from $f_D$ to $f_{D+1}$ and geometric decay beyond."
- The colon cites Corollary 1.3 as the reason the pmf "turns over at its peak". But the Corollary controls the step D→D+1, which is one step past the mode c=D−1.
- The step c→D is not controlled. The verifier's counterexample: Bin(2m+1,p) with p↑1/2.
- The statement about the peak follows from the remark at L234–238 (Vδ_c ≥ 1/5), not from the Corollary.
- Replacement: "By Corollary~\ref{cor:tail}, such a bound gives a relative drop of at least $1/(4V)$ from $f_D$ to $f_{D+1}$ and geometric decay beyond, and the remark after its proof gives $V\delta_c\geq1/5$ at the rightmost modal index itself." This also announces the remark (see 5.3).

1.4 [Low] L708–712, Lean scope: "a Lean formalization, conditional on the recurrence \eqref{eq:deficit-recurrence}, of the deficit bound and endpoint exclusion in Lemma~\ref{lem:propagation}, the bound $q_D\geq a$, and the deduction of \eqref{eq:adjacent-ratio-drop} from \eqref{eq:main-bound}."
- Core.lean confirms all three named items:
  - curvature_propagation and endpoint_exclusion, L70–125;
  - crossing_ratio_lower_bound, L128–136;
  - raw_quarter_of_effective, L177–186.
- The Lean hypothesis is narrower in form and wider in index range than (2.4). At L75 it reads `hstep : ∀ r, κ (r + 1) * (1 - κ r) ≤ κ r`: one inequality of (2.4), for an abstract nonnegative sequence, at every r ∈ ℕ. (2.4) holds only for 1≤k≤n−1.
- Two steps are not formalized: the passage from (2.4) to that hypothesis (extend the deficits by 1 beyond the support), and the conclusion D±r ∈ {1,…,n−1}.
- "Conditional on the recurrence (2.4)" is a fair summary, and the verifier also called it defensible. A referee who opens the files will still see the gap.
- Optional tightening: "conditional on the one-sided form of \eqref{eq:deficit-recurrence} for an abstract nonnegative sequence, of the deficit bound in Lemma~\ref{lem:propagation} and the exclusion of deficit one within its window, ..."

1.5 [Low] L88–90: "The known lower bounds on $\delta_k$ scale with $n$, ..."
- "The known" claims the list is exhaustive. L279–281 is hedged properly ("We are not aware").
- Replacement: "The lower bounds on $\delta_k$ that we know of depend on $n$, like ..., or on the success probabilities, like ...". "Depend on" also fits (1.7), which involves both k and n.

Checked and fine:
- "no universal constant larger than 1/3 is possible" (δ_2 → 1/3 from above at V=1);
- "for every support index D+r with r≥1" (the range of (1.6) is 1≤r≤n−D);
- the scope of the machine checks (L611, L693–695 and L706–708 correctly limit them to Prop 3.1);
- the application to real-rooted coefficient sequences (L101–106), which keeps the V≥1 hypothesis and is backed by Pitman's Proposition 1.

====================================================================
PASS 2. PROOFREAD
====================================================================

2.1 [Blocks submission; known deferral] L701: \url{ZENODO-DOI-PENDING}. Nothing else in the text says where the code is.

2.2 [Medium] L516–518 uses H four times ("on $3<H\leq4$ ... $H\geq4$ ... $3<H\leq16$ ... $H\geq16$"), but H is defined only at (4.1), L519–520. Put "Set H=... (4.1)" before the roadmap sentence. The sentence itself is item 4.1.

2.3 [Medium-low] L424: "We call the ratios $f_{c\pm r}/M$ normalized masses." The term never appears again (grep finds this line only). Delete the sentence.

2.4 [Medium-low] L354: "They extend the pmf by zero outside its support." The nearest plural antecedent is the inequalities (2.2)–(2.3), which cannot extend anything. Replacement: "Hillion and Johnson extend the pmf by zero outside its support."

2.5 [Low-medium] L661–662: "Writing $H=16+5t$ with $0\leq t\leq1$, the degree-four Bernstein coefficients of $N_5$ are" is a dangling participle. Replacement: "With $H=16+5t$ and $0\leq t\leq1$, the degree-four Bernstein coefficients of $N_5(16+5t)$ in $t$ are".

2.6 [Low-medium] L150–151: "As the smaller root of $x(1-x)=1/n$, $p_n$ lies in $(0,1/2)$, where $x\mapsto x(1-x)$ is increasing; so does $2/(n+1)$."
- "so does" must be read as "2/(n+1) also lies in (0,1/2)", but it can also attach to "is increasing".
- Replacement: "Both $p_n$, the smaller root of $x(1-x)=1/n$, and $2/(n+1)$ lie in $(0,1/2)$, where $x\mapsto x(1-x)$ is increasing."

2.7 [Low] L195–207, Corollary: "If $j^*=\max\{j:g_j>0\}$ and $\rho_Z=1-(4V_Z)^{-1}$, then [deficit and drop display]. We also have [tail display]."
- j* and ρ_Z are set up in a condition whose display does not use them.
- Replacement: "Then [deficit and drop display]. With $j^*=\max\{j:g_j>0\}$ and $\rho_Z=1-(4V_Z)^{-1}$, we also have [tail display]."

2.8 [Low] Two formulas separated only by a comma:
- L408 "Since $D\geq K+1\geq2$, $q_{D-1}$ is defined" → "Since $D\geq K+1\geq2$, the ratio $q_{D-1}$ is defined".
- L627 "For $J\geq5$, $J<J(J+1)/2\leq H$, so $J\leq K$." → "If $J\geq5$, then $J<J(J+1)/2\leq H$, so $J\leq K$."

2.9 [Low] Locator style. L345 has "\cite[(78)--(79)]", while the other locators read "equation~(1)" (L72), "equation~(41)" (L248) and "equations~(20)--(21)" (L264). Use "equations~(78)--(79)".

2.10 [Low] L297, "proves the one-variable inequality exactly", can be read as "proves precisely this inequality". Use "in exact arithmetic" (see 1.1).

2.11 [Low] L625: "Either choice is valid at a shared endpoint." → "Either choice of $J$ is valid at a shared endpoint."

2.12 [Low] L3: SHORTTITLE "at first descent" → "at the first descent", to match the title.

2.13 [Low] L570: "A computer-algebra computation (a SymPy program in the supplement) gives" → "A SymPy computation (in the supplement) gives".

2.14 [Low] L102–103 cites \cite{pitman1997} for the Bernoulli-sum representation without a locator. Add "[Proposition~1]". Pitman's Proposition 1, (i)⇔(iii), states exactly this (notes/literature/pitman_1997_coefficients_real_zeros.txt L62–70).

Checked and clean:
- Display punctuation: I checked every \] / \end{equation} / \end{align*} / \end{aligned} / \end{gathered} line. All have correct terminal punctuation; L67 and L359 end in commas where the sentence continues, and L226 continues without punctuation correctly.
- Every \ref and \eqref target exists.
- ejpecp.cls defines acks with an optional argument, supplement/\stitle/\sdescription, \orcid and \AMSSUBJSECONDARY.
- The class sets no bibliography style (grep finds no \bibliographystyle), so it does not contradict `amsplain`. ECP's author instructions are the place to confirm.
- Terminology and hyphenation are consistent (Poisson--binomial, first-descent index, maximal-mass bound, nonnegative/nonincreasing, degree-four/degree-ten/degree-$(2m+2)$).

Citation locators, checked against the texts on disk. These are arXiv or preprint versions (the Pitman text carries the preprint grant footer), so confirm the numbering against the published versions.
- Pitman, eq (1): Newton's inequality (pitman txt L107–123). Correct.
- Pitman, eqs (20)–(21): consecutive-ratio bounds (L316–365). Correct.
- Hillion–Johnson (arxiv_1303_3381.txt):
  - Definition 3.11 and (41) at L889–897. Correct.
  - Theorem A.2 = (78) at L1435–1439. Correct.
  - Corollary A.3 = (79) at L1444–1448, stated "for all k ∈ Z". Correct.
- Johnson (arxiv_1507_06268.txt): Definition 1.2 (L65–71), Lemma 5.1 (L530) and Example 5.2(2) (L549–552), with c = (Σ p_j/(1−p_j))^{−1}. The paper's display is the correct rearrangement of E(V)(k) ≥ c.
- BMM (arxiv_2007_11030.txt): Corollary 3.2, N_∞(X) ≤ 1+12Var(X) for any integer-valued X (L309–310), with N_∞ = M^{−2} at (2.1) (L187). Correct.
- Dümbgen–Wellner (arxiv_1910_03444.txt): Proposition 2, strict decrease of (x+1)b(x+1)/b(x) (L283–284). Correct.
- Marsiglietti–Melbourne (arxiv_2205_08293.txt): item 3 after Definition 1.3 defines ULC(n) as log-concavity relative to the binomial (L84–90). Correct.
- Not checked (not on disk): Darroch 1964, Baillon–Cominetti–Vaisman, Tang–Tang. Pitman's text (L293–294) supports the Darroch statement: "differs from the mean by less than 1".

====================================================================
PASS 3. DE-AI PROSE
====================================================================

Almost nothing to report. The mathematical prose is clean. I found:
- no praise-as-structure, faux-coaching, mechanism metaphors, sealed paragraph closers or unsupported evaluative participles;
- no filler intensifiers (grep for genuinely/really/truly/actually/crucial/essential/precisely/merely: none).
Agency verbs like "makes", "converts", "forces" and "transfers" are normal mathematical register.

Three low items:

3.1 [Low] L99: "Such a bound says that the pmf turns over at its peak at a rate set by the variance alone:" combines a metaphor with a colon set-up, and it is inaccurate. Replacement under 1.3.

3.2 [Low] L73–74: "its normalized slack is the normalized Tur\'an deficit" defines one new term by another carrying the same modifier. Replacement: "For $0\leq k\leq n$, we measure its slack relative to $f_k^2$ by the normalized Tur\'an deficit".

3.3 [Low] Noun stack at L468: "the lattice maximal-mass bound of Bobkov, Marsiglietti, and Melbourne" → "the maximal-mass bound of Bobkov, Marsiglietti, and Melbourne for integer-valued laws". See also 2.13.

The corrective negations at L165 and L706–708 are under pass 4.

====================================================================
PASS 4. EDITORIAL SCAR TISSUE
====================================================================

4.1 [Medium] L516–518, Section 4 opener: "We keep the asymmetric bounds on $3<H\leq4$ and symmetrize them for $H\geq4$; the range $3<H\leq16$ is handled by exact computation and $H\geq16$ symbolically."
- This was written to fix verifier item P2 ("We symmetrize the weights" was false on 3<H≤4).
- The exception now leads the section, and it is repeated twice:
  - L564–565: "On $3<H\leq4$, we instead keep the bounds $L_r$ and $R_r$ without symmetrizing."
  - L569–570: "On $3<H\leq4$ we have $K=3$, and we use $L_1,\dots,R_3$ themselves."
- Replacement: open with (4.1), then "We treat $3<H\leq16$ by exact computation and $H\geq16$ symbolically." State the exception once, at L569–570 where it is used, and cut the second sentence of L564–565.

4.2 [Medium] L165: "Proposition~\ref{prop:one-third} gives an upper bound only; it does not identify an extremal law."
- This is item 1.14 of the earlier readability report, copied word for word.
- No reader would think Prop 1.2 identifies an extremal law, and the paragraph's own last clause already says so: "we do not know which, if either, endpoint of \eqref{eq:constant-bracket} equals $\kappa_\star$".
- Delete it and start the paragraph at "The lower-bound proof combines inequalities...".
- In the same paragraph, put the parenthetical in proof order: "(the mass bounds of Section~\ref{sec:propagation}, the maximal-mass bound of Section~\ref{sec:variance}, and a symmetrization step in Section~\ref{sec:scalar})".

4.3 [Low-medium] L150–151, "...; so does $2/(n+1)$." This clause was bolted on to answer verifier nit P8. Fix as in 2.6.

4.4 [Low-medium] L261–262, "Darroch shows that a modal index lies within one unit of the mean \cite{darroch1964}." This survived the deletion of "In more detail," (verifier P1). Darroch is already cited for the same purpose at L94–96. Merge as in 5.2.

4.5 [Low] L613–615, "which uses only exact rational arithmetic from the Python standard library and does not import the generator". "Does not import the generator" answers an auditor's question, and "independently written" already covers it. Replacement: "An independently written checker, using only exact rational arithmetic from the Python standard library, recomputes them from the defining formulas."

4.6 [Low; a legitimate scope statement] L706–708, "they do not machine-verify the reduction of Theorem~\ref{thm:main} to it." Keep the content but say it positively: "the reduction of Theorem~\ref{thm:main} to Proposition~\ref{prop:scalar} (Sections~\ref{sec:propagation} and~\ref{sec:variance}) is proved in the text."

4.7 [Low] L702 "(with digests)" and L705 "a manifest gives the reference environment and replay commands" are leftovers of the lab-notebook register the earlier referee flagged (items 1.12, 4.9, 5.4). Drop "(with digests)", and write "a manifest gives the software versions and the commands to rerun the checks".

Checked; not scar tissue:
- L354 is the one boundary sentence the referee asked to keep, and it is needed for k=1 and k=n−1. Only its antecedent needs fixing (2.4).
- L529 "(so that $H-r-1>0$)" and L408 "$q_{D-1}$ is defined" are steps the proof needs.
- L105–106 "whenever the associated law has variance at least one" is a hypothesis the claim needs.
- L98 states the question once.

====================================================================
PASS 5. COHERENCE AND COHESION
====================================================================

The argument, section by section:
- Section 1 (L52–298): defines δ_k and D; Thm 1.1 (Vδ_D ≥ 1/4); the threshold is necessary; κ⋆ ∈ [1/4,1/3] via Prop 1.2; Cor 1.3 (drop, tail, X−Y); transfer to the mode; comparison with ULC(n), Johnson and the localization results; novelty; proof outline.
- Section 2.1 (L303–329): reduction to 0<p_i<1; D exists when V≥1.
- Section 2.2 (L331–451): HJ cubic inequalities give the recurrence (2.4); Lemma 2.1 propagates the deficit bound; lower bounds R_r and L_r on the masses.
- Section 3 (L453–511): V ≥ M²A(δ); the BMM bound with proof; Thm 1.1 reduced to Prop 3.1.
- Section 4 (L513–565): change of variable to H; R_r ≥ b_r; symmetrization A ≥ ST; target Q(H).
- Section 4.1 (L567–616): exact Bernstein certificates on 3<H≤16.
- Section 4.2 (L618–696): H≥16 over intervals between triangular numbers.
- Supplement and acknowledgements (L698–725).

Problems at the joins:

5.1 [Medium] L516–518: H is used before (4.1), and the 3<H≤4 exception is stated three times (2.2, 4.1).

5.2 [Medium-low] L261–266 repeats Darroch from L94–96.
- The Pitman / Dümbgen–Wellner sentence never says how those results bear on δ_D, so the paragraph sits as an unexplained list between the ULC comparison and the Johnson comparison.
- Fix: delete L261–266 and extend L94–96: "Darroch's theorem \cite{darroch1964}, Pitman's bounds on adjacent mass ratios through exponentially tilted laws \cite[equations~(20)--(21)]{pitman1997}, the strict decrease of $(k+1)f_{k+1}/f_k$ proved by D\"umbgen and Wellner \cite[Proposition~2]{duembgenwellner2020}, and the maximal-mass bounds of \cite{bailloncominettivaisman2016,bobkovmarsigliettimelbourne2022} locate the mode or bound ratios and masses, but do not control $\delta_k$."

5.3 [Medium-low] The mode-transfer remark at L234–238 appears without warning, between the Corollary's proof and the literature comparison.
- The paper never asks why δ is taken at D rather than at the mode, so the remark answers a question nobody has raised.
- Fix: announce it at L99–101 (replacement under 1.3), or open the remark with "Theorem~\ref{thm:main} is stated at $D$, one step past the mode; it also bounds the deficit at the rightmost modal index $c=D-1$."

5.4 [Low-medium] Unclear antecedents and "this" without a noun:
- L354 "They" (see 2.4).
- L703–704, "an independently written checker for it". Grammatically "it" is the generator program, but the checker checks the coefficients. → "an independently written checker that recomputes those coefficients using only the Python standard library".
- L231, "so the same argument gives the three estimates for $Z$". Which argument? → "so Theorem~\ref{thm:main} and the argument above give the three estimates for $Z$."
- L694–695, "Together with the range $3<H\leq16$, this proves Proposition~\ref{prop:scalar}." "This" has no noun, and a range does not prove anything. → "This settles $H\geq16$; with Section 4.1 it proves Proposition~\ref{prop:scalar}."

5.5 [Low] The definition at L424 is never used (see 2.3).

5.6 [Low] Repeated statements:
- "c=D−1 is the rightmost modal index" appears at L86, L234 and L406–407. Drop the "where" clause at L406–407.
- The programs are described four times: L570–571, L612–616, L693–695 and L701–706. The supplement block needs its description; L612–616 could shrink to one sentence.

5.7 [Low] Section headings:
- The title of Section 2, "The recurrence for normalized Tur\'an deficits", covers neither Section 2.1 (reductions) nor the main output of Section 2.2, the mass bounds (2.8)–(2.9).
- The title of Section 2.2, "The recurrence", repeats it. Consider "Deficit propagation and mass bounds near the mode".
- Sections 2 and 3 have no opening sentence; Section 3 starts "Define auxiliary weights by". One orienting sentence each would help; the earlier referee's item 5.2 has suggested wording.

Already known, not raised again:
- Symbol reuses the parent chose to keep (verifier header): d (L194 vs L600), m (L251 vs L584), a (L379 vs a_j at L600), and q_k defined at L145 and again at L335.
- Earlier referee SHOULD items still open. These are content questions, not editing:
  - why the threshold is V≥1. The verifier's derivation, that Vδ_D ≥ min{V,1}/4 whenever D exists, would answer it in one clause.
  - where 1/4 comes from.
  - for which laws the result is new.

====================================================================
OUTSIDE THE FIVE PASSES
====================================================================

- Acknowledgements, L717: "OpenAI's GPT-5.6 model". Project records use "GPT-5.6", "GPT-5.6 sol" and "GPT-5.6-sol Ultra" (project CLAUDE.md; DECISIONS.md L125, L281). Confirm which name the disclosure should use; the portfolio rule is to check model names rather than write them from memory.
- Rebuild before submission (see the version note at the top).
