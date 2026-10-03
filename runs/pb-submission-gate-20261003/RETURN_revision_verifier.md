<!-- Saved verbatim by the parent session from the subagent's final message (Claude Opus via responsibility-verifier, read-only, no shell; checked the revised paper/poisson_binomial/main.tex against RETURN_referee_readability.md, 2026-10-03). Subagents cannot write report files. After this report the parent applied E1–E3 and P1–P8 (except P7's d, m and a reuses) and fixed CERTIFICATE.md (P5); see DECISIONS 2026-10-03. -->

REVIEW OF THE REVISED POISSON–BINOMIAL PAPER AGAINST RETURN_referee_readability.md

No responsibility notice is needed. You need one coordination fact before anything else.

VERSION NOTE
- main.tex changed while I was reviewing it. My first read (724 lines) matched the PDF extract. By my final read it had 726 lines, with rewritten intro paragraphs (L88–105) and a rewritten literature block (L260–280). Everything else had only shifted by 2 lines.
- All line numbers below refer to the 726-line version: /Users/brettreynolds/projects/LLM-CLI-projects/papers/queue/erdos-problem-993/paper/poisson_binomial/main.tex
- The PDF extract (scratchpad/coldread2/fullpaper.txt) and main.log L687 ("Output written on main.pdf (11 pages ...)") both predate those edits. The page count needs a fresh build.
- I had no shell. Every arithmetic check below was done by hand. I ran no script and did not open the zip.

======================================================================
SUMMARY OF DEFECTS
======================================================================

Errors (wrong statements):

E1. L100–101 is false: "Corollary~\ref{cor:tail} converts it into a drop of at least $1/(4V)$ in the first step past the mode".
- L86 defines the mode as c = D−1. The first step past c is c→D, governed by q_D = f_D/f_c, and nothing in the paper keeps that away from 1. (1.5) controls the next step, D→D+1.
- Counterexample: Bin(2m+1,p) with m≥2 and p↑1/2. Then c=m, D=m+1, f_D/f_c = p/(1−p) → 1, and V → (2m+1)/4 ≥ 5/4, so V(1−f_D/f_c) → 0.
- Fix: "converts it into a drop of at least 1/(4V) from f_D to f_{D+1}, the step after the first descent, and geometric decay beyond."

E2 (minor). L91–93 overstates: "Both degenerate when summands with $p_i$ near $0$ or $1$ ... accumulate at fixed variance".
- If only near-0 summands accumulate (V fixed, D bounded), the right side of (1.7) at D tends to 1/(D+1) > 0, and Johnson's odds sum Σ p_i/(1−p_i) stays bounded. Neither bound degenerates.
- (1.7) degenerates when D and n−D both grow; Johnson's bound degenerates when near-1 summands accumulate. The supporting example (L250–258) has both kinds.
- Fix: "near 0 and near 1".

E3 (minor). L101–105 ("these conclusions also apply to the coefficient sequences of such polynomials, for example rows of the Stirling and Eulerian arrays") drops the hypothesis V ≥ 1.
- Fix: add "whenever the associated law has variance at least one".
- The examples themselves are grounded: Pitman's abstract names Stirling numbers of both kinds and Eulerian numbers (notes/literature/pitman_1997_coefficients_real_zeros.txt L17–24, L81).

Presentation issues that impair reading (fixes in brackets):

P1. L260 "In more detail, Darroch shows ..." now continues a paragraph (L94–96) about 165 lines earlier. [Delete "In more detail,", or write "As noted in the introduction,".]

P2. The Section 4 opener (L513–514) has two inconsistencies.
- "We symmetrize the weights", but L560–561 says the weights are not symmetrized on 3<H≤4.
- "H≥16 by hand", but L689–691 says a SymPy program checks it, and the intro (L296) says "symbolically".
[Suggested: "We keep the asymmetric bounds on 3<H≤4 and symmetrize for H≥4; 3<H≤16 is handled by exact computation and H≥16 symbolically."]

P3. L525–529: the claim that (4.3) is ≥ 1 needs H−r−1 > 0. This was referee SHOULD 3.11 and is still unstated. [Add "since r ≤ K−1 < H−1".]

P4. Supplement text, L704–708: in "...the deduction of (1.5) from (1.3), which takes the recurrence (2.4) as a hypothesis", the relative clause attaches grammatically to (1.3). The Lean scope is also slightly overstated (details under (h)). [Suggested: "The archive also contains a Lean formalization, conditional on the recurrence (2.4), of the deficit bound and endpoint exclusion in Lemma 2.1, the bound q_D ≥ a, and the deduction of (1.5) from (1.3)."]

P5. The archive manifest is stale: /Users/brettreynolds/projects/LLM-CLI-projects/papers/queue/erdos-problem-993/paper/poisson_binomial/CERTIFICATE.md L5.
- It says "Lemma 2.1 (the bound (2.6) and the endpoint exclusion)". In the current numbering (2.6) is the definition of K and a; the bound is (2.7).
- "discussed in the manuscript's acknowledgement" should now say "in the supplement description".

P6. Undefined at first use:
- "the generator program" (L566); it is only explained at L608–611.
- "certificate" (L698).
[Write "(a SymPy program in the supplement)" at L566; "Bernstein coefficients" for "certificate".]

P7. Residual symbol reuse, low priority:
- u: L194 (the Corollary's top of support) and L667 (J = u+6). This is the one worth fixing.
- d: L193 (first descent of Z) and L596 (degree).
- m: L250 (family size) and L580 (cell index).
- a: L378, against the coefficients a_j at L596.

P8. Nits, none affecting correctness:
- L496 "By (2.5) we may assume 0<δ<1/4" cites the display of the very inequality being assumed. ["By the reduction before (2.5)".]
- L393–395 doesn't say that the first inequality of (2.4) serves the minus sign and the second the plus sign.
- L149–156: the monotonicity step also needs 2/(n+1) ∈ (0,1/2), true for n ≥ 4 but unstated.

======================================================================
TASK 1: STATUS OF THE [MUST] ITEMS
======================================================================

1.1 FIXED. H is gone from the Introduction. L296–297: "symbolically for $\delta\leq1/17$, and by Bernstein expansions with positive rational coefficients for $1/17\leq\delta<1/4$." These ranges match H≥16 and 3<H≤16 exactly.

1.2 FIXED. One symbol is used throughout:
- L287 "the maximal mass $M$"
- L290 "$M^2\geq(1+12V)^{-1}$"
- L405 "Write $M=f_c=\max_k f_k$ for the maximal mass."
- (3.3), L470: "$V\geq\frac{M^{-2}-1}{12}$"
A grep found no p used for the maximal mass. The only bare p is Ber(p) at L118–119, which is correct.

1.3 FIXED. L178 "\label{cor:tail}". L229–230: "the same argument gives the three estimates for $Z$." "signed" occurs nowhere (grep).

1.4 FIXED for all four prescribed renames:
- L121 "$\kappa_\star$" and L278 "universal constant $\kappa>0$"
- L639–641 "$\sigma_3$, $\sigma_4$"
- L630–634 "$\tilde S_J$, $\tilde T_J$"
- maximal mass M
"Johnson's c-log-concavity" is gone. The minor reuses listed under P7 remain.

1.5 FIXED. L241–242: "ultra-log-concave of order $n$, written $\operatorname{ULC}(n)$: the sequence $f_k/\binom nk$ is log-concave." I re-derived (1.7) from this definition: δ_k ≥ (n+1)/((k+1)(n−k+1)). Correct.

1.6 FIXED, though no reference is given, which the report left optional.
- L72–74: "$f_k^2\geq f_{k-1}f_{k+1}$ for every $k$, an inequality of Tur\'an type. For $0\leq k\leq n$, its normalized slack is the normalized Tur\'an deficit"
- Abstract L29–30: "the normalized slack in the log-concavity (Tur\'an) inequality".

1.7 FIXED in the paper. L704–708: "a conditional Lean formalization of Lemma~\ref{lem:propagation}, the bound $q_D\geq a$, and the deduction of \eqref{eq:adjacent-ratio-drop} from \eqref{eq:main-bound}". L717–719: "Aristotle (Harmonic), an automated theorem prover, produced the Lean proofs ...". "first-crossing algebra", "raw-drop corollary" and "bounded proof obligations" are gone (grep). What remains is the wording in P4 and the stale CERTIFICATE.md in P5.

1.8 FIXED by deletion. "Gnedin", "mean–mode", "leading mode" and "peak-skewness" do not occur, and Gnedin is no longer in the bibliography.

1.9 FIXED. L91 "a bound of Johnson", L267 "Johnson's bound reads". The name "c-log-concavity" is gone.

2.1 FIXED. L282–297, "The proof runs as follows ...", includes the heuristic. It checks out; see (g).

2.2 PARTLY FIXED.
- The transfer to the mode is added (L233–237) and correct; see (b).
- The motivation sentence (L98–101) is present but contains error E1.
- The optional one-line origin note (Erdős 993) is absent.

2.3 FIXED. L88 and L98: "Ordinary log-concavity gives only $\delta_D\geq0$ ... We ask whether $V\delta_D$ is bounded below; the answer is yes once $V\geq1$." The unsupported "Thus" is gone. Error E2 sits in this new text.

3.1 FIXED in the text.
- N_J is written out at L648–656.
- The Bernstein relation is now an equation, L670–671: "are $\beta_i=\mu_i\pi_i(u)/2880$ for $i=0,\ldots,4$".
- L689–691: "A short SymPy program in the supplement checks every expansion in this subsection in exact arithmetic." L700–701 says the same in the supplement description.
- The script exists: /Users/brettreynolds/projects/LLM-CLI-projects/papers/queue/erdos-problem-993/scripts/verify_pb_large_h_range.py. It checks the closed forms, N_J, the J=5 coefficients, and each β_i identity in u.
- The script was untracked in git at session start. I could not confirm that poisson_binomial_certificate_supplement.zip has been rebuilt to include it; run `unzip -l` on it.

3.2 DEFERRED, as instructed. L697: "The archive, at \url{ZENODO-DOI-PENDING}, contains ...". The report also asked for the repository URL (https://github.com/BrettRey/erdos-problem-993, CERTIFICATE.md L122). It is absent from the paper, which may be deliberate.

3.3 FIXED.
- Lemma 2.1 (L382–401) states, for 0≤r≤K and each sign, "$D\pm r\in\{1,\ldots,n-1\}$ and $\delta_{D\pm r}\leq\frac{\delta}{1-r\delta}<1$", proved by one induction.
- L403: "By Lemma~\ref{lem:propagation}, $D-K\geq1$ and $D+K\leq n-1$".
- The proof is correct; see (a).

3.4 FIXED. L560–561 "On $3<H\leq4$, we instead keep the bounds ..." and L565 "On $3<H\leq4$ we have $K=3$".

3.5 FIXED. L566–567: "A computer-algebra computation (the generator program in the supplement) gives".

6.1 FIXED, pending a rebuild. The PDF extract runs to "Page 11/11", main.log L687 says 11 pages, and the \clearpage before the bibliography is gone (L723). Both predate the L88–105 and L260–280 edits; I estimate the net change at under five lines, but rebuild and check.

6.4 Partly fixed, partly deferred.
- Table 1, with its SHA-256 and 39-digit minimum columns, is removed.
- The H≥16 computation is packaged (see 3.1).
- The persistent link is deferred (see 3.2).

[SHOULD] items addressed
- 1.10, 1.11, 1.13, 1.14, 1.15.
- 1.12, mostly: "cell" is defined at L580, "degree-$(2m+2)$" is used, the subsection is titled "The range $H\geq16$", and "power-basis numerator" and "archive status" are gone.
- 2.4, mostly. Darroch now appears twice (L94 and L260).
- 2.5, 2.10.
- 3.6 (L342), 3.7 (L459–460), 3.8 (L475–482), 3.9 (L496), 3.10 (L406–407), 3.12 (L369), 3.13 (L148–156), 3.14 (L252–253), 3.15 (L227).
- 3.11, partly: K is moved to L518, and b_r ≥ 0 is stated at L523. H−r−1 > 0 is not (P3).
- 4.2–4.4, 4.6–4.8, 4.10, 4.11; most of 4.5.
- 5.2 for Section 4 only.
- 5.4: table, replay sentence and duplicate Lean passage removed.
- 6.3 wording; 6.5.

[SHOULD] items not addressed that I consider worth doing
- 2.6 (why V ≥ 1). My own derivation, for the author to test: the branch δ_D < 1/4 never uses V ≥ 1. Lemma 2.1, (3.2), (3.3) and Proposition 3.1 need only that D exists and δ < 1/4. So whenever D exists, Vδ_D ≥ min{V,1}/4. In particular V ≤ 1 forces δ_D ≥ 1/4, and the threshold is a normalization. One clause after Theorem 1.1 would answer the referee's question.
- 2.7 (where 1/4 comes from) and 2.11 (for which laws the result is new): not addressed.
- 5.2: Sections 2 and 3 still have no opening sentence. Section 3 begins "Define auxiliary weights by".
- 5.3: the Section 2 title is unchanged.
- 2.9: the X−Y extension still takes about 17 lines of the Corollary statement.

======================================================================
TASK 2: MATHEMATICS OF THE NEW OR REWRITTEN PASSAGES
======================================================================

(a) Lemma 2.1 and its proof (L374–401): NO ERROR.
- Base case. δ_D = δ < 1. D ≥ 1 by (1.2), and D ≠ n because δ_n = 1 > δ.
- Inductive step.
  - The hypothesis D±r ∈ {1,…,n−1} is exactly what (2.4) needs at k = D±r: the first inequality for the minus sign, the second for the plus sign.
  - Since δ_{D±r} < 1, it gives δ_{D±(r+1)} ≤ δ_{D±r}/(1−δ_{D±r}).
  - x↦x/(1−x) is increasing on [0,1), and φ(δ/(1−rδ)) = δ/(1−(r+1)δ).
  - δ/(1−(r+1)δ) < 1 ⟺ (r+2)δ < 1, and r < K gives (r+2)δ ≤ (K+1)δ < 1 by the definition of K.
- Endpoint. D±(r+1) ∈ {0,…,n}, and δ_0 = δ_n = 1 excludes 0 and n.
- K ≥ 3 follows from 4δ < 1, and K is finite because δ > 0.
- Range check for every later use: 0≤r≤K is right, and it is tight on the left.
  - Left product (L429–442), 1≤j≤K: uses δ_{D−s} for s ≤ K (so r = K is needed for the minus sign), and q_{D−j} with D−j ≥ D−K ≥ 1.
  - Right product (L412–421), 0≤j≤K−1: uses deficits only up to δ_{D+K−2} and q up to q_{D+K−1}, with D+K−1 ≤ n−2.
  - Masses: c+K = D+K−1 ≤ n−2 and c−K = D−K−1 ≥ 0, so f_{c±r} is in the support for r ≤ K.
  - L406–409: D ≥ K+1 ≥ 2 makes q_{D−1} defined; q_{D−1} ≥ 1 by the minimality of D; (2.1) at k = D−1 is valid; δ_{D−1} ≤ δ/(1−δ) gives q_D ≥ (1−2δ)/(1−δ) = a.
  - The telescoping products, R_r and L_r, also check.

(b) Mode transfer (L233–237): NO ERROR.
- D = 1: δ_c = δ_0 = 1, so Vδ_c = V ≥ 1 > 1/5. The text leaves this last step implicit.
- D ≥ 2: the second inequality of (2.4) at k = D−1 (valid, since 1 ≤ D−1 ≤ n−1) reads δ_D(1−δ_{D−1}) ≤ δ_{D−1}, i.e. δ_{D−1} ≥ δ_D/(1+δ_D).
- δ_D < 1/4: Vδ_c ≥ (Vδ_D)/(1+δ_D) ≥ (1/4)/(5/4) = 1/5.
- δ_D ≥ 1/4: x/(1+x) is increasing, so δ_D/(1+δ_D) ≥ 1/5, and V ≥ 1.
- c = D−1 is indeed the rightmost modal index: the ratios q_k are ≥ 1 for k < D and ≤ q_D < 1 for k ≥ D.

(c) Proof of the maximal-mass bound (L472–483): NO ERROR.
- W+U has the piecewise-constant density φ = f_k on (k−1/2, k+1/2], so φ ≤ M, and Var(W+U) = V + 1/12.
- ψ integrates to M·(1/M) = 1.
- Sign claim. On |x−x_0| ≤ 1/(2M), both φ−M ≤ 0 and (x−x_0)^2 − (2M)^{−2} ≤ 0. Outside, φ ≥ 0 and the bracket is > 0.
- Integrating (x−x_0)^2(φ−ψ) = [(x−x_0)^2 − (2M)^{−2}](φ−ψ) + (2M)^{−2}(φ−ψ) gives ≥ 0, using ∫(φ−ψ) = 0. The text leaves that last step implicit, which is acceptable.
- ∫(x−x_0)^2 ψ = M·(2/3)(2M)^{−3} = 1/(12M^2). Hence V ≥ (M^{−2}−1)/12.
- The deduction of Theorem 1.1 from Proposition 3.1 (L494–508) is correct.

(d) Proposition 1.2 (L128–162): NO ERROR, apart from nit P8.
- p_n(1−p_n) = 1/n gives p_n = 1/(n(1−p_n)) > 1/n, so np_n > 1 and q_1 = np_n/(1−p_n) > 1.
- p_n is the smaller root, so p_n ∈ (0,1/2).
- 2(n−1)/(n+1)^2 > 1/n ⟺ 2n^2−2n > n^2+2n+1 ⟺ n^2−4n−1 > 0, which holds for n ≥ 5 (roots 2±√5) and fails at n = 4.
- q_2 = ((n−1)/2)·p/(1−p) < 1 ⟺ (n+1)p < 2.
- δ_2 = 1 − q_3/q_2 = 1 − 2(n−2)/(3(n−1)) = (n+1)/(3(n−1)).

(e) "By (1.7), δ>0" (L369): NO ERROR. (1.7) is displayed at L243–246, before its use, and at k = D it gives δ_D ≥ (n+1)/((D+1)(n−D+1)) > 0. The logical order is sound, though (1.7) sits inside a comparison paragraph as a cited fact.

(f) Section 4: NO ERROR found.
- With δ = 1/(H+1), (r+1)δ < 1 ⟺ r < H, so K = max{r ≥ 1: r < H}.
- L_r = (1−δ)^{−r} Π_{j=2}^{r+1}(1−jδ) = H^{−r} Π_{s=1}^{r}(H−s) = b_r.
- C_1 = 1, and C_{r+1}/C_r = (1−2δ)(1−rδ)/(1−(r+2)δ) = (H−1)(H+1−r)/((H+1)(H−r−1)) = 1 + 2r/((H+1)(H−r−1)).
- ∂A/∂w_ℓ is as stated, A_sym = ST, and Q(H) = (3H+4)(H+1)/4.
- Compact range:
  - Degree of P_m is 2m+2, with integer coefficients.
  - Coefficient count Σ_{m=4}^{15}(2m+3) = 264, and 11 + 264 = 275.
  - From the weights, P(4) = 3486896 and A−Q = 6.81034 at H = 4, both matching the stated identity.
- Bounds for H ≥ 16:
  - (4.11): the product bound needs s/H ≤ 1, and λ_r ≥ 0 needs J(J+1)/2 ≤ H; both hold.
  - The closed forms (4.12) and (4.13), σ_3, and σ_4 are correct.
  - \tilde S_J, \tilde T_J ≥ 0, so ST ≥ \tilde S_J \tilde T_J.
  - The intervals cover [15, ∞).
- The quartic. I expanded it myself. With α = 2J+1, β = J(J+1)(J+2)/3, γ = J(J+1)(2J+1)/3 and Σ = σ_3+σ_4: H^2·\tilde S_J\tilde T_J = αγH^2 − (αΣ+βγ)H + βΣ, and H^2·Q = (3/4)H^4 + (7/4)H^3 + H^2. This matches L651–655 term by term.
- J = 5 cell:
  - Coefficients 1209, 20944, 84280.
  - I recomputed all five Bernstein coefficients on H = 16+5t: 2360, 7500, 25055/2, 254205/16, 31115/2. All match L660–661.
- J = 6, t = 0: β_0 = N_6(21) = 15557.5 = 1764·25400/2880. This also equals N_5(21), as it must, since λ_6 = 0 at H = 21.
- All five β_i at J = 6 and at J = 7 match μ_iπ_i/2880. The J = 6 values are 15557.5, 39252.0625, 65153.67, 87710 and 99568.
- I then verified all five relations β_i = μ_iπ_i(u)/2880 exactly, as polynomial identities in J, by hand expansion. For example, 2880·N_J(J(J+1)/2)/(J^2(J+1)^2) = 121J^4 − 430J^3 − 965J^2 − 510J − 736 = π_0(J−6).
- The cross-cell identity π_0(u+1) = π_4(u) holds exactly. It is forced, since β_4(J) = β_0(J+1) at the shared endpoint.
- The leading coefficient 121/2880 matches the large-J asymptotics.
- Since these are hand computations, treat `python3 scripts/verify_pb_large_h_range.py` printing `ALL CHECKS PASSED` as the authoritative replay.

(g) "The proof runs as follows" (L282–297): NOTHING FALSE OR MISLEADING.
- δ_{D±r} ≤ δ/(1−rδ) holds for r ≤ K ≈ 1/δ.
- "explicit multiples of M" refers to R_r and L_r.
- V ≥ M^2A and M^2 ≥ (1+12V)^{−1} together reduce the theorem to A ≥ (3+δ)/(4δ^2): V(1+12V) ≥ A, with equality at V = 1/(4δ).
- The heuristic is sound:
  - b_r ≈ exp(−r^2δ/2), so the lower bounds are comparable to M for r ≲ δ^{−1/2};
  - S ~ δ^{−1/2} and T ~ δ^{−3/2}, so A ≳ δ^{−2};
  - 12V^2 ≳ A then gives V ≳ 1/δ.
- One nuance: the "symbolic" range δ ≤ 1/17 also uses Bernstein expansions, with symbolic coefficients. This is not false.
- The errors in the rewritten introduction are elsewhere: E1 at L100–101, E2 at L92, E3 at L101–105.

(h) Symbols and undefined terms.
- Consistent and defined at or before first use: M (L287, L405), κ_⋆ (L121), κ (L278), n_1 (L307), n_Y (L223), σ_3 and σ_4 (L639–641), \tilde S_J and \tilde T_J (L630–634), and cor:tail (L178, referenced at L93 and L100). No stray p for the maximal mass. "Turán", "first-descent index", "normalized masses" (L422), "cell" (L580) and "complementation" (L174–175) are all defined.
- Defects: the undefined terms in P6 and the symbol reuse in P7.
- "Ultra-log-concavity" is used at L89 before its definition at L241, but a forward pointer is present. Acceptable.
- Lean scope (/Users/brettreynolds/projects/LLM-CLI-projects/papers/queue/erdos-problem-993/formalization/pb_effective_drop_aristotle/PBReserve/Core.lean L70–125):
  - `curvature_propagation` assumes an abstract sequence κ with κ_0 = d, κ ≥ 0, and κ_{r+1}(1−κ_r) ≤ κ_r for every r ∈ ℕ. It proves κ_r ≤ the bound when (r+1)d < 1.
  - `endpoint_exclusion` rules out κ_r = 1.
  - This is the kernel of Lemma 2.1. The index-membership conclusion, and the passage from (2.4) (valid only for 1≤k≤n−1) to an all-r hypothesis (extend by 1 beyond the support), are modelling steps the Lean project does not formalize.
  - `crossing_ratio_lower_bound` (L128–136) and `raw_quarter_of_effective` (L177–186) match "q_D ≥ a" and "(1.5) from (1.3)".
  - "Conditional" is defensible; P4 tightens the wording.

======================================================================
WHAT IS VERIFIED, WHAT IS INFERRED, WHAT WOULD REVERSE THIS
======================================================================

Directly verified by me, by hand, from the files:
- every proof step in (a)–(e);
- the Section 4 algebra, including the five β_i identities as polynomial identities and the J = 5 and J = 6 coefficients;
- P(4), and A−Q at H = 4 computed from the weights;
- the counterexample to E1, and the Poisson-type case behind E2;
- the [MUST] statuses, by quotation and grep;
- that the H≥16 script exists, and its contents;
- the Lean theorem statements;
- the stale numbering in CERTIFICATE.md.

Inferred, not verified:
- that the current zip contains the H≥16 script and a matching CERTIFICATE.md;
- that the current source still builds to 11 pages;
- that the Lean project builds;
- the external citations (HJ, Johnson, BMM Corollary 3.2, Darroch), except Pitman's abstract;
- the 275 compact-range coefficients, of which I checked only P(4).

Would reverse my conclusions:
- `scripts/verify_pb_large_h_range.py` failing on replay. That would point to an arithmetic error that my hand expansion and the script share; unlikely, but it is the authoritative check.
- A fresh build of the 726-line source exceeding 11 pages (reverses 6.1).
- A further edit to main.tex after my last read, which would invalidate the line numbers above.

Files consulted: main.tex (726-line version); CERTIFICATE.md (L3–5, L13, L122); scripts/verify_pb_large_h_range.py; formalization/pb_effective_drop_aristotle/PBReserve/Core.lean; notes/literature/pitman_1997_coefficients_real_zeros.txt; main.log L687; runs/pb-submission-gate-20261003/RETURN_referee_readability.md; scratchpad/main_before_s4.tex; scratchpad/coldread2/fullpaper.txt (stale relative to the current source).
