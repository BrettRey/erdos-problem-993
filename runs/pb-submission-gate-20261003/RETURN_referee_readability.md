<!-- Saved verbatim by the parent session from the subagent's final message (Claude, responsibility-reader, full-paper referee readability read of paper/poisson_binomial/main.tex, 2026-10-03). Subagents cannot write report files. -->

REFEREE-STYLE READABILITY REPORT: "A variance-scaled Turán inequality at the first descent of a Poisson–binomial mass function" (ECP)

No responsibility notice is needed. Everything below is an ordinary review finding.

SCOPE
- I read all of /Users/brettreynolds/projects/LLM-CLI-projects/papers/queue/erdos-problem-993/paper/poisson_binomial/main.tex (784 lines) and the extracted PDF text (12 pages, references included). I also consulted /Users/brettreynolds/projects/LLM-CLI-projects/papers/queue/erdos-problem-993/paper/poisson_binomial/CERTIFICATE.md and /Users/brettreynolds/projects/LLM-CLI-projects/papers/queue/erdos-problem-993/notes/poisson_binomial_venue_decision_2026-07-16.md.
- Citations checked against the source: Hillion–Johnson (78)–(79) and Definition 3.11, in notes/literature/arxiv_1303_3381.txt. Expanding the paper's (2.2) and (2.3) gives exactly HJ's C_1(k) ≥ 0 and \bar C_1(k) ≥ 0, stated for all k in Z. The citation is correct. No other citation was checked (Pitman, Dümbgen–Wellner, Gnedin, Johnson, BCV, BMM Cor. 3.2, Aravinda, Marsiglietti–Melbourne).
- Spot-checked by hand, all matching:
  - Prop. 1.2 algebra; (1.7) derived from ULC(n); L_r = b_r; (4.3); Q(H); the closed forms (4.12)–(4.13).
  - P(3) = 22272, which equals the Table 1 minimum for [3,4].
  - At H = 4, A(δ) − Q = 6.8103 computed directly from the weights, equal to P(4)/(4·4^5·5^3) = 3486896/512000.
  - J=5 cell: β_0 = 2360, β_1 = 7500, β_4 = 31115/2.
  - J=6 (u=0): β_0 = 15557.5, β_1 = 39252.0625, β_4 = 99568, each equal to multiplier × polynomial / 2880.
- Not checked: the general-u identities and the other 274 compact-range coefficients.
- Points marked "heuristic" below are my own derivations for the author to test. They are not findings.

====================================================================
1. TERMINOLOGY
====================================================================

[MUST] 1.1 H appears in the Introduction before it is defined. L308–309: "A direct symbolic argument proves that inequality for H≥16; exact Bernstein expansions ... prove it for 3<H≤16". H is defined only at (4.1), L534.
Fix: say "in terms of H = 1/δ − 1", or give the ranges in δ (δ ≤ 1/17 and 1/17 ≤ δ < 1/4).

[MUST] 1.2 The Introduction cites (3.3) using a different letter. L297: "(3.3) gives M ≥ (1+12V)^{-1/2}". But (3.3) at L486 is written with p, and p is not defined until L420.
Fix: use one symbol for the maximal mass throughout, and not p (see 1.4).

[MUST] 1.3 "signed" is never defined. The Corollary label is cor:signed (L171), and the proof ends "so the same argument gives all three signed estimates" (L225).
Fix: "the three estimates for Z".

[MUST] 1.4 Symbol collisions a probabilist will trip over:
  - p = f_c = max_k f_k (L420), used in (2.8), (2.9), (3.2), (3.3) and Section 3. It collides with the success probabilities p_i (L56), the Bernoulli parameter p (L111) and p_n (L133). A probabilist reads p as a probability.
  - c = D−1, the rightmost mode (L84). It collides with c_⋆ (L114), "Johnson's c-log-concavity" (L270) and "a universal constant c>0" (L286).
  - M is the maximal mass at L290 but M_3, M_4 are power sums at L696.
  - S_0, T_0 (L687) read as the m=0 case of S_m, T_m (L599).
  - Minor: s means two things (L320, L433) and m means three (L218, L243, L597).
  Fix: write f_max (or M) for the maximal mass, and if M, rename M_3, M_4 to σ_3, σ_4. Write κ for the generic constant at L286. Rename S_0, T_0 to \tilde S_J, \tilde T_J.

[MUST] 1.5 ULC(n) is introduced (L87–88, "written ULC(n) and displayed in (1.7) below"), but (1.7) shows a consequence for δ_k, not the definition. The definition is HJ Def. 3.11: f_k / Bin_k(n,t) is log-concave, equivalently f_k/C(n,k) is log-concave.
Fix: give the definition in one clause at L87.

[MUST] 1.6 "Turán" is in the title, the abstract, L72 and seven other places, but is never explained. Probabilists say log-concavity.
Fix: one sentence at L72, e.g. "the inequality f_k^2 ≥ f_{k−1}f_{k+1} is often called a Turán inequality; δ_k is its slack relative to f_k^2", with a reference of your choice.

[MUST] 1.7 Project-internal names leak into the supplement text. L758–759 says the Lean project "formalizes the recurrence propagation, endpoint exclusion, first-crossing algebra, and raw-drop corollary". L769–770 has "bounded proof obligations". "First-crossing algebra" and "raw-drop corollary" occur nowhere else in the paper. They don't even match the manifest's own names (CERTIFICATE.md L5: "the first-crossing ratio bound, and the raw-from-effective corollary"). This is exactly what the earlier referee complained about.
Fix: name the formalized statements by paper labels: "(2.6) and the endpoint exclusion after it, the bound q_D ≥ a, and (1.5), assuming (2.4)".

[MUST] 1.8 The Gnedin paragraph (L264–268) uses four undefined terms: "mean–mode rule", "extended Bernoulli sums", "three-mass peak-skewness statistic", "leading mode".
Fix: cut the paragraph (recommended, for page budget), or define them in one clause inside a merged literature paragraph.

[MUST] 1.9 "c-log-concavity" at L270 is undefined and collides with c.
Fix: "Johnson [8, Lemma 5.1] gives", and drop the name.

[SHOULD] 1.10 "Normalized mass" means two different things.
  - L229 uses it for f_{D+r}/f_D.
  - L439–465 use "normalized-mass lower bound" for f_{c±r}/f_c.
  - The term appears in the abstract (L33–34) and the outline (L303–304) before either use.
Fix: define it at (2.8) ("we call f_{c±r}/f_c the normalized masses"). At L229 write "the ratio f_{D+r}/f_D".

[SHOULD] 1.11 "Cubic mass inequalities" (abstract L34, L92, L303) is the author's own label. HJ write "a 'cubic' inequality".
Fix: in the abstract and introduction, write "cubic inequalities of Hillion and Johnson for Bernoulli-sum masses".

[SHOULD] 1.12 Computational jargon that is never defined:
  - "scalar-cell checks", "triangular cell", "uniform cell family" (L669–670);
  - "compact scalar cells" (L662, L751, L756);
  - "full-degree Bernstein coefficients" (L610);
  - "certificate coefficients" (L625; "certificate" is never introduced);
  - "generator program" (L657, first mention, unexplained);
  - "standard-library checker" (L661; which library? Python is named only in the supplement);
  - "power-basis numerator" (L660);
  - "reference environment" (L749, L776);
  - "archive status" (L749);
  - the subsection title "The analytic range", where the argument is algebraic and nothing is analytic.
Fix: define "cell" at L597 as the interval [m,m+1]. Replace "triangular cell" with "the intervals between consecutive triangular numbers J(J+1)/2", "full-degree" with "degree-(2m+2)", and the subsection title with "The range H ≥ 16".

[SHOULD] 1.13 "Complementation" (L168, L187) is used in the Corollary statement before the proof explains it (L220).
Fix: add "(replace each Y_j by 1−Y_j)" in the statement.

[SHOULD] 1.14 L156: "Proposition 1.2 supplies a limiting admissible family but no extremality argument." "Admissible" is undefined.
Fix: "Proposition 1.2 gives an upper bound only; it does not identify an extremal law."

[SHOULD] 1.15 Other small definition gaps:
  - "first-descent index" in the abstract (L24) is undefined there. Add "(the least k with f_k < f_{k−1})".
  - "min-entropy" (L291): cut along with the Aravinda paragraph.
  - "Aristotle (Harmonic)" (L769): readers won't know it's an automated Lean prover. Say so.
  - Acronyms are otherwise fine. pmf is expanded twice (L22, L68), which is harmless.

====================================================================
2. INTRODUCTION
====================================================================

Overall: a probabilist knows what is proved by page 2. Theorem 1.1 is at L100–108, the threshold counterexample at L110–112, the bracket at (1.4) and Prop. 1.2 follows. Novelty is argued, but in two passes with repetition. How the proof goes is only partly conveyed. The outline (L302–310) lists steps, uses an undefined symbol and gives no idea.

[MUST] 2.1 There is no proof idea. Add 2–3 sentences along these lines (heuristic built from the paper's own bounds; check the wording):
  - If δ_D is small, (2.4) spreads that smallness: δ_{D±r} ≤ δ/(1−rδ).
  - So mass ratios near the mode are at least about Π(1−s/H) ≈ exp(−r^2/(2H)).
  - The pmf therefore has a plateau of width of order δ^{−1/2} at height comparable to f_c, which contributes at least a constant times f_c^2 δ^{−2} to V.
  - Since f_c^2 ≥ (1+12V)^{−1}, this forces V ≳ 1/δ.
This is the sentence that tells a referee why δ of order 1/V is the right order.

[MUST] 2.2 Nothing explains why δ is taken at D. L97 calls this "the natural local question" but never says why the first descent rather than the mode c, or what controlling δ_D is good for.
Fix:
  - One or two sentences of motivation. The venue note (L19) keeps Erdős 993 out of the main motivation, but a one-line origin note in the acks or a footnote would still tell a referee where the question came from.
  - A remark that the result carries over to the mode (my derivation; check it). The second inequality of (2.4) at k = D−1 gives δ_{D−1} ≥ δ_D/(1+δ_D). Hence Vδ_c ≥ (1/4)/(1+δ_D) ≥ 1/8, and Vδ_c = V ≥ 1 when D = 1.

[MUST] 2.3 The motivation at L96–98 doesn't follow: "The variance of the Bernoulli summand B_i is p_i(1−p_i), which tends to zero as p_i approaches 0 or 1. Thus the natural local question is whether Vδ_D is uniformly bounded away from zero. The answer is yes." Nothing in the first sentence supports "Thus". Then "The answer is yes" is contradicted at L110 ("A positive universal bound without a lower variance threshold is impossible").
Fix: "Summands with p_i near 0 or 1 increase n but add little to V, so bounds indexed by n, such as (1.7), degenerate (see the example after (1.7)). We ask whether Vδ_D is bounded below; it is, once V ≥ 1."

[SHOULD] 2.4 The literature review is split in two and repeats itself.
  - L86–94 covers ULC, Darroch, BMM and HJ, with three citations lumped at the end ("[3,4,7]") so the reader can't tell which supports what.
  - L255–288 covers Darroch again, Pitman, DW, Gnedin, Johnson, and BCV with BMM again.
  - The ULC point is made three times (L88–89, L240–241, L252–253).
Fix: one literature block after Corollary 1.3, in this order: the ULC(n) example; one sentence covering Darroch, Pitman and DW; Johnson; BCV and BMM; then the novelty sentence (L285–288). Cut L86–94 to a single forward pointer.

[SHOULD] 2.5 The Aravinda paragraph (L290–300) is never used. Its only job is to separate a result that runs in the opposite direction.
Fix: cut it, or reduce it to one clause in the BMM sentence.

[SHOULD] 2.6 Why V ≥ 1? A threshold is shown to be necessary, but not why it is 1.
Fix: say whether the method gives some constant for V ≥ v_0 with v_0 < 1, or say that the threshold is a normalization.

[SHOULD] 2.7 The constant discussion (L156–163) is defensive. A sentence on where 1/4 comes from would do more. Heuristics for the author to test:
  - By local-CLT reasoning, log f has second difference about −1/V near the mode, so Vδ_D → 1 as V → ∞.
  - In the proof's large-H regime, S ≈ √(2πH) and T ≈ H√(2πH), so ST/Q(H) → 8π/3. The method then gives about √(π/6) ≈ 0.72.
  - So the constant 1/4 is set by small H (δ near 1/4), i.e. small variance.
If these hold, the reader learns that the theorem's content is uniformity at small and moderate V. Small-n numerics on c_⋆ would help if you have them.

[SHOULD] 2.8 Section 1 runs from L51 to L310, about 36% of the body and pages 1–5 of the 11 text pages. That's because it carries the proofs of Prop. 1.2 and Corollary 1.3.
Fix: consider moving the Corollary proof and the X−Y extension to the end of Section 3.

[SHOULD] 2.9 The X−Y generalization (L183–201) is only a translation, yet it takes about 19 lines of statement.
Fix: make it a one-sentence remark after the Corollary.

[SHOULD] 2.10 L302 says "The proof uses two external results." But BMM is reproved at L488–502, and Newton's inequalities are also external.
Fix: "The proof combines the cubic inequalities of Hillion and Johnson with the maximal-mass bound (3.3)."

[SHOULD] 2.11 Say for which laws the result is new (heuristic; check). D is within about 2 of the mean, so (1.7) gives Vδ_D ≳ Vn/(μ(n−μ)). That is already a variance-scale bound when V is comparable to μ(n−μ)/n (nearly equal p_i). The theorem is new for heterogeneous p_i. One sentence would sharpen the novelty claim.

====================================================================
3. PROOF READABILITY
====================================================================

[MUST] 3.1 The H ≥ 16 range (L672–741) can't be reproduced.
  - The J=5 Bernstein coefficients (L706–710) and the five u-polynomials (L724–731) come from a symbolic computation that is neither shown nor packaged.
  - L756 says "The programs verify only the thirteen compact scalar cells", and CERTIFICATE.md L3 agrees.
  - A referee can't check polynomials of degree 4–6 in u by hand. My endpoint checks match, but nothing in the package verifies the general-u identities.
  Fix: (a) add a short SymPy script for H ≥ 16 to the supplement and say so in the text; and/or (b) write N_J out so the computation can be redone from the page:
      N_J(H) = −(3/4)H^4 − (7/4)H^3 + (αγ − 1)H^2 − (αΣ + βγ)H + βΣ,
      where α = 2J+1, β = J(J+1)(J+2)/3, γ = J(J+1)(2J+1)/3, Σ = M_3 + M_4.
      Checked at J=5: 1209, 20944 and 84280, matching my direct expansion.
  Also state the Bernstein relation as an equation. L732–736 currently reads "After factoring out the common strictly positive scalar 1/2880, the removed multipliers are, respectively, ...", which leaves the reader to guess whether 1/2880 multiplies or divides. Write "β_i(J) = μ_i(J) π_i(J−6)/2880, with μ_0 = J^2(J+1)^2, ..." (check: at J=6, 1764·25400/2880 = 15557.5 = N_6(21)).

[MUST] 3.2 The paper doesn't say where the supplement is. main.tex has no URL or DOI. CERTIFICATE.md L112 confirms no DOI exists yet, and the GitHub URL it gives isn't in the paper. Table 1's caption ("Full coefficients and digests are in the supplement") points nowhere.
Fix: put a Zenodo DOI (and the repository URL) in the supplement block and in the Section 4.1 text.

[MUST] 3.3 (2.6) and the endpoint exclusion are a single induction, presented backwards (L400–418).
  - (2.6) is stated "provided D±r ∈ {0,…,n}".
  - But deriving it at step r applies (2.4) at k = D±(r−1), which needs 1 ≤ D±(r−1) ≤ n−1. That is the endpoint exclusion at step r−1, which is stated only afterwards (L407–409).
  - "Iterating the map x ↦ x/(1−x)" leaves the induction, the monotonicity and the identity φ^r(x) = x/(1−rx) to the reader.
Fix: a lemma, "For 1 ≤ r ≤ K and each sign, D±r ∈ {1,…,n−1} and δ_{D±r} ≤ δ/(1−rδ) < 1", proved by induction in about four lines. Then "D>K and n−D>K" (L417) follows at once, and the separate sentence at L407–409 can go.

[MUST] 3.4 The range is inconsistent: L578 says "On 3<H<4" and L583 says "On 3<H≤4".
Fix: use ≤ in both.

[MUST] 3.5 L585: "Direct simplification gives A(δ)−Q(H)=P(H)/(4H^5(H+1)^3)." "Direct simplification" suggests a hand check, but this is a degree-10 identity.
Fix: "A computer-algebra computation (generator script, supplement) gives". My spot checks at H = 3 and H = 4 agree.

[SHOULD] 3.6 L355–357 cites HJ without the hypothesis.
Fix: add "for the pmf of any finite sum of independent Bernoulli variables".

[SHOULD] 3.7 L475–480 justifies (3.2) only with "The pairwise variance identity and the preceding normalized-mass lower bounds give".
Fix: add "restrict the double sum to indices in {c−K,…,c+K}, which lie in the support, and use f_{c+i} ≥ p w_i."

[SHOULD] 3.8 The BMM proof (L488–502) is an informal bathtub argument: "Moving mass toward x_0 without exceeding the density bound lowers the integral until this density is reached."
Fix: use the two-line rigorous version. Let ψ = p·1{|x−x_0| ≤ 1/(2p)}. Then (x−x_0)^2 − 1/(4p^2) and φ−ψ have the same sign pointwise, and ∫(φ−ψ) = 0. Alternatively, cite BMM and delete the proof (saves about 12 lines).

[SHOULD] 3.9 The proof of Theorem 1.1 from Prop. 3.1 (L513–527) never recalls that δ < 1/4, which it needs before invoking Prop. 3.1.
Fix: open with "By (2.5) we may assume 0<δ<1/4."

[SHOULD] 3.10 L421–426: "Since q_{D−1}≥1>q_D, the first inequality in (2.4) gives ..." Two problems:
  - The hypothesis is used in the second display, not the first, and "1>q_D" isn't used at all.
  - q_{D−1} needs D ≥ 2, which holds because D > K ≥ 3. Say so.

[SHOULD] 3.11 Three facts about K and H are used before they are stated:
  - K = max{r : r < H} first appears at L678, inside the H ≥ 16 subsection. Move it to (4.1).
  - (4.3) is ≥ 1 because H−r−1 > 0, which follows from r < K < H. This is unstated.
  - L550, "Since every auxiliary weight is nonnegative", needs b_r ≥ 0 for r ≤ K < H. Say so.

[SHOULD] 3.12 The argument that δ > 0 (L388–392) is roundabout. Newton's inequalities give strict log-concavity directly, so δ_k > 0 for every k.
Fix: replace it with that one clause.

[SHOULD] 3.13 Gaps in the Prop. 1.2 proof:
  - It says "p_n>1/n and hence q_1>1". Spell it out: q_1 = np_n/(1−p_n) > 1 ⟺ p_n > 1/(n+1).
  - It needs p_n ≤ 1/2 for the monotonicity step. That holds because p_n is the smaller root; say so.
  - The condition n ≥ 5 comes from n^2−4n−1 > 0. Say so.

[SHOULD] 3.14 L246: "The resulting pmf is symmetric and strictly log-concave." The claim D = m+1 needs one more clause: the pmf is symmetric about m and strictly log-concave, so f_m > f_{m+1}.

[SHOULD] 3.15 The Corollary statement contains proof text (L186–188): "Complementation followed by translation to a Bernoulli sum shows that this index is well defined." Move it to the proof.

[SHOULD] 3.16 The proof of (1.5) (L205–209) should say "multiply by V and apply Theorem 1.1", and note that D = n (where f_{D+1} = 0) is trivial.

====================================================================
4. PROSE
====================================================================

Sections 2–3 are mostly clean, compressed mathematical prose. The problems cluster in two places: Section 1's comparison paragraphs and the text about computation and provenance. A referee is more likely to say "over-compressed and repetitive" than "machine-written". The Aravinda paragraph and the verification paragraphs could still draw the second comment.

[SHOULD] 4.1 The same contrast is made twice in rotating parallel clauses:
  - L89–91: "Darroch's localization theorem controls modal-index location; variance-scaled maximal-mass bounds control the maximal mass; neither controls the normalized Turán deficit."
  - L255–257: "Existing quantitative results control modal-index location, adjacent mass ratios, or individual masses rather than the variance-scaled normalized Turán deficit at the first-descent index."
  - "rather than [the deficit]" recurs at L240, L256, L267 and L284.
Fix: say it once.

[SHOULD] 4.2 Two paragraph closers just restate what came before:
  - L240–241: "This bound is indexed by the number of summands rather than by the variance V of their sum." (Already said at L88–89.)
  - L252–253: "Even at fixed variance, the ULC(n) lower bound at the first-descent index tends to zero." (Restates the display above it.)
Fix: cut both. End the example with something new: "whereas Theorem 1.1 gives δ_D ≥ 1/4 for this family."

[SHOULD] 4.3 L228–230 states the tail bound a third time (after the abstract and L165–167): "Equation (1.6) gives a variance-scaled geometric upper bound for the normalized mass ...". Cut it and start the paragraph at "The standard ULC(n) estimate".

[SHOULD] 4.4 L296–300 is padding: "Its direction is opposite to the input needed here ... Aravinda's result complements that step and reinforces the distinction between variance control of one atom and variance control of a local normalized Turán deficit." Cut (see 2.5).

[SHOULD] 4.5 Several asides answer an auditor the reader never meets:
  - L309–310: "without floating-point sampling".
  - L366–368: "Their all-integer statement understands the pmf as zero beyond its support, so these inequalities include the cases adjacent to a support boundary."
  - L376: "show the normalization explicitly".
  - L385–386: "In particular, this includes relations adjacent to the boundary, such as δ_n(1−δ_{n−1})≤δ_{n−1}."
  - L417–418: "so every index used below lies in the support and its corresponding mass is strictly positive". Under the standing assumption every mass on {0,…,n} is positive.
  - L162–163: "we make no conjecture about which, if either, endpoint equals c_⋆".
Fix: keep one sentence at L366 ("HJ state (78)–(79) for all k in Z with the pmf extended by zero") and cut the rest.

[SHOULD] 4.6 Positivity is stated three times:
  - L594–595: "has eleven strictly positive rational coefficients";
  - L610–611: "all its full-degree Bernstein coefficients are strictly positive";
  - L621–625: "... have, respectively, 11 and 264 strictly positive coefficients. Thus all 275 certificate coefficients are strictly positive."
"Strictly positive" occurs 10 times (L37, 309, 418, 595, 624, 625, 718, 732, 738, 739).
Fix: state it once, after (4.9).

[SHOULD] 4.7 L657–658 is vague and padded: "The generator program reconstructs every identity before checking the signs. This reconstruction provides an internal consistency check."
Fix: one sentence, e.g. "A SymPy script computes each P_m and its Bernstein coefficients; an independent script using only exact rational arithmetic recomputes them from the definitions."

[SHOULD] 4.8 The opener of Section 4.2 (L669–670) is all jargon: "Thirteen exact scalar-cell checks cover the compact range; one initial triangular cell and a uniform cell family cover the analytic range."
Fix: "For H ≥ 16 we use the intervals between consecutive triangular numbers. The first (J=5) is checked directly; the rest (J≥6) are handled uniformly in J."

[SHOULD] 4.9 Noun stacks:
  - "the variance-scaled normalized Turán deficit at the first-descent index" (L256–257);
  - "left and right normalized-mass lower bounds about the rightmost modal index" (abstract, L33–34);
  - "exact asymmetric normalized-mass lower bounds" (L579);
  - "coefficient-vector SHA-256 digest" (L631–632);
  - "the Bernstein-vector digest for each cell, each cell-payload digest, and the overall payload digest" (L754–755);
  - "compact summary certificate" (L748).
"Normalized Turán deficit" occurs 10 times and "normalized-mass" 11 times.
Fix: unpack these where they're needed and move the digests to the manifest.

[SHOULD] 4.10 "Exact" is misused at L579 and L583 ("exact ... normalized-mass lower bounds"). These are lower bounds, not exact values; what's meant is "not symmetrized".
Fix: "the bounds L_r, R_r themselves, without replacing R_r by b_r".

[SHOULD] 4.11 Small items:
  - L303–304 has "yield ... yields".
  - "for either sign" (L401) should be "for each sign".
  - The last sentence of the abstract (L33–37) is a 50-word chain of noun phrases. Suggested rewrite: "The proof uses cubic inequalities of Hillion and Johnson to bound f_{c±r}/f_c from below near the mode c, turns these into a lower bound on V, and closes with a maximal-mass bound of Bobkov, Marsiglietti and Melbourne; the resulting one-variable inequality is proved symbolically for large arguments and by exact Bernstein expansions on a compact range."

====================================================================
5. COHERENCE
====================================================================

[SHOULD] 5.1 The abstract matches the content: main theorem, tail bound, 1/3 ceiling, proof route. It omits the X−Y extension, which is fine. Apply 1.15 and 4.11.

[SHOULD] 5.2 Most sections don't say what they do:
  - Section 2 (L312) goes straight into 2.1.
  - Section 3 opens "Define auxiliary weights by" (L470).
  - Section 4 opens "Set H=" (L532).
  - Section 4.2 opens with jargon.
Fix: one opening sentence for each:
  - Section 2: "We show that a small δ_D spreads to neighbouring indices and deduce lower bounds for f_{c±r}/f_c."
  - Section 3: "We turn these into a lower bound on V and reduce Theorem 1.1 to Proposition 3.1."
  - Section 4: "We symmetrize, then treat 3<H≤16 by exact computation and H≥16 by hand."

[SHOULD] 5.3 The title of Section 2, "The recurrence for normalized Turán deficits", doesn't name its main output, (2.8)–(2.9). Consider "Lower bounds near the mode". Section 2.1, "Standing reductions", is setup rather than recurrence; it fits better at the end of Section 1 or under a short "Preliminaries" heading.

[SHOULD] 5.4 Material that is out of place:
  - the Aravinda paragraph;
  - the SHA-256 column in Table 1;
  - the supplement block (L744–760), which repeats L657–665;
  - the replay results in the acknowledgements (L774–777);
  - the conditional Lean project, described twice (L757–760 and L770–772), each passage pointing at the other.

[SHOULD] 5.5 There is no closing remark. The open question about c_⋆ sits in Section 1 (L156–163). That's acceptable for ECP. If you adopt the heuristics in 2.7, a 3-line closing remark would be a better home for both.

====================================================================
6. VENUE FIT
====================================================================

[MUST in effect] 6.1 Page budget. The PDF is 12/12 pages. The venue record (L8, L38) sets 12 pages all-inclusive, with 11 as the target.
  - The additions above (motivation, a 3-sentence proof idea, the induction lemma, section openers, the supplement link) need roughly 20–25 lines.
  - Available cuts:
    - Aravinda paragraph (about 10 lines);
    - the duplicate literature pass at L86–94 (about 8);
    - Gnedin (about 5);
    - the restatements at L228–230, L240–241 and L252–253 (about 6);
    - the repeated positivity statements and L657–665 (about 10);
    - the duplication in the supplement block (about 6);
    - the BMM proof, if replaced by the citation (about 12);
    - Table 1, reduced to Cell/Degree/Count or to one sentence.
  - The \clearpage at L780 leaves roughly 60% of page 11 blank, and the references would probably fit there. \clearpage is house style; whether the journal's conventions win in the submission copy is your call.
  The cuts clearly outweigh the additions, and the paper should come out at 11 pages.

[SHOULD] 6.2 The tone is appropriate overall. The defensive asides (4.5) read as replies to earlier internal reviewers.

[SHOULD] 6.3 AI disclosure. The content is honest and proportionate. Page-1 placement is portfolio policy (2026-07-16), and I wouldn't change it. Three wording points:
  - "audit status" in the title footnote (L7–8) and the byte-for-byte replay sentence in the acknowledgements (L774–777) belong in the manifest.
  - "Aristotle (Harmonic) completed bounded proof obligations" would be clearer as "an automated prover produced Lean proofs of ...".
  - The venue record documents the AI policies of CPC and ALEA but nothing for ECP. Confirm ECP's current policy before submission; I didn't check it.

[MUST] 6.4 The computer-assisted part. Exact rational Bernstein certificates for explicit polynomial inequalities are acceptable for ECP in principle. As it stands, three things need fixing:
  - there's no persistent link (3.2);
  - the H ≥ 16 computation isn't packaged (3.1);
  - Table 1 has a SHA-256 column (lab-notebook material) and a minimum-coefficient column of up to 39 digits that tells the reader nothing. CERTIFICATE.md L76 says so itself: "Its magnitude depends on the chosen normalization; strict positivity is the invariant claim." Keep Cell/Degree/Count, or replace the table with one sentence.

[SHOULD] 6.5 The supplement description (L744–760) is too long and repeats Section 4.1. A referee needs four things: what's there, what it verifies, how to run it, and where it is. Four sentences will do.

[SHOULD] 6.6 Check that \bibliographystyle{amsplain} is the style the ejpecp class expects. I didn't check.

====================================================================
OVERALL VERDICT
====================================================================

Partly yes. Taking the earlier referee's three complaints one at a time:

- Undefined non-standard terminology: present, but at the edges rather than in the mathematical core. The worst cases are H in the outline, the misreferenced ULC(n), "signed", the Gnedin and Johnson names, the computational jargon of Section 4, and above all the internal Lean names in the supplement ("first-crossing algebra", "raw-drop corollary"), which don't even match the manifest.
- An Introduction that fails to say what the paper achieves: not present. Theorem 1.1, the necessity of a threshold and the 1/4–1/3 bracket are all on page 2. What's missing is a motivation for δ at D and a proof idea: the outline is a list of steps with an undefined symbol. Novelty is argued twice, at length.
- "Research journal" register: present in the provenance material (two certificates, generator and checker, reference environment, digests in the body, byte-for-byte replay in the acknowledgements, a conditional Lean project described twice) and in the defensive asides. The proofs in Sections 2–3 read like a paper.

One substantive gap lies outside the earlier complaints: the symbolic computation behind H ≥ 16 isn't in the supplement. Every endpoint I checked is correct, but a referee has nothing to reproduce. That, the missing supplement DOI, and the backwards (2.6) and endpoint-exclusion induction are the MUST-FIXes with real consequences. The rest is editing.

I'd call it one focused revision pass from submittable, probably a day's work. All the needed additions fit within the page budget once the duplicated and lab-notebook material is cut.
