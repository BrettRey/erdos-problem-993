<!-- Saved verbatim by the parent session from the subagent's final message (Claude Opus via responsibility-reader, read-only; four prose passes over the text rewritten in the second revision, 2026-10-03). Subagents cannot write report files. The parent then applied 1.1-1.9, 2.1-2.16, 2.18, 2.20-2.35 (2.17 and 2.19 left), 3.1-3.7, J1-J5 and R1-R5. Item 1.9 was moved to Section 2.1 as suggested; 1.8b was taken (Pitman paragraph now precedes Example 1.6). -->

PROSE REVIEW: four passes over the rewritten parts of the Poisson–binomial paper (read-only)

Manuscript: /Users/brettreynolds/projects/LLM-CLI-projects/papers/queue/erdos-problem-993/paper/poisson_binomial/variance-scaled-turan-first-descent.tex
Spec: /Users/brettreynolds/projects/LLM-CLI-projects/.claude/rules/writing-style.md (contractions and house macros waived, as instructed)

No RESPONSIBILITY NOTICE. Nothing found is a material risk. The worst items are a wrongly described split in the proof outline, one missing hypothesis (m≥2), and symbol clashes.

**What I checked and how**
- **Arithmetic, by hand only.** I could not run code. I checked:
  - Prop 1.2's algebra;
  - Prop 1.4's proof, including the e^{-|t|} tilt bound, the integral, and the factorization (1-1/q_c)(1-q_{c+1})≥0;
  - the derivation of (1.7);
  - Prop 1.5: the deficit (2m+1)/(m(m+2)), the variance limit 4, and δ_m ≥ 3/4;
  - Example 1.6: the leading terms of V_N, V_{10^4}≈1.1632, 1/(4V)≈0.2149, and the (1.7) value ≈0.0004. The asymptotic formula puts V_1996≈0.99997 and V_1998≈1.00007, which is consistent with "first reaches 1 at N=1998";
  - the Pitman derivation;
  - D≥2 in Section 2.1;
  - the derivation of (2.5);
  - Lemma 2.1.
- **Build state.** The .log contains only the hyperref `pagebackref` warning. PDF pp. 1–6 match the current .tex.
- **Sources.** Darroch was checked only through Pitman's paraphrase (/Users/brettreynolds/projects/LLM-CLI-projects/papers/queue/erdos-problem-993/notes/literature/pitman_1997_coefficients_real_zeros.txt, l.293–294: "differs from the mean μ by less than 1"). Pitman (21) was read in the same file, l.364–365. I did not re-verify the cited loci in Meckes–Meckes, HKPV, Hillion–Johnson or Johnson.
- **Provenance.** Several phrases and symbols trace to /Users/brettreynolds/projects/LLM-CLI-projects/papers/queue/erdos-problem-993/runs/pb-revision2-20261003/INPUT_external_review.md (cited below as INPUT), and to DECISIONS.md l.940–941.

Line numbers are .tex lines.

---

## PASS 1: Contribution alignment

All the abstract's claims match proved statements: Thm 1.1, Prop 1.2, (1.5), Prop 1.4, Prop 1.5 and Example 1.6. Nothing is overstated, and the abstract follows the same order as Section 1. The findings are framing, accuracy and order.

**1.1 L25–26 (abstract)**
- Text: "measures how strictly log-concave $f$ is at $k$; for a normal density the corresponding curvature is $1/V$."
- Problem: the abstract never says δ_k is a transform of a curvature, so "the corresponding curvature" has no antecedent. A reader may take it that the deficit of a normal is 1/V, when it is 1-e^{-1/V}. INPUT l.146 warned about exactly this identification.
- Replace with: "is a bounded transform of the discrete curvature of $\log f$ at $k$, and for a normal density of variance $V$ that curvature is $1/V$."

**1.2 L28 and L33–36 (abstract)**
- Text: "Cubic inequalities of Hillion and Johnson then give…" / "The proof reduces the main inequality, through a maximal-mass bound of Bobkov, Marsiglietti, and Melbourne, to a one-variable inequality…"
- Problem: the proof-route sentence leaves out the cubic inequalities, which are the first step of the main proof. With "then", the abstract presents them as an add-on after the theorem.
- Replace the last sentence with: "The proof combines these cubic inequalities with a maximal-mass bound of Bobkov, Marsiglietti, and Melbourne to reduce the main inequality to a one-variable inequality, which is proved in exact arithmetic by Bernstein expansions."

**1.3 L31–33 (abstract)**
- Text: "the bound is of order $1/\log N$, against $1/N$ from the number of summands."
- Problem: "the bound" and "from the number of summands" are both compressed past clarity.
- Replace with: "the bound at $D$ is of order $1/\log N$, against order $1/N$ from ultra-log-concavity."

**1.4 L346–348 (outline) against L718–724**
- Text: "with positive rational coefficients for $1/17\leq\delta<1/4$, and for $\delta\leq1/17$ with coefficients that are polynomials…"
- Problem: this is inaccurate. The J=5 cell (16≤H≤21, so 1/22≤δ≤1/17) is certified with explicit rational coefficients (2360, 7500, …). The rational/polynomial split is at δ=1/22 (H=21). The sentence also says "coefficients" twice.
- Replace with: "with positive rational coefficients for $1/22\leq\delta<1/4$, and for $\delta\leq1/22$ with coefficients that are polynomials in an integer parameter, all of whose coefficients are positive."

**1.5 L243–244**
- Text: "Poisson--binomial pmfs are ultra-log-concave of order $n$, written ULC(n): the sequence $f_k/\binom nk$ is log-concave."
- Problem: ULC(n) is defined only through the Poisson–binomial case, but Prop 1.5 applies ULC(2m) to g, which is not a Bernoulli sum.
- Replace with: "A pmf $(f_k)_{k=0}^n$ is ultra-log-concave of order $n$, written ULC(n), if $f_k/\binom nk$ is log-concave in $k$. Poisson–binomial pmfs are ULC(n), and this gives"

**1.6 L89–91 (the motivating question)**
- Text: "We ask how much of this survives, without approximation, for every Bernoulli sum: does the variance alone force a lower bound on $\delta_k$ near the mode?"
- Problem: the question asks for "a lower bound" without naming the 1/V scale that the paragraph has just set up. Theorem 1.1 answers the sharper question. The contrast between local shape and global spread is also left implicit: the words "local", "global" and "three adjacent masses" never appear.
- Replace with: "Since $\delta_k$ depends on three adjacent masses and $V$ on the whole law, we ask whether this scale holds without approximation for every Bernoulli sum: does the variance alone force $\delta_k\geq\kappa/V$ near the mode for a universal $\kappa>0$?" Also use one term throughout: L89 says "near the centre" and L91 says "near the mode". Prefer "near the mode".

**1.7 L106–107**
- Text: "…locate the mode or bound ratios and masses, but do not control $\delta_k$."
- Problem: the proof itself uses Darroch (L199, L230, L383) and BMM (Section 3) to control δ_D, so the unqualified negative is too strong.
- Replace with: "…locate the mode, or bound mass ratios and maximal masses, but do not by themselves control $\delta_k$." This also removes the garden path in "bound ratios".

**1.8 Order of the comparisons (L174–178, L275–326)**
- (a) L174–178 names "a symmetrization step in Section 4" before the reader has the proof, and the outline at L332–348 never mentions symmetrization. Fix in 3.6.
- (b) Example 1.6 sits between the comparison with Johnson's bound (L284–292) and the comparison with Pitman's (L317–326), splitting the run of comparisons with prior bounds. Move the Pitman paragraph to directly after L292. The order then runs: comparisons, Example 1.6, priority claim, outline.

**1.9 L207–210 (differences extension)**
- Text: "Complementation (replacing each summand $Y_j$ of a Bernoulli sum by $1-Y_j$) and translation preserve the variance and the deficits, so Theorem 1.1 and Corollary 1.3 extend to differences of independent Bernoulli sums, with $D$ the first descent of their pmf."
- Problem: this is the only result in Section 1 stated loosely. What actually holds is that W-W'+m is a Bernoulli sum. The abstract does not mention it, and it interrupts the move from Cor 1.3 to Prop 1.4.
- Replace with this, moved to the end of Section 2.1 or set as a remark after Prop 1.4: "If $W'=\sum_{j=1}^mB'_j$ is a Bernoulli sum independent of $W$, then $W-W'+m=W+\sum_j(1-B'_j)$ is a Bernoulli sum with the same variance as $W-W'$. Since translation preserves the deficits, Theorem 1.1 and Corollary 1.3 apply to $W-W'$, with $D$ the first descent of its pmf."

**1.10 (1.6) and the Pitman paragraph**
- Problem: (1.6) is the only result in Section 1 that has no motivation where it is stated. Its only purpose comes eleven lines later, in the Pitman comparison. Either motivate it at L190 with one clause ("a relative drop of at least $1/(4V)$ immediately past the mode"), or shorten the Pitman paragraph as in 3.5.

**1.11** The priority claim (L328–330) states the contribution accurately. The order of the abstract matches Section 1. Nothing further.

---

## PASS 2: Proofread (grammar, LaTeX, symbols, references)

**Mathematical precision**

- **2.1 L279–280.** "Take $m$ success probabilities equal to $\varepsilon_m$ and $m$ equal to $1-\varepsilon_m$, where $2m\varepsilon_m(1-\varepsilon_m)=1$."
  - Problem: there is no real root when m=1, and the root is not specified.
  - Fix: "For $m\geq2$, take $m$ success probabilities equal to $\varepsilon_m\in(0,1/2)$ and $m$ equal to $1-\varepsilon_m$, where …"
- **2.2 L100 (and L230).** "every mode lies within one unit of the mean" / "Every mode lies within one unit of the mean \cite{darroch1964}"
  - Problem: "within one unit" reads as ≤1. The paper relies on the strict form: L199 has |c-EW|<1, L231 has open intervals, L383 has EW<1, and Prop 1.4's "<2/V" needs it. Pitman's paraphrase says "less than 1".
  - Fix at L100: "every mode lies at distance less than one from the mean". Fix at L230: "By Darroch's theorem \cite{darroch1964}, applied to the tilted laws (which are Bernoulli sums), $\mu(t_-)\in(c-1,c)$ and $\mu(t_+)\in(c,c+1)$."
- **2.3 L382.** "if $f_1<f_0$, then $0$ is a mode"
  - Problem: this follows only by log-concavity.
  - Fix: "if $f_1<f_0$, then by log-concavity $0$ is a mode". The full rewrite is in 3.3.
- **2.4 L269–270.** "increases to $1$ for each fixed $j$"
  - Problem: the limit variable is missing.
  - Fix: "increases to $1$ as $m\to\infty$, for each fixed $j$". Optionally, L271: "envelopes $2^{-|j|}$ and $j^22^{-|j|}$", since the normalizing sum also needs dominated convergence.
- **2.5 L458.** "and \eqref{eq:reciprocal} gives $1/\delta_{D\pm r}\geq1/\delta-r>1$"
  - Problem: this takes r applications of (2.5), not one.
  - Fix: "and $r$ applications of \eqref{eq:reciprocal} give …"
- **2.6 L317–318.** "gives, for $k>\mathbb EW$, $f_{k+1}/f_k<1/\theta(k)$"
  - Problem: θ(k) exists only for 0<k<n. The line also has formula, comma, formula.
  - Fix: "gives $f_{k+1}/f_k<1/\theta(k)$ for $\mathbb EW<k<n$".
- **2.7 L424–428.** "Together with $\delta_0=\delta_n=1$, which covers the two pairs at the ends of the support, this says that for $0\leq j\leq n-1$"
  - Problem: "covers the two pairs" is vague. In fact it supplies only the one direction (2.4) misses at each end.
  - Fix: "The two inequalities not given by (2.4), $1/\delta_0-1/\delta_1\leq1$ and $1/\delta_n-1/\delta_{n-1}\leq1$, hold because $1/\delta_0=1/\delta_n=1\leq1/\delta_k$ for every $k$. Hence, for $0\leq j\leq n-1$,". At L425, also write "the second inequality in \eqref{eq:deficit-recurrence}".

**Undefined at first use**

- **2.8 L181.** "The cubic inequalities behind the proof imply…"
  - Problem: this is their first mention in the body, with a definite article and no citation.
  - Fix: "Cubic inequalities of Hillion and Johnson \cite{hillionjohnson2016} for Bernoulli sums imply that $1/\delta_k$ changes by at most one per lattice step (Section~2), so the bound at $D$ extends across the whole support."
- **2.9 L302.** "$\gamma$" is undefined. Add ", where $\gamma$ is Euler's constant," after the display.
- **2.10 L318–320.** "$\mu(t)=k$", "$V(t)\leq e^tV$": μ(t) and V(t) are defined inside the proof of Prop 1.4 (L223) and leak into the main text. Fix: define the tilt family in running text just before Prop 1.4, or add here "(with $\mu(t)$, $V(t)$ the mean and variance of the tilted law $f^{(t)}$ of the proof of Proposition 1.4)". At L224 add "so $V(0)=V$".

**Notation clashes.** α, Γ and ρ are unused in the file.

- **2.11 L206.** "$|k-\mathbb EW|\leq A\sqrt V$ … for fixed $A$". A clashes with A(δ) in (3.1) and L339 (from INPUT l.86). Use $\alpha\sqrt V$ and "for fixed $\alpha$".
- **2.12 L227–234.** "put $a=-\log q_c$ … $b=-\log q_{c+1}$". These clash with a in (2.7) and b_r in (4.2) (from INPUT l.128). Use $t_-=-\log q_c\leq0$ and $t_+=-\log q_{c+1}>0$, with $f^{(t_\pm)}$, $\int_{t_-}^{t_+}$ and $2-e^{t_-}-e^{-t_+}$.
- **2.13 L300–301.** "$K_{jk}$", "$\operatorname{tr}K$". K clashes with the window K in (2.7) and throughout Sections 2–4. Use $\Gamma_{jk}$ and $\operatorname{tr}\Gamma-\operatorname{tr}\Gamma^2$.
- **2.14 L300 against L318.** θ is the angle of integration in Example 1.6 and θ(k) in the next paragraph. Rename the angle to φ ("$e^{i(j-k)\phi}\,d\phi$"), or write the tilt as $e^{t_k}$.
- **2.15 L324.** "$\operatorname{Bin}(n,p)$ with $np=k+1-\varepsilon$ … there $D=k+1$". Here k is reused as a fixed integer in the paragraph where it indexes Pitman's bound. Use ℓ: "as for $\operatorname{Bin}(n,p)$ with $np=\ell+1-\varepsilon$ and $p<\varepsilon$, where $D=\ell+1$ and $D-\mathbb EW=\varepsilon\downarrow0$."
- **2.16 L207.** "$Y_j$": the summands are B_i elsewhere. Fixed by 1.9.
- **2.17 L84.** $\mathcal C_k$ is used once and sits beside C_r=R_r/L_r in (4.3). This is a minor clash. Leave it, or rename C_r in (4.3) to $\rho_r$.
- **2.18 L260.** "$\Var(g^{(m)})\,\delta_{m+1}\to0$": Var is applied to a pmf, and δ (defined in (1.1) for f) is reused for g without notice. Fix: "Writing $V_m$ for its variance and $\delta_k$ for its deficits, $V_m\delta_{m+1}\to0$."
- **2.19 Minor, existing.** μ(t) clashes with μ_i, and π (the constant, Example 1.6) with π_i(u), both in Section 4.2. Optional: rename μ_i and π_i there.

**Grammar and punctuation.** The house convention puts a comma after a clause-initial phrase.

- **2.20 Missing commas:**
  - L87: "For a normal density of variance $V$, the …"
  - L113: "When $V\geq1$, this set …"
  - L114: "and, by log-concavity, $c:=D-1$ …"
  - L205: "In particular, …"
  - L211: "At the mode, …"
  - L221: "For $t\in\mathbb R$, let"
  - L257: "For $m\geq1$, let"
  - L269: "Finally, $\binom…$"
- **2.21 L96–97.** "like the ultra-log-concavity bound … like a bound of Johnson" → "such as … such as".
- **2.22 L127.** "A positive universal bound without a lower variance threshold is impossible." The compound "lower variance threshold" is ambiguous. Fix: "Without some lower bound on $V$, no positive constant is possible. If …"
- **2.23 L275.** "For $g^{(m)}$, $1/\delta_m\leq4/3$ while" has formula, comma, formula. Fix: "The law $g^{(m)}$ has $1/\delta_m\leq4/3$ while $1/\delta_{m+1}\to\infty$: …"
- **2.24 L456.** "For $r=0$, $D\in\{1,\ldots,n-1\}$" has formula, comma, formula. Fix: "For $r=0$, we have $D\in\{1,\ldots,n-1\}$ because …, and \eqref{eq:deficit-bound} holds with equality."
- **2.25 L457.** "If $D\pm(r-1)\in\{1,\ldots,n-1\}$ for some $1\leq r\leq K$" → "Let $1\leq r\leq K$ and suppose $D\pm(r-1)\in\{1,\ldots,n-1\}$."
- **2.26 L459.** "is not $0$ or $n$" → "is neither $0$ nor $n$".
- **2.27 L221.** "The lower bound is a case of \eqref{eq:support-bound}." → "The lower bound is \eqref{eq:support-bound} with $k=c$."
- **2.28 L237–238.** "where the last step is $(1-1/q_c)(1-q_{c+1})\geq0$. The right-hand side is $V\delta_c$." → "where the last step uses $(1-1/q_c)(1-q_{c+1})\geq0$, and $1-q_{c+1}/q_c=\delta_c$."
- **2.29 L297–298.** "By \cite[Proposition~3(1)]{meckesmeckes2017}, which rests on \cite[Theorem~7]{…}" makes a bare numeric citation the agent. Fix: "By a result of Meckes and Meckes \cite[Proposition~3(1)]{meckesmeckes2017}, which relies on Hough, Krishnapur, Peres, and Vir\'ag \cite[Theorem~7]{houghkrishnapurperesvirag2006}, …"
- **2.30 L310–312.** "The sequence $V_N$ increases and first reaches $1$ at $N=1998$. For larger even $N$" leaves the parity scope implicit, and "larger" excludes 1998 itself. Fix: "Over even $N$, the sequence $V_N$ increases and first reaches $1$ at $N=1998$. For even $N\geq1998$, …"
- **2.31 L313–314.** "these are $0.2149$ and $0.0004$". By the house precision rule, changing the fourth digit of 0.2149 alters nothing. Fix: "about $0.215$ and $0.0004$". Optional.
- **2.32 L328.** "an earlier universal constant $\kappa>0$": the constant is not "earlier"; the proof is. Fix: "We are not aware of an earlier proof that $\kappa_\star>0$." This uses the notation already defined at L131.
- **2.33 L332.** "The proof runs as follows." After five results, "the proof" is ambiguous. Fix: "The proof of Theorem~\ref{thm:main} runs as follows." At L335, add "while $(r+1)\delta<1$" after the displayed inequality.
- **2.34 L197.** "By \eqref{eq:reciprocal}" is a bare forward reference in a Section 1 proof. Fix: "By the reciprocal bound \eqref{eq:reciprocal} of Section~\ref{sec:propagation}". (2.5) does not depend on Cor 1.3, so the argument is not circular.
- **2.35 L60–91 is one LaTeX paragraph** (about 200 words, confirmed in the PDF on pp. 1–2), well over the 100-word ceiling. Break before L79 "We have $0\leq\delta_k\leq1$" or before L79 "At an index with".

**LaTeX faults.** Nothing found in the reviewed spans. Labels resolve, the log has no undefined references, and the \Bigl/\Bigr and \substack usage is correct.

---

## PASS 3: De-AI prose and editorial scar tissue

**Imported reviewer wording** (the spec's "feedback's justification imported as new content"):

- **3.1 L89.** "We ask how much of this survives". This is lifted from INPUT l.36 ("how much of that relationship survives"). "Survives" is an agency verb with an inanimate subject, and "this" is vague. Replace as in 1.6.
- **3.2 L180.** "The first descent is a starting point." This is lifted from INPUT l.7 ("the first descent can serve as the starting point") and l.68 (a referee should not "regard D as an isolated index"). The sentence exists to answer "why D?" and says nothing. Cut it, and open the paragraph with the 2.8 rewrite.
- **Abstract L30, "at every support index".** This is INPUT l.7's phrase. Replace with "at every index $0\leq k\leq n$".
- **Inherited notation.** The A√V window (INPUT l.86) and the a, b tilts (INPUT l.128) are where two of the clashes the parent asked about come from (2.11, 2.12). The reviewer's notation was carried over without adapting it to the paper's own.

**Other findings**

- **3.3 L382.** "Moreover" is on the house's hackneyed-adverb list. Replace: "If $f_1<f_0$, then by log-concavity $0$ is a mode, so $\mathbb EW<1$ by Darroch's theorem \cite{darroch1964} and $V\leq\sum_ip_i=\mathbb EW<1$. Hence $D\geq2$."
- **3.4 L326.** "The drop \eqref{eq:adjacent-ratio-drop} is uniform." This is a sealed closer, and "uniform" is vague. Replace with the content: "By contrast, \eqref{eq:adjacent-ratio-drop} gives $f_{D+1}/f_D\leq1-1/(4V)$ for every Bernoulli sum with $V\geq1$."
- **3.5 L317–326 (the Pitman paragraph as a whole).** DECISIONS.md l.940 records that it was added to answer a cold reader. It derives a consequence of Pitman's bound only in order to dismiss it, and it is the longest comparison in Section 1, serving the secondary result (1.6). Keep the content but compress and move it after L292 (see 1.8b): "Pitman's ratio bound \cite[equation~(21)]{pitman1997}, with $V(t)\leq e^tV$ for $t\geq0$, gives $f_{k+1}/f_k<V/(V+k-\mathbb EW)$ for $\mathbb EW<k<n$. This is stronger than \eqref{eq:support-bound} far from the mean, but it tends to $1$ at $k=D$ when $D-\mathbb EW\downarrow0$, as for $\operatorname{Bin}(n,p)$ with $np=\ell+1-\varepsilon$ and $p<\varepsilon$. By contrast, \eqref{eq:adjacent-ratio-drop} gives $f_{D+1}/f_D\leq1-1/(4V)$ for every such law."
- **3.6 L174–178.** "The lower-bound proof combines inequalities that need not be tight simultaneously (…a symmetrization step in Section 4), and we do not know…" This answers the referee question about the gap between 1/4 and 1/3 (DECISIONS.md l.941) before the reader knows the proof. Shorten it to "We do not know which, if either, endpoint of \eqref{eq:constant-bracket} equals $\kappa_\star$." Then end the outline (after L348) with: "Three steps can lose a constant, namely the mass bounds, the maximal-mass bound, and the symmetrization in Section~\ref{sec:scalar} that replaces each $R_r$ by $L_r$, and they need not be tight simultaneously."
- **3.7 L181.** "behind the proof" is a mild mechanism metaphor. It is replaced by the 2.8 rewrite ("Cubic inequalities of …").

**Rulings on flagged words.** Each was considered.
- L88, "natural scale": keep (mathematical sense).
- L93, "The number of summands is the wrong scale": keep. It names the rejected scale and the bounds that use it, so the opponent is real.
- L241, "The upper bound is elementary": keep. It calibrates credit, and INPUT l.144 explicitly declines a novelty claim. Optionally merge as in 4.R2.
- L297, "rests on": acceptable; "relies on" is offered in 2.29.
- L182, "propagates", and L277, "forbids": keep (standard mathematical usage).
- L205, "In particular": the phrase is fine. The sentence is a duplicate (4.R1).

**Checks that found nothing actionable**
- **Corrective negation without an opponent.** Every negation in scope targets a named position or a proved statement: L93, L107, L253, and abstract L31 ("Ultra-log-concavity alone gives no such bound", backed by Prop 1.5).
- **Colon reveals.** None is a drumroll. L243 is a definition colon, L275 explanatory, L345 elaborative, and L382 a claim followed by its proof.
- **Not found:** paraprosdokians, faux-insight setups, praise-as-structure, faux-coaching, restatement-as-revelation inside single paragraphs, "In short", trailing evaluative participles. A grep for the high-signal vocabulary list matched only "Moreover" in scope (L382) and "Indeed" outside scope (L372, L585).
- **Sealed closers.** Only L326 is a clear one, and there is no run of them. L211 and L253 are lead-in sentences to the next result at the end of a paragraph, which is normal mathematical practice; L211 is also a duplicate (4.R2).

---

## PASS 4: Coherence

**Spine, one line per prose paragraph from L79** (bracketed lines are the statement and proof blocks)

1. L79–91: δ_k is a bounded transform of log-curvature, and the normal benchmark gives 1/V. Does V alone force a lower bound near the mode?
2. L93–107: n is the wrong scale. Known bounds on δ depend on n or on the p_i and degenerate at fixed V. Results on the mode, ratios and masses do not control δ_k.
3. L109–115: D is defined. Under V≥1, 2≤D≤n and c=D-1 is the rightmost mode.
   [Thm 1.1: Vδ_D≥1/4.]
4. L127–129: some lower bound on V is necessary (the Bernoulli(p) example).
5. L131–136: κ_⋆ is defined, with 1/4≤κ_⋆≤1/3.
   [Prop 1.2 and proof: a binomial family with V=1 and δ_D→1/3.]
6. L174–178: the sources of slack. Neither endpoint is known to be sharp.
7. L180–183: the reciprocal bound extends the bound at D to the whole support.
   [Cor 1.3 and proof.]
8. L205–211: three moves in one paragraph: the δ_c bound and the √V window; the extension to differences; the lead-in to the mode result.
   [Prop 1.4 and proof: 1/(4V+1)≤δ_c<2/V.]
9. L241–253: the bracket at the mode fixes the scale. Then a pivot: the lower bounds need more than ULC; (1.7); ULC laws admit no bound scaled by the variance.
   [Prop 1.5 and proof.]
10. L275–292: why ULC fails (the jump in curvature that (2.5) forbids), and how (1.7) and Johnson's bound degenerate for Bernoulli sums at V=1.
    [Example 1.6: a count from a random unitary matrix, 1/log N against 1/N.]
11. L317–326: Pitman's bound gives Gaussian decay but no drop at D, while (1.6) does give one.
12. L328–330: the priority claim.
13. L332–348: the proof outline.

The spine is sound. Question, why n fails, theorem, sharpness, extension, scale at the mode, why ULC is not enough, a natural example, novelty, method: each result serves the question in paragraph 1. The weak points are the joins and the repeats.

**Joins**

- **J1, L107→L109 (2→3).** The text goes straight from "do not control δ_k" to "Let D:=" without saying where the answer will be stated. Open L109 with: "We state the bound at the first index where $f$ decreases,"
- **J2, L178→L180 (6→7).** "The first descent is a starting point" is a weak join. With 3.2 and 3.6 applied, Prop 1.2 leads directly into the propagation sentence.
- **J3, L211→L213 (8→Prop 1.4).** The remark on differences sits between Cor 1.3 and Prop 1.4. Move it (1.9), and break paragraph 8 so that the lead-in to Prop 1.4 stands alone.
- **J4, L239→L241 (Prop 1.4→9).** "it shows that 1/V is the right scale" repeats L211. The paragraph then changes topic to ULC. Break before L242.
- **J5, L292→L294→L317.** Example 1.6 interrupts the comparisons. Move paragraph 11 before Example 1.6 (1.8b).

**Repeated statements**

- **R1. δ_c≥1/(4V+1) appears three times in Section 1:** L205, in the statement of Prop 1.4, and at L221, plus once in the abstract. Delete "In particular $\delta_c\geq1/(4V+1)$, and" at L205 and start "On a window $|k-\mathbb EW|\leq\alpha\sqrt V$, …"
- **R2. "1/V is the right scale at the mode"** appears at L211 ("attained up to constants") and at L241. Keep one. Either make L211 "At the mode, $\delta_c$ is within constant factors of $1/V$:" and cut "it shows that 1/V is the right scale" from L241, or cut L211 and keep L241 as: "The upper bound is elementary. With the lower bound it gives $1/5\leq V\delta_c<2$."
- **R3. ULC is not enough** at L242 ("The lower bounds need more than ultra-log-concavity") and L253 ("No variance-scaled bound holds for ULC laws in general."), two sentences apart around (1.7). Cut L242. L253 introduces Prop 1.5 directly.
- **R4. Darroch's theorem is restated** at L100, L199 and L230. L230 restates the content instead of naming it. Use "By Darroch's theorem" (wording in 2.2).
- **R5. The reciprocal bound is described** at L181 and again at L332–333 in the outline. This is acceptable in an outline. Once L181 introduces the inequalities with a citation (2.8), L332 can start "The reciprocal bound \eqref{eq:reciprocal} shows that a small value of …"
- **R6. Intended forward pointer, not a defect.** The degeneration of (1.7) and Johnson's bound is promised at L98–99 and delivered at L278–292.

**Paragraph lengths over the house ceiling** (mathematics register noted): paragraph 2 (about 130 words), paragraph 10 (about 150; break after L277 "for Bernoulli sums."), paragraph 13 (about 170; break before L342 "Heuristically,"), and the fused L60–91 paragraph (2.35).

The sections outside the reviewed spans (the rest of 2.1–2.2, and Sections 3–4) were not reviewed, except where a symbol clashes with the new text.
