# Redundancy pass: PB paper (proposal stage)
<!-- SUMMARY: redundancy pass on variance-local-log-concavity-poisson-binomial.tex · status: R1 (fold variant), R2, R3 applied 2026-10-03 · updated: 2026-10-03 -->

**Outcome (2026-10-03).** Brett approved the edits. R2 and R3 were applied as
proposed. R1 was applied in its fold variant, not as a cut: the ChatGPT Pro
referee report (`runs/pb-revision2-20261003/INPUT_chatgpt_pro_referee_report.md`,
l. 283–286 and 541–544) calls $1/5\leq V\delta_c<2$ the paper's statement
about the mode, so the constants now sit in the lead-in to Proposition 1.4.
Only "The upper bound is elementary" and the restating paragraph were removed.
The resulting tex has sha256 `b5bc1d8285c5521a…`.

Manuscript: `paper/poisson_binomial/variance-local-log-concavity-poisson-binomial.tex`,
sha256 `ccdf50d820fc80db…`, repo HEAD `7a3da7e`. Criterion (registry): text
that says again what the paper has already said. Not cutting for length. The
acknowledgement is out of scope, under the credits rule.

## Section scale

No two sections make the same move. §1 states the results and outlines the
proof, §2 gives the recurrence and the mass bounds, §3 the reduction, §4 the
scalar inequality. The intro's forward pointers (l. 62, 110, 122) are
delivered later, not repeated.

## Paragraph scale

- **R1** l. 245–246, "The upper bound is elementary. With the lower bound it
  gives $1/5\leq V\delta_c<2$." This closer restates the lead-in at l. 215
  ("At the rightmost mode, δ_c is within constant factors of 1/V"), and the
  constants follow in one step from 4V+1 ≤ 5V. "Elementary" is a judgement the
  proof already makes visible. **Cut the paragraph.** (Alternative if the
  explicit constants are wanted: fold them into l. 215 and still cut l. 245–246.)

## Sentence scale

- **R2** l. 151, "This gives $\kappa_\star\leq1/3$." (inside Proposition 1.2).
  This repeats l. 140–143, "The theorem and the next proposition give
  (1.5)". **Cut from the proposition.** The paragraph before it still states
  the consequence and introduces (1.5), which l. 181 refers to.
- **R3** l. 817–821: "Thus $N_J(H)>0$ and $\tilde S_J\tilde T_J>Q(H)$ for
  $H\geq16$." and then "This settles $H\geq16$". That's a doubled closure.
  **Proposed:** "Thus … for $H\geq16$, which with Section~\ref{sec:compact}
  proves Proposition~\ref{prop:scalar}. A short SymPy program in the
  supplement checks every expansion in this subsection in exact arithmetic."
  The program sentence then ends the section. Note that this edit is in §4, so
  the "§4 verbatim" copy in
  `formalization/pb_scalar_inequality_aristotle/PROOF_CONTEXT.md` would no
  longer be verbatim. The difference is wording only, so it isn't drift; it
  should be recorded so the later Aristotle replay doesn't flag it.

## Considered and kept

- The program and independent checker are described three times (outline
  l. 396–401, §4.1 l. 735–739, supplement). Each serves a different entry
  point: the opening reader, the referee checking the proof, and the archive
  reader. The outline sentence was added on purpose after the Elicit report.
- Darroch's theorem is re-cited with its theorem number at l. 203, 235 and
  440. In a proof that helps local readability.
- l. 771, "Since $J\leq K$, … $b_r\geq\lambda_r$", restates l. 751 and (4.10).
  A referee reads this as explicitness, not repetition.
- Theorem 1.1 re-expands δ_D in its display, which keeps the statement
  self-contained.
- The outline (l. 384) restates Lemma 2.1's bound. That's the job of a proof
  sketch.
- l. 94–95, "universal κ>0, with no dependence on the number of summands or
  on the success probabilities". The elaboration sets up the next paragraph's
  list of non-variance-scaled bounds, and cold-read rounds asked for exactly
  this.
- l. 281–283, the remark after Proposition 1.5. It states the mechanism that
  separates ULC laws from Bernoulli sums; it doesn't restate the proposition.
- The section-opening roadmaps at l. 410 and 580 are conventional and one
  sentence each.

No sealed-closer runs: the paper has no two consecutive paragraphs ending in
a summarizing sentence. No restatement-as-revelation found.

## Effect

R1–R3 remove about four lines of PDF. That's a side effect, not the criterion.
