# What an ECP paper looks like, and how ours compares
<!-- SUMMARY: genre assessment of the PB paper against seven recent ECP papers (structure, phrasing, style, exemplification) · status: recommendations awaiting Brett · updated: 2026-10-03 -->

Brett asked about structure, phrasing, writing style and the amount of
exemplification, because he doesn't know what mathematics papers are usually
like.

## Evidence

**Comparators.** Seven recent ECP papers, read in their arXiv versions:
- Niles-Weed, *Reweighted information inequalities*;
- Chewi, Pooladian and Zhang, *Stability of the Kim–Milman flow map*;
- Wang, *Infinite divisibility of α-Cauchy distributions*;
- Samani, *An exponential concentration inequality for the components of a
  uniform random vector*;
- Austern and Mackey, *Bounding Hellinger distance with Stein's method*;
- Athaide et al., *A correlation bound for the one-dimensional heterogeneous
  Ising model*;
- Cai, *Rates of bulk convergence for ensembles of classical compact groups*.

**Selection.** An OpenAlex query on ECP for 2023–2026, keeping short papers on
inequalities or distributions that have open arXiv copies.

**Two independent reads.**
1. My quantitative profile (table below).
2. A Codex (GPT-6.1-sol) editorial read grounded in the same seven texts,
   saved verbatim at `runs/pb-genre-review-20261003/REPORT.md`.

The two reads agree on every main point. Seven papers show a range of
practice, not a rulebook.

## Quantitative profile (from the PDF text; word counts are rough because math symbols are excluded)

| Paper | Words | Remarks | Examples (numbered) | Tables | "we/our/us" per 1000 words | Mean sentence length | Conclusion section |
|---|---:|---:|---:|---:|---:|---:|---|
| **Ours** | 3804 | **0** | 1 | 1 | **6.8** | **32.7** | no |
| Austern–Mackey | 3789 | 3 | 0 | 0 | 16.1 | 31.7 | no |
| Athaide et al. | 3868 | 1 | 0 | 0 | 16.3 | 27.6 | no |
| Samani | 2312 | 1 | 0 | 0 | 21.2 | 30.4 | no |
| Cai | 2462 | 0 | 0 | 0 | 18.3 | 31.5 | no |
| Chewi–Pooladian–Zhang | 3662 | 2 | 0 | 0 | 22.4 | 25.6 | no |
| Wang | 2208 | 1 | 0 | 0 | 13.6 | 20.8 | no |
| Niles-Weed | 3840 | 2 | 0 | 0 | 13.3 | 25.7 | no |

Our Section 1 runs from p. 1 to p. 6 of 13, about 45% of the paper. It
contains four complete proofs (Propositions 1.2, 1.4, 1.5 and Corollary 1.3),
the balanced-family comparison with Table 1, a paragraph comparing Pitman's
ratio bound, Example 1.6, and the proof strategy.

## The conventions, for a reader from another field

1. **Introduction: recognizable moves, flexible order.** A reader should be able
   to find the setting, the question or gap, the contribution measured against
   prior results, the method, and where the details are. The order varies. What
   is shared is that each part of the opening has one visible function: Athaide
   et al. use subsections, and Austern–Mackey move most material out of the
   introduction.
2. **Results stated early, proved later.** The main theorem appears in §1 or §2.
   Secondary results are often stated near the main theorem, but their proofs
   usually go to later sections (Wang, Samani, Austern–Mackey, Athaide et al.).
3. **Remarks carry commentary.** Six of the seven papers set off interpretation,
   sharpness and extensions as numbered *Remarks* (one to three each). Ours puts
   the same kind of commentary in unlabeled paragraphs between results.
4. **Examples: no quota, and ours is not short of them.** None of the seven
   uses a numbered Example. Several have one worked application, and the
   others have none, because the model they study is already concrete. Ours
   has more exemplification than any of them: a sharpness family (Prop 1.2), a
   counterexample (Prop 1.5), the balanced family with Table 1, and the
   Haar-unitary application. Each does distinct work, so none should go, but
   more examples would be against the genre, not toward it.
5. **Voice and register.** "We" is standard even for a sole author (Wang and
   Niles-Weed are single authors). The register is direct assertion with
   procedural verbs: *let, put, define, fix, consider, hence, it remains to
   prove*. Ours is in this register. It uses "we" about half as often as the
   comparators, preferring constructions like "the theorem and the next
   proposition give". That is not wrong, only drier than usual. Our sentences
   are the longest in the set on average, and Codex points to *discourse*
   density rather than *mathematical* density: paragraphs that introduce
   notation, derive a bound, change regime and assess novelty all at once.
6. **Literature.** Cite with exact locators, say what the cited result
   supplies, and qualify novelty claims ("we are not aware of"). Ours does
   this well.
7. **Endings.** None of the seven has a conclusion section. Papers end on the
   last proof or remark. Open problems sit next to the result they concern, as
   ours does.
8. **AI acknowledgement.** One 2026 ECP paper credits ChatGPT and Claude in a
   single clause. That shows disclosure is accepted; it doesn't settle how much
   detail is right when the AI contributions were substantive, as ours were.

## Recommendations

**Likely to be raised by a referee or editor:**

| # | Change | Where | Why |
|---|---|---|---|
| A | **Restructure §1.** Keep the setting, gap, Theorem 1.1, the *statements* of Prop 1.2, Cor 1.3, Prop 1.4 and Prop 1.5 (each with one line of significance), a short version of Example 1.6, and the proof strategy. Move the four proofs, the balanced-family comparison with Table 1, the Pitman ratio paragraph and the Example 1.6 calculation to a new section after the main proof ("Proofs of the results of Section 1" or "Consequences, sharpness and examples"). | §1 → new §5 | §1 drops from about 5½ pages to about 3, and the proof strategy moves from p. 6 to about p. 3. This matches conventions 1–2. Both reads give it top priority. |
| B | **Turn the unlabeled commentary paragraphs into numbered Remarks** (the √V window, the reciprocal-bound explanation after Prop 1.5, the real-rooted-polynomial application, the Pitman ratio comparison). | §1 | Convention 3. Readers can see what is commentary and can skip it. |
| C | **Split the balanced-family/Pitman paragraph** (construction, Johnson's bound, the Pitman-derived bound and its asymptotics, the Skellam limit). Display the Pitman-derived bound and define θ where it is used. | §1 | Codex calls it the least readable passage, because it changes task five times. |
| D | **Trim the proof roadmap.** Keep the strategy and the heuristic. Drop the list of where constants may be lost, which uses H, R_r, L_r, b_r and λ_r before they are defined, and move the parameter ranges to §4. | end of §1 | It introduces notation before its definition. |
| E | **Split the supplement description into three parts:** what the proof needs, extra checks, and formalization scope. | supplement box | One paragraph currently does four jobs. |
| F | **Abstract:** make the normal-curvature comparison its own clause, so that "which" can only refer to the curvature, and write "the proof of the main theorem" (the other results are not formalized). | abstract | Removes an ambiguity, and stops the abstract from implying every result is formalized. |

**Matters of taste (Codex has a strong preference on G):**

| # | Change | Why |
|---|---|---|
| G | **Acknowledgement:** keep every specific credit (your rule), but move the chronology (dates, "Pro effort", the sequence of reports) to the archive's provenance file, and call the ChatGPT and Elicit outputs "AI-generated reviews", not "referee reports". | It currently reads like an audit log, and "referee report" could be confused with the journal's own peer review. |
| H | **Example 1.6:** keep the Bernoulli representation, the growth of the variance, the change of order from 1/N to 1/log N, and one number. Move the exact first-applicable N = 1998 to the supplement. | The threshold draws attention away from the point of the example. |
| I | **Smaller points:** restate the setting inside Theorem 1.1 so it can be quoted alone; drop the forward reference to S and T at the start of §4.1; optionally give only the minimum of the 12 printed margins. | Codex suggests these; I'd keep the margins, which you wanted for readability. |

**Leave alone:** the mathematical voice and register; the curvature
motivation; every example and counterexample; Table 1; the heuristic and the
explicit mass calculations; the visible positivity certificates; exact source
locators and qualified novelty claims; the open-constant question beside its
bracket; the absence of a conclusion. Don't add elementary examples, a
conclusion, or generic "applications" rhetoric to look more conventional.

## Cost

A and C move a lot of text. B, D, E and F are local. None of them changes the
mathematics, but A changes pages 1–2, which is what a cold reader sees.
