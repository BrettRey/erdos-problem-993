# Pre-flight re-check of the last two manuscript commits
<!-- SUMMARY: diff-scoped re-check of 052db5a and 40a356c for the five prose passes last run at 20:49; one editorial scar found and fixed (sentence order after demoting the window remark) · status: clean build (errors checked); 13 pp · updated: 2026-10-03 -->

**Why.** Five prose passes were last recorded at 20:49: lakoffian-metaphor,
paragraph-opening-audit, rhetoric-and-humour, redundancy and density-leavening.
Two manuscript commits followed: 052db5a (21:04, the four review items) and
40a356c (21:23, readability edits). de-ai-prose and editorial-scar-tissue were
re-checked after each, but these five were not. Following Brett's instruction,
no new cold read was run. The check covers only the word diff `5efde19..40a356c`.

**Findings.**

- **Metaphor.** None introduced. "It arises in reliability analysis…" and
  "found a boundary case" are literal.
- **Paragraph openings.** Two remarks were turned into prose, so two
  paragraphs now have new first sentences. "A polynomial with nonnegative
  coefficients and only real zeros is…" is OK. "On a window…" needed the fix
  below.
- **Rhetoric.** Nothing in the diff.
- **Redundancy.** §1 now repeats "independent" within one sentence ("the
  summands are independent, so that $W$ counts the successes in independent
  trials…"). The second occurrence belongs to the plain-words gloss Brett
  asked for, so it was left. The abstract and §1 share the same gloss, which
  is normal.
- **Density.** The edits make the text less dense. There is no new
  compression.
- **Editorial scar (fixed).** Demoting the window remark joined it with ", and"
  to an unrelated sentence: that Theorem 1.1 and Corollary 1.3 extend to
  differences of Bernoulli sums. The open question that follows ("Whether a
  bound that varies over the support in this way…") then had to reach back
  past that clause to the bound it refers to. Fix: the sentences were
  reordered with no new words. The scope sentence now comes first, then the
  window bound, then the open question.

  Before: "On a window …, the bound (1.5) stays on the 1/V scale for fixed α,
  and Theorem 1.1 and Corollary 1.3 also hold for differences of independent
  Bernoulli sums (Section 2.1). Whether a bound …"

  After: "Theorem 1.1 and Corollary 1.3 also hold for differences of
  independent Bernoulli sums (Section 2.1). On a window …, the bound (1.5)
  stays on the 1/V scale for fixed α. Whether a bound …"

**Checks after the fix.**
- Clean rebuild: 0 lines beginning `!`, 0 undefined, 0 overfull, 0
  underfull, 13 pages.
- The source bundle was rebuilt (SHA-256 `e1abd93e…07bd16`). It compiles alone,
  and its text is identical to the manuscript PDF.
- The Zenodo archive holds no manuscript file, so the deposit is unaffected.

**Passes made stale by the fix itself.** The reorder moved §1's prose, so eight
more passes went stale: adversarial-cold-read, level-category-audit,
terminological-hygiene, negative-claims-audit, reader-pass, coherence-cohesion,
external-review-triage and source-reread. Each was judged against this diff,
which reorders three existing sentences without adding words. So the change
adds no term, level, claim, citation or review item. It is on page 3, outside
the first two pages a cold read covers. It also improves the cohesion and
reader dimensions, which were checked by reading the result in the source and
on the rendered page. Each was recorded with its reason in the pass ledger.
