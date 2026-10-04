# Re-verification after the four review items (Brett: "do all four")
<!-- SUMMARY: acknowledgement condensed, two remarks demoted, three references added, further-work sentences added · status: sources read; clean build (errors checked); 13 pp · updated: 2026-10-03 -->

**Changes.**
1. **Acknowledgement condensed.**
   - The substantive credits stay: GPT-5.6's proof search; the first ChatGPT
     review's proposed results; the second review's Pitman bound and
     reconstruction; Aristotle's proofs and monotonicity route.
   - The review-level checks (ChatGPT, Elicit, Codex/GPT-6.1-sol) are now one
     sentence that names every system, plus a pointer to the archive's
     provenance table.
   - That table now has a row for the third ChatGPT review: the D = n
     boundary case, plus the suggestions adopted in this round.
2. **Remarks 1.4 and 1.8 are now prose.** Remark 1.7 is kept, renumbered 1.6.
   The rest renumbers as follows: Prop 1.5 → 1.4 (mode), Prop 1.6 → 1.5
   (ULC), Example 1.9 → 1.7. Current files that cite results by number were
   updated: CERTIFICATE.md, the CUE program's docstring, the portal fields and
   both replay records. A final grep finds no stale numbers.
3. **References.**
   - Saumard and Wellner 2014, cited as a review in the "why it matters"
     paragraph. Its abstract covers log-concavity and strong log-concavity "in
     both discrete and continuous settings", including concentration of
     measure and log-Sobolev inequalities.
   - Brändén 2015, cited for the real-rootedness of the Stirling and Eulerian
     rows. Its Example 1.3 proves the Stirling (surjection) polynomials
     real-rooted, Eulerian real-rootedness is discussed, and §2 covers
     probabilistic consequences.
   - Liggett 1997, cited only for ultra-logconcave sequences (its title, and
     the citation Saumard–Wellner give: "see e.g. Liggett [1997]"). The paper
     itself is behind a ScienceDirect CAPTCHA, which I did not attempt, so its
     convolution theorem is not stated.
   - MR numbers are from AMS MR Lookup: 3290441, 1462561, 3409348. Brändén's
     pages follow MR (437–483); Crossref gives 461–508.
4. **Further work.**
   - After (1.4): whether the binomial laws of Prop 1.2 are asymptotically
     extremal, which would give κ⋆ = 1/3.
   - After Cor 1.3: whether a bound that varies over the support yields a
     functional inequality. This is posed as an open question, not a claim.

**Sources.** The Saumard–Wellner and Brändén PDFs (arXiv) were added to
`literature/`.

**Checks.**
- Build: 0 errors, 0 undefined references, 0 overfull boxes, 13 pages, 14
  references.
- Supplement: rebuilt (77 members), manifest verified, all four programs pass
  from the extracted archive.
