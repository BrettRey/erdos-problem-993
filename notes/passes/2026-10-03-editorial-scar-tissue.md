# Editorial scar tissue: PB paper (full pass)
<!-- SUMMARY: full editorial-scar-tissue pass after the day's revisions · status: 3 hits fixed; acknowledgement left as is by genre · updated: 2026-10-03 -->

Manuscript: `paper/poisson_binomial/variance-local-log-concavity-poisson-binomial.tex`,
after the genre restructure (commit `ff4d15a`). This is a full read. The ledger
had listed the pass as current, but only through re-checks scoped to each diff,
and the text had been through many correction rounds today:
- ChatGPT Pro (the Pitman bound);
- the verifier's fixes;
- Codex (zero extension, formalization scope);
- the density-leavening additions;
- the restructure.

Scar tissue collects at exactly those places, so I read them most closely.
Method: the five tells in `.claude/rules/writing-style.md`, "Revise for the
end-reader, not the conversation", plus a grep for corrective-negation markers
(not, no, rather than, instead, itself, by themselves, indeed, earlier, only,
alone, whereas, but). Every match was read in context.

## Hits (applied)

| # | Location | Tell | Old text | Fresh draft |
|---|---|---|---|---|
| 1 | Supplement description, "Additional check" | 4: a clause written so a misreading can't arise (from the verifier's finding that "these programs" was ambiguous) | "… recomputed by an independently written checker. The proof does not need it." | Sentence deleted. The label "Additional check" and the preceding paragraph, which names the two computations the proof needs, already say this. |
| 2 | §4 preamble, §4.1 opening, §4.2 | 3: residue of two relocations (the density pass, then the verifier's forward-reference fix). The same justification appeared three times. | §4.1: "Here we keep the weights as they stand; the symmetrized bound of Section 4.2 is too weak near H=3." §4.2: "This bound is too weak near H=3, where ST<Q(H) (at H=7/2, for example), which is why Section 4.1 does not use it. It remains to prove ST≥Q(H) for H≥16." | §4.1 opening deleted; the §4 roadmap already says §4.1 keeps the weights as they stand. §4.2: "This bound is too weak near H=3 (ST<Q(H) at H=7/2, for example), so we use it only for H≥16, where it remains to prove ST≥Q(H)." |
| 3 | §1, paragraph on existing bounds | 2/3: the correction appended as its own sentence. The earlier draft wrongly said Pitman gives no lower bound, and the fix added "Pitman's two-sided bounds … also give a positive lower bound on δ_k" after the list. | "… proved by Dümbgen and Wellner [4, Proposition 2]. Pitman's two-sided bounds on adjacent mass ratios through exponentially tilted laws [10, equation (20)] also give a positive lower bound on δ_k. All four degenerate …" | Pitman now sits in the list where it belongs: "… on the success probabilities, such as a bound of Johnson [7] and the bound that follows from Pitman's two-sided estimates of adjacent mass ratios through exponentially tilted laws [10, equation (20)]; or on k alone, …". Pitman's bound depends on the law through the tilt function θ, so it fits the "success probabilities" slot. "All four" still counts correctly. |

## Considered and kept

- **"… but by themselves give no lower bound on δ_k"** (§1, Darroch and the
  maximal-mass bounds). This is positioning, not a correction: the paper uses
  both results, and the venue record asks pages 1–2 to separate the new result
  from them.
- **"Hillion and Johnson extend the pmf by zero outside its support"** (§2.2).
  Since Codex's fix, §1 extends f by zero everywhere, so this sentence might
  look like residue. It is a source-fidelity note: it records that the cited
  inequalities are stated "for every k ∈ ℤ" under the same convention. Kept.
- **"Ultra-log-concavity alone gives no variance-scaled lower bound on δ_D"**
  (abstract) and **"No variance-scaled lower bound on δ_D holds for ULC laws in
  general"** (§1). Both state a result (Proposition 1.6); neither repairs a
  draft.
- **"The other results of Section 1 are not formalized"** (supplement). A
  needed scope statement, not a defence against a misreading.
- **Remark 5.1 (Pitman's ratio bound).** It was first added to answer a cold
  reader's question, but in the current paper it supports the novelty of
  (1.6), which is how a referee would read it.
- **Acknowledgement.** Phrases such as "prompted corrections to …" and
  "prompted the reorganization of Section 1" refer to drafts the reader never
  saw. That is the job of an acknowledgement: it reports the history of
  contributions, and Brett's rule is to credit each one specifically. No body
  text depends on it. Left as is.

No contrastive negations of earlier wording and no opponent-free corrective
negations remain in the body.
