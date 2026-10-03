# Re-verification after the prose-pass edits: PB paper
<!-- SUMMARY: re-checks of the passes the 2026-10-03 prose edits made stale, scoped to the diff · status: all re-checked, no new findings · updated: 2026-10-03 -->

Manuscript: `paper/poisson_binomial/variance-local-log-concavity-poisson-binomial.tex`.
Before the edits: sha256 `ccdf50d820fc80db…` (HEAD `7a3da7e`). After:
`b5bc1d8285c5521a…`. PDF after: `aceeb103b81235c5…`, 12 pages. Page 12 now
opens with the acknowledgement; the references follow.

## Why scoped to the diff

Brett, 2026-10-03: the pass rules exist for a purpose and aren't
box-ticking. Each stale pass was therefore re-checked against what it exists to
catch, applied to the text that actually changed. Passes weren't re-run from
scratch because the edits were small. Where a pass's input didn't change in a
way that bears on its purpose, the earlier full run still stands, and this note
says so.

## The diff (seven hunks)

1. Prop 1.2: cut "This gives $\kappa_\star\leq1/3$." (R2)
2. After Prop 1.2: "We are not aware of an earlier proof that $\kappa_\star>0$." moved from after Example 1.6 to follow "We do not know which, if either, endpoint … equals $\kappa_\star$." (P1)
3. l. 186: "extends across the whole support" → "gives a bound at every point of the support". (M1)
4. Lead-in to Prop 1.4 now ends "; since $V\geq1$, the next proposition gives $1/5\leq V\delta_c<2$". Cut "The upper bound is elementary. With the lower bound it gives $1/5\leq V\delta_c<2$." (R1, fold variant)
5. Removed the old position of hunk 2.
6. Paragraph break before "Section~\ref{sec:scalar} proves the one-variable inequality". (P2)
7. §4.2 close: "Thus … for $H\geq16$, which with Section~\ref{sec:compact} proves Proposition~\ref{prop:scalar}. A short SymPy program … exact arithmetic." Replaces "… for $H\geq16$. A short SymPy program … This settles $H\geq16$; with Section~\ref{sec:compact}, it proves …". (R3)

No citation, displayed equation, theorem statement, label or number was
added. The only number removed is the inline $1/3$ in hunk 1; $(1.4)$ still
states $\kappa_\star\leq1/3$.

## Per pass

| Pass | Purpose | Check on the diff | Result |
|---|---|---|---|
| adversarial-cold-read | Is the problem, the gap and the contribution legible on pages 1–2? | `coldread.py extract --prompt-only` before and after: hunks 1–4 fall in the opening. None touches the problem paragraph (l. 82–115) or Theorem 1.1. Hunk 2 adds an explicit novelty claim to page 2. | The round-9 verdict (4/4 advance) stands. No new readers. |
| coherence-cohesion | Does each paragraph follow from the last? | Read the junctions in the PDF: Prop 1.2 proof → κ⋆ paragraph → Cor 1.3 lead-in; Prop 1.4 proof → ULC definition; Example 1.6 → proof outline; the split outline paragraph. | Better: the κ⋆ statements are together now. No breaks. |
| reader-pass | Pointers, forward references, things a reader trips on. | No pointer targets removed. Prop 1.2 no longer states its consequence, but the sentence before it ("The theorem and the next proposition give (1.4)") still does, and the proof ends at $\to1/3$. | No issue. |
| level-category-audit | Category errors (philosophical sense). | New predications: "the bound at $D$ gives a bound" (one inequality implying another); "which with Section 4.1 proves Proposition 3.1" (the same section-proves metonymy as before, standard in mathematics). | No issue. |
| terminological-hygiene | Consistent terms; no undefined terms. | No new terms. "rightmost mode", "first-descent", κ⋆ are used as before. | No issue. |
| negative-claims-audit | Is every negative claim backed by a stated search? | The one negative claim in the diff ("not aware of an earlier proof that $\kappa_\star>0$") moved with its wording unchanged. Its search record is `notes/literature/poisson_binomial_novelty_update_2026-10-03.md` (five databases plus citation chains; MathSciNet not run). | No issue. |
| numbers-audit | Is every number right? | Removed inline $1/3$ (still in (1.4)). The moved $1/5\leq V\delta_c<2$ follows from Prop 1.4: $1/(4V+1)\geq1/(5V)$ iff $V\geq1$, and $\delta_c<2/V$. | Correct. |
| source-reread | Do attributed claims match their sources? | No hunk changes a cited or attributed claim. Hunk 3 changes our own inference in a sentence citing Hillion–Johnson, and the attributed part ("$1/\delta_k$ changes by at most one per lattice step", via Section 2) is unchanged. | No issue. |
| external-review-triage | Do the edits undo an accepted review point? | The ChatGPT Pro referee report calls $1/5\leq V\delta_c<2$ the statement about the mode. R1 was changed from a cut to a fold so that statement stays. No other hunk touches a review point. | Kept. |
| build-integrity | Clean build. | pdflatex/bibtex/pdflatex ×2: 0 undefined references, 0 overfull boxes, 12 pages. The one hyperref warning (`pagebackref` already used) involves the option that `ejpecp.cls` sets; no manuscript text is involved. | Clean. |
| proofread | Errors in the rendered text. | Read every hunk in `pdftotext -layout` output. | No errors. |
| de-ai-prose | AI tics. | The new wording has no listed tics. The semicolon clause in hunk 4 is plain apposition with no colon reveal. | No issue. |
| editorial-scar-tissue | Does the revision show its history? | No hunk reads as a correction of earlier wording, and nothing refers to cut text. | No issue. |
| lakoffian-metaphor, paragraph-opening-audit, rhetoric-and-humour, redundancy | The passes that produced these edits. | Their proposals are applied as reported (R1 as the fold variant). | Re-recorded against the new text. |

## Side effect outside the manuscript

R3 changes §4.2's closing sentences, so
`formalization/pb_scalar_inequality_aristotle/PROOF_CONTEXT.md` (sent to
Aristotle as "§4 verbatim") now differs from the paper in wording only. The
mathematics is identical. The Aristotle replay should not count this as drift.
