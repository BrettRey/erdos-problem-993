# Paragraph opening audit: PB paper
<!-- SUMMARY: paragraph-opening audit of variance-local-log-concavity-poisson-binomial.tex · status: 1 cut (as redundancy R1) and 2 cohesion moves, applied 2026-10-03 · updated: 2026-10-03 -->

Manuscript: `paper/poisson_binomial/variance-local-log-concavity-poisson-binomial.tex`,
sha256 `ccdf50d820fc80db…`, repo HEAD `7a3da7e`. Openings came from
`extract_openings.py` and were checked against the source. Lines that only
continue a paragraph after a display (80, 488, 599, 817) are not openings and
are left out. Proposition, lemma and proof openings ("Let…", "Set…", "For…",
"We induct…") are standard mathematical openings, and all are OK.

```text
line 53  OK: "Let W = …" (standard setup).
line 61  OK: deals with degenerate summands before the definitions need them.
line 82  OK: begins the curvature reading on which the question depends.
line 97  OK: "The number of summands is the wrong scale for this question." Picks up "no dependence on the number of summands" from the paragraph before and makes the claim the paragraph then argues. It's the only "The X is…" opener in §1.
line 117 OK: introduces D, which the theorem needs.
line 135 OK: explains why the theorem needs V≥1.
line 138 OK: defines κ⋆.
line 181 OK, but see 379: "We do not know which, if either, endpoint … equals κ⋆."
line 184 OK: names the tool (Hillion–Johnson cubics) and the consequence the corollary states.
line 210 OK.
line 215 OK: introduces Proposition 1.4.
line 245 CUT: "The upper bound is elementary." This is an evaluative opener, and the rest of the paragraph ("With the lower bound it gives 1/5≤Vδ_c<2") restates line 215. Handled as R1 in the redundancy report.
line 248 OK: definition of ULC(n).
line 281 OK: says why the ULC example is not a Bernoulli sum.
line 285 OK: delivers the intro's forward pointer (l. 108–110).
line 326 OK: long, because the derivation follows the colon. Splitting it at the colon would add words without making it clearer.
line 337 OK: changes topic to scope of application, and its last sentence leads into Example 1.6.
line 379 COHESION CHECK: "We are not aware of an earlier proof that κ⋆>0." This one-sentence paragraph sits between Example 1.6 and the proof outline, but κ⋆ was last discussed at l. 138–182. Proposed (P1): move it to follow l. 181, making a two-sentence paragraph on what is known about κ⋆.
line 381 OK: "The proof of Theorem 1.1 runs as follows." Conventional signpost for a proof outline.
line 393 COHESION CHECK: the paragraph (128 words, above the house maximum of 100) covers three things: the heuristic scale, which ranges depend on machine computation, and where constants are lost. Proposed (P2): start a new paragraph at "Section~\ref{sec:scalar} proves the one-variable inequality…" (l. 396). No words are added or removed.
line 410 OK: one-sentence section roadmap, conventional in a mathematics paper. It overlaps the outline (see the redundancy report, kept).
line 415, 428, 443 OK.
line 452, 462 OK.
line 473 OK: the second paragraph in a row to open with "Hillion and Johnson". That's mild repetition and a different subject would be no clearer.
line 499, 530, 539, 555 OK: "On the left, let…" mirrors the right-hand argument.
line 580, 593 OK.
line 612 OK: "It remains to prove one scalar estimate." A conventional lead-in to the proposition.
line 641, 660, 693, 707, 723, 731, 744, 776, 785, 793 OK.
```

Cadence across the paper: four openings start with "The" (97, 285, 381, 731)
and seven with "We" (117, 181, 379, 410, 593, 707, 744). That's about 45
openings in all, with no run of three. No "This…" openings with an unclear
referent, no throat-clearers, no mic-drop runs.

## Proposed changes (P1 and P2 applied 2026-10-03 with Brett's approval; new tex sha256 `b5bc1d8285c5521a…`)

- **P1** Move l. 379 to follow l. 181.
- **P2** Paragraph break before l. 396.
- The cut at l. 245 is R1 in `2026-10-03-redundancy.md`.
