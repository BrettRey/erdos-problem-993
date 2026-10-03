# Lakoff-style metaphor audit: PB paper
<!-- SUMMARY: metaphor audit of variance-local-log-concavity-poisson-binomial.tex · status: one minor steering case (M1), applied 2026-10-03 · updated: 2026-10-03 -->

Manuscript: `paper/poisson_binomial/variance-local-log-concavity-poisson-binomial.tex`,
sha256 `ccdf50d820fc80db…`, repo HEAD `7a3da7e`. Whole manuscript read, including
abstract, Table 1, supplement description and acknowledgement.

## 1. Source domains

Most candidate frames are mathematical terms of art, so they're not imported
metaphors: *curvature* (C_k really is the negative second difference of
log f), *slack* (l. 74), *window* (l. 210, 842), *cell* (l. 645, 707),
*cover* (l. 744), *degenerate* (l. 108, 285), *tilt* (l. 318), *tight*
(l. 405). No section is organized by one of these the way an argument paper is
organized by a frame.

Four frames carry non-technical entailments:

| Frame | Passages | Entailment |
|---|---|---|
| SCALE / NATURAL UNIT | l. 91 "the natural scale for δ_k near the mode"; l. 97 "the wrong scale"; l. 211 "stays on the 1/V scale" | there is a canonical unit, and the quantity is of that order in both directions |
| FORCE (logical implication as compulsion) | l. 382, 395, 410 "forces" | the premise compels the conclusion |
| FOUNDATION | l. 350 "on which this rests"; l. 397 "the proof rests on 275 rational coefficients" | the upper structure fails if the base fails |
| SPATIAL EXTENSION | l. 186 "the bound at D extends across the whole support" | the thing extended keeps its size as it spreads |

## 2–3. Do the entailments fit?

- **Natural scale (l. 91): fits.** The frame implies the claim is two-sided, of
  order 1/V, while the main theorem only gives a lower bound at D. The sentence
  says "suggests" and the paragraph goes on to ask a question, so nothing is
  granted. The two-sided claim also holds near the mode. Proposition 1.4 gives
  1/(4V+1) ≤ δ_c < 2/V, and with the reciprocal bound (2.5) this yields
  1/δ_k > V/2 − |k−c|, so δ_k = O(1/V) on |k−c| ≤ α√V. The frame is earned,
  although the paper states the upper bound only at c. **No change.**
- **Force: fits.** This is the standard mathematical sense of "forces" =
  implies, and each use is followed by the explicit inequality. **No change.**
- **Foundation: fits, and it helps the disclosure.** "The proof rests on 275
  rational coefficients computed by a program" says exactly what depends on the
  computation, which is the point the Elicit report asked to have made.
  **No change.**
- **Spatial extension (l. 186): mild steering.** "The bound at D extends
  across the whole support" implies the bound keeps its strength. In fact it
  weakens by one unit of 1/δ per step (Corollary 1.3: δ_k ≥ 1/(4V+|k−D|)).
  The corollary corrects this four lines later, and l. 210 limits the 1/V-scale
  claim to a √V window. Still, the sentence briefly suggests more than is
  proved. **Proposed (M1):** "so the bound at D gives a bound at every point of
  the support." Same length, and no implication about strength.

## 4. Mixed frames

None found. The two FORCE uses in the outline (l. 382, 395) sit next to the
heuristic "comparable to M" without crossing into another domain.

## Figures and diagrams

The paper has no figures. Table 1 lists bounds and has no layout metaphor.

## Load-bearing cap

`grep -i load-bearing` finds no matches. The cap of two uses is satisfied.

## Proposed change (applied 2026-10-03 with Brett's approval; new tex sha256 `b5bc1d8285c5521a…`)

- **M1** l. 186: "so the bound at $D$ extends across the whole support" →
  "so the bound at $D$ gives a bound at every point of the support".
