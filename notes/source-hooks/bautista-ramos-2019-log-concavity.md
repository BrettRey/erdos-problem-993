# Source hook: Bautista-Ramos et al. (2019), log-concavity via a partial ordering
<!-- SUMMARY: Discrete Math paper proving log-concavity for some tree families via a new partial order; directly on the Alavi question · status: unread against project state · updated: 2026-08-18 -->

**Source.** Bautista-Ramos, Guillén-Galván & Gómez-Salgado (2019), "Log-concavity
of some independence polynomials via a partial ordering," *Discrete Mathematics*
342(1), 18–28. Central note:
`literature/bautista-ramos_etal_2019_log_concavity_independence_polynomials-note.md`.

**Why this project cares.** It states the target question explicitly: Schwenk
(1981) proved the edge independence polynomial is unimodal; Alavi et al. showed
the vertex independence polynomial need not be, and asked whether a *tree's* is.
That is the question behind #993. The paper then proves log-concavity (hence
unimodality) for some families of trees and related graphs.

**What it may support.** A method rather than an instance: a partial order on
non-negative polynomials giving sufficient conditions for log-concavity, plus an
iterative construction of self-similar graphs from n-paths with pendant graphs
that preserves the property. If the construction composes, it extends the class
of trees settled rather than adding one more example.

**What must be checked before it is used.**
1. Whether these tree families overlap what the project already covers, or extend
   it. A method covering a new family is worth much more than a repeat.
2. Whether the partial order interacts with the mean-bounds and
   structural-reduction machinery already in the manuscript.
3. Whether a 2019 result has been superseded. This area has moved (Lorentzian
   polynomials, Mason's conjecture proved 2018, Rota disproved 2026), so check
   forward citations before building on it.
4. The markdown mangles some formulas (subscripts, relation symbols). Verify
   anything quantitative against the PDF, not the .md.

**Metadata is settled**, so no bibliographic work is owed: *Discrete Mathematics*
342(1), 18-28, DOI `10.1016/j.disc.2018.09.010`, verified against both the PDF
and Crossref. Suggested key `bautistaramos2019logconcavity`.
