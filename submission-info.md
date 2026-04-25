# Electronic Journal of Combinatorics -- Submission Info

## Title

Mean bounds, structural reductions, and exhaustive verification for tree independence polynomial unimodality

## Author

Brett Reynolds
ORCID: 0000-0003-2407-9448
Humber College, Toronto, Canada
brett.reynolds@humber.ca

## Abstract

Alavi, Malde, Schwenk, and Erdős (1987) conjectured that the independence polynomial of every tree is unimodal. I prove $\mu(T) < n/3$ for every $d_{\mathrm{leaf}} \le 1$ tree on $n \ge 3$ vertices via Steiner peeling on a compensation function. Hub Exclusion and Transfer give a conditional reduction for 1-Private maximal independent sets: Conjecture A on the $d_{\mathrm{leaf}} \le 1$ regime plus a Case B hub bound would prove that every such set has size at least the mode. Conditional on $\mathrm{mode} \le \lceil\mu\rceil$, the mean bound settles the Conjecture A part. A subdivision-contraction identity shows that, conditional on a conjectured Edge Contraction Mode Stability property and a tail condition verified through $n \le 19$, any minimal counterexample is homeomorphically irreducible. Neither conditional route closes the full conjecture. Computationally, I verify unimodality for all 8,691,747,673 trees on $n \le 29$ vertices. Multi-arm stars are the empirically extremal family in the tested regimes, with near-miss ratio $\mathrm{nm}(s) = 1 - C/s + O(1/s^2)$, $C \in [4, 8]$.

## Keywords

independence polynomial, unimodality, trees, Erdős Problem 993, exhaustive verification, graph polynomials

## MSC Codes

05C31 (Graph polynomials), 05C05 (Trees), 05-04 (Explicit machine computation and programs)

## Disclosure Statement

The author reports no competing interests.

## Data Availability Statement

All code, data, and reproduction scripts are available at https://github.com/BrettRey/erdos-problem-993. The frozen artifact for this paper is archived at Zenodo: concept DOI https://doi.org/10.5281/zenodo.18745546, version DOI https://doi.org/10.5281/zenodo.19100781.

## AI Disclosure

Claude Opus 4.6 (Anthropic) assisted with computational exploration, Python code development, and prose drafting. ChatGPT 5.2 Pro and Codex 5.3 (OpenAI) assisted with proof-strategy exploration, Lean formalization support, and editorial revision suggestions. Gemini 3 Pro (Google) assisted with computational exploration. Aristotle (Harmonic) was used for bounded Lean proof-completion tasks. The author manually reviewed all machine-generated code, proofs, and prose, and takes responsibility for all claims and any errors.

## Submission URL

https://www.combinatorics.org/ojs/index.php/eljc/about/submissions

Use the E-JC OJS login/register flow from that page.

## Files to Upload

- Manuscript PDF: `paper/main_v2.pdf`

E-JC asks for the article itself as a PDF for the initial submission, and explicitly says not to upload source files at this stage.

## Notes

- E-JC checklist requires that the paper isn't under consideration elsewhere.
- The OJS form asks for an HTML abstract in addition to the manuscript abstract. Use `paper/abstract_submission.txt`; it avoids local macros such as `\dleaf` and `\nm`.
- E-JC has an AI policy. The manuscript's AI disclosure is in the reproducibility section.
- If accepted, prepare the final source with the E-JC style file and consolidate the LaTeX source and bibliography into one file, as requested by the journal.
