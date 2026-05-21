# Literature — minimal citation chain for Erdős #993

Annotated list of the papers a research mathematician would need to cite or rely on when attacking the tree independence-polynomial unimodality conjecture. Dates reflect preprint or publication, whichever is earlier. If a preprint was later published, both are noted.

## The conjecture itself

- **Alavi, Malde, Schwenk, Erdős (1987).** "The vertex independence sequence of a graph is not constrained." *Congressus Numerantium* 58, 15–23. The original conjecture: every tree has a unimodal independence sequence. They also showed that general graphs can realize any prescribed sequence shape, so tree structure is load-bearing.

## Log-concavity failures (the conjecture's stronger form, falsified)

- **Kadrawi, Levit, Yosef, Mizrachi (2023).** "The independence polynomial of trees is not always log-concave starting from order 26." Exhaustive computation: exactly two trees on 26 vertices fail log-concavity. Now published as *Ars Math. Contemp.* 25(4), P4.03. DOI: 10.26493/1855-3974.3207.2ad.
- **Ramos, Sun (2025).** "An AI enhanced approach to the tree unimodality conjecture." arXiv:2510.18826. PatternBoost heuristic search finds tens of thousands of log-concavity counterexamples in the range 27–101 vertices; their publicly shared experimental output is concentrated at $n = 60$. Code: https://github.com/ericgramos/TreeUnimodalityPatternBoost.
- **Galvin (2025).** Infinite families of subdivided stars with log-concavity failures arbitrarily far from the ends of the coefficient sequence.
- **Bautista-Ramos (2025).** Shows that log-concavity can break at an arbitrary number of indices simultaneously.

## Positive results (unimodality proved for sub-classes)

- **Levit (2006).** "On the strong descendant property of tree independence polynomials." $i_k(T)$ is strictly decreasing for $k \geq \lceil (2\alpha(T) - 1)/3 \rceil$. This is the descending-tail result that pairs with Conjecture A to close #993.
- **Basit, Galvin (2020).** Random tree shape asymptotics: almost surely the initial ~49.5% of the sequence is increasing, the terminal ~38.8% is decreasing.
- **Heilman (2025).** "Four-fifths true": ~46.8% increasing-prefix complement to Basit–Galvin.
- **Li et al. (2025).** arXiv:2501.04245. "A symmetric function approach to log-concavity of independence polynomials." Proves log-concavity for all spiders; extends to brooms as a special case.
- **Bendjeddou, Hardiman (2024/2025).** "Lorentzian polynomials and the independence sequences of graphs." arXiv:2405.00511; now in *Bull. LMS* 57 (2025), 1305–1323, DOI: 10.1112/blms.70031. A subclass of trees (those in the image of an edge-replacement operator) have pre-Lorentzian multivariate independence polynomials, hence log-concave sequences.
- **Li (2026).** arXiv:2603.03025. "Unimodality of independence polynomials of two family of trees." Proves unimodality for two specific families including the Kadrawi–Levit LC-failing trees. Directly relevant post-log-concavity work.

## Root-location techniques

- **Chudnovsky, Seymour (2007).** "The roots of the independence polynomial of a clawfree graph." Real-rootedness implies log-concavity. The two LC-failing trees at $n = 26$ have complex roots, so the Chudnovsky–Seymour framework doesn't apply directly to close the conjecture, but the technique is load-bearing when real-rootedness can be established.

## AI-assisted combinatorics precedents (for framing)

- **Nagda, Raghavan, Thakurta (2026).** "Reinforced Generation of Combinatorial Structures: Ramsey Numbers." arXiv:2603.09172. AlphaEvolve-based search improves lower bounds for several classical Ramsey numbers. Cited in Reynolds (2026) as broader context.
- **OpenAI / Sawhney / Alexeev (2026).** "On the ratio of $R(k, \ell)$ and $R(k, \ell+1)$." Internal GPT-5.5 proof that $\lim_{\ell \to \infty} R(k, \ell+1) / R(k, \ell) = 1$ for all fixed $k$. Erdős Problem #1014. Uses dependent random choice. Lean-verified by Alexeev. Proof PDF: https://cdn.openai.com/pdf/6dc7175d-d9e7-4b8d-96b8-48fe5798cd5b/Ramsey.pdf. Lean repo: https://github.com/plby/lean-proofs. Announced 2026-04-23.

## External tools referenced

- **Dependent random choice** — Fox, Sudakov (2011). "Dependent random choice." *Random Structures Algorithms* 38, 68–99. Lemma 2.1 is the standard statement. Textbook reference: Zhao (2023), "Graph theory and additive combinatorics," Theorem 1.7.5.
- **Erdős–Szekeres (1935).** "A combinatorial problem in geometry." *Compositio Math.* 2, 463–470. $R(k, \ell) \leq \binom{k + \ell - 2}{k - 1}$. Baseline upper bound for off-diagonal Ramsey.
- **Bohman, Keevash (2010).** "The early evolution of the H-free process." *Invent. Math.* 181, 291–336. Best known lower bound for off-diagonal Ramsey via the random process.

## Infrastructure / reproducibility

- **Reynolds (2026).** Full manuscript: https://doi.org/10.5281/zenodo.19100781. Github: https://github.com/BrettRey/erdos-problem-993. Reproducibility appendix in `main_v2.pdf`. Concept DOI: 10.5281/zenodo.18745546; version DOI: 10.5281/zenodo.19100781.
- **Tao's AI-contributions wiki.** https://github.com/teorth/erdosproblems/wiki/AI-contributions-to-Erdős-problems. Reynolds (2026) is listed in Section 2(e).
- **Erdős problems website.** https://www.erdosproblems.com/993 (for #993) and https://www.erdosproblems.com/1014 (for the analogous Ramsey ratio result).

## Source grounding — what you can and cannot trust

- Anything with a preprint link above: read the preprint before citing.
- Anything by Reynolds (2026): read `main_v2.pdf`.
- Anything stated as "Conjectured" in `conjecture_and_state.md`: do not cite as proved.
- Bibliography metadata in Reynolds (2026) may be slightly out of date — the review of the manuscript flagged that Kadrawi–Levit and Bendjeddou–Hardiman both have updated published references. Use the DOIs above.
