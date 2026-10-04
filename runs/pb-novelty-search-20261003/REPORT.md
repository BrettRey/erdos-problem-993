# Free novelty and delta search (zbMATH Open, arXiv), 2026-10-03
<!-- SUMMARY: unauthenticated novelty sweep for Theorem 1.1 (V·δ_D ≥ 1/4); no document found that states or implies it; MathSciNet pass still not done · status: bounded negative · updated: 2026-10-03 -->

**Why this was run.** The venue decision record (2026-07-16) requires an
authenticated MathSciNet pass and a submission-date delta search before
submission. A side agent pointed out on 2026-10-03 that neither had been done,
and that the portal record wrongly showed MathSciNet as settled. MathSciNet
needs an institutional login. This sweep uses the free indexes instead: zbMATH
Open, which indexes the same journal literature with its own reviews, and
arXiv, for preprints and for the delta since 1 July 2026. It does not
replace the MathSciNet pass.

**Scripts and raw data (this folder).**
- `zbmath_search.py`: 17 queries to `api.zbmath.org/v1/document/_search`, with 100 results per query. Raw output in `zbmath_raw.json`, deduplicated hits in `zbmath_hits.json`.
- `arxiv_search.py`: 13 general queries over all dates, plus 7 delta queries filtered to items submitted on or after 2026-07-01. Up to 100 results per query, newest first. Output in `arxiv_results.json`.

## Queries and hit counts

**zbMATH Open** (free-text search_string):

| Hits | Query |
|---:|---|
| 18 | Poisson binomial log-concave |
| 2 | Poisson binomial Turan |
| 9 | Poisson binomial variance mode |
| 3 | Poisson binomial ultra log-concave |
| 7 | Bernoulli sums log-concave |
| 1 | sums of independent Bernoulli log-concavity |
| 0 | independent Bernoulli ultra log-concave |
| 1 | strongly log-concave discrete Bernoulli |
| 0 | Newton inequalities improvement real-rooted |
| 3 | Newton inequalities refinement coefficients |
| 0 | real-rooted polynomial coefficients log-concavity variance |
| 1 | Turan inequalities real-rooted polynomial coefficients |
| 5 | Shepp-Olkin conjecture |
| 1 | discrete log-Sobolev Bakry-Emery log-concave |
| 0 | c-log-concave discrete curvature |
| 0 | Poisson binomial mode ratio bound |
| 0 | log-concave sequence quantitative curvature variance |

That is 45 distinct documents.

**arXiv** (export API). General queries, over all dates:

| Hits | Query |
|---:|---|
| 8 | `all:"Poisson binomial" AND all:"log-concave"` |
| 0 | `all:"Poisson binomial" AND all:Turan` |
| 8 | `all:"Poisson-binomial" AND all:"log-concavity"` |
| 0 | `all:"sums of independent Bernoulli" AND all:"log-concave"` |
| 0 | `all:"Bernoulli sums" AND all:"log-concave"` |
| 0 | `all:"real-rooted" AND all:Turan AND all:variance` |
| 1 | `all:"Newton inequalities" AND all:variance` |
| 1 | `all:"Newton's inequalities" AND all:improvement` |
| 0 | `all:"ultra log-concave" AND all:variance AND all:mode` |
| 2 | `all:"strongly log-concave" AND all:Bernoulli` |
| 0 | `all:"discrete curvature" AND all:"log-concave"` |
| 0 | `all:"c-log-concave"` |
| 7 | `all:"Shepp-Olkin"` |

Delta queries, submitted on or after 2026-07-01:

| Hits | Query |
|---:|---|
| 2 | `all:"Poisson binomial"` |
| 2 | `all:"Poisson-binomial"` |
| 4 | `all:"log-concave" AND all:Bernoulli` |
| 6 | `all:"ultra-log-concave"` |
| 6 | `all:"ultra log-concave"` |
| 6 | `all:"real-rooted" AND all:probability` |
| 1 | `all:Turan AND all:"log-concave"` |

## Screening

All titles were read. Abstracts or reviews were read for the candidates
below; every other hit is off-topic or already cited.

- **Guo, "An inequality for coefficients of the real-rooted polynomials", J. Number Theory (2021), pp. 294–309; arXiv 2012.03530.**
  - For $f(x)=\sum\binom nk a_kx^k$ with real zeros only, Guo bounds the ratio $(a_{k+1}^2-a_ka_{k+2})/(a_k^2-a_{k-1}a_{k+1})$ above and below in terms of $c_k=a_ka_{k+2}/a_{k+1}^2$.
  - That controls *adjacent* Turán differences, which is the role the Hillion–Johnson cubics play in the paper. It gives no variance-scaled lower bound on a single deficit, so it does not imply Theorem 1.1.
  - It is citable as related work if Brett wants it. It is not needed.
- **Chen and Wang, "High-Dimensional Ultra-Log-Concave Distributions", arXiv 2609.23994 (21 Sep 2026).**
  - They define δ-ultra-log-concavity on $\mathbb N^d$ and prove Poincaré, Brascamp–Lieb and modified log-Sobolev inequalities, with concentration.
  - The paper shares the motivation (functional inequalities from quantitative ULC), but the summary states no variance-scaled curvature bound at a descent or mode for Bernoulli sums. Proposition 1.5 of the paper concerns classical ULC.
  - It is a recent neighbour and might be cited beside Johnson's c-log-concavity. That is Brett's choice; it doesn't affect novelty.
- **Jaramillo and Melbourne (2022), "Quantitative limit theorems via relative log-concavity".** This is the same research line as the cited Marsiglietti–Melbourne paper. It gives limit theorems, not a deficit lower bound.
- **Wang, arXiv 2609.27433 (23 Sep 2026).** Rényi and Tsallis thresholds in the Shepp–Olkin problem: entropy concavity in the parameters, not deficits.
- **Kunisky–Spielman–Yu, arXiv 2608.28520.** Permanent inequalities and finite free convolutions; not relevant.
- **Banerjee–Bringmann–Rolen, arXiv 2609.14295.** Turán and Laguerre inequalities for partition-type sequences; not relevant.
- **Sharma–Sharma–Saini, arXiv 1702.04665.** Refinements of AM–GM through Newton's inequalities in terms of variance, for positive reals; not about coefficient deficits.
- **Rosenman, arXiv 1907.09053.** Poisson–binomial likelihoods in ecological inference; not relevant.

## Conclusion (bounded)

In these 37 queries, no document states or implies a lower bound of order
$1/V$ on the normalized Turán deficit of a Poisson–binomial law at its first
descent or its mode. The statement "we are not aware of an earlier proof" is
consistent with this sweep.

This is not an exhaustive negative. zbMATH's free-text search and arXiv's
phrase search are coarser than MathSciNet's reviewed index, and Google Scholar
was not searched (it is bot-walled for agents). **The authenticated MathSciNet
pass the venue record requires is still not done.**
