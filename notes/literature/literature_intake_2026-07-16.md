# Literature intake for Erdős Problem 993

**Date:** 2026-07-16  
**Scope:** tree independence-polynomial unimodality, current structural and spectral methods, and the finite Poisson--binomial effective-drop theorem  
**Cutoff:** arXiv and publisher records checked through 2026-07-16  
**Status:** reviewed intake; selected actions integrated after review (see `literature_followthrough_2026-07-16.md`)

## Executive judgment

The search found no solution of Erdős Problem 993 and no non-unimodal tree. It did, however, change the project in four material ways.

1. Four sources should be integrated promptly: Bautista-Ramos--Guillén-Galván--Gómez-Salgado's now-published recurrence paper, Hibi--Kara--Vien's symmetric-unimodal tree paper, Levit--Kadrawi's adjacent-mode argument, and Jerrum--Patel's February 2026 revision on tree zero sets.
  
2. Jerrum--Patel rule out a universal fixed positive-axis sector for the zeros of all trees. The sector observed in the July computation remains valid as a finite-search fact, and the near-dominant angle--modulus cusp remains open, but the unrestricted sector conjecture does not survive.
  
3. The finite Poisson--binomial theorem
   \[
   V\left(1-\frac{f_{D-1}f_{D+1}}{f_D^2}\right)\ge \frac14
   \]
   at the first strict descent appears new after the documented search, with moderate--high confidence. Two of its ingredients are exact prior results, however: Hillion--Johnson's cubic inequalities and Bobkov--Marsiglietti--Melbourne's max-atom/variance inequality. The novelty lies in their endpoint-aware propagation, the forced modal mass windows, and the variance synthesis.
  
4. The abstract of Ramos--Sun reports examples from 27 to 101 vertices, while their Section 4.1 says early searches found none at 27 or 29 and their public artifacts contain no 27-vertex certificate. Our exact $n=27$ audit finds no log-concavity failure. This should remain explicitly unresolved pending a certificate or author confirmation.
  

My recommended priority order is:

1. repair the probability citations and novelty language;
  
2. update the manuscript's direct 2026 literature and the Ramos--Sun discrepancy paragraph;
  
3. test the Jerrum--Patel and Bautista recurrence families against the certified root Pareto frontier;
  
4. develop the Poisson--binomial theorem as a potentially standalone result;
  
5. only then broaden the background bibliography.
  
## Confidence labels

- **[paper]**: statement checked in the full paper or version of record.
  
- **[source abstract]**: checked only against the publisher/arXiv abstract or metadata.
  
- **[inference]**: a project implication inferred from a source; not claimed by that source.
  
- **[open]**: unresolved, unverified, or dependent on missing artifacts.
  
## Priority 1: sources that change the project now

| Source | Source-grounded result | Project consequence | Action |
| --- | --- | --- | --- |
| [Bautista-Ramos, Guillén-Galván, and Gómez-Salgado, arXiv:2603.14204](https://arxiv.org/abs/2603.14204); [_Graphs and Combinatorics_ 42, article 59](https://doi.org/10.1007/s00373-026-03054-4), published 2026-06-22 | **[paper]** Rooted pattern graphs unify the main non-log-concave tree families. Linear recurrences yield infinite families with one, two, and three consecutive log-concavity breaks, finite families with four and five, and the Beraha--Kahane--Weiss limit circle $\lvert z+1/3\rvert=1/3$ for the basic $P_2$-based reservoirs. | The introduction currently cites arbitrary _nonconsecutive_ breaks from 2025 but misses the stronger 2026 consecutive-break result and the publication. Any heuristic that log-concavity defects must be isolated is now false. The recurrence families are also the best explicit bridge from tree composition to limiting root geometry. | Add the journal article to the bibliography and introduction. Add its families to the exact stress corpus and root-frontier overlay. |
| [Hibi, Kara, and Vien, arXiv:2604.18824](https://arxiv.org/abs/2604.18824) | **[paper]** For every order except $2,4,5,7,10$, there exists a tree with symmetric unimodal independence polynomial; for every degree except 3, such a tree exists. A Bridge Lemma preserves rooted γ-admissibility under a coefficientwise condition. The paper cites Brett's preprint for exhaustive verification through $n\le29$. | This is direct 2026 tree-unimodality literature and a reciprocal citation is warranted. The bridge construction is a real operation-preservation result, though only for a restrictive symmetric/γ-positive safe class. | Add to the introduction and structural discussion. Treat its finite Macaulay2 counts as source-reported unless independently rerun. |
| [Jerrum and Patel, arXiv:2510.01466v2](https://arxiv.org/abs/2510.01466), revised 2026-02-06 | **[paper, with proof-scope caveat]** Lemma 17 states positive-axis accumulation for complete binary trees with every edge subdivided $2k$ times. A post-intake line audit found that the cited normal-family step directly yields phase-truncated periodic balanced trees and appears to shift Buys's $H_\Delta$ index; this proof-supported family is still maximum-degree~3 and suffices for the zero-region conclusion. | This rules out any universal fixed sector $\lvert\arg z\rvert\ge\theta>0$ for all zeros of all trees. It does **not** rule out the narrower claim that a root close in modulus to the minimum root must have small angle from the negative axis. The exact-family formulation merits clarification. | Amend the private root notes and test both the stated and proof-supported families with certified modulus ranking. |
| [Hillion and Johnson, arXiv:1303.3381](https://arxiv.org/abs/1303.3381), _Annals of Probability_ 44 (2016), 276--306 | **[paper]** Theorem A.2, equation (78), and Corollary A.3, equation (79), are exactly the two cubic inequalities that become $\delta_{k\pm1}(1-\delta_k)\le\delta_k$. | The current Poisson--binomial note points to the later [arXiv:1503.01570](https://arxiv.org/abs/1503.01570), which restates the inequalities and attributes them to the 2016 paper. The primary citation should be the 2016 article. | Correct the attribution. It is fine to retain the 2017 source as a secondary restatement. Do not call the cubic inequalities new. |
| [Bobkov, Marsiglietti, and Melbourne, arXiv:2007.11030](https://arxiv.org/abs/2007.11030), [_Combinatorics, Probability and Computing_ 31 (2022), 54--72](https://doi.org/10.1017/S096354832100016X) | **[paper]** Theorem 1.1 and Corollary 3.2 give exactly $M(X)\ge(1+12\operatorname{Var}X)^{-1/2}$, equivalently $V\ge(M^{-2}-1)/12$, using the same uniform-smoothing proof. No log-concavity is needed for this direction. | The max-atom lemma in the project is known, including its proof and sharp constant. This does not subsume the effective-drop theorem. | Cite it as a known input and remove novelty language from that lemma. A short self-contained proof can remain. |
## Direct 2026 tree-unimodality landscape

### The three genuinely direct papers

| Source | Core claim and method | Relationship to the manuscript | Confidence |
| --- | --- | --- | --- |
| [Grace Li, arXiv:2603.03025](https://arxiv.org/abs/2603.03025), _Unimodality of independence polynomials of two family of trees_ | Proves unimodality of both Kadrawi--Levit families $T_{3,m,n}$ and $T^*_{3,m,n}$. The proof establishes a long log-concave prefix using symmetric-function/2-Schur machinery and finishes with universal tail monotonicity. | Already cited accurately. This is the strongest direct template for working below global log-concavity: partial shape control plus an independent tail theorem. It should not be generalized to all known non-log-concave families. | **[paper]** |
| [Bautista-Ramos, Guillén-Galván, and Gómez-Salgado, arXiv:2603.14204](https://arxiv.org/abs/2603.14204) | Gives recurrences, consecutive log-concavity breaks, and limiting zero curves for known pattern families. It does not prove these new families unimodal. | Missing and essential. Its explicit families should replace generic random search as the first structural stress test. | **[paper]** |
| [Hibi, Kara, and Vien, arXiv:2604.18824](https://arxiv.org/abs/2604.18824) | Constructs at least one symmetric unimodal tree in almost every order and degree using γ-admissible rooted bridges. | Missing and essential, but existential rather than universal. It supports a restricted gluing lane, not Erdős 993 itself. | **[paper]** |
### Strongly adjacent 2026 sources

| Source | Core result | Use / limit | Confidence |
| --- | --- | --- | --- |
| [Levit and Kadrawi, arXiv:2603.17114](https://arxiv.org/abs/2603.17114), _Closing Trees into Unicyclic Counterexamples_ | Explicit infinite unicyclic families remain unimodal while the penultimate log-concavity inequality fails. The proof uses a dominant convolution, a real-rooted correction, Ibragimov strong unimodality, Darroch localization, and an adjacent-mode bridge. | The graph class is unicyclic, but its Lemma 2.16 is essentially the same elementary adjacent-mode lemma as the one in `main_v2.tex`. Cite it there as independent nearby literature. The convolution-plus-correction architecture is very close to the subdivision sum. | **[paper]** |
| [Bhardwaj, Chau, Ikram, Lather, Nanjangud, and Venugopal, arXiv:2607.08480](https://arxiv.org/abs/2607.08480), submitted 2026-07-09 | If $\operatorname{mult}_{-1}P_G\ge\alpha(G)-2$, then $P_G$ is log-concave. The paper also characterizes $P_T(-1)=0$ trees using grafting on 3-subdivisions. | This is the freshest relevant preprint and has a real subdivision connection, but the multiplicity condition is far too restrictive for general trees. Add to a positive-class/structural survey, not as a proof route. | **[paper]** |
| [Zhang and Xu, arXiv:2604.01717v2](https://arxiv.org/abs/2604.01717), _On expectations and variances in the hard-core model_ | Gives occupancy and variance bounds, including $V_G(\lambda)\ge \lambda/(1+n\lambda)^2$ and a maximum-degree refinement, where $V_G$ is variance per vertex. | Relevant context for the mean/variance program, but its occupancy upper bound cannot prove $\mu<n/3$ for trees: at $\lambda=1$ and $\alpha\ge n/2$, the bound is at least $n/3$. Contextual citation only. | **[paper]** |
| [Hibi, Kara, and Vien, arXiv:2603.16695](https://arxiv.org/abs/2603.16695), _Independence polynomials of graphs_ | Gives an explicit leaf-attachment formula. If every base vertex receives at least one leaf, symmetry holds exactly when every vertex receives two leaves; it also treats big stars, caterpillars, and $P_G(-1)$. | Neighbors the manuscript's leaf-attachment asymptotics and provides useful formula/symmetry context, but does not subsume them. | **[paper]** |
| [Pham and Vu, arXiv:2604.10269](https://arxiv.org/abs/2604.10269), _Contractible independence complexes of trees_ | Characterizes $I(T;-1)\in{-1,0,1}$ by branch-truncation reductions. | Directly about trees but only about evaluation at one fixed point, not coefficient shape. Screened background; optional alongside Bhardwaj et al. | **[paper]** |
| [Pereyra, arXiv:2605.14076](https://arxiv.org/abs/2605.14076) | Proves the 2-quasi-regularizability conjecture for $W_2$ graphs and gives coefficient regions for log-concavity/unimodality. | The connected-tree intersection with the relevant $W_2$ class is essentially trivial beyond $K_2$. Do not spend manuscript space here. | **[paper]** |
## Late-2025 direct literature that should not be lost

| Source | Contribution | Decision |
| --- | --- | --- |
| [Ramos and Sun, arXiv:2510.18826](https://arxiv.org/abs/2510.18826) | PatternBoost found tens of thousands of non-log-concave trees and no non-unimodal example; the linked public corpus is at order 60. | Already central. Revise only the $n=27$ discrepancy wording, below. |
| [Bautista-Ramos, arXiv:2511.00334](https://arxiv.org/abs/2511.00334) | Constructs trees with arbitrarily many log-concavity breaks, not necessarily consecutive. | Already cited. Place the 2026 consecutive-break sequel beside it. |
| [Galvin, arXiv:2502.10654v2](https://arxiv.org/abs/2502.10654) | Constructs late-tail log-concavity failures; v2 was revised in January 2026. | Already cited. Check bibliography/version metadata when updating. |
| [Li, Li, Yang, and Zhang, arXiv:2501.04245](https://arxiv.org/abs/2501.04245) | Equates global 2-Schur positivity with log-concavity and proves spiders and pineapple graphs log-concave. | Already cited for spiders. Its main strategic lesson is that global 2-Schur positivity is not a weaker route around failed log-concavity; Li 2026 succeeds by proving only the needed prefix. |
| [Liu, Tang, and Zhao, DOI 10.1007/s10255-025-0082-x](https://doi.org/10.1007/s10255-025-0082-x), _Trees with Independence Polynomials Having Only Real Zeros_ | Constructs additional infinite tree/composite families with real-rooted independence polynomials. | Add only if the introduction aims for a fuller positive-family survey. Full publisher text was not acquired in this intake. **[source abstract]** |
| [Xie and Zhang, DOI 10.3969/j.issn.1006-8074.2025.02.004](https://mta.csu.edu.cn/EN/10.3969/j.issn.1006-8074.2025.02.004), _Unimodality of Independence Polynomials of Rooted Products of a Kind of Tree_ | Proves real-rootedness preservation for specified rooted products, hence log-concavity and unimodality for resulting tree families. | Narrow but direct. Add only in a fuller family survey. The publisher PDF endpoint returned HTTP 403 during intake; the abstract and metadata were accessible. **[source abstract]** |
## Root geometry and the July spectral program

### Ranked sources

1. **Prakash--Sharma:** [arXiv:2510.09197](https://arxiv.org/abs/2510.09197). **[paper]**  
  Their Theorem 1.1 quantitatively isolates the unique minimum-modulus root of a connected graph, but only by an exponentially small additive gap. More promising for this project are the local-univalence argument and recursively monotone angular majorants for occupation ratios. Paths show that a uniform radial collar cannot hold: the second-root ratio can approach 1 at polynomial rate while all relevant roots remain real. This is consistent with, but does not prove, the surviving cusp claim that near-equal modulus forces vanishing non-real angle.
  
2. **Jerrum--Patel:** [arXiv:2510.01466v2](https://arxiv.org/abs/2510.01466). **[paper, with proof-scope caveat]**  
  Lemma 17 states a bounded-degree tree family with zeros tending to the positive real axis; the proof directly supplies phase-truncated periodic balanced trees, which are enough for that conclusion but are not plainly the exact uniformly subdivided family in the lemma statement. For the unsplit binary recurrence the positive neutral point is $+4$, while the binary-tree negative threshold is $-4/27$; thus the positive-axis mechanism is remote from minimum modulus. This kills the global sector but leaves the near-dominant tradeoff open.
  
3. **Bautista-Ramos--Guillén-Galván--Gómez-Salgado:** [arXiv:2603.14204](https://arxiv.org/abs/2603.14204). **[paper]**  
  The circle $\lvert z+1/3\rvert=1/3$ is a limiting zero locus, not a ranking of roots by modulus. Its families need an explicit ranked-root overlay before they can support or refute the empirical cusp envelope.
  
4. **Foundational occupation-ratio and bounded-degree work:** [Scott--Sokal](https://arxiv.org/abs/math/0309352), [Bencs--Buys--Peters](https://arxiv.org/abs/2111.06451), [Bencs--Csikvári--Srivastava--Vondrák](https://arxiv.org/abs/2204.04868), and [Hlushchanka--Peters](https://arxiv.org/abs/2411.14791). **[paper]**  
  These papers rigorously connect recursive gluing, occupation ratios, rational dynamics, and absolute zero-free regions. They do not compare the arguments of roots at nearly equal moduli within one tree polynomial. They supply the recursion technology, not the needed relative spectral theorem.
  
5. **Special algebraic safe classes:** [Bendjeddou--Hardiman](https://arxiv.org/abs/2405.00511) and [Hibi--Kara--Vien](https://arxiv.org/abs/2604.18824). **[paper]**  
  Pre-Lorentzian gluing and γ-admissible bridging are the best existing operation-preservation technologies to test on restricted subdivision blocks. Neither proves that ordinary single-edge subdivision preserves unimodality.
  
### Root-space contradiction table

| Project reading | Literature verdict | Revised position |
| --- | --- | --- |
| Every tree root avoids a fixed positive-axis sector | False by Jerrum--Patel Lemma 17. | Keep only the observed $0.90284$ minimum over the certified finite sample; any theorem needs degree/modulus/critical-activity restrictions. |
| Every tree has a fixed modulus collar above its smallest root | False already for paths; also killed computationally in the project. | Couple modulus ratio to angle. A cusp envelope remains viable. |
| Existing bounded-degree zero-free regions prove the cusp | No. They locate absolute regions in the activity plane, not relative root order within one polynomial. | Try to specialize Prakash--Sharma's angular majorants to tree occupation messages. |
| The Bautista limit circle identifies dangerous roots | Not without ordering the roots by modulus. | Compute the root ranking for its explicit recurrence families. |
| Known gluing theorems imply arbitrary subdivision preservation | No. They cover special symmetric, pre-Lorentzian, or fixed-block classes. | Use them as certificate templates for restricted subclasses only. |
## Poisson--binomial novelty audit

### Verdict

No exact theorem overlap was found for the finite first-descent bound

\[
V\delta_D
=
V\left(1-\frac{f_{D-1}f_{D+1}}{f_D^2}\right)
\ge \frac14.
\]

No screened source appears to imply it by routine substitutions. The safe novelty claim is:

> We prove a variance-normalized modal curvature bound for finite Poisson--binomial laws. The argument combines Hillion--Johnson's cubic coefficient inequalities with a known lattice max-atom/variance inequality and a new propagation-to-variance synthesis.

Avoid claiming a new cubic inequality, a new max-atom bound, or an unqualified historical first. A reasonable sentence is: "We are not aware of a previous finite bound at the first strict descent."
### Exact prior ingredients

| Source | Exact overlap | What remains new |
| --- | --- | --- |
| [Hillion--Johnson 2016](https://arxiv.org/abs/1303.3381), Theorem A.2 and Corollary A.3, equations (78)--(79) | The two cubic inequalities become exactly $\delta_{k-1}(1-\delta_k)\le\delta_k$ and $\delta_{k+1}(1-\delta_k)\le\delta_k$. | Iterating them with endpoint values, turning them into two-sided modal mass windows, and combining the windows with variance. |
| [Bobkov--Marsiglietti--Melbourne 2022](https://arxiv.org/abs/2007.11030), Theorem 1.1 / Corollary 3.2 | Exactly $V\ge(M^{-2}-1)/12$, including the uniform-smoothing proof. | Its use together with the propagated mass window and the exact scalar optimization. |
| [Darroch 1964](https://doi.org/10.1214/aoms/1177703287) | One or two adjacent modes and $\lvert m-\mu\rvert<1$. | Any variance-scale lower bound on the local descent. |
| [Pitman 1997](https://statistics.berkeley.edu/tech-reports/453) | Tilted-mean bounds for consecutive coefficient ratios of real-rooted polynomials. | A bound on the _difference_ of adjacent ratios using the original variance. |
| [Johnson, arXiv:1507.06268](https://arxiv.org/abs/1507.06268) | $c$-log-concavity and discrete curvature; connects monotonic curvature to Hillion--Johnson's $C_1(k)$. | The odds-based curvature constant can be arbitrarily smaller than $1/V$; it does not yield the target. |
| [Dümbgen--Wellner, arXiv:1910.03444](https://arxiv.org/abs/1910.03444), Proposition 2 | Strict decrease of $(k+1)f_{k+1}/f_k$. | No quantitative variance-normalized size of the decrease. |
### Nearest 2026 probability papers

| Source | Why it is nearby | Why it does not overlap |
| --- | --- | --- |
| [Kontorovich, arXiv:2601.04079](https://arxiv.org/abs/2601.04079) | Poisson--binomial laws, variance, unimodality, max-atom bounds, and total variation. | No neighboring-ratio or modal-curvature estimate. |
| [Broadie--Petkova, arXiv:2603.09019](https://arxiv.org/abs/2603.09019) | Splits Poisson-trinomial laws into two interleaved conditional Poisson--binomial parts; uses Hermite--Biehler, real-rootedness, and Darroch. | Locates modes and proves log-concavity of parts but gives no drop magnitude. |
| [Kurauskas, arXiv:2603.11043](https://arxiv.org/abs/2603.11043) | Asymptotically optimal concentration-function comparison for sums of independent integer variables. | Asymptotic max-atom result, no adjacent ratios or finite curvature. |
| [Kovačević, arXiv:2605.11831](https://arxiv.org/abs/2605.11831) | Ternary sums, parity-conditioned ultra-log-concavity, Hermite--Biehler, Newton inequalities. | Entropy optimization rather than local modal curvature. |
| [Atminas--Kurauskas, arXiv:2510.07899v2](https://arxiv.org/abs/2510.07899) | 2026 revision proving variance-minimizing rearrangements on the integers. | No Poisson--binomial adjacent-ratio result. |
| [Gaxiola--Melbourne--Pigno--Pollard, arXiv:2505.05793](https://arxiv.org/abs/2505.05793) | Very recent sharp mode-atom/variance inequalities for log-concave integer laws. | Still controls atoms rather than the normalized three-mass deficit. |

The current local-limit literature is too coarse: central adjacent-mass cancellations occur at the same order as the available approximation error, and local limits do not identify the exact first strict descent.
## Ramos--Sun $n=27$ discrepancy

The present manuscript is right to flag the conflict, but the alternatives can now be narrowed.

- **[paper]** Ramos--Sun define $N$ as the number of vertices and use the standard independence polynomial, so a different size convention or polynomial is not a well-supported explanation.
  
- **[paper]** Their abstract reports counterexamples with vertex counts from 27 to 101.
  
- **[paper]** Their Section 4.1 says the early attempts found examples at $N=26,28,30$, found none at 27 or 29, and still had no success at odd $N$ for the $N/2-1$ target.
  
- **[paper]** Their appendix and linked public corpus do not supply a 27-vertex certificate; the public corpus is at 60 vertices.
  
- **[project computation]** The retained exact $n=27$ artifact reports 751,065,460 trees and zero log-concavity failures.
  

Recommended wording:

> Ramos and Sun report in their abstract that their non-log-concave examples range from 27 to 101 vertices. Our exhaustive enumeration finds no log-concavity failure on 27 vertices. Their Section 4.1 states that early searches did not find failures at $N=27$ or $29$, and neither the paper's appendix nor its linked 60-vertex corpus supplies a 27-vertex certificate. We therefore leave the abstract's lower endpoint unresolved: it may refer to a later example not included in the public artifacts, or it may be a typographical error. A 27-vertex certificate or confirmation from the authors would settle the discrepancy.

Do not call 27 a typo without author confirmation.
## Other active work, not formal literature

The [Erdős Problems #993 discussion](https://www.erdosproblems.com/forum/thread/993) now includes a June 2026 report by Will Blair, with [public code and witnesses](https://github.com/willblair0708/verified-combinatorics/tree/main/erdos-993). It reports 4,445 non-log-concave rooted-bush tree polynomials up to 60 vertices and 253,695 products/powers/path-products of the strongest seeds, all unimodal. The upstream verifier and all 25 publicly materialized rooted bushes have now been independently replayed; see `blair_rooted_bush_comparison_2026-07-16.md`. The remaining 4,420 claimed bush polynomials are not serialized in the public result. This is not peer-reviewed literature, but it is a relevant parallel forest/product search. **[source report; partial independent replay]**

The corpus comparison is now complete. Blair's V2 grammar is wholly contained in the local spider-bouquet grammar. Of the 25 public bushes, one is an exact Galvin duplicate and 24 are new exact coefficient/tree records relative to the bounded index; five improve the available same-order pressure comparator, but none extends the order-at-most Pareto frontier. The duplicate-enumeration gate is therefore closed. The bounded next step, only if the missing 4,420 records become available, is archival ingestion rather than another broad search.
## Search method and coverage

This was not a flat keyword search. The intake used four passes:

1. **Direct 2026 arXiv sweep.** An arXiv API query for exact phrase `"independence polynomial"` in submissions from 2026-01-01 through 2026-07-16 returned 27 records. Every mathematically plausible graph-theory result was screened; full texts were acquired for the direct and structurally adjacent papers.
  
2. **Variant searches.** Queries covered `independent set sequence`, tree unimodality/log-concavity, independence complexes of trees, hard-core expectation/variance, zero-free regions, occupation ratios, Poisson--binomial/Bernoulli sums, ultra-log-concavity, adjacent ratios, Turán deficits, max atoms, and local limits.
  
3. **Citation chaining.** Citations were followed backward and forward from Li 2026, Bautista et al. 2026, Hibi et al. 2026, Prakash--Sharma, Jerrum--Patel, Hillion--Johnson, Johnson, Pitman, and Bobkov--Marsiglietti--Melbourne.
  
4. **Project collision check.** The strongest papers were compared against `paper/main_v2.tex`, the current proof-status notes, the July root-conjecture notes, and the finite Poisson--binomial theorem note.
  

Limitations:

- The subsequent multi-database audit includes a reproducible zbMATH query layer, but no authenticated MathSciNet or proprietary-index session was available.
  
- Search engines can miss differently worded coefficient inequalities and non-indexed journals.
  
- The cutoff is 2026-07-16; the July arXiv feed is still moving.
  
- Source-reported computations were not treated as independently verified project computations.
  
## Local acquisition inventory

The full-text corpus is in `notes/literature/`. New PDF/text pairs acquired during this intake include:

- `arxiv_2510_01466` -- Jerrum--Patel;
  
- `arxiv_2510_18826` -- Ramos--Sun;
  
- `arxiv_2603_14204` -- Bautista-Ramos--Guillén-Galván--Gómez-Salgado;
  
- `arxiv_2603_16695` -- Hibi--Kara--Vien, _Independence polynomials of graphs_;
  
- `arxiv_2603_17114` -- Levit--Kadrawi;
  
- `arxiv_2604_01717` -- Zhang--Xu;
  
- `arxiv_2604_10269` -- Pham--Vu;
  
- `arxiv_2605_14076` -- Pereyra;
  
- `arxiv_2607_08480` -- Bhardwaj et al.;
  
- `arxiv_1303_3381` and `arxiv_1503_01570` -- Hillion--Johnson source pair;
  
- `arxiv_2007_11030` -- Bobkov--Marsiglietti--Melbourne;
  
- `arxiv_2601_04079`, `arxiv_2603_09019`, `arxiv_2603_11043`, and `arxiv_2605_11831` -- nearest 2026 probability papers;
  
- `arxiv_1507_06268`, `arxiv_2408_06477`, `arxiv_2505_05793`, and `pitman_1997_coefficients_real_zeros` -- probability curvature, modal, and adjacent-ratio sources;
  
- `arxiv_2411_14791`, `arxiv_2510_09197`, `arxiv_2111_06451`, `arxiv_2204_04868`, and Scott--Sokal -- spectral/occupation-ratio sources.
  

PDFs are ignored by the repository's Git rules. The retained text conversions and audit notes provide the reviewable source layer; unrelated project files were not changed by this follow-up.
### Access still useful

- Liu--Tang--Zhao 2025: publisher full text was not acquired.
  
- Xie--Zhang 2025: the publisher PDF endpoint returned HTTP 403, though the abstract and bibliographic metadata were accessible.
  
- The publication-grade public-database and citation-chain audit is recorded in `poisson_binomial_novelty_database_audit_2026-07-16.md`. Only an authenticated MathSciNet pass remains useful immediately before submission, especially if wording stronger than “we are not aware” is contemplated.
  
## Recommended next actions after review

1. **Bibliographic repair, small and high value.** Add Bautista et al. 2026, Hibi et al. 2604.18824, Levit--Kadrawi, and the relevant fixed-root/subdivision source; update Bautista et al. to the journal version; cite the adjacent-mode lemma.
  
2. **Probability provenance repair.** Change the cubic primary source to Hillion--Johnson 2016 and attribute the max-atom inequality to Bobkov--Marsiglietti--Melbourne. Reframe the theorem's novelty around the synthesis.
  
3. **Root-note correction.** Add Jerrum--Patel as a decisive counterexample to an unrestricted positive-axis sector; retain the near-dominant cusp as the live conjecture.
  
4. **Two discriminating computations.** Run the Jerrum--Patel subdivided-binary family and Bautista recurrence families through certified root isolation and record modulus ranking plus angle.
  
5. **Manuscript wording.** Replace the Ramos--Sun discrepancy paragraph with the source-specific version above.
  
6. **Publication triage.** Consider a standalone Poisson--binomial note before burying the universal $1/4$ theorem inside the tree paper. A cautious expert check with Hillion, Johnson, or Melbourne would be high leverage before submission.
  
7. **Optional breadth.** Add the 2025 real-rooted tree-family papers only if the introduction is being expanded into a fuller survey; they are not needed for the main argumentative spine.
