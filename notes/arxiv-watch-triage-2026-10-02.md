# arXiv watch triage, 2026-10-02
<!-- SUMMARY: The 17 papers in the arxiv-watch digest for 29 Sep to 2 Oct, read for #993 relevance: 1 relevant (Liu & Tang 2609.37553, a second Lorentzian edge-replacement tree family, replayed: no LC failures), 8 peripheral, 8 vocabulary coincidences; found the watch's 2-day window dropping Thu-afternoon to early-Sunday submissions (Fang et al. 2609.20961 among them), since fixed (5-day default, weekly sweep) with the backlog triaged in arxiv-watch-backlog-triage-2026-10-02.md · status: done · updated: 2026-10-02 -->

Covers the four digest entries in `Project-Management/arxiv-watch/digest.md`
dated 29 Sep, 30 Sep, 1 Oct and 2 Oct. All 17 PDFs were fetched from
arxiv.org. The two papers about independence polynomials (2609.37553 and
2609.35102) were read in full; every other paper was read through its first
pages, with a full-text grep for `independence polynomial`, `independent set`,
`tree`, `unimodal` and `log-concav`. No verdict below comes from the abstract
alone.

The `erdos993-depth3` profile matched none of the 17. Every hit came through
`portfolio-generic` terms.

## Relevant: 1

- **arXiv:2609.37553v1, Liu & Tang, "Lorentzian polynomials and
  log-concavity of the independence polynomials of graphs"** (29 Sep).
  This extends Bendjeddou–Hardiman. They define coloured graphs
  F(l,m,t,s): a centre c adjacent to a free vertex x_i, with cliques hung off
  both (Definition 3.1, Fig. 5). The operator E_{G4(l,m,t,s)} replaces
  every edge uv of a graph G by a copy of F at c_u and one at c_v, with the
  two free vertices joined. Using Bendjeddou–Hardiman's pre-Lorentzian gluing
  lemma (their Lemma 2.1), they prove the images have log-concave
  independence polynomials (Theorem 3.3), with no condition at all when G
  has no vertex of degree 2 (Theorem 3.2, via Remark 3.3).
  - **Tree cases.** E_{G4(1,0,0,0)} is Bendjeddou–Hardiman's W4 operator
    (their Corollary 3.2). The new tree family is E_{G4(2,1,0,0)} applied
    to forests (Corollary 3.3 and the closing paragraph). In that family,
    for each incident edge, every centre gets one pendant leaf and one
    pendant P2, and so does every free vertex. Definition 3.1 forces these
    two to be the only tree cases. If s or t is positive, an edge from c
    into a clique that x_i also reaches closes a cycle through the edge
    c–x_i, and a clique of size 3 or more contains a triangle. F(1,1,0,0)
    is F(1,0,0,0) plus two isolated vertices.
  - **Replay** (`scripts/verify_liu_tang_2609_37553_20261002.py`, exact
    integer arithmetic). The construction was checked against their
    closed form for I(F_n), and against brute force on small cases. I then
    took every tree T on 2 to 12 vertices (986 trees, from geng) and
    computed I(E_{G4}(T)) for both families. Results: **0 LC failures and
    0 unimodality failures** in either family. Image orders were 8 to 78
    for (1,0,0,0) and 16 to 166 for (2,1,0,0). The smallest
    a_k^2/(a_{k-1}a_{k+1}) was 1.0992 for (1,0,0,0) and 1.0490 for
    (2,1,0,0), both at |T| = 12.
  - **Not real-rooted, so LC is not inherited.** An exact Sturm count gave
    real-rooted images only when T was a path (max degree at most 2): 11 of
    986 in each family. Every other image has non-real zeros. On the
    K_{1,3} input, certified Arb isolation (python-flint) agrees, giving 4
    certified non-real roots in each family. So LC for these trees doesn't
    follow from the real-rootedness theorems for trees (Zhu et al.; Liu,
    Tang and Zhao). I didn't check whether the other routes they list, Zhu
    and Chen's LC results (their [29]) and the symmetric-function method of
    Li et al., already cover either family.
  - **Fit with the manuscript.** It is consistent with the Schweitzer
    boundary paragraph in `paper/main_v2.tex` (lines 96–99). The family is
    another special edge-replacement subclass, not a general tree theorem.
    Image orders are 7v−6 and 15v−14 for |T| = v, so neither family has a
    tree of order 26, and the Kadrawi–Levit LC failures lie outside both
    trivially. They cite Bendjeddou–Hardiman, Bencs 2018, Galvin–Hilyard and
    Kadrawi–Levit. They do not cite Schweitzer 2608.23262: `grep -c -i
    schweitzer` on the full text returns 0.
  - **Not checked.** Their Hessian computations (Lemma 3.1 and the n > 2
    case of Theorem 3.1) were not replayed. The non-tree graphs in Theorem
    3.2 were not tested. The replay tests their conclusion on trees, not
    their proof.
  - **If the paper is reactivated**, the Bendjeddou–Hardiman sentence at
    `main_v2.tex:96` could add Liu–Tang as a second edge-replacement family.
    Parked, so no manuscript change now. Filed to
    `literature/liu_tang_2026_lorentzian_independence_polynomials.{pdf,md}`
    (PDF SHA-256
    `8db929f6c0810672bdef13e750805ad8086d5f8982c746e5d68d741ab47e34c9`).

## Peripheral: method-ecosystem awareness, no action

- **arXiv:2609.35102, Hlushchanka & Peters, "Zeros of the independence
  polynomial on recursive sequences of graphs."** Take a recursion with
  uniformly bounded degrees and diverging distances between labelled
  vertices. Their Theorem 1.1 shows the zeros of Z_{G_n} avoid a uniform
  neighbourhood of [0, ∞), so the hard-core model has no phase transition.
  Their Theorem 1.4 adds bounded zeros and a zero-free cone when G_0 is
  "maximally independent". They use complex dynamics of a renormalization
  map on P^{2^k−1}. Both theorems constrain zeros near the positive axis,
  where the free energy lives. Neither constrains the shape of the
  coefficient sequence. Their examples are Sierpiński tripods and
  hierarchical lattices. The d-ary rooted trees fall outside their framework:
  §1.5.2 says the recursion that produces them isn't of the kind Theorem 1.1
  covers. The
  introduction also points to Jerrum–Patel [JP26] on zero-free regions for
  H-free bounded-degree graphs, with H a subdivided claw. Inward check:
  `grep -i -e peters -e bencs -e hlushchanka gpt_attack/literature.md`
  returns nothing, and `bencs2018` is already in `paper/references.bib`.
- **arXiv:2609.39917, Leake & Mohammadi Yekta, "Log-concavity and
  approximate counting for totally unimodular polytopes."** Their Theorem
  1.4: for unimodular M, the fibre counts #{y ≥ 0 : My = b} are log-concave
  along every line in b. Corollary 1.5 gives log-concave Ehrhart
  evaluations for unimodular polytopes. A tree's stable-set polytope is TU
  (bipartite incidence), so the theorem makes E_{STAB(T)}(k) log-concave in
  the dilation k. That is not the sequence i_k(T). Getting i_k needs an
  extra all-ones size row. If that augmented system met their hypothesis
  for every tree, Theorem 1.4 would make every tree LC, which is false at
  n = 26. So for the Kadrawi–Levit trees it fails the hypothesis. This is
  an inference from their theorem; I didn't compute it. Their AI
  declaration (§1.4) says Theorem 1.4 was "the main AI input" and came from
  ChatGPT 6 Astra.
- **arXiv:2609.37439, Fu**, rank-two matroid Ehrhart h*-polynomials.
  Disproves Ferroni's real-rootedness conjecture with parallel-class sizes
  (4, 561, 600), while proving ULC for every matroid of rank or corank two.
  It's the "real-rootedness fails, ULC survives" pattern again. The objects
  are matroid base polytopes, so nothing transfers.
- **arXiv:2609.40084, Xie & Zhang**, unimodality of Kazhdan–Lusztig
  polynomials of sparse paving matroids. Matroid-side, nothing transfers.
- **arXiv:2609.33651, Gao & Yuan**, log-concavity of external
  semi-activity polynomials of real flat arrangements, via mixed volumes and
  Alexandrov–Fenchel. Applications: spanning-tree polynomials of Eulerian
  digraphs and Alexander polynomials. It's the Hodge/mixed-volume pattern
  applied to graphic and cographic objects, not independence sequences.
- **arXiv:2609.36818, F. Liu, Tang, Tao & Zhang**, Ehrhart coefficients of
  Hermite normal form simplices. Unimodality, plus a classification of the
  LC and real-rooted cases. Lattice-polytope side, nothing transfers.
- **arXiv:2609.38084, Alexandersson & Leite**, peak and peak-nesting
  polynomials over unit interval graphs. These are distributions of a
  Dyck-path statistic over a graph class, not graph polynomials of single
  graphs. The full text never mentions independence polynomials (grep count
  0).
- **arXiv:2609.33732, Chern & Shi**, strong q-log-convexity of d-Hoggatt
  polynomials via Schur positivity. Symmetric-function positivity, same
  standing as the Schur papers in the 26 Aug triage.

## Vocabulary coincidences: 8

- **2609.33587, Guo et al.** (Hopf algebra of the quantum Magnus expansion).
  The trees are directed trees of the Magnus expansion, and the "chromatic
  polynomial" hit has no coefficient-shape content.
- **2609.38567, Hylock, Lettington & Schmidt** (c-polynomial
  factorisations). The chromatic polynomial of a path only appears as a
  counting formula.
- **2609.38295, Mikulincer & Zadik** (KLS for unconditional log-concave
  measures). Log-concave means a density on R^n. The abstract says the core
  ideas "were generated by AI generative tools".
- **2609.38904, Friedland, Guédon & Souli** (Gaussian partition functions).
  Log-concave profiles on R.
- **2609.39252, Gozlan, Malamut & Waldspurger** (functional inequalities
  along Wasserstein geodesics). Strongly log-concave measures.
- **2609.39505, Qu & Yang** (algebraic integral points on curves). "Totally
  positive" refers to algebraic integers.
- **2609.35404, Roy Choudhury & Shuddhodan** (Artin vanishing on abelian
  varieties). Hard Lefschetz is only mentioned in passing.
- **2610.01049, Sutherland** (Watkins's conjecture for infinite groups).
  Cayley graphs and GRRs, no independence-sequence content.

## Watch coverage hole: Thursday 14:00 ET to early Sunday submissions are dropped

`arxiv_watch.py` keeps an entry only if its `published` time (v1
submission) falls inside the last `--days` of the run time (`_parse_feed`,
lines 349–365). The LaunchAgent runs daily at 06:15 with `--days 2`
(`~/Library/LaunchAgents/com.brettreynolds.arxiv-watch.plist`). arXiv's
schedule (info.arxiv.org/help/availability.html, fetched 2026-10-02) is:
Thu 14:00 to Fri 14:00 ET is announced Sun 20:00, and Fri 14:00 to Mon
14:00 is announced Mon 20:00. A paper only reaches the API after its
announcement. The Monday run therefore sees Thursday-afternoon and Friday
submissions only once they are already past its Saturday-morning cutoff.
The Tuesday run sees Friday-afternoon to early-Sunday submissions only once
they are past its Sunday-morning cutoff. Neither run keeps them. The hole
runs from Thursday 14:00 ET to about Sunday 06:00 ET, roughly 2.7 of every
7 days of submission time, for every profile, since the watch started on 14
Aug. The seen-set dedup would make a wider window harmless.

The digest agrees. It has never had a Monday entry: 17, 24 and 31 Aug and
14, 21 and 28 Sep are all absent, though the agent runs daily. That includes
a Saturday entry on 29 Aug, and the log shows two zero-retrieval runs
between 25 and 29 Sep. 2609.30257 shows where the boundary falls: submitted
Thu 24 Sep 17:59 UTC (13:59 EDT, one minute inside the Thursday-evening
batch), it was caught.

Evidence: five arXiv API queries run today (`abs:"independence polynomial"`;
`abs:"independent sets" AND abs:tree AND cat:math.CO`; `abs:unimodal AND
cat:math.CO`; `abs:"log-concave" AND cat:math.CO`; `abs:"hard-core
model"`), newest 40 results each since 10 Sep. Each result was compared
against `papers.json`, `evaluations.json` and `seen.json` (111 known IDs).
Misses whose timing fits the hole:

| arXiv | submitted (UTC) | title |
|---|---|---|
| 2609.20961 | Thu 17 Sep 18:16 | Unimodality of Independence Polynomials for Sufficiently Large Forests |
| 2609.21083 | Thu 17 Sep 20:59 | Unique Minimizers for Permanents, Mixed Discriminants, and Log-concave Polynomials |
| 2609.22078 | Fri 18 Sep 17:58 | Stellahedral geometry of partially ordered sets |
| 2609.25100 | Sat 19 Sep 22:25 | Latin Eulerian Numbers |
| 2609.31592 | Fri 25 Sep 17:48 | Riesz kernels of hyperbolic polynomials |
| 2609.13888 | Sat 12 Sep 11:24 | Independence polynomials and the weak Lefschetz property for tadpole graphs |
| 2609.13417 | Fri 11 Sep 18:28 | Prescribed Turán sign patterns for inverse Kazhdan-Lusztig polynomials of matroids |

Every paper in the sample submitted between Thursday 14:00 ET and Saturday
was missed.

Possible fix, not made: raise the default `--days` from 2 to 5. Four days
covers Thu 14:00 to the Monday-morning run (3.7 days); five also absorbs
deferred mailings and API index lag. Any backfill replay since 14 Aug has
to be chunked by week or have a higher cap. `MAX_RESULTS = 100` with a
newest-first break means a long window on a busy term like `log-concave`
truncates without warning. A mocked weekend-lag test belongs in
`test_arxiv_watch.py`.

A separate gap: `launchagent.err` holds 22 UNSCREENED and 26
partial-coverage warnings from arXiv 429s. They carry no dates, so I can't
say which days were lost. A 2-day window recovers a failed day only if
someone reruns it by hand. Sunday submissions (2609.33651, 2609.33732) were caught, as the
hypothesis predicts. Two other misses come from term or category coverage,
not the window: 2609.35332 (math.PR, "hard-core model") and 2609.18611
(Wed 16 Sep, independent sets of matroids). Out-of-scope categories (e.g.
2609.34952, SAT) are misses by design.

The two misses submitted inside this week's announcement window were read
through their first pages. Both are peripheral:

- **2609.31592, Wei**, Riesz kernels of hyperbolic polynomials. It answers
  a Scott–Sokal question on half-plane-property polynomials. Stable
  polynomials in several variables, with no independence-sequence content.
- **2609.35332, Briceño, Kosloff & Paviato**, CLTs for finitary factors of
  iid processes on Z^d, with the hard-core model in the uniqueness regime as
  one application. These are lattice random fields, not a CLT uniform over
  trees of the Fang et al. kind.

2609.13888 (tadpole graphs, weak Lefschetz property) is the earlier miss
most likely to bear on #993. It is untriaged; I read only the title.

## Operational notes

- **The watch has gone untriaged since 4 Sep.** None of the 22 papers in
  the 22–25 Sep entries, nor three sampled from 9–18 Sep (2609.07325,
  2609.14386, 2609.18118), appears in any project note, STATUS or
  DECISIONS. The check was a grep for each arXiv ID over
  `papers/queue/erdos-problem-993` and `Project-Management/harvest-review.md`.
  The last triage records are this note's predecessor (26 Aug) and
  `notes/arxiv_astra_transfer_2026-09-04.md`. arXiv:2609.20961 (Fang et al.)
  reached the project through erdosproblems.com. It appears in no watch
  file because of the coverage hole above.
- **Two of this week's papers declare AI-generated core arguments**:
  Leake–Mohammadi Yekta (ChatGPT 6 Astra) and Mikulincer–Zadik.
