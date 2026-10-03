# arXiv watch backlog triage, 28 Aug to 2 Oct 2026
<!-- SUMMARY: 124 untriaged papers (68 from digest entries 28 Aug-25 Sep, 56 found by the 52-day backfill after the watch fix), read for #993 relevance: 2 relevant (Zhang-Tu spider stability, certified on 66,272 spiders; Vatter's subsequence LC, whose mechanism gives an LC certificate found on every tree n<=13), 62 peripheral, 60 coincidences · status: done · updated: 2026-10-02 -->

Companion to `notes/arxiv-watch-triage-2026-10-02.md`, which covered the
29 Sep to 2 Oct entries. This note covers everything else the watch had
queued but nobody had read, plus what the repaired watch found when replayed.

## Scope

- **68 digest papers with no project record.** These are all IDs in the digest
  entries dated 28 Aug to 25 Sep that appear nowhere in `notes/`, STATUS,
  DECISIONS or `paper/` (grep by arXiv ID). The 26 Aug triage covered entries
  up to 26 Aug. The 4 Sep audit discussed the 1 to 4 Sep entries, which are
  excluded here.
- **56 backfill papers.** The 52-day replay after the watch fix (DECISIONS
  2026-10-02) queued 58 papers the old 2-day window had never kept. Two were
  already read (2609.20961 Fang et al.; 2609.31592 in this week's note),
  leaving 56.
- **Total: 124.**

## Method

- All 124 PDFs were fetched from arxiv.org and converted with `pdftotext`.
  None failed, and none came out empty.
- **I read 2 myself in full:** 2609.04694 and 2609.13888, the two
  independence-polynomial papers that scored 8 in the watch.
- **Seven read-only Claude Sonnet agents read the other 122,** 17 or 18 each,
  under one brief (`runs/arxiv-triage-20261002/BRIEF.md`). Every verdict had to
  quote the paper and report keyword counts. Their returns are saved verbatim
  as `RETURN_batch_{1..7}.md`.
- **I read every paper an agent marked "relevant?"** at the cited passages, and
  checked where needed.
- **Second-family panel:** qwen3.8:27b, local, excerpt-only. It ran first on
  all 122 and was stopped after 14 because it took about 70 s per paper. It
  was then run on the 40 papers whose text mentions independent sets, the
  hard-core model, trees or forests. Results are in
  `runs/arxiv-triage-20261002/local_panel_qwen3.8-27b.jsonl`; the comparison is
  below.
- **Paid route unused:** OpenRouter returned HTTP 402 (insufficient credits),
  so the paid second-opinion route was not used.

## Relevant: 2

- **arXiv:2609.04694, Zhang & Tu, "Stability of independence polynomials of
  spiders"** (4 Sep).
  - **What it proves.** Theorem 4: every independence root of every spider
    S(l_1, ..., l_d) has strictly negative real part. This answers part of
    Brown and Cameron's question about which trees are stable. Brown and
    Cameron proved that stars are stable, built trees with roots arbitrarily
    far into the right half-plane, and found all trees of order at most 20
    computationally stable. The proof is analytic: a homotopy J_t, the
    argument principle, and a technical estimate on path-polynomial ratios.
  - **My check.** `scripts/verify_spider_stability_2609_04694_20261002.py`
    builds I(S) exactly and isolates every root with Arb (python-flint). Over
    all 66,272 spiders on 2 to 35 vertices it found no counterexample and no
    undecided ball. The largest certified upper bound on Re(root) was
    −0.0639, for the star K_{1,34}. The spider polynomial was cross-checked
    against the tree DP on S(1,2,3).
  - **What it means here.** It is a roots-lane result. Left-half-plane
    stability does not imply unimodality, since a quadratic factor
    z^2 + bz + c with 0 < b^2 < c is stable and not log-concave. Spider
    log-concavity is already known (Li et al., symmetric functions), so it
    adds nothing to the unimodality question. It does sharpen what is known
    about where tree independence roots can sit.
  - Filed as `literature/zhang_tu_2026_stability_independence_polynomials_spiders.{pdf,md}`
    (PDF SHA-256 `fb1940a9…6c4d82`).
- **arXiv:2608.22147, Vatter, "Log-concavity of subsequence counts of words"**
  (23 Aug; one page; flagged by the batch-3 agent).
  - **What it proves.** A short proof of Chase's theorem. Classify
    subsequences by first letter, so ρ_{k+1}(w) is a weighted average of the
    ρ_k(tails), with ρ_k = count_k / count_{k−1}. Claim 2 shows that a suffix
    never has a larger ratio, ρ_k(suffix) ≤ ρ_k(word).
  - **The transfer to trees** (my construction, not in the paper). Fix a
    vertex order and classify independent sets by least vertex. Then
    i_{k+1}(T) = Σ_v i_k(G_v), where G_v is induced on the later
    non-neighbours of v. So T is log-concave whenever some order has
    ρ_k(G_v) ≤ ρ_k(T) for every v and k with positive weight. Call that a
    *Vatter order*.
  - **Probe.** `scripts/probe_vatter_order_certificate_20261002.py` runs an
    exact DFS over the sets of vertices still to come. **Every tree on 2 to 13
    vertices has a Vatter order** (2,287 trees).
  - **Sanity checks.** The decomposition identity held exactly on 600 random
    orders of the two order-26 LC-failing trees. None of those orders passed,
    as it must not, since a Vatter order would imply log-concavity.
  - **Status: a lead, not a result.** The certificate is strictly stronger
    than log-concavity, so it must fail somewhere at or below order 26. What
    matters for the open moderate-n range is a *window* version, ρ_k(G_v) ≤
    ρ_k(T) only for k in the central window [⌈n/4⌉, q], which would certify
    central log-concavity. That is untested. The project's ratio-dominance
    work (`check_ratio_dominance.py` and related notes) compares E and J at
    support vertices through Karlin TP2. It is a neighbour of this certificate,
    not the same thing. The LP discharging certificates of 30 Sep redistribute
    defect between neighbours, also a different mechanism.
  - Filed as `literature/vatter_2026_log_concavity_subsequence_counts.{pdf,md}`
    (PDF SHA-256 `4efdcc9b…7304d9`).

## Read in full or at the cited passages, downgraded to peripheral: 6

- **2609.13888, Hoa, Phuoc & Son, tadpole graphs** (12 Sep). Their unimodality
  criterion (Theorem 3.2) says that if G and Tail(G,v,1) are unimodal with
  modes 0 or 1 apart, every pendant path at v keeps unimodality. It is the
  adjacent-mode sum lemma already in `main_v2.tex:865`, applied repeatedly
  along a path, and it duplicates `notes/degree2_chain_peeling.md`. Nothing new
  for #993. The rest of the paper (WLP of A(T_{m,n})) is commutative algebra.
- **2609.11589, Liu & Mao, log-concavity under Hadamard products.** If
  W(p) and W(q) are LC-NIZ, so is W(pq), via TP2/RR2 kernels; this answers
  Brändén–Ferroni–Jochemko Question 6.1. The operation is the Hadamard product
  of polynomial sequences (Ehrhart series of Cartesian products, Segre
  products). No tree or forest operation corresponds to it.
- **2609.23994, Chen & Wang, high-dimensional ULC.** §3.6 bounds Var|I| for
  the hard-core model on Δ-regular graphs under 1 + λ^{-1} ≥ −λ_min(A_G).
  Remark 3.21 gives δ-ULC under a low-fugacity spectral condition. For a star
  K_{1,n} that condition means λ ≲ 1/√n, far below the fugacity where the
  central window lives. It says nothing about the size sequence i_k.
- **2609.04654, Davies, random independent sets and local sparsity.** Lower
  bounds on Z_G(λ) and on vertex marginals for locally sparse graphs, trees
  included (triangle-free). These are free-energy bounds with no
  coefficient-shape content.
- **2609.28983, antichain polynomials of [k] × P.** Theorem 2.1 expands
  N_P(x) over antichains B outside a maximum antichain L as
  x^{|B|}(1 + x)^{d − |Γ_L(B)|}. For a tree (a comparability graph) this is the
  defect expansion around a maximum independent set, the bookkeeping of the
  project's depth-3 defect classes. I did not check exact equivalence. The
  non-unimodal family of Prop. 5.2 is a clique plus isolated vertices, not a
  tree.
- **2609.23779, Mao, Shi & Zhu, Almkvist's conjecture.** Unimodality of
  ∏_k (1 + q^k + … + q^{(r−1)k}). The method is exact rational inequalities for
  n ≤ 11, then Fourier inversion with local estimates near roots of unity for
  n ≥ 12. That two-regime pattern is generic, and the paper is about a product
  of palindromic factors, not independence sequences.

## Tally

| verdict | agents (122) | after parent reads | final (124) |
|---|---|---|---|
| relevant | 0 | 1 (Vatter) | 2 |
| relevant? | 6 | 0 | 0 |
| peripheral | 56 | 61 | 62 |
| coincidence | 60 | 60 | 60 |

The agents' nearest misses were 2609.25038 (polymer cluster expansion with a
Scott–Sokal independent-set link; zero-free polydisc only), 2609.07728 (loopy
polynomial; determines independence polynomials, no shape content),
2609.10238 (dually Lorentzian ⇔ LC-NIZ, which excludes exactly the hard
cases) and 2609.23694 (Mu–Welker counterexamples). I left them peripheral as
the agents recommended.

## Second-family panel

The qwen3.8:27b panel classified 54 papers from excerpts, with no errors: the
first 14 of all 122, plus all 40 that mention independent sets, the hard-core
model, trees or forests. It agreed with the final verdicts on 45. **No
disagreement surfaced a relevant paper the reading agents and I had missed.**

- **Panel more generous, overridden (3).**
  - 2609.04654 (Davies): panel "relevant", final peripheral. Read; free-energy
    lower bounds only.
  - 2609.23994 (high-dimensional ULC): panel "relevant?", final peripheral.
    Read §3.6; it holds only at low fugacity.
  - 2609.05341 (bounded ratios for Lorentzian polynomials): panel "relevant?",
    final peripheral. Its hypothesis is that the polynomial is Lorentzian.
    Tree independence polynomials are not Lorentzian in general (the
    Schweitzer boundary), and nothing in the paper connects to independence
    systems.
- **Panel more cautious, overridden (1).** 2608.22147 (Vatter): panel
  "coincidence", final relevant. The excerpt shows only words and
  subsequences. The relevance comes from the transplanted mechanism and the
  probe above, which an excerpt cannot show.
- **Peripheral vs coincidence (5),** where the two families draw the boundary
  differently: 2608.19780, 2608.28484, 2609.11323 and 2609.25038 (panel
  coincidence, agents peripheral), and 2609.29796 (panel peripheral, agent
  coincidence). Left as the agents had them. Neither verdict changes anything
  for #993.

Coverage limit: 68 of the 122 agent-read papers were not seen by the panel.
They are the ones with no independent-set, hard-core, tree or forest mention
anywhere in the text, so the risk of a missed relevant paper there is low.

## Notes for the record

- **Five papers declare AI-generated mathematics**: 2608.21507 (an LLM found
  the non-unimodal h*-vector), 2608.13836 (Brenti counterexample, AI tools),
  2609.07964 (Nakano counterexamples, with AI), 2609.24542 (generalized Lax
  conjecture, AI-generated proof) and 2609.39917 (Astra, this week's batch).
- **Log-concavity fails while unimodality survives** in several neighbouring
  objects (Bruhat up-degree 2609.13764; genus distributions 2609.27628;
  inverse KL polynomials 2609.13417; Brenti's R̃-polynomials 2608.13836). It is
  the same pattern as Kadrawi–Levit for trees, and none of these papers offers
  a mechanism that transfers.
