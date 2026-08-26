# arXiv watch triage, 2026-08-26
<!-- SUMMARY: All 10 papers harvested by arxiv-watch since 2026-08-14 collected and read for #993 relevance; 1 relevant (Schweitzer, matroid-intersection ULC), 5 peripheral, 4 vocabulary coincidences · status: done · updated: 2026-08-26 -->

Covers every paper in `Project-Management/arxiv-watch/digest.md` from the
watch's first run (2026-08-14) through today. The 06:15 run this morning
failed wholesale (arXiv API returned 429/503 on all 21 queries); a manual
catch-up at 09:54 with `--days 3` closed the gap and found one more paper. All
ten PDFs were collected into `~/pdf-inbox/` with `.md` sidecars under
descriptive stems. Each verdict below is from reading the paper (first pages
at minimum), not the abstract alone.

## Relevant: 1

- **arXiv:2608.23262, Schweitzer, "A note on the ultra log-concavity of
  matroid intersection."** ULC for common independent sets of any two
  matroids; fails at three partition matroids. Maps the sharp boundary of the
  Lorentzian lane relative to graph/tree independence systems. Full hook with
  the graph reading of his counterexample and pre-use checks:
  `source-hooks/schweitzer-2026-matroid-intersection-ulc.md`. Also surfaces a
  companion worth one skim: arXiv:2601.02547 (Ardila-Mantilla, Cristancho,
  Denham, Eur, Huh, Wang), "Tree metrics and log-concavity for matroids,"
  absent from `gpt_attack/literature.md`.

## Peripheral: method-ecosystem awareness, no action

- **arXiv:2608.18519, Sinclair, "The radial derivative on the graded Möbius
  algebra."** Post-Rota-disproof response: a new matroid invariant $H_\beta$
  on the lattice of flats, conjectured log-concave and top-heavy where flat
  counts fail (cites the Divoux-Larson-Lowen-Wang non-unimodality result that
  motivated building this watch). Objects are matroid flats, not independence
  sequences. The method (a canonical degree-lowering operator certifying
  shape) is the Hodge-style pattern, but nothing transfers directly.
- **arXiv:2608.12266, Zhang, "Normalized skew Schur polynomials are
  Lorentzian"** and **arXiv:2608.13544, Le & Nguyen, "Skew Hives, Skew Skeps,
  Skew Schur Log-Concavity."** Symmetric-function positivity (skew Schur,
  Kostka, Littlewood-Richardson). Only thread to here: Li et al. 2025
  (arXiv:2501.04245) proved spider LC via symmetric functions, so growth in
  Schur log-concavity technology is worth knowing exists. No import now.
- **arXiv:2608.15682, Yan, "Variations of colored multiset Eulerian
  polynomials."** Real-rootedness via interlacing recurrences for
  Eulerian/Ehrhart $h^*$-polynomials. Tree independence polynomials are not
  real-rooted in general (the two $n=26$ LC-failing trees have complex
  roots), so the technique lives where #993 provably cannot.
- **arXiv:2608.13247, Menon, "Gamma-positivity for octopuses."** Despite the
  tree-shaped vocabulary (arbors, octopuses), the polynomial counts lattice
  points of Chapoton's arbor polytopes and is palindromic; gamma-positivity
  needs palindromicity, which independence polynomials lack. No transfer.

## Vocabulary coincidences: 4

- **arXiv:2608.15906, Gwóźdź.** Optimal transport; "log-concave" is a density
  on $\mathbb{R}^d$. The exact false-positive class the watch's ranking
  predicts (scored [3], non-core category).
- **arXiv:2608.17702, Ayyer & Prasad.** "Totally positive" matrices over
  finite fields; no order structure, so nothing for the TP2/Polya-frequency
  machinery used here.
- **arXiv:2608.19773, Zhang & Dong.** Chromatic polynomial vs list-color
  function equality persistence; coloring enumeration, no coefficient-shape
  content.
- **arXiv:2608.24552, Sanchez Galan.** Secretary-problem thresholds;
  "unimodality" is a hypothesis on payoff increments in optimal stopping.
  Scored [0] by the watch, correctly.

## Operational note

The watch's first hard failure was this morning: every query 429/503, likely
arXiv API rate-limiting or an outage at 06:15. The two-day window plus
seen-set dedup absorbed it exactly as designed once rerun by hand. If morning
429s recur, consider adding retry-with-backoff to `arxiv_watch.py` rather
than widening the window.

## Follow-up, same day (acting session)

- **Retry-with-backoff added to `arxiv_watch.py`** (delays 15 s / 60 s,
  Retry-After honored up to 120 s, retryable codes 429/500/502/503), plus the
  fix the morning actually needed: a wholesale failure now exits 1, warns on
  stderr, and posts a macOS notification under `--notify` instead of printing
  "nothing new". Partial failure warns that coverage is partial. Mocked tests
  cover all retry paths; a live `--replay --days 1` run afterwards completed
  all 22 queries cleanly (API recovered; everything already seen).
- **Schweitzer integrated into the manuscript.** Pre-use checks all run
  (derivation chain replayed by hand, witness replayed exactly via
  `scripts/verify_schweitzer_witness_20260826.py`, still v1 as of 14:35 UTC);
  two-sentence paragraph added to the `main_v2.tex` intro after the
  Bendjeddou-Hardiman paragraph, key `schweitzer2026`. Details in the source
  hook.
- **Companion arXiv:2601.02547 skimmed and filed** (PDF + md to
  `literature/`). Matroid-side confirmed: "tree metrics" = ultrametric
  distance matrices (Graham-Pollak refinement as PSD base case). Nothing for
  #993 beyond being the result Schweitzer strengthens. No further action.
