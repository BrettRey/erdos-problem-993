# Poisson–binomial novelty search: update to 2026-10-03
<!-- SUMMARY: Reran the July novelty audit with the cutoff moved to 2026-10-03; 31 titles new since July, none overlaps the first-descent theorem; the closest two were checked by abstract; authenticated MathSciNet still not run · status: done · updated: 2026-10-03 -->

Companion to `poisson_binomial_novelty_database_audit_2026-07-16.md`. That
audit's verdict, method and exclusions stand. This note extends its window
only.

- **Procedure.** `scripts/audit_poisson_binomial_databases.py` was run
  unchanged except for one line, `CUTOFF = "2026-10-03"`; the copy was kept
  in the session scratchpad. Output:
  `results/poisson_binomial_novelty_queries_20261003.json`, with exit 0 and
  the same queries, databases (zbMATH Open, OpenAlex, Semantic Scholar,
  arXiv, Crossref) and citation chains as in July.
- **Diff.** The July record has 692 distinct titles and the October record
  711. 31 titles appear in October that were absent in July. Almost all are
  off-topic: regression, item response theory, KV-cache eviction,
  changepoint detection, permutation families, Ehrhart theory and similar.
- **Closest two.** Each was read at the abstract:
  - *High-dimensional ultra-log-concave distributions* (Chen–Wang,
    arXiv:2609.23994). It concerns δ-ULC measures on ℕ^d and hard-core
    variance bounds at low fugacity. It gives no variance-scaled bound on the
    normalized Turán deficit at the first descent. It was already read in
    full in the arXiv-watch backlog triage of 2026-10-02.
  - *The sharp Rényi and Tsallis threshold in the Shepp–Olkin concavity
    problem* (arXiv:2609.27433, 23 Sep 2026). It concerns joint concavity of
    Rényi and Tsallis entropies of a Bernoulli sum in the parameter vector,
    using the Hillion–Johnson transport inequality. It is about entropy
    concavity, not local pmf curvature at the first descent. No overlap.
- **Verdict.** Unchanged. The manuscript's cautious sentence, "We are not
  aware of an earlier universal constant c > 0 …", remains the correct form.
  Do not strengthen it.
- **Still not done.** An authenticated MathSciNet search needs an
  institutional login (see the July note, which records the LibLynx
  redirect).
