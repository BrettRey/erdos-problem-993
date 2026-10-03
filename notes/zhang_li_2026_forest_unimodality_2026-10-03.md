# Zhang & Li (2026): a claimed proof of #993 for all forests
<!-- SUMMARY: Reading note and local checks for Zhang-Li, Zenodo 10.5281/zenodo.22999166 (27 Sep 2026), which claims #993 for every finite forest (computer-assisted; arXiv submission withdrawn before announcement), and for Vallier's Lean formalization v1.0-claim; local Lean replay 2026-10-03 builds (8,709 jobs, no sorry), headline axioms exactly as claimed (standard three plus Lean.ofReduceBool/trustCompiler from native_decide); n<=60 part kernel-checked with standard axioms only · status: under local review · updated: 2026-10-03 -->

## The record

- **Paper.** Tong Zhang and Wei Li (Northwestern Polytechnical University),
  *Unimodality of Forest Independence Polynomials*, Zenodo record 22999166,
  DOI 10.5281/zenodo.22999166, dated 27 September 2026, 28,229 words of
  extracted text. Brett sent the link on 2026-10-03.
- **arXiv status.** The arXiv submission (`submit/8138169`) was withdrawn on 28
  September, before announcement, "to make its exposition clearer". So no
  arXiv watch could have seen it.
- **AI credit.** "OpenAI's GPT 6.0 Astra assisted with the computational work;
  the authors provided the overall proof strategy."
- **Local copies.** `literature/zhang_li_2026_unimodality_forest_independence_polynomials.{pdf,md}`.
  PDF SHA-256 `bd25c409…5b2f8`; the archived LaTeX source zip has SHA-256
  `649e88dd…452f206`. Both were fetched from the Zenodo API on 2026-10-03.

## What it claims (Theorem 1.1)

Every finite forest has a unimodal independence sequence. The proof has two
parts that overlap at orders 59 and 60.

- **Proposition 1.2, n ≤ 60.** Fix a maximum independent set I with complement
  Y. Every independent set is a choice of J ⊆ Y plus any subset of the free
  part of I, so I_G(z) = (1+z)^a + Σ t_{j,m} z^j (1+z)^m (§2.1). Structural
  counting inequalities place the counts t_{j,m} in a polyhedron. Exact
  rational certificates then separate every one of the 17,100 parameter
  triples (n, a, k) from a valley.
- **Proposition 1.3, n ≥ 59.** The ends of the sequence are monotone (Lemma
  3.1): increasing to ⌈n/4⌉, by an injection, and decreasing from
  ⌊(n+2a)/6⌋+1, by ratio domination against a maximum matching's polynomial
  (1+2x)^v (1+x)^{a−v}. Each intermediate rank is centred at an activity
  λ_k ∈ [1/3, 9/4) through the mean bound 3µ(2) ≥ a + 29n/60 (Lemma 3.2,
  Theorem 6.1). The remaining work is a conditional-binomial mixture
  K = Y + Bin(M, q) over a maximum-weight independent set B (eq. (41)),
  with exact kernel certificates. The verification code replays an inherited
  chain of thresholds, 900 → 500 → 350 → 170 → 130 → 99 → 80 → 59.

The 30 Sep attack found Fang's CLT architecture floored at N0 of 10^66 or
more even with idealized constants. That finding does not bear on this paper,
whose large-order route is different: it controls the conditional-binomial
mixture with exact certificates, not a Kolmogorov-rate CLT.

## Verification status: three separate things

| Part | What checks it | Status |
|---|---|---|
| n ≤ 60 (Zhang–Li Prop. 1.2) | Vallier's Lean, kernel-checked by `decide +kernel`, standard axioms; the authors' Python checker | **Replayed locally 2026-10-03**: all 24 kernel-check modules built; `forest_unimodal_of_card_le_sixty_kernel` and `certificatesSound_kernel` depend on `[propext, Classical.choice, Quot.sound]` only |
| n ≥ 61, Lean route | Vallier's own argument (O1–O6 on Zhang–Li's mixture), 271 `native_decide` certificate checks; "has not been independently reviewed" (README) | **Replayed locally 2026-10-03**: builds; `erdos993` depends on `[propext, Classical.choice, Lean.ofReduceBool, Lean.trustCompiler, Quot.sound]`. Accepted by Lean under compiler trust; the argument itself is unreviewed (by us as well) |
| n ≥ 59, the manuscript's own argument (§§3–8) | The authors' Python checker only; not formalized (the Lean project replaced it) | **Authors' checker run locally 2026-10-03**: `all_mathematical_acceptance_checks_passed: true`, `full_inherited_chain_replayed: true`. This is certificate acceptance only; the written reductions of §§3–8 are not checked by it and are not reviewed by us |

A green Lean build means Lean accepts Vallier's large-order argument under
compiler trust (`Lean.ofReduceBool`, `Lean.trustCompiler`). It is not a
review of that argument, and it says nothing about Zhang–Li's §§3–8.

## Lean formalization (Kevin Vallier, `github.com/selfreferencing/erdos993-lean`)

- **Release.** Tag `v1.0-claim` at commit
  `865e81498ecacda3de0d47647927353a7fadbba5` (29 Sep 2026). Lean v4.28.0,
  Mathlib `v4.28.0`. 641 `.lean` files, 124,970 lines.
- **Headline.** `Erdos993Lean.Analytic.erdos993 : Erdos993Statement`, with no
  hypotheses (`Analytic/Erdos993.lean`). It is the only definition of
  `Erdos993Statement`.
- **Statement faithfulness.** I read `Statement.lean` and the Mathlib sources
  directly. Forests are `SimpleGraph (Fin n)` with Mathlib's `IsAcyclic` (no
  cycle walks). `independenceCount F k` counts Mathlib's `indepSetFinset k`,
  the finsets of card k that are pairwise non-adjacent (`IsNIndepSet`). The
  claim is `UnimodalUpTo (indepNum) (independenceCount)`, a single peak
  through α, which suffices because i_k = 0 beyond α. The statement is
  faithful to #993 for all forests.
- **Pre-build safety scan.** The lakefile is declarative TOML with no
  scripts. In the non-`.lake` sources, `#eval`, `unsafe` and `implemented_by`
  appear only in comments. There are no `run_cmd`, `initialize`, `IO.`,
  `Process`, `@[extern]`, `opaque` or `partial def` hits.
- **Local replay (2026-10-03, this machine, Lean v4.28.0, Mathlib cache).**
  Driver `scratchpad/zhangli/replay.sh`.
  - `lake exe cache get` fetched the cache; `lake build Erdos993Lean` built
    the default library (8,063 jobs).
  - The 24 modules of `Erdos993Lean.ZhangKernel.Checks` were built one at a
    time. Every module succeeded. Each runs one `decide +kernel` theorem per
    (n, a); the kernel time for C45, recompiled directly with the profiler,
    was 10.6 s, and modules took up to about 66 s at n = 60.
  - `lake build Erdos993LeanTheorem` reported "Build completed successfully
    (8709 jobs)", with zero `declaration uses 'sorry'` warnings and zero
    error lines in the log.
  - `#print axioms`:
    - `Erdos993Lean.Analytic.erdos993`: `[propext, Classical.choice,
      Lean.ofReduceBool, Lean.trustCompiler, Quot.sound]`.
    - `Erdos993Lean.Zhang.forest_unimodal_of_card_le_sixty_kernel`:
      `[propext, Classical.choice, Quot.sound]`.
    - `Erdos993Lean.Zhang.certificatesSound_kernel`: `[propext,
      Classical.choice, Quot.sound]`.

    All three match the README exactly.

## The authors' verification code

- **Repository and contents.** `github.com/zhangzenozhang-jpg/forest-unimodality-arxiv-verification`
  at commit `806789352f7ed8e64594907a2ed334652ce11160`, directory
  `companion/2026-09-27/`. The entry point `anc/proof_reproduce.py` uses only
  the Python standard library and replays both the finite certificates and
  the large-order chain from zipped inputs.
- **Pre-run scan.** No network or deletion calls appear in any extracted
  script, two levels of nested zips included.
- **Run.** Python 3.12.9 (the authors recorded 3.12.14), `--workers 2`.
  - **Outcome.** Every stage passed:
    - the finite certificates (`forest_n60_extension_and_n100_gap`: verify,
      grid, countermodels);
    - `forest_threshold_59_release`: check_release, and the replay of the
      whole chain 900 → 500 → 350 → 170 → 130 → 99 → 80 → 59, 1,535 s.

    The summary reports `all_mathematical_acceptance_checks_passed: true`
    and `full_inherited_chain_replayed: true` in 1,565 s. The log ends "ALL
    PROOF-ONLY EXACT CERTIFICATE CHECKS PASSED".
  - **Caveat.** The entry point runs the original audit scripts with 13
    authenticated patches (`original_unmodified_full_entrypoint_used:
    false`). Each patch omits supplementary diagnostics or named-host
    cross-checks. One omits 378 named-host attachment cross-checks; others
    omit 16, 16 and 24. Each patch states that the certificate, coverage and
    kernel checks are kept. That is the authors' description; I did not
    audit the omitted code paths.
  - **Scope, in the authors' words.** "Exact certificate acceptance plus the
    written universal mathematical reductions; no forest counterexample
    search", and explicitly not a proof-assistant formalization.

## My own exact checks of the paper's lemmas

`/private/tmp/…/zhangli/killtest.py`, mirrored as
`scripts/killtest_zhang_li_lemmas_20261003.py`. The checks below cover every
tree on 2 to 17 vertices plus paths, spiders, caterpillars and the hub
killers (H(9,2), H(75,5), MSH(4×9,4×10)). Forests follow, because each
statement is additive or componentwise.

- **Lemma 3.1, ratio domination.** c_{k+1} q_k ≤ c_k q_{k+1} against
  Q = (1+2x)^v (1+x)^{a−v}: holds everywhere. The decreasing tail from
  U = ⌊(n+2a)/6⌋+1 also holds everywhere.
- **Theorem 6.1, (47).** µ(2) ≥ n/3 + 2c/15: holds, with equality exactly
  for K_2, as stated.
- **Theorem 6.1, (48).** 3µ(2) ≥ α + 29n/60: holds. Paths are the tight
  case, with slack about n/60, odd paths tighter.

This is consistency with the paper, not verification of it.

## What it means for the project

- **If it holds, #993 is solved.** Our work then becomes context: the manuscript's mean bound, the census to n ≤ 32, the depth-3 and central-window lanes, and the Vatter certificate.
- **Positioning.** The "if reactivated" positioning in STATUS now has to account for two papers: Fang et al. (large forests) and Zhang–Li (all forests).
- **No public action.** Nothing here is posted or claimed by us.
