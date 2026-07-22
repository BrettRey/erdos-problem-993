# Erdős Problem #993 — an open-problem work packet

You are being handed a famous open problem and a small, exact toolkit. The
framing below is deliberately light: pick your own angle. Read `known_results.md`
before you start so you build past what is already done rather than redo it.

## The problem

Let `T` be a tree on `n` vertices and let `i_k(T)` be the number of independent
sets of size `k` (an independent set is a vertex set with no two adjacent; the
empty set counts, so `i_0 = 1`). The **independence polynomial** is
`I(T; x) = sum_k i_k x^k`, and `alpha` is the largest `k` with `i_k > 0`.

> **Conjecture (Alavi, Malde, Schwenk, Erdős 1987).** For every tree `T`,
> the sequence `(i_0, i_1, ..., i_alpha)` is **unimodal**: it rises to a peak,
> then falls (ties allowed). Catalogued as Erdős problem #993. Open since 1987.

General graphs are not constrained this way — Alavi–Malde–Schwenk–Erdős showed
a general graph's independence sequence can be any positive sequence you like
(e.g. a 26-vertex graph realizes `1, 26, 15, 20, 15, 6, 1`, which dips and
rises). The conjecture is that *trees* are rigid enough to forbid that.

## Two ways to win — refutation is an explicit success mode

- **PROVE.** Unimodality for all trees; or a new nontrivial tree class; or a
  new reduction that provably shrinks the open case (e.g. discharges one of the
  conditional hypotheses in `known_results.md`).
- **REFUTE.** Exhibit one tree whose sequence dips then rises again (a
  **valley**). A single valid edge list ends a 39-year-old conjecture. This is
  not a consolation prize; a verified counterexample is the maximal outcome.

## Graded sub-targets (claim the highest you can actually reach)

- **T1.** Reproduce the state of the art: run `python3 erdos993_handoff.py`,
  understand the near-miss ratio, and state precisely where a counterexample
  must live (see the frontier facts in the embedded `DATA` dict).
- **T2.** A new *sufficient* structural condition under which a tree is unimodal,
  with proof, covering trees not already handled by the known families.
- **T3.** Prove or refute the tree-specific **mode ≤ ⌈mean⌉** localization, or
  the **Edge Contraction Mode Stability** property, either of which is a named
  conditional hinge in the current partial results (`known_results.md`).
- **T4.** Push the near-miss frontier: a tree (any size) with near-miss ratio
  strictly above the best known for its vertex count, *or* a proof that a whole
  family cannot exceed `1 - c/n` for an explicit `c` — pinning where the wall is.
- **T5.** Settle it: a full proof, or a verified valley.

## Binding verification rules

1. **Exact integer arithmetic only.** Independence-polynomial coefficients are
   integers; compute them exactly (the kit uses Python big ints). **Never** use
   floating-point root-finding or float coefficient arithmetic to decide
   unimodality or log-concavity. (In this project a claw-free control once
   produced a fabricated "result" from `numpy.roots` rounding noise; float roots
   are banned for any zero/shape claim. Use exact Sturm/`sympy` or certified Arb
   isolation if you go near polynomial zeros at all.)
2. **Theorem vs. evidence.** Label every claim `THEOREM` (proved, with the proof
   written out) or `EVIDENCE` (computational, with the exact range checked, e.g.
   "all trees on n ≤ N"). Do not blur the two. "Verified for n ≤ 20" is evidence,
   not a theorem.
3. **Effective constants.** If you assert a bound, give the constant.
4. **A counterexample must ship a certificate**: the exact edge list *and* the
   exact integer sequence it produces, reproducible by
   `check_counterexample(n, edges)` in `erdos993_handoff.py`. Nothing else counts.

## Self-grade your reply (put this line first)

Begin your response with one of:
`SOLVED` / `COUNTEREXAMPLE` / `PARTIAL` / `NO-PROGRESS`, then grade yourself
honestly. **An honest `PARTIAL` with one genuinely proved lemma is worth more
than a grand claim that dissolves when checked.** If you refute, lead with
`COUNTEREXAMPLE` and the certificate. If you only reproduced known facts, that
is `NO-PROGRESS` (T1) — say so; it is still useful.

## The kit

| File | What it is |
|------|------------|
| `erdos993_handoff.py` | **Everything in one file, std-lib only, no external data.** Exact engine — `indpoly(n, edges)`, `analyze(poly)`, `check_counterexample`, tree generators (`path`, `star`, `spider`, `broom`, `multiarm`, `caterpillar`, `random_tree`), `graph6_decode` — plus the embedded data (the two n=26 log-concavity breakers, the near-miss champions, the frontier facts) as a `DATA` dict, and this brief as its module docstring. Run `python3 erdos993_handoff.py` to self-test. |
| `README.md` | This brief. |
| `known_results.md` | Compressed map of what is proved, what is conditional, and what has been searched — read this first. |

If your interface only lets you upload one file, upload `erdos993_handoff.py`:
it is self-contained (the brief is its docstring, the data is embedded, the
self-test needs nothing else).

Full paper (methods, proofs, the conditional reductions): Brett Reynolds,
"Mean bounds, structural reductions, and exhaustive verification for tree
independence polynomial unimodality" (2026), doi:10.5281/zenodo.19100781.
Code repository: https://github.com/BrettRey/erdos-problem-993

Minimal example:

```python
from erdos993_handoff import multiarm, indpoly, analyze, check_counterexample
n, edges = multiarm(66, [6, 2])          # the n=75 near-miss champion
print(analyze(indpoly(n, edges)))         # near_miss ~ 0.9437, unimodal True
is_valley, poly, info = check_counterexample(n, edges)
```
