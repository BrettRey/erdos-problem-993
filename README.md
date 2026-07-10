# Erdos Problem #993

Computational and structural investigation of the conjecture that the independence polynomial of every tree is unimodal.

## The problem

Given a tree *T* on *n* vertices, let *i_k(T)* = the number of independent sets of size *k*. The conjecture (Alavi, Malde, Schwenk, and Erdos 1987) states that the sequence *(i_0, i_1, ..., i_alpha)* is always unimodal: it increases to a peak and then decreases.

This is known to be false for general graphs. The tree case remains open.

- [erdosproblems.com #993](https://erdosproblems.com/993)

## Current public targets, July 2026

This repository does not contain a proof of Erdos Problem #993.

The current active targets are:

1. Signed-reserve/hub-bouquet route. The one-sided low-probability
   Poisson-binomial bound has now been lifted to every finite
   Poisson-binomial law: at a supported first descent,
   `V * Delta_eff >= 1/4`, hence the raw reserve has the same lower bound.
   This closes the general signed bridge and the product-term reserve for
   `A=(1+x)^s Q`. The active target is now perturbation by the hub-included
   term `xR`, especially for growing broom arms. This is not a general-tree
   reduction.

2. STP2/tree-DP route. `Formal/STP2Closure.lean` now guards the LC/STP2
   inequalities at `k >= 1` and records two abstract counterexamples showing
   that coefficient-shape hypotheses are too weak, even with contiguous
   support. The remaining target is to identify a genuine tree-DP realizability
   invariant.

3. Fixed-r certificate route. The abstract Route-2 Lean bridge is packaged; the
   next step is to emit a concrete `Route2SplitCertificateFor` or
   `Route2FamilyCertificate` instance from exact rational data.

4. Forest/product search route. Products of known non-log-concave tree
   polynomials should be searched using a direct valley/non-unimodality score
   rather than log-concavity defect.

5. Computation frontier. Exhaustive tree unimodality is verified through
   `n <= 29`; the analogous `n = 29` log-concavity / near-miss audit has not
   been completed.

## Results

The manuscript (`paper/main_v2.tex`) reports:

**Proved theorems:**
- **Subdivision-contraction identity:** I(T_e) = I(T) + x I(T/e) for any tree edge e
- **Conditional subdivision lemma:** If edge contraction shifts the mode by at most 1 (ECMS), then subdivision preserves unimodality, making any minimal counterexample homeomorphically irreducible
- **Mean bound:** mu(T) < n/3 for every d_leaf <= 1 tree on n >= 3 vertices
- **Conditional PNP framework:** Hub Exclusion + Transfer reduce the 1-Private maximal-IS bound to Conjecture A together with a separate Case-B hub bound; PNP itself does not prove unimodality
- **Edge bound:** P(u) + P(v) < 2/3 for every tree edge (hard-core model)
- **Leaf-attachment asymptotics:** nm(s) = 1 - C/s + O(1/s^2) with C in [4, 8)

**Computational verification:**
- Exhaustive: all 8,691,747,673 trees on n <= 29 are unimodal (0 violations)
- n = 29: 5,469,566,585 trees, 0 unimodality failures
- n = 28: 2,023,443,032 trees, 0 unimodality failures, 19 log-concavity failures (all at k = 14), best near-miss ratio 0.8565666
- n = 27: 751,065,460 trees, 0 unimodality failures, 0 log-concavity failures, best near-miss ratio 0.8571425
- ECMS verified for 24.7M edges (n <= 20), 0 violations
- Conjecture A verified for 931,596 trees (n <= 23), 0 violations
- Multi-arm stars identified as the extremal family (surpassing brooms)

## Setup

```bash
pip install networkx numpy
```

The July 2026 proof-audit harnesses additionally use:

```bash
pip install mpmath sympy
```

Optional: install [nauty](https://pallini.di.uniroma1.it/) for fast tree enumeration (`brew install nauty` on macOS).

## Reproduce

```bash
# Unit tests (55 tests)
python3 -m unittest test_all.py -v

# Exhaustive search, n <= 26 locally (8 workers, requires geng)
python3 search.py --max-n 26 --workers 8

# Modal cloud runs for larger n
modal run search_modal_exhaustive.py::dispatch --n 28 --workers 1024
modal run analyze_modal_lc_nm_n28.py::dispatch --n 28 --workers 1024 --top-k 200 --lc-top-k 200
python3 scripts/collect_modal_results.py collect --kind unimodality --n 28 --workers 1024 --expected 2023443032 --out results/analysis_n28_modal_unimodality.json
python3 scripts/collect_modal_results.py collect --kind lc_nm --n 28 --workers 1024 --top-k 200 --lc-top-k 200 --out results/analysis_n28_modal_lc_nm.json
modal run search_modal_exhaustive_n29.py::dispatch --n 29 --workers 1024
python3 scripts/collect_modal_results.py collect --kind unimodality --n 29 --workers 1024 --expected 5469566585 --out results/analysis_n29_modal_unimodality.json

# Targeted family search (n up to 500)
python3 targeted.py --max-n 500 --random-count 5000

# Evolutionary nm optimizer (multi-arm star search)
python3 nm_optimizer.py --min-n 50 --max-n 200

# Signed-reserve theorem audit harnesses
python3 scripts/verify_skellam_effective_drop.py
python3 scripts/verify_one_reflected_bernoulli.py
python3 scripts/verify_two_reflected_bernoulli.py
python3 scripts/verify_universal_pb_effective_drop.py
python3 scripts/verify_universal_pb_finite_bernstein.py
```

## Key files

| File | Description |
|------|-------------|
| `paper/main_v2.tex` | Current manuscript (XeLaTeX + biber) |
| `paper/main.tex` | Previous version |
| `indpoly.py` | Core DP + analysis functions |
| `search.py` | Exhaustive parallel search |
| `targeted.py` | Structured family search |
| `nm_optimizer.py` | Evolutionary near-miss optimizer |
| `test_all.py` | Unit + integration tests |
| `notes/` | Detailed analysis notes |
| `results/` | JSON result snapshots |

## Build the paper

```bash
cd paper
xelatex main_v2.tex && biber main_v2 && xelatex main_v2.tex && xelatex main_v2.tex
```

## License

MIT
