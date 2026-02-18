# Erdos Problem #993

Computational and structural investigation of the conjecture that the independence polynomial of every tree is unimodal.

## The problem

Given a tree *T* on *n* vertices, let *i_k(T)* = the number of independent sets of size *k*. The conjecture (Alavi, Malde, Schwenk, and Erdos 1987) states that the sequence *(i_0, i_1, ..., i_alpha)* is always unimodal: it increases to a peak and then decreases.

This is known to be false for general graphs. The tree case remains open.

- [erdosproblems.com #993](https://erdosproblems.com/993)

## Results

The manuscript (`paper/main_v2.tex`) reports:

**Proved theorems:**
- **Subdivision-contraction identity:** I(T_e) = I(T) + x I(T/e) for any tree edge e
- **Conditional subdivision lemma:** If edge contraction shifts the mode by at most 1 (ECMS), then subdivision preserves unimodality, making any minimal counterexample homeomorphically irreducible
- **PNP reduction:** Hub Exclusion + Transfer Lemmas reduce unimodality for all trees to Conjecture A about trees with d_leaf <= 1
- **Edge bound:** P(u) + P(v) < 2/3 for every tree edge (hard-core model)
- **Leaf-attachment asymptotics:** nm(s) = 1 - C/s + O(1/s^2) with C in [4, 8)

**Computational verification:**
- Exhaustive: all 447,672,596 trees on n <= 26 are unimodal (0 violations)
- ECMS verified for 24.7M edges (n <= 20), 0 violations
- Conjecture A verified for 528K trees (n <= 23), 0 violations
- Multi-arm stars identified as the extremal family (surpassing brooms)

## Setup

```bash
pip install networkx numpy
```

Optional: install [nauty](https://pallini.di.uniroma1.it/) for fast tree enumeration (`brew install nauty` on macOS).

## Reproduce

```bash
# Unit tests (37 tests)
python3 -m unittest test_all.py -v

# Exhaustive search, n <= 26 (8 workers, requires geng)
python3 search.py --max-n 26 --workers 8

# Targeted family search (n up to 500)
python3 targeted.py --max-n 500 --random-count 5000

# Evolutionary nm optimizer (multi-arm star search)
python3 nm_optimizer.py --min-n 50 --max-n 200
```

## Key files

| File | Description |
|------|-------------|
| `paper/main_v2.tex` | Current manuscript (11pp, XeLaTeX + biber) |
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
