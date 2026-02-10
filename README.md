# Erdős Problem #993

Computational search for a counterexample to the conjecture that the independent set sequence of every tree is unimodal.

## The problem

Given a tree *T* on *n* vertices, let *i_k(T)* = the number of independent sets of size *k* in *T*. The conjecture (Alavi, Malde, Schwenk, and Erdős 1987) states that the sequence *(i_0, i_1, ..., i_α)* is always unimodal: it increases to a peak and then decreases.

This is known to be false for general graphs. The tree case remains open.

- [erdosproblems.com #993](https://erdosproblems.com/993)

## Approach

1. Enumerate non-isomorphic trees on *n* vertices (using [nauty](https://pallini.di.uniroma1.it/)'s `geng` when available, otherwise custom generation)
2. Compute the independent set polynomial for each tree via dynamic programming on the rooted tree
3. Check whether the coefficient sequence is unimodal
4. Report any counterexample

## Status (manuscript)

The manuscript in `paper/main.tex` reports:
- Exhaustive verification for all trees with `n <= 26` (447,399,080 trees), with no unimodality violations and exactly two log-concavity failures at `n = 26`.
- Targeted search of 145,362 trees across five structured families up to `n = 500`, with no unimodality violations and best near-miss ratio `~0.9917` from brooms.
- Broom asymptotics consistent with `nm(s) = 1 - C/s + O(1/s^2)` and `C ~ 4.12` for `p = 13`.

Local result snapshots used by the manuscript:
- `results/analysis_n26.json` (n=26 log-concavity + near-miss analysis)
- `results/targeted_n500.json` (targeted search summary + top near-misses)

## Setup

```bash
pip install networkx numpy
```

Optional: install [nauty](https://pallini.di.uniroma1.it/) for fast tree enumeration (`brew install nauty` on macOS).

## Reproduce

```bash
python test_all.py
python search.py --max-n 26 --workers 8
python analyze.py 26 --workers 8
python targeted.py --max-n 500 --random-count 5000
python broom_asymptotic.py
```

## License

MIT
